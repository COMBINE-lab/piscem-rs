#!/usr/bin/env python3
"""Gate scATAC's fine calibration publisher against the generic cadence.

Both arms use one release binary and the same fixed decode split, so the broker
does not run and allocation cannot confound the result.  The only difference is
``PISCEM_SCATAC_PROGRESS_FLUSH_EVERY=64|256``. Runs are position-balanced within
paired blocks, canonical RAD content and read counts must agree, and the formal
gate is a one-sided bootstrap 95% upper bound of 1% for mapping wall, process
wall, and aggregate process CPU.
"""

import argparse
import hashlib
import json
import math
import os
import random
import statistics
import subprocess
import sys
import time
from pathlib import Path


FLUSH_ENV = "PISCEM_SCATAC_PROGRESS_FLUSH_EVERY"
SPLIT_ENV = "PISCEM_DECODE_SLOTS"
ARMS = {"fine-64": "64", "generic-256": "256"}


def canonical_digest(binary, rad):
    digest = hashlib.sha256()
    process = subprocess.Popen([str(binary), "atac", str(rad)], stdout=subprocess.PIPE)
    while True:
        chunk = process.stdout.read(1024 * 1024)
        if not chunk:
            break
        digest.update(chunk)
    code = process.wait()
    if code:
        raise subprocess.CalledProcessError(code, process.args)
    return digest.hexdigest()


def percentile(values, fraction):
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(len(ordered) * fraction))]


def stratified_ratio_estimate(ratios_by_position):
    """Combine equal crossover positions without reintroducing order bias."""
    log_medians = [
        statistics.median(math.log(value) for value in values)
        for values in ratios_by_position.values()
    ]
    return math.exp(statistics.mean(log_medians))


def stratified_ratio_upper(ratios_by_position, seed, iterations=10000):
    """Bootstrap within each randomized position stratum, preserving 15/15."""
    rng = random.Random(seed)
    estimates = []
    for _ in range(iterations):
        sampled = {
            position: [rng.choice(values) for _ in values]
            for position, values in ratios_by_position.items()
        }
        estimates.append(stratified_ratio_estimate(sampled))
    return percentile(estimates, 0.95)


def run_one(args, root, arm, repetition, position):
    output = root / arm / ("rep-%02d" % repetition)
    output.parent.mkdir(parents=True, exist_ok=True)
    log_path = output.parent / (output.name + ".run.log")
    time_path = output.parent / (output.name + ".time.json")
    argv = [
        str(args.binary),
        "map-scatac",
        "-i",
        str(args.index),
        "-1",
        str(args.read1),
        "-b",
        str(args.barcode),
        "-2",
        str(args.read2),
        "--quiet",
        "--dict",
        args.dictionary,
        "-t",
        str(args.threads),
        "--decoder",
        "auto",
        "-o",
        str(output),
    ]
    env = os.environ.copy()
    for name in ("PISCEM_RAPIDGZIP_THREADS", SPLIT_ENV, FLUSH_ENV):
        env.pop(name, None)
    env[SPLIT_ENV] = str(args.decode_slots)
    env[FLUSH_ENV] = ARMS[arm]
    timed = [
        "/usr/bin/time",
        "-f",
        '{"user_seconds":%U,"system_seconds":%S,"peak_rss_kib":%M,"wall_seconds":%e}',
        "-o",
        str(time_path),
        *argv,
    ]
    started = time.monotonic()
    with log_path.open("wb") as log:
        result = subprocess.run(timed, env=env, stdout=log, stderr=subprocess.STDOUT)
    record = {
        "arm": arm,
        "flush_every": int(ARMS[arm]),
        "repetition": repetition,
        "block_position": position,
        "command": argv,
        "exit_code": result.returncode,
        "runner_wall_seconds": time.monotonic() - started,
    }
    if result.returncode:
        record["error"] = "command failed; see %s" % log_path
        return record
    try:
        info = json.loads((output / "map_info.json").read_text())
        process = json.loads(time_path.read_text())
        plan = info["execution_plan"]
        if plan["allocation"]["kind"] != "pinned_aggregate":
            raise ValueError("run was not aggregate-pinned: %r" % plan)
        if plan["decode_slots"] != args.decode_slots or info.get("thread_broker") is not None:
            raise ValueError("fixed split was not honored: %r" % plan)
        record.update(
            {
                "mapping_seconds": info["mapping_seconds"],
                "process_wall_seconds": process["wall_seconds"],
                "process_cpu_seconds": process["user_seconds"] + process["system_seconds"],
                "peak_rss_kib": process["peak_rss_kib"],
                "num_reads": info["num_reads"],
                "num_mapped": info["num_mapped"],
                "canonical_output_sha256": canonical_digest(
                    args.canonical_rad, output / "map.rad"
                ),
            }
        )
        return record
    finally:
        # Keep compact timing, logs, and map-info; discard multi-megabyte output.
        for name in ("map.rad", "unmapped_bc_count.bin"):
            path = output / name
            if path.exists():
                path.unlink()


def summarize(records, args):
    valid = [record for record in records if record.get("exit_code") == 0]
    cells = {
        arm: {record["repetition"]: record for record in valid if record["arm"] == arm}
        for arm in ARMS
    }
    complete = all(len(cell) == args.repetitions for cell in cells.values())
    digests = {record["canonical_output_sha256"] for record in valid}
    counts = {(record["num_reads"], record["num_mapped"]) for record in valid}
    first_position_counts = {
        arm: sum(record.get("block_position") == 0 for record in cells[arm].values())
        for arm in ARMS
    }
    summary = {
        "repetitions": args.repetitions,
        "complete": complete,
        "content_equal": len(digests) == 1 and len(counts) == 1,
        "first_position_counts": first_position_counts,
        "comparisons": {},
        "failures": [],
    }
    metrics = ("mapping_seconds", "process_wall_seconds", "process_cpu_seconds")
    if complete:
        for index, metric in enumerate(metrics):
            pairs = [
                (
                    cells["generic-256"][repetition][metric],
                    cells["fine-64"][repetition][metric],
                )
                for repetition in range(args.repetitions)
            ]
            ratios = [fine / generic for generic, fine in pairs]
            ratios_by_position = {position: [] for position in range(len(ARMS))}
            for repetition, ratio in enumerate(ratios):
                position = cells["fine-64"][repetition]["block_position"]
                ratios_by_position[position].append(ratio)
            summary["comparisons"][metric] = {
                "generic_median": statistics.median(left for left, _ in pairs),
                "fine_median": statistics.median(right for _, right in pairs),
                "paired_delta_median": statistics.median(right - left for left, right in pairs),
                "paired_ratio_median": statistics.median(ratios),
                "position_stratified_ratio": stratified_ratio_estimate(
                    ratios_by_position
                ),
                "paired_ratio_upper_95": stratified_ratio_upper(
                    ratios_by_position, args.seed + index * 101
                ),
            }
    if args.repetitions < 30:
        summary["failures"].append("formal gate requires at least 30 repetitions")
    if not complete:
        summary["failures"].append("one or more paired runs are missing")
    if not summary["content_equal"]:
        summary["failures"].append("canonical output or read counts differ")
    if max(first_position_counts.values(), default=0) - min(
        first_position_counts.values(), default=0
    ) > 1:
        summary["failures"].append("block positions are not counterbalanced")
    for metric, comparison in summary["comparisons"].items():
        if comparison["paired_ratio_upper_95"] > 1.01:
            summary["failures"].append(
                "%s upper ratio %.6f exceeded 1.010000"
                % (metric, comparison["paired_ratio_upper_95"])
            )
    summary["gate_passed"] = not summary["failures"]
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=Path("target/release/piscem-rs"))
    parser.add_argument(
        "--canonical-rad", type=Path, default=Path("target/release/examples/canonical_rad")
    )
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--read1", type=Path, required=True)
    parser.add_argument("--barcode", type=Path, required=True)
    parser.add_argument("--read2", type=Path, required=True)
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--decode-slots", type=int, default=6)
    parser.add_argument("--dictionary", default="auto")
    parser.add_argument("--repetitions", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20260807)
    parser.add_argument("--summarize-only", action="store_true")
    args = parser.parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    records_path = args.results_dir / "records.jsonl"
    if args.summarize_only:
        records = [json.loads(line) for line in records_path.read_text().splitlines()]
    else:
        rng = random.Random(args.seed)
        arm_names = list(ARMS)
        first_arms = [arm_names[i % len(arm_names)] for i in range(args.repetitions)]
        rng.shuffle(first_arms)
        records = []
        with records_path.open("w") as stream:
            for repetition, first_arm in enumerate(first_arms):
                order = [first_arm, *(arm for arm in arm_names if arm != first_arm)]
                for position, arm in enumerate(order):
                    print(
                        "[%d/%d] %s rep %d"
                        % (len(records) + 1, args.repetitions * 2, arm, repetition),
                        flush=True,
                    )
                    try:
                        record = run_one(args, args.results_dir, arm, repetition, position)
                    except Exception as error:
                        record = {
                            "arm": arm,
                            "repetition": repetition,
                            "block_position": position,
                            "exit_code": -1,
                            "error": repr(error),
                        }
                    records.append(record)
                    stream.write(json.dumps(record, sort_keys=True) + "\n")
                    stream.flush()
    summary = summarize(records, args)
    (args.results_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["gate_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
