#!/usr/bin/env python3
"""Gate policy-scoped scATAC reader/progress geometry on a long real run.

Serial and aggregate-pinned policies compare their coarse scoped defaults
(16384/256) with the former adaptive-only settings (1024/64). Adaptive compares
the inverse. Every paired block is position balanced and must preserve canonical
RAD content, counts, the requested policy, and serialized tuning telemetry.
"""

import argparse
import hashlib
import json
import math
import os
import random
import statistics
import subprocess
import time
from pathlib import Path


BATCH_ENV = "PISCEM_SCATAC_READER_BATCH_SIZE"
PROGRESS_ENV = "PISCEM_SCATAC_PROGRESS_FLUSH_EVERY"
SPLIT_ENV = "PISCEM_DECODE_SLOTS"
POLICY_ENV = "PISCEM_THREAD_BROKER_POLICY"
PROBE_ENV = "PISCEM_THREAD_BROKER_PROBE_INTERVAL_MS"
POLICIES = ("serial", "pinned", "adaptive")
ARMS = ("scoped", "opposite")


def canonical_digest(binary, rad):
    digest = hashlib.sha256()
    process = subprocess.Popen([str(binary), "atac", str(rad)], stdout=subprocess.PIPE)
    for block in iter(lambda: process.stdout.read(1024 * 1024), b""):
        digest.update(block)
    code = process.wait()
    if code:
        raise subprocess.CalledProcessError(code, process.args)
    return digest.hexdigest()


def percentile(values, fraction):
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(len(ordered) * fraction))]


def stratified_estimate(by_position):
    return math.exp(
        statistics.mean(
            statistics.median(math.log(value) for value in values)
            for values in by_position.values()
        )
    )


def stratified_upper(by_position, seed, iterations=10000):
    rng = random.Random(seed)
    estimates = []
    for _ in range(iterations):
        sample = {
            position: [rng.choice(values) for _ in values]
            for position, values in by_position.items()
        }
        estimates.append(stratified_estimate(sample))
    return percentile(estimates, 0.95)


def expected_geometry(policy, arm, fixed_batch):
    scoped = {
        "serial": (1024, 256),
        "pinned": (1024, 256),
        "adaptive": (1024, 64),
    }[policy]
    opposite = {
        "serial": (fixed_batch, 256),
        "pinned": (fixed_batch, 256),
        "adaptive": (16384, 256),
    }[policy]
    return scoped if arm == "scoped" else opposite


def duration_seconds(value):
    if value is None:
        return None
    if isinstance(value, dict):
        return value.get("secs", 0) + value.get("nanos", 0) / 1e9
    return float(value)


def run_one(args, policy, arm, repetition, position):
    output = args.results_dir / policy / arm / ("rep-%02d" % repetition)
    output.parent.mkdir(parents=True, exist_ok=True)
    log_path = output.parent / (output.name + ".run.log")
    time_path = output.parent / (output.name + ".time.json")
    decoder = "serial" if policy == "serial" else "auto"
    command = [
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
        decoder,
        "-o",
        str(output),
    ]
    environment = os.environ.copy()
    for name in (
        "PISCEM_RAPIDGZIP_THREADS",
        BATCH_ENV,
        PROGRESS_ENV,
        SPLIT_ENV,
        POLICY_ENV,
        PROBE_ENV,
    ):
        environment.pop(name, None)
    if policy == "pinned":
        environment[SPLIT_ENV] = str(args.decode_slots)
    if arm == "opposite":
        batch, progress = expected_geometry(policy, arm, args.coarse_batch)
        environment[BATCH_ENV] = str(batch)
        environment[PROGRESS_ENV] = str(progress)

    timed = [
        "/usr/bin/time",
        "-f",
        '{"user_seconds":%U,"system_seconds":%S,"peak_rss_kib":%M,"wall_seconds":%e}',
        "-o",
        str(time_path),
    ] + command
    started = time.monotonic()
    with log_path.open("wb") as log:
        result = subprocess.run(timed, env=environment, stdout=log, stderr=subprocess.STDOUT)
    record = {
        "policy": policy,
        "arm": arm,
        "repetition": repetition,
        "block_position": position,
        "command": command,
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
        tuning = info["pipeline_tuning"]
        expected_batch, expected_progress = expected_geometry(
            policy, arm, args.coarse_batch
        )
        if (tuning["reader_batch_size"], tuning["progress_flush_every"]) != (
            expected_batch,
            expected_progress,
        ):
            raise ValueError("unexpected geometry: %r" % tuning)
        expected_allocation = {
            "serial": "serial",
            "pinned": "pinned_aggregate",
            "adaptive": "adaptive",
        }[policy]
        if plan["allocation"]["kind"] != expected_allocation:
            raise ValueError("unexpected plan: %r" % plan)
        report = info.get("thread_broker")
        failure = info.get("thread_broker_failure")
        if policy == "adaptive" and arm == "scoped":
            if report is None or failure is not None:
                raise ValueError("adaptive broker did not complete cleanly: %r" % failure)
        elif policy != "adaptive" and (report is not None or failure is not None):
            raise ValueError("fixed policy unexpectedly ran the broker")
        record.update(
            {
                "mapping_seconds": info["mapping_seconds"],
                "process_wall_seconds": process["wall_seconds"],
                "process_cpu_seconds": process["user_seconds"] + process["system_seconds"],
                "peak_rss_kib": process["peak_rss_kib"],
                "num_reads": info["num_reads"],
                "num_mapped": info["num_mapped"],
                "pipeline_tuning": tuning,
                "execution_plan": plan,
                "broker": report,
                "broker_failure": failure,
                "canonical_output_sha256": canonical_digest(
                    args.canonical_rad, output / "map.rad"
                ),
            }
        )
        return record
    finally:
        for name in ("map.rad", "unmapped_bc_count.bin"):
            path = output / name
            if path.exists():
                path.unlink()


def summarize(records, args, policies):
    valid = [record for record in records if record.get("exit_code") == 0]
    summary = {
        "repetitions": args.repetitions,
        "complete": len(valid) == args.repetitions * len(policies) * len(ARMS),
        "policies": {},
        "failures": [],
    }
    digests = {record["canonical_output_sha256"] for record in valid}
    counts = {(record["num_reads"], record["num_mapped"]) for record in valid}
    summary["content_equal"] = len(digests) == 1 and len(counts) == 1
    metrics = (
        "mapping_seconds",
        "process_wall_seconds",
        "process_cpu_seconds",
        "peak_rss_kib",
    )
    for policy_index, policy in enumerate(policies):
        cells = {
            arm: {
                record["repetition"]: record
                for record in valid
                if record["policy"] == policy and record["arm"] == arm
            }
            for arm in ARMS
        }
        policy_summary = {"comparisons": {}, "failures": []}
        complete = all(len(cell) == args.repetitions for cell in cells.values())
        policy_summary["complete"] = complete
        policy_summary["scoped_first"] = sum(
            record["block_position"] == 0 for record in cells["scoped"].values()
        )
        if complete:
            for metric_index, metric in enumerate(metrics):
                ratios = []
                by_position = {0: [], 1: []}
                for repetition in range(args.repetitions):
                    scoped = cells["scoped"][repetition]
                    opposite = cells["opposite"][repetition]
                    ratio = scoped[metric] / opposite[metric]
                    ratios.append(ratio)
                    by_position[scoped["block_position"]].append(ratio)
                comparison = {
                    "scoped_median": statistics.median(
                        record[metric] for record in cells["scoped"].values()
                    ),
                    "opposite_median": statistics.median(
                        record[metric] for record in cells["opposite"].values()
                    ),
                    "paired_ratio_median": statistics.median(ratios),
                    "position_stratified_ratio": stratified_estimate(by_position),
                    "paired_ratio_upper_95": stratified_upper(
                        by_position, args.seed + policy_index * 101 + metric_index
                    ),
                }
                policy_summary["comparisons"][metric] = comparison
                upper_limit = (
                    (1.10 if policy == "serial" else 1.05)
                    if metric == "peak_rss_kib"
                    else (
                        1.25
                        if policy == "serial" and args.short_run
                        else (
                            1.15
                            if policy == "serial"
                            else (1.10 if policy == "adaptive" else 1.05)
                        )
                    )
                )
                if comparison["paired_ratio_upper_95"] > upper_limit:
                    policy_summary["failures"].append(
                        "%s scoped/opposite upper %.6f exceeded %.6f"
                        % (metric, comparison["paired_ratio_upper_95"], upper_limit)
                    )
                if (
                    policy == "serial"
                    and metric != "peak_rss_kib"
                    and comparison["paired_ratio_median"] > (1.20 if args.short_run else 1.05)
                ):
                    policy_summary["failures"].append(
                        "%s serial scoped/opposite median %.6f exceeded %.6f"
                        % (
                            metric,
                            comparison["paired_ratio_median"],
                            1.20 if args.short_run else 1.05,
                        )
                    )
            if policy == "adaptive":
                for record in cells["scoped"].values():
                    report = record["broker"]
                    converge = duration_seconds(report.get("time_to_converge"))
                    if converge is None or (
                        converge > 5.0 and converge > 0.20 * record["mapping_seconds"]
                    ):
                        policy_summary["failures"].append(
                            "scoped adaptive run exceeded both 5 s absolute and 20% fractional convergence"
                        )
                        break
            if policy == "serial":
                scoped_rss = statistics.median(
                    record["peak_rss_kib"] for record in cells["scoped"].values()
                )
                opposite_rss = statistics.median(
                    record["peak_rss_kib"] for record in cells["opposite"].values()
                )
                if scoped_rss - opposite_rss > 128 * 1024:
                    policy_summary["failures"].append(
                        "serial scoped median RSS increase exceeded 128 MiB"
                    )
        else:
            policy_summary["failures"].append("one or more paired runs are missing")
        if abs(policy_summary["scoped_first"] * 2 - args.repetitions) > 1:
            policy_summary["failures"].append("block positions are not balanced")
        summary["policies"][policy] = policy_summary
        summary["failures"].extend(
            "%s: %s" % (policy, failure) for failure in policy_summary["failures"]
        )
    if not summary["complete"]:
        summary["failures"].append("matrix is incomplete")
    if not summary["content_equal"]:
        summary["failures"].append("canonical output or counts differ")
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
    parser.add_argument("--coarse-batch", type=int, default=16384)
    parser.add_argument("--policies", default=",".join(POLICIES))
    parser.add_argument("--dictionary", default="sshash")
    parser.add_argument("--repetitions", type=int, default=6)
    parser.add_argument("--seed", type=int, default=20260807)
    parser.add_argument("--summarize-only", action="store_true")
    parser.add_argument(
        "--short-run",
        action="store_true",
        help="use the documented <=20%% serial startup trade gate; long runs require no median regression",
    )
    args = parser.parse_args()
    policies = tuple(value.strip() for value in args.policies.split(",") if value.strip())
    if not policies or any(policy not in POLICIES for policy in policies):
        parser.error("--policies must be a comma-separated subset of %s" % ",".join(POLICIES))
    args.results_dir.mkdir(parents=True, exist_ok=True)
    records_path = args.results_dir / "records.jsonl"
    if args.summarize_only:
        records = [json.loads(line) for line in records_path.read_text().splitlines()]
    else:
        rng = random.Random(args.seed)
        records = []
        with records_path.open("w") as stream:
            for policy in policies:
                first_arms = [ARMS[i % 2] for i in range(args.repetitions)]
                rng.shuffle(first_arms)
                for repetition, first in enumerate(first_arms):
                    order = (first, ARMS[1] if first == ARMS[0] else ARMS[0])
                    for position, arm in enumerate(order):
                        print(
                            "[%d/%d] %s/%s rep %d"
                            % (
                                len(records) + 1,
                                args.repetitions * len(policies) * len(ARMS),
                                policy,
                                arm,
                                repetition,
                            ),
                            flush=True,
                        )
                        try:
                            record = run_one(args, policy, arm, repetition, position)
                        except Exception as error:
                            record = {
                                "policy": policy,
                                "arm": arm,
                                "repetition": repetition,
                                "block_position": position,
                                "exit_code": -1,
                                "error": repr(error),
                            }
                        records.append(record)
                        stream.write(json.dumps(record, sort_keys=True) + "\n")
                        stream.flush()
    summary = summarize(records, args, policies)
    (args.results_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["gate_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
