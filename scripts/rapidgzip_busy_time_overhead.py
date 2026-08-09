#!/usr/bin/env python3
"""Paired no-measurable-overhead gate for rapidgzip busy-time accounting."""

import argparse
import hashlib
import json
import os
import random
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def percentile(values, fraction):
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(len(ordered) * fraction))]


def bootstrap_upper(ratios, seed, iterations=10000):
    rng = random.Random(seed)
    estimates = []
    for _ in range(iterations):
        estimates.append(statistics.median(rng.choice(ratios) for _ in ratios))
    return percentile(estimates, 0.95)


def run(binary, argv, repetition, variant, order, output_dir):
    with tempfile.NamedTemporaryFile(
        prefix="rapidgzip-busy-time-", suffix=".json", dir=str(output_dir), delete=False
    ) as timing:
        timing_path = Path(timing.name)
    command = [
        "/usr/bin/time",
        "-f",
        '{"user_seconds":%U,"system_seconds":%S,"wall_seconds":%e,"peak_rss_kib":%M}',
        "-o",
        str(timing_path),
        str(binary),
        *argv,
    ]
    started = time.monotonic()
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    record = {
        "variant": variant,
        "repetition": repetition,
        "block_order": order,
        "command": command,
        "exit_code": result.returncode,
        "runner_wall_seconds": time.monotonic() - started,
        "stdout": result.stdout.decode("utf-8", errors="replace"),
        "stderr": result.stderr.decode("utf-8", errors="replace"),
    }
    if result.returncode == 0:
        record["process"] = json.loads(timing_path.read_text())
    timing_path.unlink()
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--accounting", required=True, type=Path)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--decoded-bytes", required=True, type=int)
    parser.add_argument("--members", default=1, type=int)
    parser.add_argument("--iterations", default=8, type=int)
    parser.add_argument("--repetitions", default=30, type=int)
    parser.add_argument("--workers", default=8, type=int)
    parser.add_argument("--seed", default=973, type=int)
    parser.add_argument("--margin", default=0.01, type=float)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summarize-existing", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    argv = [
        str(args.input), "1", str(args.workers), str(args.workers), str(1024 * 1024),
        str(args.decoded_bytes), str(args.members), str(args.iterations), "private", "read",
    ]
    binaries = {"baseline": args.baseline, "accounting": args.accounting}
    provenance = {
        "binary_sha256": {key: digest(path) for key, path in binaries.items()},
        "input_sha256": digest(args.input),
        "argv": argv,
    }
    rng = random.Random(args.seed)
    runs_path = args.output / "runs.jsonl"
    if args.summarize_existing:
        records = [json.loads(line) for line in runs_path.read_text().splitlines()]
    else:
        records = []
        for repetition in range(1, args.repetitions + 1):
            order = ["baseline", "accounting"]
            rng.shuffle(order)
            for variant in order:
                record = run(
                    binaries[variant], argv, repetition, variant, order, args.output
                )
                record.update(provenance)
                records.append(record)
                with runs_path.open("a") as stream:
                    stream.write(json.dumps(record, sort_keys=True) + "\n")
                if record["exit_code"] != 0:
                    raise SystemExit("failed run: %s" % json.dumps(record, indent=2))

    # The benchmark intentionally prints elapsed time and throughput. The exact
    # decoded byte/member lines are the canonical semantic output.
    outputs = {
        tuple(
            line
            for line in record["stdout"].splitlines()
            if line.startswith(("decoded_bytes\t", "member_count_per_decoder\t"))
        )
        for record in records
    }
    by_variant = {
        variant: {
            record["repetition"]: record
            for record in records
            if record["variant"] == variant
        }
        for variant in binaries
    }
    summary = {
        "repetitions": args.repetitions,
        "margin": args.margin,
        "stdout_equal": len(outputs) == 1,
        "metrics": {},
    }
    for metric in ("wall_seconds", "cpu_seconds"):
        ratios = []
        for repetition in range(1, args.repetitions + 1):
            baseline = by_variant["baseline"][repetition]["process"]
            accounting = by_variant["accounting"][repetition]["process"]
            if metric == "cpu_seconds":
                baseline_value = baseline["user_seconds"] + baseline["system_seconds"]
                accounting_value = accounting["user_seconds"] + accounting["system_seconds"]
            else:
                baseline_value = baseline[metric]
                accounting_value = accounting[metric]
            ratios.append(accounting_value / baseline_value)
        median = statistics.median(ratios)
        upper = bootstrap_upper(ratios, args.seed + len(summary["metrics"]))
        summary["metrics"][metric] = {
            "median_ratio": median,
            "ratio_upper_95": upper,
            "passed": upper <= 1.0 + args.margin,
        }
    summary["gate_passed"] = summary["stdout_equal"] and all(
        metric["passed"] for metric in summary["metrics"].values()
    )
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    if not summary["gate_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
