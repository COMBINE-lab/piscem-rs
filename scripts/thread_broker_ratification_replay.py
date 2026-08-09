#!/usr/bin/env python3
"""Replay production ratification against variance from real broker windows."""

import argparse
import json
import math
import random
import statistics
from pathlib import Path


def estimate(values):
    mean = sum(values) / len(values)
    if len(values) > 1:
        variance = sum((value - mean) ** 2 for value in values) / (len(values) - 1)
        half_width = 1.96 * math.sqrt(variance) / math.sqrt(len(values))
    else:
        half_width = float("inf")
    return mean, half_width


def compare(baseline_values, achieved_values, tolerance=0.05):
    baseline, baseline_width = estimate(baseline_values)
    achieved, achieved_width = estimate(achieved_values)
    if baseline <= 0:
        return "kept"
    if len(achieved_values) > 1 and achieved_values.count(0.0) == 1:
        return "inconclusive"
    boundary = baseline * (1.0 - tolerance)
    difference_width = (1.645 / 1.96) * math.hypot(
        baseline_width * (1.0 - tolerance), achieved_width
    )
    if achieved + difference_width < boundary:
        return "regressed"
    if achieved - difference_width < boundary:
        return "inconclusive"
    return "kept"


def real_blocks(paths, samples):
    traces = []
    sources = []
    for path in paths:
        for line in path.read_text().splitlines():
            record = json.loads(line)
            if record.get("exit_code") != 0 or "map_info" not in record:
                continue
            rates = record["map_info"].get("thread_broker", {}).get("throughput_trace", [])
            positive = [rate for rate in rates if math.isfinite(rate) and rate > 0]
            if len(positive) < 4:
                continue
            center = statistics.median(positive)
            normalized = [rate / center for rate in positive]
            if len(normalized) >= samples * 2:
                traces.append(normalized)
            sources.append({
                "path": str(path), "mode": record.get("mode"),
                "cadence_ms": record.get("cadence_ms"), "windows": len(positive),
            })
    residuals = [value for trace in traces for value in trace]
    if len(residuals) < 16 or not traces:
        raise ValueError("fewer than 16 positive real throughput windows")
    return traces, residuals, sources


def centered(values, target):
    mean = sum(values) / len(values)
    return [target * value / mean for value in values]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("runs", nargs="+", type=Path)
    parser.add_argument("--traces", default=1000, type=int)
    parser.add_argument("--samples", default=4, type=int)
    parser.add_argument("--seed", default=20260807, type=int)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    traces, residuals, sources = real_blocks(args.runs, args.samples)
    false_equal = 0
    detected_regression = 0
    false_zero = 0
    for seed in range(args.seed, args.seed + args.traces):
        rng = random.Random(seed)
        trace = rng.choice(traces)
        start = rng.randrange(0, len(trace) - 2 * args.samples + 1)
        baseline_shape = trace[start:start + args.samples]
        achieved_shape = trace[start + args.samples:start + 2 * args.samples]
        baseline = centered(baseline_shape, 100.0)
        equal = centered(achieved_shape, 100.0)
        regressed = centered(achieved_shape, 90.0)
        isolated_zero = equal[:]
        isolated_zero[seed % args.samples] = 0.0
        false_equal += compare(baseline, equal) == "regressed"
        detected_regression += compare(baseline, regressed) == "regressed"
        false_zero += compare(baseline, isolated_zero) == "regressed"
    mean = sum(residuals) / len(residuals)
    coefficient_of_variation = statistics.stdev(residuals) / mean
    summary = {
        "traces": args.traces,
        "ratification_samples": args.samples,
        "real_windows": len(residuals),
        "real_coefficient_of_variation": coefficient_of_variation,
        "sources": sources,
        "false_confident_equal_rate": false_equal / args.traces,
        "true_ten_percent_regression_detection": detected_regression / args.traces,
        "isolated_zero_false_rejection": false_zero / args.traces,
    }
    summary["gate_passed"] = (
        summary["false_confident_equal_rate"] < 0.01
        and summary["true_ten_percent_regression_detection"] >= 0.95
        and summary["isolated_zero_false_rejection"] == 0
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    if not summary["gate_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
