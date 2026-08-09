#!/usr/bin/env python3
"""Run and analyze real-input Gate E producer-measurement traces."""

import argparse
import bisect
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path


CADENCES_MS = (25, 50, 73, 100, 137, 250)


def digest_command(command, output):
    if not command:
        return None
    argv = [part.replace("{output}", str(output)) for part in command]
    value = hashlib.sha256()
    process = subprocess.Popen(argv, stdout=subprocess.PIPE)
    for chunk in iter(lambda: process.stdout.read(1024 * 1024), b""):
        value.update(chunk)
    if process.wait() != 0:
        raise subprocess.CalledProcessError(process.returncode, argv)
    return value.hexdigest()


def output_paths(root, job, repetition):
    output = root / job["name"] / ("rep-%d" % repetition)
    if job["output_kind"] == "stem":
        output.parent.mkdir(parents=True, exist_ok=True)
        return output, Path("%s.map_info.json" % output)
    output.mkdir(parents=True, exist_ok=True)
    return output, output / "map_info.json"


def execute(manifest, job, repetition):
    root = Path(manifest["results_dir"])
    output, map_info = output_paths(root, job, repetition)
    trace = Path("%s.busy-trace.json" % map_info)
    timing = Path("%s.time.json" % map_info)
    log = Path("%s.run.log" % map_info)
    argv = [
        manifest["binary"], *job["args"], "-t", str(job["budget"]),
        "--decoder", job.get("decoder", "auto"),
        job.get("output_argument", "-o"), str(output),
    ]
    env = os.environ.copy()
    for name in (
        "PISCEM_RAPIDGZIP_THREADS", "PISCEM_THREAD_BROKER_POLICY",
        "PISCEM_THREAD_BROKER_PROBE_INTERVAL_MS",
        "PISCEM_FIXED_DECODE_MEASUREMENT",
        "PISCEM_FIXED_DECODE_MEASUREMENT_TRACE",
        "PISCEM_DECODE_SLOTS",
    ):
        env.pop(name, None)
    if not job.get("sequential", False):
        env["PISCEM_DECODE_SLOTS"] = str(job["pin"])
        env["PISCEM_FIXED_DECODE_MEASUREMENT"] = "calibration"
        env["PISCEM_FIXED_DECODE_MEASUREMENT_TRACE"] = str(trace)
    timed = [
        "/usr/bin/time", "-f",
        '{"user_seconds":%U,"system_seconds":%S,"wall_seconds":%e,"peak_rss_kib":%M}',
        "-o", str(timing), *argv,
    ]
    started = time.monotonic()
    with log.open("wb") as stream:
        result = subprocess.run(timed, env=env, stdout=stream, stderr=subprocess.STDOUT)
    record = {
        "job": job["name"], "repetition": repetition, "command": argv,
        "exit_code": result.returncode,
        "runner_wall_seconds": time.monotonic() - started,
    }
    if result.returncode != 0:
        record["error"] = "command failed; see %s" % log
        return record
    info = json.loads(map_info.read_text())
    plan = info["execution_plan"]
    if plan["requested_budget"] != job["budget"]:
        raise ValueError("budget mismatch for %s" % job["name"])
    if job.get("sequential", False):
        if plan["decode_slots"] != 0 or info.get("producer_measurement") is not None:
            raise ValueError("sequential cell emitted producer work")
    else:
        if plan["allocation"]["kind"] != "pinned_aggregate" or plan["decode_slots"] != job["pin"]:
            raise ValueError("measurement cell was not pinned: %r" % plan)
        record["trace"] = json.loads(trace.read_text())
    record["map_info"] = info
    record["process"] = json.loads(timing.read_text())
    record["canonical_output_sha256"] = digest_command(
        job.get("canonical_digest_command"), output
    )
    if job.get("cleanup_output", False):
        payload = Path("%s.rad" % output) if job["output_kind"] == "stem" else output / "map.rad"
        payload.unlink()
    return record


def percentile(values, fraction):
    values = sorted(values)
    return values[min(len(values) - 1, int(len(values) * fraction))]


def trace_windows(trace, cadence_ms, signal_field, reference_field):
    points = trace["points"]
    elapsed = [point["elapsed_nanos"] for point in points]
    cadence = cadence_ms * 1000000
    windows = []
    # Five deterministic phase offsets prevent a favorable alignment from
    # hiding cadence aliasing. Each result spans exactly three controller
    # windows, matching the production cost estimator.
    for phase_fraction in range(5):
        cursor = cadence * phase_fraction // 5
        while cursor + 3 * cadence <= elapsed[-1]:
            left = bisect.bisect_left(elapsed, cursor)
            right = bisect.bisect_left(elapsed, cursor + 3 * cadence)
            if right >= len(points):
                break
            sampled = points[right][signal_field] - points[left][signal_field]
            reference = points[right][reference_field] - points[left][reference_field]
            if reference > 0:
                windows.append((sampled, reference))
            cursor += cadence
    return windows


def analyze_record(record, job):
    if job.get("sequential", False):
        return {"sequential_zero_signal": True, "passed": True}
    trace = record["trace"]
    points = trace["points"]
    if len(points) < 4:
        return {"passed": False, "failure": "fewer than four trace points"}
    signal_field = (
        "accounted_busy_nanos"
        if job.get("signal", "sampled") == "accounted"
        else "sampled_busy_nanos"
    )
    reference_field = (
        "sampled_busy_nanos"
        if job.get("reference", "accounted") == "sampled"
        else "accounted_busy_nanos"
    )
    sampled_total = points[-1][signal_field] - points[0][signal_field]
    reference_total = points[-1][reference_field] - points[0][reference_field]
    whole_error = abs(sampled_total - reference_total) / reference_total
    measurement = record["map_info"]["producer_measurement"]
    if measurement.get("accounted_busy_nanos") != points[-1]["accounted_busy_nanos"]:
        raise ValueError("trace/map-info exact counter mismatch")

    consumer = record["map_info"]["consumer_measurement"]
    consumer_total = sum(
        consumer.get(name, 0)
        for name in ("busy_nanos", "callback_setup_nanos", "output_flush_nanos")
    )
    consumer_ratio = consumer_total / max(1, reference_total)
    cadence = {}
    all_errors = []
    median_shares = []
    for cadence_ms in CADENCES_MS:
        windows = trace_windows(trace, cadence_ms, signal_field, reference_field)
        errors = [abs(sampled - reference) / reference for sampled, reference in windows]
        shares = [
            sampled / (sampled + reference * consumer_ratio)
            for sampled, reference in windows
            if sampled + reference * consumer_ratio > 0
        ]
        exact_share = 1.0 / (1.0 + consumer_ratio)
        median_share = percentile(shares, 0.5) if shares else None
        cadence[str(cadence_ms)] = {
            "windows": len(windows),
            "p95_three_window_relative_error": percentile(errors, 0.95) if errors else None,
            "median_solved_producer_share": median_share,
            "exact_solved_producer_share": exact_share,
        }
        all_errors.extend(errors)
        if median_share is not None:
            median_shares.append(median_share)

    worker_cpu = measurement.get("completed_worker_cpu_nanos")
    auxiliary_cpu = measurement.get("completed_auxiliary_cpu_nanos")
    producer_coverage = None
    if worker_cpu is not None and auxiliary_cpu is not None:
        producer_coverage = (
            measurement["accounted_busy_nanos"] + auxiliary_cpu
        ) / max(1, worker_cpu + auxiliary_cpu)
    excluded_consumer_fraction = (
        consumer.get("callback_setup_nanos", 0) + consumer.get("output_flush_nanos", 0)
    ) / max(1, consumer_total)
    process_cpu = (
        record["process"]["user_seconds"] + record["process"]["system_seconds"]
    ) * 1000000000
    accounted_process = consumer_total
    for name in (
        "completed_worker_cpu_nanos", "completed_auxiliary_cpu_nanos",
        "sampler_cpu_nanos",
    ):
        accounted_process += measurement.get(name) or 0
    process_coverage = accounted_process / max(1, process_cpu)
    p95 = percentile(all_errors, 0.95) if all_errors else float("inf")
    share_drift = max(median_shares) - min(median_shares) if median_shares else float("inf")
    expected_paths = set(job.get("expected_decoder_paths", []))
    observed_paths = set(trace.get("decoder_paths", []))
    path_passed = not expected_paths or bool(expected_paths & observed_paths)
    passed = all((
        whole_error <= 0.03,
        p95 <= 0.10,
        share_drift <= 0.02,
        producer_coverage is not None and producer_coverage >= 0.95,
        excluded_consumer_fraction < 0.02 or job.get("model_excluded_consumer", False),
        path_passed,
    ))
    return {
        "whole_run_relative_error": whole_error,
        "production_signal": job.get("signal", "sampled"),
        "reference_signal": job.get("reference", "accounted"),
        "p95_three_window_relative_error": p95,
        "cadence_solved_share_drift": share_drift,
        "producer_component_coverage": producer_coverage,
        "whole_process_component_coverage": process_coverage,
        "excluded_consumer_fraction": excluded_consumer_fraction,
        "decoder_paths": sorted(observed_paths),
        "path_passed": path_passed,
        "cadences": cadence,
        "passed": passed,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--summarize-existing", action="store_true")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    root = Path(manifest["results_dir"])
    root.mkdir(parents=True, exist_ok=True)
    runs_path = root / "runs.jsonl"
    if args.summarize_existing:
        records = [json.loads(line) for line in runs_path.read_text().splitlines()]
    else:
        records = []
        for repetition in range(1, manifest.get("repetitions", 1) + 1):
            for job in manifest["jobs"]:
                record = execute(manifest, job, repetition)
                records.append(record)
                with runs_path.open("a") as stream:
                    stream.write(json.dumps(record, sort_keys=True) + "\n")
                if record["exit_code"] != 0:
                    raise SystemExit(record["error"])
    jobs = {job["name"]: job for job in manifest["jobs"]}
    analyses = [
        {"job": record["job"], "repetition": record["repetition"],
         **analyze_record(record, jobs[record["job"]])}
        for record in records if record.get("exit_code") == 0
    ]
    digest_sets = {
        name: {record.get("canonical_output_sha256") for record in records if record["job"] == name}
        for name in jobs
    }
    summary = {
        "gate_passed": all(item["passed"] for item in analyses)
        and all(len(values) == 1 for values in digest_sets.values()),
        "cadences_ms": CADENCES_MS,
        "analyses": analyses,
        "canonical_digest_sets": {name: sorted(values) for name, values in digest_sets.items()},
    }
    (root / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    if not summary["gate_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
