#!/usr/bin/env python3
"""Real alternating-compressibility/lull validation for Gate G."""

import argparse
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path


def seconds(duration):
    return duration.get("secs", 0) + duration.get("nanos", 0) / 1e9


def canonical_digest(command, output):
    argv = [part.replace("{output}", str(output)) for part in command]
    value = hashlib.sha256()
    process = subprocess.Popen(argv, stdout=subprocess.PIPE)
    for chunk in iter(lambda: process.stdout.read(1024 * 1024), b""):
        value.update(chunk)
    if process.wait() != 0:
        raise subprocess.CalledProcessError(process.returncode, argv)
    return value.hexdigest()


def execute(manifest, cadence, lulls, repetition):
    mode = "lulls" if lulls else "control"
    root = Path(manifest["results_dir"]) / ("%dms" % cadence) / mode
    root.mkdir(parents=True, exist_ok=True)
    output = root / ("rep-%d" % repetition)
    directory_output = manifest.get("output_kind", "stem") == "directory"
    info_path = output / "map_info.json" if directory_output else Path("%s.map_info.json" % output)
    sidecar = root / ("rep-%d" % repetition)
    time_path = Path("%s.time.json" % sidecar)
    log_path = Path("%s.run.log" % sidecar)
    lull_path = Path("%s.lulls.jsonl" % sidecar)
    argv = [
        manifest["binary"], *manifest["args"], "-t", str(manifest["budget"]),
        "--decoder", "auto", "-o", str(output),
    ]
    env = os.environ.copy()
    for name in (
        "PISCEM_DECODE_SLOTS", "PISCEM_FIXED_DECODE_MEASUREMENT",
        "PISCEM_FIXED_DECODE_MEASUREMENT_TRACE", "PISCEM_THREAD_BROKER_POLICY",
        "PISCEM_THREAD_BROKER_PROBE_INTERVAL_MS", "PISCEM_THREAD_BROKER_SAMPLE_INTERVAL_MS",
        "PISCEM_VALIDATION_READ_LULL", "PISCEM_VALIDATION_READ_LULL_TRACE",
        "PISCEM_VALIDATION_READ_LULL_MATCH",
    ):
        env.pop(name, None)
    env["PISCEM_THREAD_BROKER_SAMPLE_INTERVAL_MS"] = str(cadence)
    if lulls:
        env["PISCEM_VALIDATION_READ_LULL"] = manifest["lull_schedule"]
        env["PISCEM_VALIDATION_READ_LULL_TRACE"] = str(lull_path)
        if manifest.get("lull_path_match"):
            env["PISCEM_VALIDATION_READ_LULL_MATCH"] = manifest["lull_path_match"]
    timed = [
        "/usr/bin/time", "-f",
        '{"user_seconds":%U,"system_seconds":%S,"wall_seconds":%e,"peak_rss_kib":%M}',
        "-o", str(time_path), *argv,
    ]
    started = time.monotonic()
    with log_path.open("wb") as log:
        result = subprocess.run(timed, env=env, stdout=log, stderr=subprocess.STDOUT)
    record = {
        "cadence_ms": cadence, "mode": mode, "repetition": repetition,
        "command": argv, "exit_code": result.returncode,
        "runner_wall_seconds": time.monotonic() - started,
        "binary_sha256": manifest["_binary_sha256"],
    }
    if result.returncode != 0:
        record["error"] = "command failed; see %s" % log_path
        return record
    info = json.loads(info_path.read_text())
    report = info.get("thread_broker")
    if report is None or info.get("thread_broker_failure") is not None:
        raise ValueError("broker did not complete cleanly: %r" % info.get("thread_broker_failure"))
    plan = info["execution_plan"]
    if plan["allocation"]["kind"] != "adaptive" or plan["effective_budget"] != manifest["budget"]:
        raise ValueError("unexpected execution plan: %r" % plan)
    measurement = report["producer_measurement"]
    if (
        measurement.get("accounted_busy_nanos") != measurement.get("busy_nanos")
        or measurement.get("calibration_samples") != 0
        or measurement.get("monitoring_samples") != 0
        or measurement.get("observation_nanos") != 0
        or measurement.get("sampler_cpu_nanos") is not None
    ):
        raise ValueError("adaptive run did not use the sampler-free native signal")
    record["map_info"] = info
    record["process"] = json.loads(time_path.read_text())
    record["lulls"] = (
        [json.loads(line) for line in lull_path.read_text().splitlines()] if lulls else []
    )
    record["canonical_output_sha256"] = canonical_digest(
        manifest["canonical_digest_command"], output
    )
    rad_path = output / "map.rad" if directory_output else Path("%s.rad" % output)
    rad_path.unlink()
    return record


def trajectory_state(report, at_seconds):
    elapsed = [seconds(value) for value in report["trajectory_elapsed"]]
    producer = report["producer_trajectory"]
    current = producer[0]
    for when, value in zip(elapsed, producer):
        if when > at_seconds:
            break
        current = value
    return current


def exact_whole_run_share(record):
    info = record["map_info"]
    producer = info["producer_measurement"]
    producer_busy = producer.get("accounted_busy_nanos")
    consumer_busy = info["consumer_measurement"]["busy_nanos"]
    if producer_busy is None or producer_busy + consumer_busy == 0:
        raise ValueError("missing exact nonzero whole-run stage counters")
    return producer_busy / (producer_busy + consumer_busy)


def recovery_delay(report, lull_end, required):
    if trajectory_state(report, lull_end) >= required:
        return 0.0
    for when, producer in zip(
        [seconds(value) for value in report["trajectory_elapsed"]],
        report["producer_trajectory"],
    ):
        if when >= lull_end and producer >= required:
            return when - lull_end
    return float("inf")


def first_reaches(report, required):
    for when, producer in zip(
        [seconds(value) for value in report["trajectory_elapsed"]],
        report["producer_trajectory"],
    ):
        if producer >= required:
            return when
    return float("inf")


def cap_delay(report, lull_start, source_ceiling):
    observations = report.get("cap_trajectory", [])
    after = [
        seconds(observation["elapsed"]) - lull_start
        for observation in observations
        if seconds(observation["elapsed"]) >= lull_start
        and observation["reason"] != "none"
        and observation["useful_cap"] <= source_ceiling + 1
    ]
    return min(after) if after else float("inf")


def attributed_lull(lulls, onset):
    candidates = [
        event for event in lulls
        if event["started_nanos"] / 1e9 <= onset
        <= event["ended_nanos"] / 1e9 + 1.0
    ]
    return candidates[-1] if candidates else None


def cap_retention_intervals(report, lulls):
    observations = report.get("cap_trajectory", [])
    active_since = None
    intervals = []
    for observation in observations:
        when = seconds(observation["elapsed"])
        if observation["reason"] == "none":
            if active_since is not None:
                event = attributed_lull(lulls, active_since)
                if event is not None:
                    intervals.append(when - active_since)
                active_since = None
        elif active_since is None:
            active_since = when
    if active_since is not None and attributed_lull(lulls, active_since) is not None:
        intervals.append(float("inf"))
    return intervals


def harmful_cap_clear_delays(report, lulls, required):
    active_since = None
    delays = []
    for observation in report.get("cap_trajectory", []):
        when = seconds(observation["elapsed"])
        harmful = (
            observation["reason"] != "none"
            and observation["useful_cap"] < required
        )
        if harmful and active_since is None:
            active_since = when
        elif not harmful and active_since is not None:
            event = attributed_lull(lulls, active_since)
            if event is not None:
                delays.append(max(0.0, when - event["ended_nanos"] / 1e9))
            active_since = None
    if active_since is not None and attributed_lull(lulls, active_since) is not None:
        delays.append(float("inf"))
    return delays


def summarize(records, manifest):
    analyses = []
    persistent_source = manifest.get("persistent_source", False)
    require_lull_cap = manifest.get("require_lull_cap", not persistent_source)
    for cadence in manifest["cadences_ms"]:
        for repetition in range(1, manifest["repetitions"] + 1):
            control = next(
                record for record in records
                if record["cadence_ms"] == cadence and record["mode"] == "control"
                and record["repetition"] == repetition
            )
            lull = next(
                record for record in records
                if record["cadence_ms"] == cadence and record["mode"] == "lulls"
                and record["repetition"] == repetition
            )
            control_report = control["map_info"]["thread_broker"]
            report = lull["map_info"]["thread_broker"]
            # Recovery means returning to the allocation retained by the clean
            # workload. Using the perturbed run's final split asks when the
            # control reached a state it may correctly never choose.
            required = control_report["final_producer_limit"]
            # Exact cumulative busy counters integrate the same complete
            # logical work and exclude the injected source sleeps. Their paired
            # share is the ground-truth bias check; one final smoothed model
            # snapshot is intentionally not used as a workload-wide statistic.
            control_share = exact_whole_run_share(control)
            lull_share = exact_whole_run_share(lull)
            control_ready = first_reaches(control_report, required)
            recovery_lulls = [
                event for event in lull["lulls"]
                if persistent_source
                or (
                    event["started_nanos"] / 1e9 >= control_ready
                    and trajectory_state(
                        report, event["started_nanos"] / 1e9
                    ) >= required
                )
            ]
            recoveries = [
                recovery_delay(report, event["ended_nanos"] / 1e9, required)
                for event in recovery_lulls
            ]
            cap_delays = [
                cap_delay(
                    report,
                    event["started_nanos"] / 1e9,
                    manifest.get("source_ceiling", required),
                )
                for event in lull["lulls"]
                if event["duration_nanos"] / 1e9 >= manifest.get("cap_persistence_seconds", 0.3)
            ] if require_lull_cap else []
            retentions = cap_retention_intervals(report, lull["lulls"])
            cap_clear_delays = harmful_cap_clear_delays(
                report, lull["lulls"], required
            )
            source_ceiling = manifest.get("source_ceiling", required)
            control_caps = [
                observation for observation in control_report.get("cap_trajectory", [])
                if observation["reason"] != "none"
                and observation["useful_cap"] <= source_ceiling + 1
            ]
            source_identification = (
                min(seconds(observation["elapsed"]) for observation in control_caps)
                if control_caps
                else float("inf")
            )
            analysis = {
                "cadence_ms": cadence,
                "repetition": repetition,
                "required_producer_slots": required,
                "exact_whole_run_share_bias_points": abs(lull_share - control_share),
                "max_cap_identification_seconds": (
                    source_identification
                    if persistent_source
                    else max(cap_delays)
                    if cap_delays
                    else (float("inf") if require_lull_cap else 0.0)
                ),
                "max_post_lull_recovery_seconds": max(recoveries) if recoveries else float("inf"),
                "max_post_lull_cap_clear_seconds": (
                    max(cap_clear_delays) if cap_clear_delays else 0.0
                ),
                "cap_retention_seconds": retentions,
                "flow_transient_windows": report["flow_transient_windows"],
                "peak_controlled_slots": report["peak_controlled_slots"],
                "lulls": len(lull["lulls"]),
                "recovery_lulls": len(recovery_lulls),
                "paired_control_ready_seconds": control_ready,
                "cap_observations": len(report.get("cap_trajectory", [])),
                "persistent_source_observed": bool(control_caps),
            }
            criteria = [
                analysis["exact_whole_run_share_bias_points"] <= 0.05,
                analysis["max_cap_identification_seconds"] <= 1.0,
                analysis["max_post_lull_recovery_seconds"] <= 1.0,
                analysis["peak_controlled_slots"] <= manifest["budget"],
                analysis["lulls"] > 0,
                analysis["recovery_lulls"] > 0,
            ]
            if persistent_source:
                criteria.append(analysis["persistent_source_observed"])
            else:
                criteria.extend((
                    analysis["max_post_lull_cap_clear_seconds"] <= 1.0,
                ))
            analysis["passed"] = all(criteria)
            analyses.append(analysis)
    digests = {record["canonical_output_sha256"] for record in records}
    retention_by_cadence = {}
    for item in analyses:
        finite = [value for value in item["cap_retention_seconds"] if value != float("inf")]
        if finite:
            retention_by_cadence.setdefault(item["cadence_ms"], []).extend(finite)
    retention_medians = {
        cadence: sorted(values)[len(values) // 2]
        for cadence, values in retention_by_cadence.items()
    }
    positive = [value for value in retention_medians.values() if value > 0]
    retention_ratio = max(positive) / min(positive) if len(positive) > 1 else 1.0
    return {
        "analyses": analyses,
        "canonical_equal": len(digests) == 1,
        "cap_retention_ratio": retention_ratio,
        "cap_retention_median_seconds_by_cadence": retention_medians,
        "gate_passed": all(item["passed"] for item in analyses)
        and len(digests) == 1
        and (persistent_source or retention_ratio <= 1.20),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--summarize-existing", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    manifest["_binary_sha256"] = hashlib.sha256(
        Path(manifest["binary"]).read_bytes()
    ).hexdigest()
    root = Path(manifest["results_dir"])
    root.mkdir(parents=True, exist_ok=True)
    runs_path = root / "runs.jsonl"
    if args.summarize_existing:
        records = [json.loads(line) for line in runs_path.read_text().splitlines()]
    else:
        records = []
        completed = set()
        if args.resume:
            if not runs_path.exists():
                raise ValueError("--resume requested but runs.jsonl does not exist")
            previous = [
                json.loads(line) for line in runs_path.read_text().splitlines() if line
            ]
            records = [
                record
                for record in previous
                if record.get("exit_code") == 0 and "map_info" in record
            ]
            if any(
                record.get("binary_sha256") != manifest["_binary_sha256"]
                for record in records
            ):
                raise ValueError("cannot resume with a different or unrecorded binary")
            completed = {
                (record["cadence_ms"], record["mode"], record["repetition"])
                for record in records
            }
            runs_path.write_text(
                "".join(json.dumps(record, sort_keys=True) + "\n" for record in records)
            )
        elif manifest.get("reuse_control_records_from"):
            if runs_path.exists():
                raise ValueError("control reuse requires a new results_dir")
            source = Path(manifest["reuse_control_records_from"])
            previous = [
                json.loads(line) for line in source.read_text().splitlines() if line
            ]
            expected_command_prefix = [
                manifest["binary"], *manifest["args"], "-t", str(manifest["budget"]),
                "--decoder", "auto", "-o",
            ]
            records = [record for record in previous if record.get("mode") == "control"]
            expected_keys = {
                (cadence, "control", repetition)
                for repetition in range(1, manifest["repetitions"] + 1)
                for cadence in manifest["cadences_ms"]
            }
            completed = {
                (record.get("cadence_ms"), record.get("mode"), record.get("repetition"))
                for record in records
            }
            if completed != expected_keys:
                raise ValueError("reused controls do not exactly cover the requested matrix")
            if any(
                record.get("exit_code") != 0
                or "map_info" not in record
                or "canonical_output_sha256" not in record
                or record.get("binary_sha256") != manifest["_binary_sha256"]
                or record.get("command", [])[:-1] != expected_command_prefix
                for record in records
            ):
                raise ValueError("reused controls do not match this binary and command")
            runs_path.write_text(
                "".join(json.dumps(record, sort_keys=True) + "\n" for record in records)
            )
        elif runs_path.exists():
            raise ValueError("runs.jsonl already exists; use --resume or a new results_dir")
        for repetition in range(1, manifest["repetitions"] + 1):
            for cadence in manifest["cadences_ms"]:
                for lulls in (False, True):
                    key = (cadence, "lulls" if lulls else "control", repetition)
                    if key in completed:
                        continue
                    record = execute(manifest, cadence, lulls, repetition)
                    records.append(record)
                    with runs_path.open("a") as stream:
                        stream.write(json.dumps(record, sort_keys=True) + "\n")
                    if record["exit_code"] != 0:
                        raise SystemExit(record["error"])
    summary = summarize(records, manifest)
    (root / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    if not summary["gate_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
