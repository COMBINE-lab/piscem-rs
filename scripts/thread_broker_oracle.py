#!/usr/bin/env python3
"""Run randomized, same-binary thread-broker oracle sweeps.

The JSON manifest format is documented by `thread_broker_oracle.example.json`.
Every invocation writes one immutable JSON record per run plus a summary. A run
is rejected before timing is accepted if `map_info.json` does not prove the
requested/effective budget and fixed/adaptive allocation semantics.
"""

import argparse
import hashlib
import json
import os
import random
import statistics
import subprocess
import sys
import time
from pathlib import Path


def run_text(argv):
    return subprocess.run(
        argv, check=True, universal_newlines=True, stdout=subprocess.PIPE
    ).stdout.strip()


def dependency_revisions():
    metadata = json.loads(
        run_text(
            [
                "cargo",
                "metadata",
                "--offline",
                "--format-version",
                "1",
                "--features",
                "rapidgzip",
            ]
        )
    )
    return {
        package["name"]: package.get("source") or "workspace"
        for package in metadata["packages"]
        if package["name"] in {"paraseq", "rapidgzip-core", "thread-broker", "piscem-rs"}
    }


def expand_pins(budget, fractions, minimum_points=5):
    pins = {max(1, min(budget - 1, round(budget * fraction))) for fraction in fractions}
    if budget >= 6 and minimum_points >= 5 and len(pins) < 5:
        pins.update({1, budget - 1, budget // 4, budget // 2, 3 * budget // 4})
    return sorted(pins)


def oracle_winner_is_bracketed(fixed_modes, winner):
    """Whether a fixed-split winner has a measured point on both sides."""
    pins = sorted(int(mode[4:]) for mode in fixed_modes)
    winner_pin = int(winner[4:])
    return len(pins) >= 3 and pins[0] < winner_pin < pins[-1]


def output_paths(root, job, budget, mode, repetition):
    output = root / job["name"] / f"t{budget}" / mode / f"rep-{repetition}"
    if job["output_kind"] == "stem":
        output.parent.mkdir(parents=True, exist_ok=True)
        return output, Path(f"{output}.map_info.json")
    if job["output_kind"] == "directory":
        output.mkdir(parents=True, exist_ok=True)
        return output, output / "map_info.json"
    raise ValueError(f"unknown output_kind for {job['name']}: {job['output_kind']}")


def validate_telemetry(info, requested_budget, pin, require_components=False):
    plan = info["execution_plan"]
    effective = plan["effective_budget"]
    if plan["requested_budget"] != requested_budget or info["num_threads"] != effective:
        raise ValueError(f"budget telemetry mismatch: {plan}")
    if not 1 <= effective <= requested_budget:
        raise ValueError(f"invalid effective budget: {plan}")
    allocation = plan["allocation"]["kind"]
    if allocation == "serial":
        if plan["map_threads"] != effective or plan["decode_slots"] != 0:
            raise ValueError(f"serial plan did not return the whole budget: {plan}")
    elif plan["map_threads"] + plan["decode_slots"] != effective:
        raise ValueError(f"allocation does not equal effective budget: {plan}")

    if pin is not None:
        if allocation != "pinned_aggregate" or plan["requested_decode_slots"] != pin:
            raise ValueError(f"aggregate pin was not honored as a fixed request: {plan}")
        if "thread_broker" in info:
            raise ValueError("fixed aggregate run unexpectedly started the broker")
    elif allocation == "adaptive":
        report = info.get("thread_broker")
        if report is None:
            raise ValueError("adaptive run emitted no broker report")
        requested_consumers = report["consumer_trajectory"]
        requested_producers = report["producer_trajectory"]
        actual_consumers = report["actual_consumer_trajectory"]
        actual_producers = report["actual_producer_trajectory"]
        if not (
            len(requested_consumers)
            == len(requested_producers)
            == len(actual_consumers)
            == len(actual_producers)
        ):
            raise ValueError("broker requested/actual trajectories have different lengths")
        if not requested_consumers:
            raise ValueError("adaptive run emitted no requested/actual occupancy samples")
        # A two-phase move first shrinks one side and waits for acknowledgement,
        # then grows the other.  The requested total may therefore be below the
        # budget while a move is in flight; exceeding it is the invariant.
        if any(
            c < 0 or p < 0 or c + p > effective
            for c, p in zip(requested_consumers, requested_producers)
        ):
            raise ValueError(f"broker requested trajectory exceeded budget: {report}")
        if any(c + p > effective for c, p in zip(actual_consumers, actual_producers)):
            raise ValueError(f"broker actual trajectory exceeded budget: {report}")
        if report["final_consumer_threads"] + report["final_producer_limit"] != effective:
            raise ValueError(f"broker final request violated budget: {report}")
        if report["final_consumer_live"] + report["final_producer_active"] > effective:
            raise ValueError(f"broker final actual occupancy exceeded budget: {report}")
        if report["peak_controlled_slots"] > effective:
            raise ValueError(f"controlled occupancy exceeded budget: {report}")
        if report.get("terminal_error") is not None:
            raise ValueError(f"broker terminated with an error: {report['terminal_error']}")

    if not isinstance(info.get("mapping_seconds"), (int, float)) or info["mapping_seconds"] <= 0:
        raise ValueError("missing positive numeric mapping_seconds")
    if require_components:
        measurement = info.get("producer_measurement")
        if measurement is None:
            raise ValueError("component-enabled run emitted no producer measurement")
        if measurement.get("completed_worker_cpu_nanos") is None:
            raise ValueError("component-enabled run emitted no completed worker CPU time")
        if measurement.get("completed_auxiliary_cpu_nanos") is None:
            raise ValueError("component-enabled run emitted no completed auxiliary CPU time")
        if measurement.get("cpu_accounting_failures") != 0:
            raise ValueError(f"producer CPU accounting failed: {measurement}")
        if measurement.get("accounted_busy_nanos") != measurement.get("busy_nanos"):
            raise ValueError(f"fixed run did not use exact producer accounting: {measurement}")
        if any(
            measurement.get(name) != 0
            for name in ("calibration_samples", "monitoring_samples", "observation_nanos")
        ) or measurement.get("sampler_cpu_nanos") is not None:
            raise ValueError(f"fixed exact accounting unexpectedly ran a sampler: {measurement}")


def canonical_digest(command, output):
    if not command:
        return None
    argv = [part.replace("{output}", str(output)) for part in command]
    digest = hashlib.sha256()
    process = subprocess.Popen(argv, stdout=subprocess.PIPE)
    while True:
        chunk = process.stdout.read(1024 * 1024)
        if not chunk:
            break
        digest.update(chunk)
    return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, argv)
    return digest.hexdigest()


def canonical_digest_expected(manifest, job, budget, mode):
    budgets = job.get(
        "canonical_digest_budgets", manifest.get("canonical_digest_budgets")
    )
    modes = job.get("canonical_digest_modes", manifest.get("canonical_digest_modes"))
    return (budgets is None or budget in budgets) and (modes is None or mode in modes)


def manifest_jobs(manifest):
    included = manifest.get("include_jobs")
    return [
        job
        for job in manifest["jobs"]
        if included is None or job["name"] in included
    ]


def cleanup(job, output):
    for template in job.get("cleanup_paths", []):
        path = Path(template.replace("{output}", str(output)))
        if path.is_file():
            path.unlink()


def percentile(values, fraction):
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(len(ordered) * fraction))]


def paired_ratio_upper(baseline, treatment, seed, iterations=10000):
    """Upper one-sided 95% bound for the median of paired treatment ratios."""
    baseline_by_rep = {item["repetition"]: item for item in baseline}
    treatment_by_rep = {item["repetition"]: item for item in treatment}
    repetitions = sorted(baseline_by_rep.keys() & treatment_by_rep.keys())
    ratios = [
        treatment_by_rep[repetition]["map_info"]["mapping_seconds"]
        / baseline_by_rep[repetition]["map_info"]["mapping_seconds"]
        for repetition in repetitions
    ]
    if len(ratios) < 3:
        raise ValueError("at least three complete paired repetitions are required")
    rng = random.Random(seed)
    bootstrap_medians = []
    for _ in range(iterations):
        bootstrap_medians.append(statistics.median(rng.choice(ratios) for _ in ratios))
    return percentile(bootstrap_medians, 0.95), ratios


def duration_seconds(value):
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return value
    return value.get("secs", 0) + value.get("nanos", 0) / 1_000_000_000


def convergence_result(item, job):
    report = item["map_info"].get("thread_broker", {})
    converged_in = duration_seconds(report.get("time_to_converge"))
    final_drift = duration_seconds(report.get("current_drift"))
    mapping_seconds = item["map_info"]["mapping_seconds"]
    # A finite job ending before the five-second gate horizon cannot demonstrate
    # five-second convergence. It must still remain live and bounded, without a
    # terminal error, resurvey, or excessive moves. Dedicated short-run jobs opt
    # back into their stricter duration-relative convergence requirement.
    evaluable = mapping_seconds >= 5.0 or job.get("short_run", False)
    passed = report.get("moves", 6) <= 5 and report.get("terminal_error") is None
    if evaluable:
        passed &= (
            converged_in is not None
            and converged_in <= 5.0
            and report.get("final_phase") == "steady"
            and final_drift == 0
        )
    if job.get("short_run", False):
        passed &= converged_in is not None and converged_in <= 0.20 * mapping_seconds
    if job.get("stable_workload", True):
        passed &= report.get("resurveys", 1) == 0
    return {
        "repetition": item["repetition"],
        "seconds": converged_in,
        "moves": report.get("moves"),
        "resurveys": report.get("resurveys"),
        "evaluable": evaluable,
        "passed": passed,
    }


def linear_prediction(points, field, target):
    ordered = sorted(points, key=lambda point: point["producer_slots"])
    if len(ordered) < 2:
        raise ValueError("allocation-aware prediction requires at least two training points")
    left, right = ordered[0], ordered[1]
    for first, second in zip(ordered, ordered[1:]):
        left, right = first, second
        if first["producer_slots"] <= target <= second["producer_slots"]:
            break
        if target > second["producer_slots"]:
            continue
        break
    x0, x1 = left["producer_slots"], right["producer_slots"]
    weight = (target - x0) / (x1 - x0)
    return left[field] + weight * (right[field] - left[field])


def fixed_cost_analysis(cells, budget):
    points = []
    for mode, runs in cells.items():
        if not mode.startswith("pin-"):
            continue
        producer = [run["map_info"]["producer_measurement"] for run in runs]
        consumer = [run["map_info"]["consumer_measurement"] for run in runs]
        points.append({
            "mode": mode,
            "producer_slots": int(mode[len("pin-"):]),
            "mapping_seconds": statistics.median(
                run["map_info"]["mapping_seconds"] for run in runs
            ),
            "producer_busy_nanos": statistics.median(
                measurement["busy_nanos"] for measurement in producer
            ),
            "consumer_busy_nanos": statistics.median(
                measurement["busy_nanos"] for measurement in consumer
            ),
        })
    points.sort(key=lambda point: point["producer_slots"])
    if len(points) < 5:
        return {"gate_evaluable": False, "points": points}

    def span(field):
        values = [point[field] for point in points]
        return max(values) / min(values) - 1.0

    constant_errors = []
    allocation_aware_errors = []
    predictions = []
    for held_out in points:
        training = [point for point in points if point is not held_out]

        def calibrated_prediction(producer_busy, consumer_busy):
            raw = max(
                producer_busy / held_out["producer_slots"],
                consumer_busy / (budget - held_out["producer_slots"]),
            ) / 1e9
            calibration = statistics.median(
                point["mapping_seconds"]
                / (
                    max(
                        point["producer_busy_nanos"] / point["producer_slots"],
                        point["consumer_busy_nanos"]
                        / (budget - point["producer_slots"]),
                    )
                    / 1e9
                )
                for point in training
            )
            return raw * calibration

        constant_prediction = calibrated_prediction(
            statistics.median(point["producer_busy_nanos"] for point in training),
            statistics.median(point["consumer_busy_nanos"] for point in training),
        )
        allocation_prediction = calibrated_prediction(
            linear_prediction(training, "producer_busy_nanos", held_out["producer_slots"]),
            linear_prediction(training, "consumer_busy_nanos", held_out["producer_slots"]),
        )
        # Throughput is inverse wall time, so predicted/observed throughput is
        # measured_time/predicted_time.
        constant_error = abs(held_out["mapping_seconds"] / constant_prediction - 1.0)
        allocation_error = abs(held_out["mapping_seconds"] / allocation_prediction - 1.0)
        constant_errors.append(constant_error)
        allocation_aware_errors.append(allocation_error)
        predictions.append({
            "mode": held_out["mode"],
            "measured_mapping_seconds": held_out["mapping_seconds"],
            "constant_cost_predicted_seconds": constant_prediction,
            "constant_cost_throughput_error": constant_error,
            "allocation_aware_predicted_seconds": allocation_prediction,
            "allocation_aware_throughput_error": allocation_error,
        })
    producer_span = span("producer_busy_nanos")
    consumer_span = span("consumer_busy_nanos")
    return {
        "gate_evaluable": True,
        "points": points,
        "producer_busy_relative_span": producer_span,
        "consumer_busy_relative_span": consumer_span,
        "allocation_dependence_exceeds_ten_percent": max(producer_span, consumer_span) > 0.10,
        "constant_cost_max_held_out_throughput_error": max(constant_errors),
        "allocation_aware_max_held_out_throughput_error": max(allocation_aware_errors),
        "allocation_aware_prediction_gate_passed": max(allocation_aware_errors) <= 0.10,
        "held_out_predictions": predictions,
    }


def execute(spec, manifest, provenance, dry_run):
    output, map_info = output_paths(
        Path(manifest["results_dir"]), spec["job"], spec["budget"], spec["mode"], spec["repetition"]
    )
    argv = [manifest["binary"], *spec["job"]["args"], "-t", str(spec["budget"]), "--decoder", "auto"]
    argv.extend([spec["job"].get("output_argument", "-o"), str(output)])
    env = os.environ.copy()
    env.pop("PISCEM_DECODE_SLOTS", None)
    env.pop("PISCEM_RAPIDGZIP_THREADS", None)
    env.pop("PISCEM_THREAD_BROKER_POLICY", None)
    env.pop("PISCEM_THREAD_BROKER_PROBE_INTERVAL_MS", None)
    env.pop("PISCEM_FIXED_DECODE_MEASUREMENT", None)
    if spec["pin"] is not None:
        env["PISCEM_DECODE_SLOTS"] = str(spec["pin"])
        if spec["job"].get("measure_fixed_components", False):
            env["PISCEM_FIXED_DECODE_MEASUREMENT"] = "native"

    record = {
        **provenance,
        **{
            key: spec[key]
            for key in ["budget", "mode", "pin", "repetition", "block_order", "block_position"]
        },
    }
    record["job"] = spec["job"]["name"]
    record["command"] = argv
    if dry_run:
        return record

    log_path = Path(f"{map_info}.run.log")
    time_path = Path(f"{map_info}.time.json")
    timed = [
        "/usr/bin/time", "-f", '{"user_seconds":%U,"system_seconds":%S,"peak_rss_kib":%M,"wall_seconds":%e}',
        "-o", str(time_path), *argv,
    ]
    started = time.monotonic()
    with log_path.open("wb") as log:
        result = subprocess.run(timed, env=env, stdout=log, stderr=subprocess.STDOUT)
    record["runner_wall_seconds"] = time.monotonic() - started
    record["exit_code"] = result.returncode
    try:
        if result.returncode != 0:
            record["error"] = f"command failed; see {log_path}"
            return record

        info = json.loads(map_info.read_text())
        validate_telemetry(
            info,
            spec["budget"],
            spec["pin"],
            spec["job"].get("measure_fixed_components", False),
        )
        record["map_info"] = info
        record["process"] = json.loads(time_path.read_text())
        digest_repetitions = spec["job"].get(
            "canonical_digest_repetitions",
            manifest.get("canonical_digest_repetitions"),
        )
        record["canonical_output_sha256"] = (
            canonical_digest(spec["job"].get("canonical_digest_command"), output)
            if canonical_digest_expected(
                manifest, spec["job"], spec["budget"], spec["mode"]
            )
            and (digest_repetitions is None or spec["repetition"] in digest_repetitions)
            else None
        )
        return record
    finally:
        cleanup(spec["job"], output)


def summarize(records, manifest):
    summary = {
        "tier": manifest.get("tier", "custom"),
        "jobs": {},
        "gate_passed": True,
        "failures": [],
    }
    successful = [record for record in records if record.get("exit_code", 0) == 0 and "map_info" in record]
    for job_index, job in enumerate(manifest_jobs(manifest)):
        job_summary = summary["jobs"].setdefault(job["name"], {})
        for budget_index, budget in enumerate(job["threads"]):
            auto_only_threads = set(
                job.get("auto_only_threads", manifest.get("auto_only_threads", []))
            )
            cells = {}
            for record in successful:
                if record["job"] == job["name"] and record["budget"] == budget:
                    cells.setdefault(record["mode"], []).append(record)
            medians = {
                mode: statistics.median(item["map_info"]["mapping_seconds"] for item in runs)
                for mode, runs in cells.items()
            }
            fixed = {mode: value for mode, value in medians.items() if mode.startswith("pin-")}
            if "auto" not in medians:
                summary["gate_passed"] = False
                summary["failures"].append(f"{job['name']} t{budget}: incomplete cells")
                continue
            if not fixed:
                if budget not in auto_only_threads:
                    summary["gate_passed"] = False
                    summary["failures"].append(
                        f"{job['name']} t{budget}: fixed oracle cells are absent"
                    )
                    continue
                convergence = []
                convergence_passed = True
                for item in cells["auto"]:
                    result = convergence_result(item, job)
                    convergence.append(result)
                    convergence_passed &= result["passed"]
                counts = {
                    (
                        item["map_info"]["num_reads"],
                        item["map_info"]["num_mapped"],
                        item["map_info"]["num_poisoned"],
                    )
                    for item in cells["auto"]
                }
                digests = {
                    item["canonical_output_sha256"] for item in cells["auto"]
                    if item.get("canonical_output_sha256") is not None
                }
                digest_coverage = any(
                    item.get("canonical_output_sha256") is not None
                    for item in cells["auto"]
                )
                digest_required = canonical_digest_expected(
                    manifest, job, budget, "auto"
                )
                correctness = len(counts) == 1 and (
                    (len(digests) == 1 and digest_coverage)
                    if digest_required
                    else (len(digests) == 0 and not digest_coverage)
                )
                passed = convergence_passed and correctness
                job_summary[f"t{budget}"] = {
                    "evidence_level": "adaptive-only",
                    "median_mapping_seconds": medians,
                    "convergence": convergence,
                    "convergence_gate_passed": convergence_passed,
                    "correctness": correctness,
                    "gate_passed": passed,
                }
                if not passed:
                    summary["gate_passed"] = False
                    summary["failures"].append(
                        f"{job['name']} t{budget}: adaptive-only gate failed"
                    )
                continue
            oracle_candidate_mode = min(fixed, key=fixed.get)
            oracle_bracketed = oracle_winner_is_bracketed(fixed, oracle_candidate_mode)
            oracle_mode = oracle_candidate_mode
            positions = range(1, len(cells) + 1)
            block_position_counts = {
                mode: {
                    str(position): sum(
                        item["block_position"] == position for item in runs
                    )
                    for position in positions
                }
                for mode, runs in cells.items()
            }
            positions_balanced = all(
                max(counts.values(), default=0) - min(counts.values(), default=0) <= 1
                for counts in block_position_counts.values()
            )
            oracle = fixed[oracle_mode]
            flat_mode = f"pin-{max(1, min(budget - 1, round(budget / 3)))}"
            auto = medians["auto"]
            if flat_mode not in medians:
                summary["gate_passed"] = False
                summary["failures"].append(f"{job['name']} t{budget}: flat split is absent")
                continue

            auto_oracle_upper, auto_oracle_pairs = paired_ratio_upper(
                cells[oracle_mode],
                cells["auto"],
                20260806 + job_index * 100 + budget_index * 10,
            )
            auto_flat_upper, auto_flat_pairs = paired_ratio_upper(
                cells[flat_mode],
                cells["auto"],
                20260807 + job_index * 100 + budget_index * 10,
            )
            short_safe_mode = None
            short_safe_upper = None
            short_safe_pairs = None
            short_run_passed = True
            if "short_run_safe_pin" in job:
                short_safe_mode = f"pin-{job['short_run_safe_pin']}"
                if short_safe_mode not in cells:
                    summary["gate_passed"] = False
                    summary["failures"].append(
                        f"{job['name']} t{budget}: short-run safe split is absent"
                    )
                    continue
                short_safe_upper, short_safe_pairs = paired_ratio_upper(
                    cells[short_safe_mode],
                    cells["auto"],
                    20260808 + job_index * 100 + budget_index * 10,
                )
                short_run_passed = short_safe_upper <= 1.05
            near_boundary = auto / oracle >= 1.07 or auto / medians[flat_mode] >= 0.97
            flat_regression_margin = job.get(
                "flat_regression_margin",
                manifest.get("flat_regression_margin", 0.0),
            )
            enough_repetitions = (
                not near_boundary
                or manifest["repetitions"]
                >= manifest.get("boundary_repetitions_required", 5)
            )
            wall_passed = (
                oracle_bracketed
                and auto_oracle_upper <= 1.10
                and auto / medians[flat_mode] <= 1.0 + flat_regression_margin
                and enough_repetitions
                and positions_balanced
            )

            cpu_medians = {
                mode: statistics.median(
                    item["process"]["user_seconds"] + item["process"]["system_seconds"]
                    for item in runs
                )
                for mode, runs in cells.items()
            }
            rss_medians = {
                mode: statistics.median(item["process"]["peak_rss_kib"] for item in runs)
                for mode, runs in cells.items()
            }
            component_medians = {}
            for mode, runs in cells.items():
                producer = [
                    item["map_info"].get("producer_measurement")
                    for item in runs
                    if item["map_info"].get("producer_measurement") is not None
                ]
                consumer = [
                    item["map_info"].get("consumer_measurement")
                    for item in runs
                    if item["map_info"].get("consumer_measurement") is not None
                ]
                component_medians[mode] = {
                    "producer_busy_nanos": (
                        statistics.median(item["busy_nanos"] for item in producer)
                        if producer
                        else None
                    ),
                    "completed_worker_cpu_nanos": (
                        statistics.median(
                            item["completed_worker_cpu_nanos"]
                            for item in producer
                            if item.get("completed_worker_cpu_nanos") is not None
                        )
                        if any(
                            item.get("completed_worker_cpu_nanos") is not None
                            for item in producer
                        )
                        else None
                    ),
                    "completed_auxiliary_cpu_nanos": (
                        statistics.median(
                            item["completed_auxiliary_cpu_nanos"]
                            for item in producer
                            if item.get("completed_auxiliary_cpu_nanos") is not None
                        )
                        if any(
                            item.get("completed_auxiliary_cpu_nanos") is not None
                            for item in producer
                        )
                        else None
                    ),
                    "consumer_busy_nanos": (
                        statistics.median(item["busy_nanos"] for item in consumer)
                        if consumer
                        else None
                    ),
                    "consumer_setup_nanos": (
                        statistics.median(item["callback_setup_nanos"] for item in consumer)
                        if consumer
                        else None
                    ),
                    "consumer_flush_nanos": (
                        statistics.median(item["output_flush_nanos"] for item in consumer)
                        if consumer
                        else None
                    ),
                }
            # "Best fixed split" is the wall-time oracle established above,
            # not a different pathological point selected independently for
            # each resource. A split that is 3x slower because it gives almost
            # every slot to decode may use slightly less mapping-state memory;
            # that does not make it the comparison configuration in Gate H.
            oracle_cpu = cpu_medians[oracle_mode]
            oracle_rss = rss_medians[oracle_mode]
            resource_passed = (
                cpu_medians["auto"] <= 1.03 * oracle_cpu
                and rss_medians["auto"] <= 1.02 * oracle_rss
            )
            if job.get("allow_resource_tradeoff", False):
                resource_passed = True
            cost_analysis = fixed_cost_analysis(cells, budget)

            convergence = []
            convergence_passed = True
            if job.get("require_convergence", True):
                for item in cells["auto"]:
                    result = convergence_result(item, job)
                    convergence.append(result)
                    convergence_passed &= result["passed"]

            passed = wall_passed and resource_passed and convergence_passed and short_run_passed
            counts = {
                (record["map_info"]["num_reads"], record["map_info"]["num_mapped"], record["map_info"]["num_poisoned"])
                for runs in cells.values() for record in runs
            }
            digests = {
                record["canonical_output_sha256"] for runs in cells.values() for record in runs
                if record.get("canonical_output_sha256") is not None
            }
            digested_modes = {
                mode for mode, runs in cells.items()
                if any(record.get("canonical_output_sha256") is not None for record in runs)
            }
            expected_digested_modes = {
                mode
                for mode in cells
                if canonical_digest_expected(manifest, job, budget, mode)
            }
            # Counts alone cannot detect changed mappings. Representative cells
            # must also provide a canonical, order-independent digest. Normal
            # and light may select budgets/modes; comprehensive leaves both
            # selectors unset and therefore requires every mode in every cell.
            correctness = (
                len(counts) == 1
                and (len(digests) == 1 if expected_digested_modes else len(digests) == 0)
                and digested_modes == expected_digested_modes
            )
            passed &= correctness
            job_summary[f"t{budget}"] = {
                "median_mapping_seconds": medians,
                "oracle_mode": oracle_mode if oracle_bracketed else None,
                "oracle_candidate_mode": oracle_candidate_mode,
                "oracle_bracketed": oracle_bracketed,
                "auto_to_oracle": auto / oracle,
                "auto_to_oracle_upper_95": auto_oracle_upper,
                "auto_to_oracle_paired_ratios": auto_oracle_pairs,
                "flat_mode": flat_mode,
                "auto_to_flat": auto / medians[flat_mode],
                "auto_to_flat_upper_95": auto_flat_upper,
                "flat_regression_margin": flat_regression_margin,
                "auto_to_flat_paired_ratios": auto_flat_pairs,
                "short_run_safe_mode": short_safe_mode,
                "auto_to_short_run_safe_upper_95": short_safe_upper,
                "auto_to_short_run_safe_paired_ratios": short_safe_pairs,
                "short_run_gate_passed": short_run_passed,
                "enough_boundary_repetitions": enough_repetitions,
                "block_position_counts": block_position_counts,
                "block_positions_balanced": positions_balanced,
                "median_cpu_seconds": cpu_medians,
                "auto_to_oracle_cpu": cpu_medians["auto"] / oracle_cpu,
                "median_peak_rss_kib": rss_medians,
                "median_components": component_medians,
                "allocation_cost_analysis": cost_analysis,
                "auto_to_oracle_rss": rss_medians["auto"] / oracle_rss,
                "resource_gate_passed": resource_passed,
                "convergence": convergence,
                "convergence_gate_passed": convergence_passed,
                "correctness": correctness,
                "gate_passed": passed,
            }
            if not passed:
                summary["gate_passed"] = False
                reason = (
                    f"unbracketed fixed winner {oracle_candidate_mode}"
                    if not oracle_bracketed
                    else "oracle/correctness gate failed"
                )
                summary["failures"].append(f"{job['name']} t{budget}: {reason}")
    if any(record.get("exit_code", 0) != 0 or "error" in record for record in records):
        summary["gate_passed"] = False
        summary["failures"].append("one or more commands failed")
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--summarize-only",
        action="store_true",
        help="recompute summary.json from an existing runs.jsonl",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        help="override the manifest results directory (useful for retained evidence)",
    )
    parser.add_argument(
        "--include-job",
        action="append",
        help="override the manifest job selection; may be supplied more than once",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="retain successful cells from an interrupted matrix and run the rest",
    )
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    if args.results_dir is not None:
        manifest["results_dir"] = str(args.results_dir)
    if args.include_job is not None:
        manifest["include_jobs"] = args.include_job
    if manifest.get("repetitions", 0) < 3:
        raise ValueError("at least three repetitions are required")

    results_dir = Path(manifest["results_dir"])
    if args.summarize_only:
        records = [
            json.loads(line)
            for line in (results_dir / "runs.jsonl").read_text().splitlines()
            if line
        ]
        summary = summarize(records, manifest)
        (results_dir / "summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n"
        )
        return 0 if summary["gate_passed"] else 1

    provenance = {
        "commit": run_text(["git", "rev-parse", "HEAD"]),
        "dirty": bool(run_text(["git", "status", "--porcelain"])),
        "dependencies": dependency_revisions(),
        "binary_sha256": hashlib.sha256(Path(manifest["binary"]).read_bytes()).hexdigest(),
    }
    blocks = []
    fractions = manifest.get("pin_fractions", [0.125, 0.25, 1 / 3, 0.5, 0.75])
    fractions_by_budget = manifest.get("pin_fractions_by_budget", {})
    minimum_pin_points = manifest.get("minimum_pin_points", 5)
    minimum_pin_points_by_budget = manifest.get("minimum_pin_points_by_budget", {})
    seed = manifest.get("seed", 20260806)
    for job_index, job in enumerate(manifest_jobs(manifest)):
        for budget_index, budget in enumerate(job["threads"]):
            budget_fractions = fractions_by_budget.get(str(budget), fractions)
            budget_minimum = minimum_pin_points_by_budget.get(
                str(budget), minimum_pin_points
            )
            auto_only_threads = set(
                job.get("auto_only_threads", manifest.get("auto_only_threads", []))
            )
            modes = [("auto", None)]
            if budget not in auto_only_threads:
                modes.extend(
                    (f"pin-{pin}", pin)
                    for pin in expand_pins(budget, budget_fractions, budget_minimum)
                )
            base_order = modes[:]
            random.Random(seed + job_index * 1000 + budget_index).shuffle(base_order)
            for repetition in range(1, manifest["repetitions"] + 1):
                offset = (repetition - 1) % len(base_order)
                order = base_order[offset:] + base_order[:offset]
                blocks.append(
                    [
                        {
                            "job": job,
                            "budget": budget,
                            "mode": mode,
                            "pin": pin,
                            "repetition": repetition,
                            "block_position": position,
                        }
                        for position, (mode, pin) in enumerate(order, 1)
                    ]
                )
    random.Random(seed).shuffle(blocks)
    specs = []
    for block_order, block in enumerate(blocks, 1):
        for spec in block:
            spec["block_order"] = block_order
            specs.append(spec)

    results_dir.mkdir(parents=True, exist_ok=True)
    raw_path = results_dir / "runs.jsonl"
    records = []
    completed = set()
    if args.resume:
        if not raw_path.exists():
            raise ValueError("--resume requested but runs.jsonl does not exist")
        previous = [json.loads(line) for line in raw_path.read_text().splitlines() if line]
        records = [
            record
            for record in previous
            if record.get("exit_code") == 0 and "map_info" in record
        ]
        expected_binary = provenance["binary_sha256"]
        if any(record.get("binary_sha256") != expected_binary for record in records):
            raise ValueError("cannot resume with a different release binary")
        completed = {
            (record["job"], record["budget"], record["mode"], record["repetition"])
            for record in records
        }
        print(
            "[resume] retained %d successful cells; failed/incomplete cells will rerun"
            % len(records),
            flush=True,
        )
    elif raw_path.exists():
        raise ValueError("runs.jsonl already exists; use --resume or a new results_dir")
    with raw_path.open("w") as raw:
        for record in records:
            raw.write(json.dumps(record, sort_keys=True) + "\n")
        raw.flush()
        for index, spec in enumerate(specs, 1):
            key = (spec["job"]["name"], spec["budget"], spec["mode"], spec["repetition"])
            if key in completed:
                continue
            print(f"[{index}/{len(specs)}] {spec['job']['name']} t{spec['budget']} {spec['mode']} rep {spec['repetition']}", flush=True)
            try:
                record = execute(spec, manifest, provenance, args.dry_run)
            except Exception as error:  # retain evidence and continue the randomized matrix
                record = {
                    **provenance,
                    "job": spec["job"]["name"],
                    **{
                        key: spec[key]
                        for key in [
                            "budget",
                            "mode",
                            "pin",
                            "repetition",
                            "block_order",
                            "block_position",
                        ]
                    },
                    "error": repr(error),
                }
            records.append(record)
            raw.write(json.dumps(record, sort_keys=True) + "\n")
            raw.flush()
    if args.dry_run:
        return 0
    summary = summarize(records, manifest)
    (results_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return 0 if summary["gate_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
