#!/usr/bin/env python3
"""Measure whole-broker cost against an oracle-pinned, no-broker baseline.

Each randomized block runs the same release binary in the modes selected by
the manifest (the original four remain the default):

* ``fixed``: the decode split is pinned at the expected settled allocation;
* ``responsive``: the full broker keeps monitoring after convergence;
* ``normal-responsive``: responsive with an explicit 25 ms steady interval;
* ``sparse-responsive``: the same policy at the configured low probe rate;
* ``freeze``: the full broker terminates after its first honest convergence.
* ``full-calibration-freeze``: complete bounded nonlinear calibration, then
  terminate recurring broker work.

The runner reports paired absolute deltas and fractional ratios for mapping wall
time, process wall time, and aggregate process CPU time. That one-percent-style
comparison is a broad regression backstop, not the primary administrative-cost
gate: controller and any compatibility-sampler lifetime thread CPU are also
measured directly using two clock reads per auxiliary thread. Canonical output, allocation,
policy, convergence, and CPU-accounting telemetry are validated before a timing
sample is admitted. See ``thread_broker_policy_overhead.example.json``.
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


POLICY_ENV = "PISCEM_THREAD_BROKER_POLICY"
PROBE_ENV = "PISCEM_THREAD_BROKER_PROBE_INTERVAL_MS"
SPLIT_ENV = "PISCEM_DECODE_SLOTS"
MEASUREMENT_ENV = "PISCEM_FIXED_DECODE_MEASUREMENT"
ALL_VARIANTS = (
    {"name": "fixed", "allocation": "fixed"},
    {"name": "responsive", "allocation": "adaptive", "policy": "responsive"},
    {
        "name": "normal-responsive",
        "allocation": "adaptive",
        "policy": "responsive",
        "probe_interval_ms": 25,
    },
    {
        "name": "sparse-responsive",
        "allocation": "adaptive",
        "policy": "responsive",
        "sparse": True,
    },
    {
        "name": "freeze",
        "allocation": "adaptive",
        "policy": "freeze-after-convergence",
    },
    {
        "name": "full-calibration-freeze",
        "allocation": "adaptive",
        "policy": "freeze-after-full-calibration",
    },
)
DEFAULT_VARIANT_NAMES = ("fixed", "responsive", "sparse-responsive", "freeze")


def selected_variants(manifest):
    by_name = {variant["name"]: variant for variant in ALL_VARIANTS}
    names = manifest.get("variants", DEFAULT_VARIANT_NAMES)
    if len(names) != len(set(names)):
        raise ValueError("manifest variants must be unique")
    baseline = manifest.get("baseline_variant", "fixed")
    if baseline not in names:
        raise ValueError("manifest variants must include baseline_variant %r" % baseline)
    unknown = [name for name in names if name not in by_name]
    if unknown:
        raise ValueError("unknown manifest variants: %r" % unknown)
    return tuple(by_name[name] for name in names)


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
                "rapidgzip-cpu-accounting",
            ]
        )
    )
    wanted = {"paraseq", "rapidgzip-core", "thread-broker", "piscem-rs"}
    return {
        package["name"]: package.get("source") or "workspace"
        for package in metadata["packages"]
        if package["name"] in wanted
    }


def output_paths(root, job, mode, repetition):
    output = root / job["name"] / mode / ("rep-%d" % repetition)
    if job["output_kind"] == "stem":
        output.parent.mkdir(parents=True, exist_ok=True)
        return output, Path("%s.map_info.json" % output)
    if job["output_kind"] == "directory":
        output.mkdir(parents=True, exist_ok=True)
        return output, output / "map_info.json"
    raise ValueError("unknown output_kind: %s" % job["output_kind"])


def canonical_digest(command, output):
    if not command:
        raise ValueError("canonical_digest_command is required")
    argv = [part.replace("{output}", str(output)) for part in command]
    digest = hashlib.sha256()
    process = subprocess.Popen(argv, stdout=subprocess.PIPE)
    while True:
        chunk = process.stdout.read(1024 * 1024)
        if not chunk:
            break
        digest.update(chunk)
    code = process.wait()
    if code:
        raise subprocess.CalledProcessError(code, argv)
    return digest.hexdigest()


def duration_seconds(value):
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    return float(value.get("secs", 0)) + float(value.get("nanos", 0)) / 1e9


def validate_telemetry(info, job, variant):
    plan = info["execution_plan"]
    effective = plan["effective_budget"]
    if plan["requested_budget"] != job["budget"] or info["num_threads"] != effective:
        raise ValueError("budget telemetry mismatch: %r" % plan)
    if plan["map_threads"] + plan["decode_slots"] != effective:
        raise ValueError("opening allocation violated budget: %r" % plan)

    report = info.get("thread_broker")
    if variant["allocation"] == "fixed":
        if plan["allocation"]["kind"] != "pinned_aggregate":
            raise ValueError("baseline was not aggregate-pinned: %r" % plan)
        if plan["decode_slots"] != job["pin"] or report is not None:
            raise ValueError("baseline did not honor pin or unexpectedly ran broker")
        return

    if plan["allocation"]["kind"] != "adaptive" or report is None:
        raise ValueError("adaptive policy emitted no broker report")
    if report["final_consumer_threads"] + report["final_producer_limit"] != effective:
        raise ValueError("final broker allocation violated budget")
    tolerance = job.get("final_split_tolerance", 0)
    if abs(report["final_producer_limit"] - job["pin"]) > tolerance:
        raise ValueError(
            "policy settled at %d, expected %d +/- %d"
            % (report["final_producer_limit"], job["pin"], tolerance)
        )
    serialized_policy = variant["policy"].replace("-", "_")
    if report["steady_state_policy"] != serialized_policy:
        raise ValueError("policy telemetry mismatch: %r" % report)
    if not report.get("time_to_converge") or report["final_phase"] != "steady":
        raise ValueError("policy did not converge: %r" % report)
    frozen = variant["policy"].startswith("freeze-")
    if report["monitoring_stopped_after_convergence"] != frozen:
        raise ValueError("policy shutdown telemetry mismatch")
    if report["controller_samples"] <= 0 or duration_seconds(report["controller_elapsed"]) <= 0:
        raise ValueError("missing controller activity telemetry")
    if report.get("controller_cpu_nanos") is None:
        raise ValueError("controller CPU accounting was not enabled")
    if report.get("controller_cpu_accounting_failures") != 0:
        raise ValueError("controller CPU accounting reported failures")
    measurement = report.get("producer_measurement")
    if measurement is None:
        raise ValueError("adaptive policy emitted no producer measurement")
    native = measurement.get("accounted_busy_nanos") == measurement.get("busy_nanos")
    if native:
        if any(
            measurement.get(name) != 0
            for name in ("calibration_samples", "monitoring_samples", "observation_nanos")
        ) or measurement.get("sampler_cpu_nanos") is not None:
            raise ValueError("native cumulative signal unexpectedly ran a sampler")
    elif measurement["observation_nanos"] <= 0:
        raise ValueError("sampled producer measurement reported no observations")
    if job.get("require_producer_cpu_accounting", True):
        if measurement.get("completed_worker_cpu_nanos") is None:
            raise ValueError("worker CPU accounting was not enabled")
        if measurement.get("completed_auxiliary_cpu_nanos") is None:
            raise ValueError("auxiliary CPU accounting was not enabled")
        if measurement.get("cpu_accounting_failures") != 0:
            raise ValueError("producer CPU accounting reported failures")
    if not native and measurement.get("sampler_cpu_nanos") is None:
        raise ValueError("sampled producer CPU accounting was not completed")
    if measurement.get("sampler_cpu_accounting_failures") != 0:
        raise ValueError("sampler CPU accounting reported failures")


def cleanup(job, output):
    for template in job.get("cleanup_paths", []):
        path = Path(template.replace("{output}", str(output)))
        if path.is_file():
            path.unlink()


def execute(spec, manifest, provenance):
    job = spec["job"]
    variant = spec["variant"]
    output, map_info = output_paths(
        Path(manifest["results_dir"]), job, variant["name"], spec["repetition"]
    )
    argv = [
        manifest["binary"],
        *job["args"],
        "-t",
        str(job["budget"]),
        "--decoder",
        "auto",
        job.get("output_argument", "-o"),
        str(output),
    ]
    env = os.environ.copy()
    for name in (
        "PISCEM_RAPIDGZIP_THREADS",
        SPLIT_ENV,
        POLICY_ENV,
        PROBE_ENV,
        MEASUREMENT_ENV,
    ):
        env.pop(name, None)
    if variant["allocation"] == "fixed":
        env[SPLIT_ENV] = str(job["pin"])
    else:
        env[POLICY_ENV] = variant["policy"]
        if "probe_interval_ms" in variant:
            env[PROBE_ENV] = str(variant["probe_interval_ms"])
        elif variant.get("sparse"):
            env[PROBE_ENV] = str(job["sparse_probe_interval_ms"])
    for name, value in job.get("environment", {}).items():
        env[name] = str(value)

    record = {
        **provenance,
        "job": job["name"],
        "mode": variant["name"],
        "repetition": spec["repetition"],
        "block_order": spec["block_order"],
        "block_position": spec["block_position"],
        "command": argv,
        "environment": job.get("environment", {}),
    }
    log_path = Path("%s.run.log" % map_info)
    time_path = Path("%s.time.json" % map_info)
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
    record["runner_wall_seconds"] = time.monotonic() - started
    record["exit_code"] = result.returncode
    if result.returncode:
        record["error"] = "command failed; see %s" % log_path
        return record

    try:
        info = json.loads(map_info.read_text())
        validate_telemetry(info, job, variant)
        record["map_info"] = info
        record["process"] = json.loads(time_path.read_text())
        record["canonical_output_sha256"] = canonical_digest(
            job["canonical_digest_command"], output
        )
        return record
    finally:
        # Benchmark outputs can be many gigabytes. Keep logs, map-info, raw
        # timing, and digests, but honor manifest cleanup even when telemetry
        # validation rejects a run before its digest is admitted.
        cleanup(job, output)


def percentile(values, fraction):
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(len(ordered) * fraction))]


def paired_values(cells, baseline_mode, treatment, metric):
    baseline = {r["repetition"]: metric(r) for r in cells[baseline_mode]}
    compared = {r["repetition"]: metric(r) for r in cells[treatment]}
    repetitions = sorted(set(baseline) & set(compared))
    pairs = [(baseline[r], compared[r]) for r in repetitions]
    deltas = [right - left for left, right in pairs]
    return {
        "pairs": len(pairs),
        "baseline_median": statistics.median(left for left, _ in pairs),
        "treatment_median": statistics.median(right for _, right in pairs),
        "paired_absolute_delta_median": statistics.median(deltas),
        "paired_ratio_median": statistics.median(right / left for left, right in pairs),
        "deltas": deltas,
        "ratios": [right / left for left, right in pairs],
    }


def bootstrap_median_upper(ratios, seed, iterations=10000):
    rng = random.Random(seed)
    medians = [
        statistics.median(rng.choice(ratios) for _ in ratios) for _ in range(iterations)
    ]
    return percentile(medians, 0.95)


def administrative_cpu_nanos(record):
    report = record["map_info"]["thread_broker"]
    return (
        report["controller_cpu_nanos"]
        + (report["producer_measurement"].get("sampler_cpu_nanos") or 0)
    )


def summarize(records, manifest, variants):
    summary = {
        "tier": manifest.get("tier", "custom"),
        "jobs": {},
        "gate_passed": True,
        "failures": [],
    }
    minimum_repetitions = manifest.get("minimum_repetitions", 30)
    formal_gate = manifest.get("formal_gate", True)
    summary["formal_gate"] = formal_gate
    summary["minimum_repetitions"] = minimum_repetitions
    metrics = {
        "mapping_wall_seconds": lambda r: r["map_info"]["mapping_seconds"],
        "process_wall_seconds": lambda r: r["process"]["wall_seconds"],
        "process_cpu_seconds": lambda r: r["process"]["user_seconds"]
        + r["process"]["system_seconds"],
    }
    valid = [r for r in records if r.get("exit_code") == 0 and "map_info" in r]
    baseline_mode = manifest.get("baseline_variant", "fixed")
    for job_index, job in enumerate(manifest["jobs"]):
        repetitions = job.get("repetitions", manifest["repetitions"])
        statistical_gate = job.get("statistical_gate", formal_gate)
        fractional_gate = job.get("fractional_gate", statistical_gate)
        job_minimum = job.get("minimum_repetitions", minimum_repetitions)
        cells = {
            variant["name"]: [
                r
                for r in valid
                if r["job"] == job["name"] and r["mode"] == variant["name"]
            ]
            for variant in variants
        }
        result = {"modes": {name: len(cell) for name, cell in cells.items()}}
        result["block_position_counts"] = {
            name: {
                str(position): sum(r["block_position"] == position for r in cell)
                for position in range(len(variants))
            }
            for name, cell in cells.items()
        }
        complete = all(len(cell) == repetitions for cell in cells.values())
        digests = {r["canonical_output_sha256"] for cell in cells.values() for r in cell}
        counts = {
            (r["map_info"]["num_reads"], r["map_info"]["num_mapped"]) 
            for cell in cells.values()
            for r in cell
        }
        result["content_equal"] = len(digests) == 1 and len(counts) == 1
        result["evidence_level"] = "statistical" if statistical_gate else "descriptive"
        result["fractional_gate"] = fractional_gate
        result["repetitions"] = repetitions
        result["minimum_repetitions"] = job_minimum
        result["comparisons"] = {}
        margin = job.get("overhead_margin", manifest.get("overhead_margin", 0.01))
        treatments = [
            variant["name"] for variant in variants if variant["name"] != baseline_mode
        ]
        for treatment_index, treatment in enumerate(treatments):
            comparison = {}
            if complete:
                for metric_index, (name, metric) in enumerate(metrics.items()):
                    values = paired_values(cells, baseline_mode, treatment, metric)
                    values["paired_ratio_upper_95"] = bootstrap_median_upper(
                        values.pop("ratios"),
                        manifest["seed"] + job_index * 101 + treatment_index * 17 + metric_index,
                    )
                    values["paired_absolute_delta_upper_95"] = bootstrap_median_upper(
                        values.pop("deltas"),
                        manifest["seed"] + 10_000 + job_index * 103
                        + treatment_index * 19 + metric_index,
                    )
                    comparison[name] = values
            result["comparisons"][treatment] = comparison

        result["policy_telemetry"] = {}
        adaptive_modes = [
            variant["name"]
            for variant in variants
            if variant["allocation"] == "adaptive"
        ]
        for mode_index, mode in enumerate(adaptive_modes):
            reports = [r["map_info"]["thread_broker"] for r in cells[mode]]
            if not reports:
                continue
            admin_cpu = [administrative_cpu_nanos(r) for r in cells[mode]]
            admin_rate = [
                nanos / 1e9 / r["map_info"]["mapping_seconds"]
                for nanos, r in zip(admin_cpu, cells[mode])
            ]
            admin_process_fraction = [
                nanos
                / 1e9
                / (r["process"]["user_seconds"] + r["process"]["system_seconds"])
                for nanos, r in zip(admin_cpu, cells[mode])
            ]
            result["policy_telemetry"][mode] = {
                "controller_samples_median": statistics.median(
                    report["controller_samples"] for report in reports
                ),
                "controller_elapsed_median": statistics.median(
                    duration_seconds(report["controller_elapsed"]) for report in reports
                ),
                "producer_calibration_samples_median": statistics.median(
                    report["producer_measurement"]["calibration_samples"]
                    for report in reports
                ),
                "producer_monitoring_samples_median": statistics.median(
                    report["producer_measurement"]["monitoring_samples"]
                    for report in reports
                ),
                "producer_observation_nanos_median": statistics.median(
                    report["producer_measurement"]["observation_nanos"]
                    for report in reports
                ),
                "administrative_cpu_nanos_median": statistics.median(admin_cpu),
                "administrative_cpu_nanos_upper_95": bootstrap_median_upper(
                    admin_cpu, manifest["seed"] + job_index * 131 + mode_index * 19
                ),
                # Fraction of one fully occupied core, not of all process CPU.
                # 0.001 is one millisecond of broker CPU per wall second.
                "administrative_cpu_per_mapping_wall_median": statistics.median(
                    admin_rate
                ),
                "administrative_cpu_per_mapping_wall_upper_95": bootstrap_median_upper(
                    admin_rate,
                    manifest["seed"] + job_index * 137 + mode_index * 23,
                ),
                "administrative_fraction_of_process_cpu_median": statistics.median(
                    admin_process_fraction
                ),
            }
        freeze_reports = [r["map_info"]["thread_broker"] for r in cells.get("freeze", [])]
        if freeze_reports:
            result["freeze_controller_samples_median"] = statistics.median(
                r["controller_samples"] for r in freeze_reports
            )
            result["freeze_controller_elapsed_median"] = statistics.median(
                duration_seconds(r["controller_elapsed"]) for r in freeze_reports
            )

        position_balanced = all(
            max(counts.values(), default=0) - min(counts.values(), default=0) <= 1
            for counts in result["block_position_counts"].values()
        )
        result["block_positions_balanced"] = position_balanced
        passed = complete and result["content_equal"] and position_balanced
        if statistical_gate and repetitions < job_minimum:
            passed = False
            summary["failures"].append(
                "%s statistical gate requires at least %d repetitions"
                % (job["name"], job_minimum)
            )
        if complete:
            for treatment, comparison in result["comparisons"].items():
                if fractional_gate:
                    for metric in metrics:
                        if comparison[metric]["paired_ratio_upper_95"] > 1.0 + margin:
                            passed = False
                            summary["failures"].append(
                                "%s/%s/%s upper ratio %.6f exceeded %.6f"
                                % (
                                    job["name"],
                                    treatment,
                                    metric,
                                    comparison[metric]["paired_ratio_upper_95"],
                                    1.0 + margin,
                                )
                            )
                absolute_limit = job.get("absolute_wall_overhead_seconds")
                if statistical_gate and absolute_limit is not None:
                    for metric in ("mapping_wall_seconds", "process_wall_seconds"):
                        observed = comparison[metric]["paired_absolute_delta_upper_95"]
                        if observed > absolute_limit:
                            passed = False
                            summary["failures"].append(
                                "%s/%s/%s absolute upper %.6f s exceeded %.6f s"
                                % (job["name"], treatment, metric, observed, absolute_limit)
                            )
            responsive_limit = manifest.get(
                "responsive_administrative_cpu_per_wall_second", 0.001
            )
            responsive_modes = [
                variant["name"]
                for variant in variants
                if variant.get("policy") == "responsive"
            ]
            for mode in responsive_modes:
                observed = result["policy_telemetry"][mode][
                    "administrative_cpu_per_mapping_wall_upper_95"
                ]
                if observed > responsive_limit:
                    passed = False
                    summary["failures"].append(
                        "%s/%s administrative CPU rate %.9f cores exceeded %.9f"
                        % (job["name"], mode, observed, responsive_limit)
                    )
            freeze_limit = manifest.get("freeze_administrative_cpu_nanos", 5_000_000)
            freeze_modes = [
                variant["name"]
                for variant in variants
                if variant.get("policy", "").startswith("freeze-")
            ]
            for mode in freeze_modes:
                observed = result["policy_telemetry"][mode][
                    "administrative_cpu_nanos_upper_95"
                ]
                if observed > freeze_limit:
                    passed = False
                    summary["failures"].append(
                        "%s/%s administrative CPU %.0f ns exceeded %.0f ns"
                        % (job["name"], mode, observed, freeze_limit)
                    )
        if not complete:
            summary["failures"].append("%s has incomplete cells" % job["name"])
        if not result["content_equal"]:
            summary["failures"].append("%s canonical output/counts differ" % job["name"])
        if not position_balanced:
            summary["failures"].append("%s block positions are not balanced" % job["name"])
        result["gate_passed"] = passed
        summary["jobs"][job["name"]] = result
        summary["gate_passed"] &= passed

    # Freeze is intended to turn a scaling cost into a fixed startup cost.
    for group in sorted(
        {j.get("scale_group") for j in manifest["jobs"] if j.get("scale_group")}
    ):
        if baseline_mode != "fixed":
            continue
        jobs = [j for j in manifest["jobs"] if j.get("scale_group") == group]
        if any(
            not summary["jobs"][job["name"]].get("comparisons", {}).get("freeze")
            for job in jobs
        ):
            continue
        jobs.sort(
            key=lambda j: summary["jobs"][j["name"]]
            .get("comparisons", {})
            .get("freeze", {})
            .get("mapping_wall_seconds", {})
            .get("baseline_median", 0)
        )
        samples = [
            summary["jobs"][j["name"]].get("freeze_controller_samples_median") for j in jobs
        ]
        if any(value is None for value in samples):
            continue
        ratio = max(samples) / max(1, min(samples))
        allowed = manifest.get("freeze_sample_scaling_margin", 1.25)
        if ratio > allowed:
            summary["gate_passed"] = False
            summary["failures"].append(
                "%s freeze sample count scaled %.3fx (allowed %.3fx)" % (group, ratio, allowed)
            )
        if len(jobs) >= 2:
            shortest = summary["jobs"][jobs[0]["name"]]["comparisons"]["freeze"]
            longest = summary["jobs"][jobs[-1]["name"]]["comparisons"]["freeze"]
            ratio_tolerance = manifest.get("freeze_fraction_scaling_tolerance", 0.002)
            absolute_tolerance = manifest.get("freeze_absolute_scaling_seconds", 0.10)
            for metric in metrics:
                short_ratio = shortest[metric]["paired_ratio_median"]
                long_ratio = longest[metric]["paired_ratio_median"]
                if long_ratio > short_ratio + ratio_tolerance:
                    summary["gate_passed"] = False
                    summary["failures"].append(
                        "%s/%s freeze fractional overhead grew %.6f -> %.6f"
                        % (group, metric, short_ratio, long_ratio)
                    )
                short_delta = max(0.0, shortest[metric]["paired_absolute_delta_median"])
                long_delta = max(0.0, longest[metric]["paired_absolute_delta_median"])
                if long_delta > short_delta + absolute_tolerance:
                    summary["gate_passed"] = False
                    summary["failures"].append(
                        "%s/%s freeze absolute overhead grew %.6f -> %.6f s"
                        % (group, metric, short_delta, long_delta)
                    )
    return summary


def schedule(manifest, variants):
    rng = random.Random(manifest["seed"])
    specs = []
    variants = list(variants)
    # Build crossover-balanced variant orders per job. In every complete group
    # of four repetitions each variant occupies every block position exactly
    # once; random shuffling alone does not guarantee that balance.
    variant_orders = {}
    for job in manifest["jobs"]:
        repetitions = job.get("repetitions", manifest["repetitions"])
        orders = []
        while len(orders) < repetitions:
            base = list(variants)
            rng.shuffle(base)
            for offset in range(len(base)):
                orders.append(base[offset:] + base[:offset])
                if len(orders) == repetitions:
                    break
        variant_orders[job["name"]] = orders
    maximum_repetitions = max(
        job.get("repetitions", manifest["repetitions"])
        for job in manifest["jobs"]
    )
    for repetition in range(maximum_repetitions):
        jobs = [
            job for job in manifest["jobs"]
            if repetition < job.get("repetitions", manifest["repetitions"])
        ]
        rng.shuffle(jobs)
        for job in jobs:
            order = variant_orders[job["name"]][repetition]
            names = [variant["name"] for variant in order]
            for position, variant in enumerate(order):
                specs.append(
                    {
                        "job": job,
                        "variant": variant,
                        "repetition": repetition,
                        "block_order": names,
                        "block_position": position,
                    }
                )
    return specs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--summarize-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    root = Path(manifest["results_dir"])
    root.mkdir(parents=True, exist_ok=True)
    records_path = root / "records.jsonl"
    provenance = {
        "binary_sha256": hashlib.sha256(Path(manifest["binary"]).read_bytes()).hexdigest(),
        "dependency_revisions": dependency_revisions(),
        "manifest": str(args.manifest.resolve()),
    }
    records = []
    if args.summarize_only:
        records = [json.loads(line) for line in records_path.read_text().splitlines()]
    else:
        variants = selected_variants(manifest)
        specs = schedule(manifest, variants)
        completed = set()
        if args.resume:
            if not records_path.exists():
                raise ValueError("--resume requested but records.jsonl does not exist")
            previous = [
                json.loads(line) for line in records_path.read_text().splitlines() if line
            ]
            records = [
                record
                for record in previous
                if record.get("exit_code") == 0 and "map_info" in record
            ]
            if any(
                record.get("binary_sha256") != provenance["binary_sha256"]
                for record in records
            ):
                raise ValueError("cannot resume with a different release binary")
            completed = {
                (record["job"], record["mode"], record["repetition"])
                for record in records
            }
            print(
                "[resume] retained %d successful cells; failed/incomplete cells will rerun"
                % len(records),
                flush=True,
            )
        elif records_path.exists():
            raise ValueError(
                "records.jsonl already exists; use --resume or a new results_dir"
            )
        with records_path.open("w") as stream:
            for record in records:
                stream.write(json.dumps(record, sort_keys=True) + "\n")
            stream.flush()
            for index, spec in enumerate(specs, 1):
                key = (
                    spec["job"]["name"],
                    spec["variant"]["name"],
                    spec["repetition"],
                )
                if key in completed:
                    continue
                print(
                    "[%d/%d] %s %s rep %d"
                    % (
                        index,
                        len(specs),
                        spec["job"]["name"],
                        spec["variant"]["name"],
                        spec["repetition"],
                    ),
                    flush=True,
                )
                try:
                    record = execute(spec, manifest, provenance)
                except Exception as error:
                    record = {
                        **provenance,
                        "job": spec["job"]["name"],
                        "mode": spec["variant"]["name"],
                        "repetition": spec["repetition"],
                        "exit_code": -1,
                        "error": repr(error),
                    }
                records.append(record)
                stream.write(json.dumps(record, sort_keys=True) + "\n")
                stream.flush()
    variants = selected_variants(manifest)
    summary = summarize(records, manifest, variants)
    (root / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["gate_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
