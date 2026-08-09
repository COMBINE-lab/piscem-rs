#!/usr/bin/env python3
"""Measure fixed-split decoder-observation overhead from one release binary.

The mapping/decode split is pinned for every run. The only changed input is
PISCEM_FIXED_DECODE_MEASUREMENT={off,monitoring,calibration}. A formal pass
requires at least 30 counterbalanced randomized repetitions, canonical output equivalence, and
a one-sided 95% bootstrap upper bound no greater than the configured margin for
both mapping wall time and process CPU time.
"""

import argparse
import hashlib
import itertools
import json
import os
import random
import statistics
import subprocess
import sys
import time
from pathlib import Path


MEASUREMENT_ENV = "PISCEM_FIXED_DECODE_MEASUREMENT"


def variants(manifest):
    return manifest.get(
        "variants",
        [
            {"name": "off", "measurement": "off"},
            {"name": "monitoring", "measurement": "monitoring"},
            {"name": "calibration", "measurement": "calibration"},
        ],
    )


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


def output_paths(root, job, mode, repetition):
    output = root / job["name"] / mode / ("rep-%d" % repetition)
    if job["output_kind"] == "stem":
        output.parent.mkdir(parents=True, exist_ok=True)
        return output, Path("%s.map_info.json" % output)
    if job["output_kind"] == "directory":
        output.mkdir(parents=True, exist_ok=True)
        return output, output / "map_info.json"
    raise ValueError("unknown output_kind for %s: %s" % (job["name"], job["output_kind"]))


def canonical_digest(command, output):
    if not command:
        raise ValueError("overhead gates require canonical_digest_command")
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


def validate_telemetry(info, job, variant):
    plan = info["execution_plan"]
    if plan["requested_budget"] != job["budget"] or info["num_threads"] != plan["effective_budget"]:
        raise ValueError("budget telemetry mismatch: %r" % plan)
    if plan["allocation"]["kind"] != "pinned_aggregate":
        raise ValueError("overhead run was not aggregate-pinned: %r" % plan)
    if plan["requested_decode_slots"] != job["pin"] or plan["decode_slots"] != job["pin"]:
        raise ValueError("fixed decode allocation was not honored: %r" % plan)
    if plan["map_threads"] + plan["decode_slots"] != plan["effective_budget"]:
        raise ValueError("fixed allocation violated budget: %r" % plan)
    if "thread_broker" in info:
        raise ValueError("fixed overhead run unexpectedly started the broker")

    measurement_mode = variant.get("measurement", "off")
    measurement = info.get("producer_measurement")
    if measurement_mode == "off":
        if measurement is not None:
            raise ValueError("off run unexpectedly emitted producer measurement")
        return
    if measurement is None or measurement["final_mode"] != measurement_mode:
        raise ValueError("measurement mode telemetry mismatch: %r" % measurement)
    if measurement_mode == "calibration" and measurement["calibration_samples"] == 0:
        raise ValueError("calibration run took no calibration samples")
    if measurement_mode == "monitoring" and measurement["monitoring_samples"] == 0:
        raise ValueError("monitoring run took no monitoring samples")
    if measurement["observation_nanos"] <= 0:
        raise ValueError("measurement run reported no observation cost")
    if variant.get("expect_cpu_accounting") is True:
        if measurement.get("completed_worker_cpu_nanos") is None:
            raise ValueError("CPU-accounting variant emitted no completed worker CPU")
        if measurement.get("completed_auxiliary_cpu_nanos") is None:
            raise ValueError("CPU-accounting variant emitted no completed auxiliary CPU")
        if measurement.get("cpu_accounting_failures") != 0:
            raise ValueError("CPU-accounting variant reported clock failures")
    if variant.get("expect_cpu_accounting") is False:
        if measurement.get("completed_worker_cpu_nanos") is not None:
            raise ValueError("baseline unexpectedly emitted completed worker CPU")
        if measurement.get("completed_auxiliary_cpu_nanos") is not None:
            raise ValueError("baseline unexpectedly emitted completed auxiliary CPU")


def cleanup(job, output):
    for template in job.get("cleanup_paths", []):
        path = Path(template.replace("{output}", str(output)))
        if path.is_file():
            path.unlink()


def execute(spec, manifest, provenance):
    job = spec["job"]
    variant = spec["variant"]
    mode = variant["name"]
    output, map_info = output_paths(
        Path(manifest["results_dir"]), job, mode, spec["repetition"]
    )
    argv = [
        variant.get("binary", manifest["binary"]),
        *job["args"],
        "-t",
        str(job["budget"]),
        "--decoder",
        "auto",
        job.get("output_argument", "-o"),
        str(output),
    ]
    env = os.environ.copy()
    env.pop("PISCEM_RAPIDGZIP_THREADS", None)
    env.pop("PISCEM_THREAD_BROKER_POLICY", None)
    env.pop("PISCEM_THREAD_BROKER_PROBE_INTERVAL_MS", None)
    env["PISCEM_DECODE_SLOTS"] = str(job["pin"])
    measurement_mode = variant.get("measurement", "off")
    if measurement_mode == "off":
        env.pop(MEASUREMENT_ENV, None)
    else:
        env[MEASUREMENT_ENV] = measurement_mode

    record = {
        **provenance,
        "job": job["name"],
        "mode": mode,
        "repetition": spec["repetition"],
        "budget": job["budget"],
        "pin": job["pin"],
        "block_order": spec["block_order"],
        "block_position": spec["block_position"],
        "command": argv,
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
    if result.returncode != 0:
        record["error"] = "command failed; see %s" % log_path
        return record

    info = json.loads(map_info.read_text())
    validate_telemetry(info, job, variant)
    record["map_info"] = info
    record["process"] = json.loads(time_path.read_text())
    record["canonical_output_sha256"] = canonical_digest(
        job.get("canonical_digest_command"), output
    )
    cleanup(job, output)
    return record


def percentile(values, fraction):
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(len(ordered) * fraction))]


def paired_ratios(baseline_records, treatment_records, metric):
    baseline = {record["repetition"]: metric(record) for record in baseline_records}
    treatment = {record["repetition"]: metric(record) for record in treatment_records}
    repetitions = sorted(set(baseline) & set(treatment))
    return [treatment[repetition] / baseline[repetition] for repetition in repetitions]


def bootstrap_median_upper(ratios, seed, iterations=10000):
    rng = random.Random(seed)
    medians = []
    for _ in range(iterations):
        medians.append(statistics.median(rng.choice(ratios) for _ in ratios))
    return percentile(medians, 0.95)


def summarize(records, manifest):
    margin = manifest.get("equivalence_margin", 0.01)
    configured_variants = variants(manifest)
    mode_names = [variant["name"] for variant in configured_variants]
    baseline_mode = mode_names[0]
    summary = {
        "required_repetitions": 30,
        "equivalence_margin": margin,
        "jobs": {},
        "gate_passed": True,
        "failures": [],
    }
    for job_index, job in enumerate(manifest["jobs"]):
        cells = {
            mode: [
                record
                for record in records
                if record.get("job") == job["name"]
                and record.get("mode") == mode
                and record.get("exit_code") == 0
                and "map_info" in record
            ]
            for mode in mode_names
        }
        complete = all(len(cells[mode]) == manifest["repetitions"] for mode in mode_names)
        digests = {
            record["canonical_output_sha256"]
            for mode in mode_names
            for record in cells[mode]
        }
        counts = {
            (
                record["map_info"]["num_reads"],
                record["map_info"]["num_mapped"],
                record["map_info"]["num_poisoned"],
            )
            for mode in mode_names
            for record in cells[mode]
        }
        correctness = len(digests) == 1 and len(counts) == 1
        formal_repetitions = manifest["repetitions"] >= 30
        job_summary = {
            "complete": complete,
            "canonical_correctness": correctness,
            "formal_repetitions": formal_repetitions,
            "modes": {},
            "gate_passed": complete and correctness and formal_repetitions,
        }
        if complete:
            for mode_index, mode in enumerate(mode_names[1:], 1):
                wall_ratios = paired_ratios(
                    cells[baseline_mode],
                    cells[mode],
                    lambda record: record["map_info"]["mapping_seconds"],
                )
                cpu_ratios = paired_ratios(
                    cells[baseline_mode],
                    cells[mode],
                    lambda record: record["process"]["user_seconds"]
                    + record["process"]["system_seconds"],
                )
                wall_ratio = statistics.median(wall_ratios)
                cpu_ratio = statistics.median(cpu_ratios)
                wall_upper = bootstrap_median_upper(
                    wall_ratios, 20260807 + job_index * 10 + mode_index
                )
                cpu_upper = bootstrap_median_upper(
                    cpu_ratios, 20260870 + job_index * 10 + mode_index
                )
                passed = wall_upper <= 1.0 + margin and cpu_upper <= 1.0 + margin
                job_summary["modes"][mode] = {
                    "median_mapping_ratio": wall_ratio,
                    "mapping_ratio_upper_95": wall_upper,
                    "median_cpu_ratio": cpu_ratio,
                    "cpu_ratio_upper_95": cpu_upper,
                    "equivalent": passed,
                }
                job_summary["gate_passed"] &= passed
        summary["jobs"][job["name"]] = job_summary
        if not job_summary["gate_passed"]:
            summary["gate_passed"] = False
            summary["failures"].append("%s: overhead/correctness gate incomplete or failed" % job["name"])
    if any(record.get("exit_code") != 0 or "error" in record for record in records):
        summary["gate_passed"] = False
        summary["failures"].append("one or more commands failed")
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    if manifest.get("repetitions", 0) < 3:
        raise ValueError("at least three repetitions are required")

    configured_variants = variants(manifest)
    binary_hashes = {
        variant["name"]: hashlib.sha256(
            Path(variant.get("binary", manifest["binary"])).read_bytes()
        ).hexdigest()
        for variant in configured_variants
    }
    provenance = {
        "commit": run_text(["git", "rev-parse", "HEAD"]),
        "dirty": bool(run_text(["git", "status", "--porcelain"])),
        "dependencies": dependency_revisions(),
        "binary_sha256_by_variant": binary_hashes,
    }
    rng = random.Random(manifest.get("seed", 20260807))
    blocks = []
    for job in manifest["jobs"]:
        # Counterbalance the six possible within-block mode orders before
        # randomizing block order. This prevents warm-cache or thermal order
        # effects from being mistaken for observation overhead.
        orders = list(itertools.permutations(configured_variants))
        mode_orders = [orders[index % len(orders)] for index in range(manifest["repetitions"])]
        rng.shuffle(mode_orders)
        blocks.extend(
            (job, repetition, mode_orders[repetition - 1])
            for repetition in range(1, manifest["repetitions"] + 1)
        )
    rng.shuffle(blocks)
    specs = []
    for job, repetition, ordered_variants in blocks:
        specs.extend(
            {
                "job": job,
                "mode": variant["name"],
                "variant": variant,
                "repetition": repetition,
                "block_order": [item["name"] for item in ordered_variants],
                "block_position": position,
            }
            for position, variant in enumerate(ordered_variants)
        )

    results_dir = Path(manifest["results_dir"])
    results_dir.mkdir(parents=True, exist_ok=True)
    records = []
    with (results_dir / "runs.jsonl").open("w") as raw:
        for index, spec in enumerate(specs, 1):
            print(
                "[%d/%d] %s %s rep %d"
                % (index, len(specs), spec["job"]["name"], spec["mode"], spec["repetition"]),
                flush=True,
            )
            try:
                record = execute(spec, manifest, provenance)
            except Exception as error:
                record = {
                    **provenance,
                    "job": spec["job"]["name"],
                    "mode": spec["mode"],
                    "repetition": spec["repetition"],
                    "error": repr(error),
                }
            records.append(record)
            raw.write(json.dumps(record, sort_keys=True) + "\n")
            raw.flush()
    summary = summarize(records, manifest)
    (results_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return 0 if summary["gate_passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
