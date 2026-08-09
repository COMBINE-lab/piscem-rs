#!/usr/bin/env python3
"""End-to-end FIFO/mixed-input correctness matrix for map-bulk.

This is a machine-local qualification harness, not a portable repository test.
It feeds every FIFO from a one-shot writer: any speculative open consumes that
only stream and makes the real open hang, which the per-run timeout catches.
"""

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import threading
import time


def make_prefix(source, destination, records):
    with gzip.open(source, "rb") as reader, destination.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=1, mtime=0) as writer:
            for record in range(records):
                for line in range(4):
                    value = reader.readline()
                    if not value:
                        raise RuntimeError(
                            f"{source} ended at record {record}, FASTQ line {line}"
                        )
                    writer.write(value)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_mapping(
    *, binary, index, output, read1, read2, decoder, threads, fifo_payloads, timeout
):
    for fifo in fifo_payloads:
        os.mkfifo(fifo)

    writer_errors = []

    def write_once(fifo, payload):
        try:
            with fifo.open("wb", buffering=0) as destination, payload.open("rb") as source:
                shutil.copyfileobj(source, destination, length=1024 * 1024)
        except BaseException as error:  # retained in the result, including BrokenPipeError
            writer_errors.append(f"{fifo}: {type(error).__name__}: {error}")

    writers = [
        threading.Thread(target=write_once, args=(fifo, payload), daemon=True)
        for fifo, payload in fifo_payloads.items()
    ]
    for writer in writers:
        writer.start()

    command = [
        str(binary),
        "map-bulk",
        "-i",
        str(index),
        "-1",
        ",".join(map(str, read1)),
        "-2",
        ",".join(map(str, read2)),
        "-o",
        str(output),
        "-t",
        str(threads),
        "--decoder",
        decoder,
        "--dict",
        "sshash",
        "--quiet",
    ]
    environment = os.environ.copy()
    environment["RUST_LOG"] = "info"
    started = time.monotonic()
    timed_out = False
    try:
        completed = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=environment,
            timeout=timeout,
            check=False,
        )
        exit_code = completed.returncode
        stdout = completed.stdout.decode(errors="replace")
        stderr = completed.stderr.decode(errors="replace")
    except subprocess.TimeoutExpired as error:
        timed_out = True
        exit_code = None
        stdout = (error.stdout or b"").decode(errors="replace")
        stderr = (error.stderr or b"").decode(errors="replace")

    for writer in writers:
        writer.join(timeout=2)
    lingering_writers = sum(writer.is_alive() for writer in writers)
    for fifo in fifo_payloads:
        try:
            fifo.unlink()
        except FileNotFoundError:
            pass

    map_info_path = Path(f"{output}.map_info.json")
    rad_path = Path(f"{output}.rad")
    info = json.loads(map_info_path.read_text()) if map_info_path.exists() else None
    return {
        "command": command,
        "seconds": time.monotonic() - started,
        "exit_code": exit_code,
        "timed_out": timed_out,
        "stdout": stdout,
        "stderr": stderr,
        "writer_errors": writer_errors,
        "lingering_writers": lingering_writers,
        "map_info": info,
        "rad_sha256": sha256(rad_path) if rad_path.exists() else None,
        "rad_path": str(rad_path),
    }


def assert_success(label, result):
    assert not result["timed_out"], f"{label}: timed out; possible speculative FIFO open"
    assert result["exit_code"] == 0, f"{label}: exit {result['exit_code']}: {result['stderr']}"
    assert not result["writer_errors"], f"{label}: one-shot writer failed: {result['writer_errors']}"
    assert result["lingering_writers"] == 0, f"{label}: FIFO writer remains blocked"
    assert result["map_info"] is not None, f"{label}: missing map_info"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=Path("target/release/piscem-rs"))
    parser.add_argument(
        "--index",
        type=Path,
        required=True,
    )
    parser.add_argument("--r1", type=Path, default=Path("/tmp/tb-bulk-pe-R1.fq.gz"))
    parser.add_argument("--r2", type=Path, default=Path("/tmp/tb-bulk-pe-R2.fq.gz"))
    parser.add_argument("--records", type=int, default=20_000)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--timeout", type=float, default=120.0)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("/tmp/thread-broker-fifo-matrix")
    )
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    prefix1 = args.output_dir / "prefix-R1.fq.gz"
    prefix2 = args.output_dir / "prefix-R2.fq.gz"
    make_prefix(args.r1, prefix1, args.records)
    make_prefix(args.r2, prefix2, args.records)

    scenarios = {
        "fifo_only": {
            "baseline": ([prefix1], [prefix2]),
            "treatment": lambda d: (
                [d / "fifo-only-R1.fq.gz"],
                [d / "fifo-only-R2.fq.gz"],
            ),
        },
        "mixed": {
            "baseline": ([prefix1, prefix1], [prefix2, prefix2]),
            "treatment": lambda d: (
                [prefix1, d / "mixed-R1.fq.gz"],
                [prefix2, d / "mixed-R2.fq.gz"],
            ),
        },
        "split_group": {
            "baseline": ([prefix1], [prefix2]),
            "treatment": lambda d: ([prefix1], [d / "split-R2.fq.gz"]),
        },
    }

    evidence = {
        "records_per_group": args.records,
        "threads": args.threads,
        "cells": {},
    }
    for scenario, spec in scenarios.items():
        evidence["cells"][scenario] = {}
        for decoder in ("auto", "parallel"):
            baseline_r1, baseline_r2 = spec["baseline"]
            baseline = run_mapping(
                binary=args.binary.resolve(),
                index=args.index,
                output=args.output_dir / f"{scenario}-{decoder}-regular",
                read1=baseline_r1,
                read2=baseline_r2,
                decoder=decoder,
                threads=args.threads,
                fifo_payloads={},
                timeout=args.timeout,
            )
            assert_success(f"{scenario}/{decoder}/regular", baseline)

            treatment_r1, treatment_r2 = spec["treatment"](args.output_dir)
            fifo_payloads = {
                path: prefix1 if "R1" in path.name else prefix2
                for path in treatment_r1 + treatment_r2
                if path != prefix1 and path != prefix2
            }
            treatment = run_mapping(
                binary=args.binary.resolve(),
                index=args.index,
                output=args.output_dir / f"{scenario}-{decoder}-fifo",
                read1=treatment_r1,
                read2=treatment_r2,
                decoder=decoder,
                threads=args.threads,
                fifo_payloads=fifo_payloads,
                timeout=args.timeout,
            )
            assert_success(f"{scenario}/{decoder}/fifo", treatment)

            comparable_fields = ("num_reads", "num_mapped", "num_poisoned")
            for field in comparable_fields:
                assert treatment["map_info"][field] == baseline["map_info"][field], (
                    f"{scenario}/{decoder}: {field} differs"
                )
            assert treatment["rad_sha256"] == baseline["rad_sha256"], (
                f"{scenario}/{decoder}: RAD bytes differ from regular-file baseline"
            )

            logs = treatment["stdout"] + treatment["stderr"]
            fifo_names = [path.name for path in fifo_payloads]
            assert all(name in logs for name in fifo_names), (
                f"{scenario}/{decoder}: message does not name every stream: {logs}"
            )
            if decoder == "auto":
                assert " INFO " in logs and "not regular files" in logs, logs
                assert " WARN " not in logs, logs
            else:
                assert " WARN " in logs and "--decoder parallel was requested" in logs, logs

            plan = treatment["map_info"]["execution_plan"]
            if scenario == "fifo_only":
                assert plan["allocation"]["kind"] == "serial", plan
                assert plan["decode_slots"] == 0, plan
                assert plan["map_threads"] == plan["effective_budget"], plan
            else:
                assert plan["allocation"]["kind"] == "adaptive", plan
                assert plan["decode_slots"] > 0, plan
                assert plan["map_threads"] + plan["decode_slots"] == plan["effective_budget"], plan

            evidence["cells"][scenario][decoder] = {
                "baseline": baseline,
                "treatment": treatment,
                "byte_identical": True,
            }
            (args.output_dir / "results.json").write_text(
                json.dumps(evidence, indent=2, sort_keys=True) + "\n"
            )

    print(json.dumps({"passed": True, "results": str(args.output_dir / "results.json")}))


if __name__ == "__main__":
    main()
