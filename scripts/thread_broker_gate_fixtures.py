#!/usr/bin/env python3
"""Create deterministic local Gate E/G FASTQ compression archetypes."""

import argparse
import gzip
import random
from pathlib import Path


def concatenate(source, destination, repetitions):
    with destination.open("wb") as stream:
        for _ in range(repetitions):
            with source.open("rb") as input_stream:
                for chunk in iter(lambda: input_stream.read(4 * 1024 * 1024), b""):
                    stream.write(chunk)


def store_gzip(source, destination, repetitions):
    with destination.open("wb") as raw_output:
        with gzip.GzipFile(
            filename="", mode="wb", compresslevel=0, fileobj=raw_output, mtime=0
        ) as output_stream:
            for _ in range(repetitions):
                with source.open("rb") as input_stream:
                    for chunk in iter(lambda: input_stream.read(4 * 1024 * 1024), b""):
                        output_stream.write(chunk)


def alternating_fastq(destination, regions, records_per_region, read_length):
    rng = random.Random(20260807)
    bases = b"ACGT"
    random_pool = bytes(bytearray(rng.getrandbits(8) for _ in range(4 * 1024 * 1024)))
    base_table = bytes(bytearray(bases[value & 3] for value in range(256)))
    quality_table = bytes(bytearray(33 + value % 41 for value in range(256)))
    with destination.open("wb") as raw_output:
        record = 0
        for region in range(regions):
            with gzip.GzipFile(
                filename="", mode="wb", compresslevel=6, fileobj=raw_output, mtime=0
            ) as member:
                compressible = region % 2 == 0
                for _ in range(records_per_region):
                    if compressible:
                        sequence = b"A" * read_length
                        quality = b"I" * read_length
                    else:
                        offset = (record * read_length * 2) % (len(random_pool) - 2 * read_length)
                        random_bytes = random_pool[offset:offset + 2 * read_length]
                        sequence = random_bytes[:read_length].translate(base_table)
                        quality = random_bytes[read_length:].translate(quality_table)
                    member.write(
                        b"@gate-g-" + str(record).encode() + b"\n"
                        + sequence + b"\n+\n" + quality + b"\n"
                    )
                    record += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat-source", type=Path)
    parser.add_argument("--repeat-output", type=Path)
    parser.add_argument("--repeat-count", default=1, type=int)
    parser.add_argument("--bursty-source", type=Path)
    parser.add_argument("--bursty-output", type=Path)
    parser.add_argument("--bursty-repetitions", default=100, type=int)
    parser.add_argument("--plain-source", type=Path)
    parser.add_argument("--stored-output", type=Path)
    parser.add_argument("--stored-repetitions", default=10, type=int)
    parser.add_argument("--alternating-output", type=Path)
    parser.add_argument("--alternating-regions", default=24, type=int)
    parser.add_argument("--alternating-records", default=25000, type=int)
    parser.add_argument("--alternating-read-length", default=100, type=int)
    args = parser.parse_args()
    if args.repeat_source is not None or args.repeat_output is not None:
        if args.repeat_source is None or args.repeat_output is None:
            parser.error("--repeat-source and --repeat-output must be supplied together")
        if args.repeat_count <= 0:
            parser.error("--repeat-count must be positive")
        concatenate(args.repeat_source, args.repeat_output, args.repeat_count)
        print("repeated_bytes\t%d" % args.repeat_output.stat().st_size)
        return
    if any(
        value is None
        for value in (
            args.bursty_source,
            args.bursty_output,
            args.plain_source,
            args.stored_output,
        )
    ):
        parser.error(
            "the archetype mode requires --bursty-source, --bursty-output, "
            "--plain-source, and --stored-output"
        )
    concatenate(args.bursty_source, args.bursty_output, args.bursty_repetitions)
    store_gzip(args.plain_source, args.stored_output, args.stored_repetitions)
    if args.alternating_output is not None:
        alternating_fastq(
            args.alternating_output,
            args.alternating_regions,
            args.alternating_records,
            args.alternating_read_length,
        )
    print("bursty_bytes\t%d" % args.bursty_output.stat().st_size)
    print("stored_bytes\t%d" % args.stored_output.stat().st_size)
    if args.alternating_output is not None:
        print("alternating_bytes\t%d" % args.alternating_output.stat().st_size)


if __name__ == "__main__":
    main()
