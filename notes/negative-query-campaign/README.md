# Untimed SSHash lookup trace

Build with `--features query-trace` and use `PISCEM_QUERY_TRACE=/path/prefix`
with `--dict sshash -t 1`. Each worker writes `prefix.N.trace`. This captures the
byte-oriented lookup path only; it is not a timing build. `R` marks an explicit
query reset. `Q` records read position, ASCII window, and the eight LookupResult
fields in declaration order, with minimizer_found encoded as 0/1. Nonconsecutive
positions imply the same automatic engine reset as normal Piscem lookup_at.

The companion SSHash negative_pipeline harness reconstructs the queried portions
of each read, preserves request order, and checks every captured result against
the corrected scalar stream before testing scheduler variants. Empty reset groups
are ignored. Unqueried sequence regions stay N and are never queried in replay.
The existing mapper algorithm and default builds are unchanged.

Trace capture was checked against byte-identical RAD output for 10,000 pairs each
of simulated, random and real reads. Capture cost is excluded from performance
claims. No batched mapping throughput is claimed from replay alone.

Use a fresh prefix: opening a worker's trace truncates an existing file of the
same name. A missing environment variable disables file output even in a trace
build. Files are buffered per worker and closed at worker exit; I/O failures
while writing abort capture. Trace files may be large and are not committed.
The attached `validation.json` contains the RAD SHA-256 values. Dataset labels
retain their source-file names (`*100k`), but each capture used only the first
10,000 pairs, not the entire 100,000-pair source.

The checked baseline is Piscem `9d3bd92` with SSHash main `644fbdb` (the merged
streaming-offset correction) plus its experimental query code. The local Cargo
patch used for capture is a command-line build override; it is not committed to
this repository or its lockfile. These changes add diagnostic capture only.
