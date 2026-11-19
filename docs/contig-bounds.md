# Contig bounds during mapping

Mapping commands can replace Elias–Fano contig offsets with a fixed-width packed
array. The entries themselves are moved, not copied. The on-disk index format,
index construction, ordinary `ReferenceIndex::load`, and saving remain unchanged.

`--dense-bounds-budget-mib` sets the maximum new allocation for this conversion.
The default is 64 MiB; zero retains EF. Tiny indexes use their existing table and
ignore this setting. The budget is an allowance for this one allocation, not a
whole-process memory limit or an estimate of free RAM.

For N boundaries of width w bits on a W-bit machine, the decision uses
`max(1, ceil(N*w/W)) * sizeof(usize)` bytes, with checked arithmetic. It requires
only metadata and is made before allocating or iterating through the offsets.
The original EF remains live while converting, then is dropped. Consequently the
cap bounds the temporary memory increment as well as the retained increase.

There is no percentage-of-contig-table condition. That condition excluded useful
cases merely because their occurrence entries were small. In the review data,
the dense allocation is 9.23 MiB for full GENCODE, 2.38 MiB for one million
singleton contigs, 30.99 MiB for ten million mixed contigs, and 321.87 MiB for one
hundred million singleton contigs. The default admits the first three. Users who
can afford the larger allocation can raise the budget. These are measured sizes,
not evidence that 64 MiB is an optimal runtime threshold.

## Dispatch and library callers

`ContigTable` retains the original EF representation and lookup implementation.
`DenseContigTable` is a separate type implementing `ContigTableLike`.
`ReferenceIndex::try_into_dense(max_bytes)` consumes an EF index and returns either
the dense index or the original EF index. Match this result **outside** the
mapping pipeline and call a function generic over `ContigTableLike` in each arm.
The mapping commands do this before starting workers, after selecting any Tiny
representation. There is no representation test per lookup or per read, and no
indirect call to select an implementation.

A library caller that continues using `ReferenceIndex::load` gets exactly its
existing EF representation. Adopting automatic selection requires a one-time
match around its generic mapping entry point; it must not put an enum or trait
object back into the hot loop. This is relevant to downstream Salmon integration.

## Limits of the performance claim

The fallback performs the original lookup work and retains the original table
layout. The load-time decision still takes constant time. This is a structural
property, not a guarantee of identical wall-clock time: binary layout, scheduling,
and cache effects can change measurements.

Dense conversion has a startup cost and may increase cache pressure. An unknown
workload with few or no successful queries cannot be guaranteed to amortize that
cost. No size-only automatic choice can guarantee lower total runtime for every
workload. The budget bounds resource use; workload measurements determine whether
acceleration is worthwhile. Zero is available when avoiding conversion matters
more than possible throughput gains.
