# Parallel Batch-Dynamic k-Core Decomposition
--------

This repository contains the implementations for our
approximate k-core algorithms.

`benchmarks/EdgeOrientation/` contains the implementations for our batch-dynamic
approximate k-core algorithm. Specifically, `ParallelLDS` (PLDS) contains the
implementation for our parallel batch-dynamic
![equation](https://latex.codecogs.com/gif.latex?%282&plus;%5Cvarepsilon%29)-approximate
k-core decomposition algorithm, while `LDS` contains our implementation of the sequential 
level data structure of Henzinger, Neumann, and Wiese but adapted to output core numbers.

`benchmarks/KCore/ApproximateKCore` contains the implementation for our static
![equation](https://latex.codecogs.com/gif.latex?%282&plus;%5Cvarepsilon%29)-approximate
k-core decomposition algorithm.

`scripts/test_approx_kcore.py` is a script to run tests for ParallelLDS, 
LDS, and KCore.

Compilation
--------

Compiler:
* g++ &gt;= 7.4.0 with support for Cilk Plus
* g++ &gt;= 7.4.0 with pthread support (Parlay Scheduler)

Build system:
* Make

To build:
```sh
$ cd benchmarks/EdgeOrientation/ParallelLDS  # go to a benchmark
$ make
```

The following commands cleans the directory:
```sh
$ make clean  # removes executables for the current directory
```

Most optionality from the [Graph Based Benchmark Suite (GBBS)](https://github.com/ParAlg/gbbs) and
[ParlayLib](https://github.com/cmuparlay/parlaylib) apply. In particular, to compile benchmarks for graphs with
more than 2^32 edges, the `LONG` command-line parameter should be set.

### Using Testing Script

To run the script `scripts/test_approx_kcore.py`, Python 3.5 or higher is
required. Also, `scripts/approx_kcore_setup.txt` should be configured
as follows.

The entries `Input graph directory` and `Output directory` contain
the directory where the dynamic graph files are stored and the desired output 
directory respectively. `Dynamic graphs` contains a comma-separated
list of the dynamic graph filenames.

The `Benchmarks` entry is a comma-separated list of the desired programs 
to run: ParallelLDS, LDS, or KCore.

The `Numbers of workers`, `Epsilons`, `Lambdas/Deltas`, and `Batch sizes`
entries are comma-separated lists of the desired input parameters
on which to run the programs. The `Output stats` entry is either True or
False, depending on if comparisons to exact k-core values is desired.

Note that the script runs KCore given a dynamic graph, although KCore
can be independently run using a static graph as well.

### Using Standard Library Allocator

The default compilation uses a scalable memory allocator developed at CMU
(Parlay). However, it is also possible to use the C++ standard library,
by adding the compile definition `-DPARLAY_USE_STD_ALLOC`.

### Using Cilk, OpenMP, or TBB

The default compilation uses a lightweight scheduler developed at CMU (Parlay)
for parallelism, which results in comparable performance to Cilk Plus.
However, it is also possible to use a Cilk, OpenMP, or Thread Building
Blocks with our implementations instead of the Parlay scheduler; simply add
the appropriate compile definition as below.

```
-DPARLAY_CILK
-DPARLAY_OPENMP
-DPARLAY_TBB
```

### Graph Format

The applications take dynamic graphs as input (in all three benchmarks) or 
static graphs as input (in the static benchmark).

For the static graph format, we support the adjacency graph format used by
[Graph Based Benchmark Suite (GBBS)](https://github.com/ParAlg/gbbs).

For the dynamic graph format, the file should contain one dynamic edge per line
(in order), with a `+` or a `-` indicating an insertion or a deletion
respectively, followed by the two vertex IDs. We assume all
vertices are in the range [0, n]. For instance, a single line of the file
could be `+ 0 1`, to insert edge (0, 1). The file should be represented
as plain text.
