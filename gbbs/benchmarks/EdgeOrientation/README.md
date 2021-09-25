# Parallel Batch-Dynamic k-Core Decomposition
--------

This repository contains the implementations for our
approximate k-core algorithms.

`benchmarks/EdgeOrientation/` contains the implementations for our batch-dynamic
approximate k-core algorithm. Specifically, `ParallelLDS` (PLDS) contains the
implementation for our parallel batch-dynamic
![equation](https://latex.codecogs.com/gif.latex?%282&plus;%5Cvarepsilon%29)-approximate
k-core decomposition algorithm, while `LDS` contains our implementation of the sequential 
level data structure of [Henzinger, Neumann, and Wiese](https://arxiv.org/abs/2002.10142) 
but adapted to output core numbers.

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

`Input graph directory:` The directory used to store your input graph files.

`Dynamic graphs:` Comma separated list of file names for your batches of
updates, e.g., `dblp_insertion_edges, livejournal_insertion_edges`.

`Output directory:` Directory to write the outputs of the experiments.

`Benchmarks:` Comma separated list of benchmarks you want to test. The available benchmarks
are `PLDS`, `LDS`, `EKCore` (ExactKCore), and `AKCore` (ApproxKCore).

`Numbers of workers:` Comma separated list of number of hyper-threads you want
to test, e.g., `30, 60`.

`Deltas:` Comma separated list of deltas you want to test, e.g., `0.4, 0.8,
1.6`. (Can only accomodate positive values.)

`Lambdas:` Comma separated list of lambdas you want to test, e.g., `3, 6, 12`.
(Can only accomodate positive values.)

`Batch sizes:` Comma separated list of batch sizes you want to test, e.g.,
`1000000, 10000000`.

`Output stats:` `True` if you want to output the error ratios; `False` if not.

`Output sizes:` `True` if you want to output the sizes of the data structures;
`False` if not.

`Initial Graph File:` Comma separated list of file names for experiments where
an initial file is loaded and then batches of updates are applied. This is
relevant for the deletion and mixed experiments. For the deletion experiments,
the `[graph]_insertion_edges` file is an input to this field. For the mixed
experiments, the `[graph]_initial` file is an input to this field. The comma
separated list must be in the *same order* as the list under `Dynamic graphs:` for the graphs 
you want to test.

`Number of Levels Divisor:` The divisor for the number of levels per group. The
default, 50, is used in all of our experiments for PLDSOpt given in our paper.

`Opt:` `True` if you want to use the heuristic where (2 + 3/\lambda) is set to
1.1 (so \lambda = -10/3); `False` if not. 
Note that this setting is not theoretically time efficient.

Note that the script runs ApproxKCore and ExactKCore given a dynamic graph, 
although both can be independently run using a static graph as well.

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
