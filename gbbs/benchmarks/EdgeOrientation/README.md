# Concurrent Batch-Dynamic k-Core Decomposition
--------

This repository contains the implementations for our concurrent LDS (CLDS),
synchronized reads (SynchronizedReads), and nonlinearizable reads (NonlinearizableReads)
 k-core algorithms.

`benchmarks/EdgeOrientation/` contains the implementations for our concurrent batch-dynamic
approximate k-core algorithm. Specifically, `ConcurrentLDS` (CLDS) contains the
implementation for our concurrent batch-dynamic
![equation](https://latex.codecogs.com/gif.latex?%282&plus;%5Cvarepsilon%29)-approximate
k-core decomposition algorithm, while `NonlinearizableReads` (Nonlin)
contains our implementation of our nonlinearizable read algorithm and
`SynchronizedReads` (Delayed)
contains the implementation for our synchronized reads algorithm where
reads occur at the end of the batch.

`scripts/test_approx_kcore.py` is a script to run tests for CLDS, Delayed, and Nonlin.

Compilation
--------

Compiler:
* g++ &gt;= 7.4.0 with support for Cilk Plus
* g++ &gt;= 7.4.0 with pthread support (Parlay Scheduler)

Build system:
* Make

To build:
```sh
$ cd benchmarks/EdgeOrientation/ConcurrentLDS  # go to a benchmark
$ make
```

The following commands cleans the directory:
```sh
$ make clean  # removes executables for the current directory
```

Most optionality from the [Graph Based Benchmark Suite (GBBS)](https://github.com/ParAlg/gbbs) and
[ParlayLib](https://github.com/cmuparlay/parlaylib) apply. In particular, to compile benchmarks for graphs with
more than 2^32 edges, the `LONG` command-line parameter should be set.
