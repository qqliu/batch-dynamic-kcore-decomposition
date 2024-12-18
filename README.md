# Batch-Dynamic K-Core Decomposition Algorithms
--------

[![DOI](https://zenodo.org/badge/372586195.svg)](https://zenodo.org/doi/10.5281/zenodo.10253798)

First run

```
sh setup.sh
```
to setup the environment for running our experimental code.

Directory Contents
--------

This directory contains the source code for two different papers (citations in bibtex are included): 

1) Liu, Quanquan C., Jessica Shi, Shangdi Yu, Laxman Dhulipala, and Julian Shun. "Parallel batch-dynamic algorithms for k-core decomposition and related graph problems." In Proceedings of the 34th ACM Symposium on Parallelism in Algorithms and Architectures, pp. 191-204. 2022.

```
@inproceedings{LSYDS2022parallel,
  title={Parallel batch-dynamic algorithms for k-core decomposition and related graph problems},
  author={Liu, Quanquan C and Shi, Jessica and Yu, Shangdi and Dhulipala, Laxman and Shun, Julian},
  booktitle={Proceedings of the 34th ACM Symposium on Parallelism in Algorithms and Architectures},
  pages={191--204},
  year={2022}
}
```

2) Liu, Quanquan C., Julian Shun, and Igor Zablotchi. "Parallel k-Core Decomposition with Batched Updates and Asynchronous Reads." In Proceedings of the 29th ACM SIGPLAN Annual Symposium on Principles and Practice of Parallel Programming, pp. 286-300. 2024.

```
@inproceedings{LSZ2024parallel,
  title={Parallel k-Core Decomposition with Batched Updates and Asynchronous Reads},
  author={Liu, Quanquan C and Shun, Julian and Zablotchi, Igor},
  booktitle={Proceedings of the 29th ACM SIGPLAN Annual Symposium on Principles and Practice of Parallel Programming},
  pages={286--300},
  year={2024}
}
```

Code for the parallel batch-dynamic k-core decomposition paper are given in the [ParallelLDS directory](https://github.com/qqliu/batch-dynamic-kcore-decomposition/tree/master/gbbs/benchmarks/EdgeOrientation/) and code for the concurrent-reads CLDS code are given in the
[ConcurrentLDS directory](https://github.com/qqliu/batch-dynamic-kcore-decomposition/tree/master/gbbs/benchmarks/EdgeOrientation/ConcurrentPLDS).

PLDS Algorithm and Experiments
--------

This repository contains code for our parallel batch-dynamic k-core
decomposition algorithms. Our code for our parallel batch-dynamic data
structures uses the framework of the [Graph-Based Benchmark Suite (GBBS)](https://github.com/ParAlg/gbbs).
The `gbbs/benchmarks/EdgeOrientation/` directory ([EdgeOrientation Directory Link](./gbbs/benchmarks/EdgeOrientation/README.md)) contains all relevant information
to our k-core decomposition algorithms, including how to run the experiments
from our paper, while the README within `gbbs` explains details about GBBS overall.

To run our script within [`gbbs/scripts/`](./gbbs/scripts):

```
python3 test_approx_kcore.py
```

The script outputs three different output files for each of the three different
programs: `PLDS_dblp_insertion_edges_0.4_3_1000000_60_.out`,
`LDS_dblp_insertion_edges_0.4_3_1000000_60_.out`,
`AKCore_dblp_insertion_edges_0.4_3_1000000_60_.out`, and an additional
`EKCore_dblp_insertion_edges_0.4_3_1000000_60_.out` for the ExactKCore
benchmark.

Updating Submodules
--------

Our code uses the [parlaylib](https://github.com/cmuparlay/parlaylib) library as a submodule. Please run the below command to update this submodule:

```
git submodule update --init --recursive
```

Customizing Experiments
--------

You can change the parameters of the experiments by changing the
[`approx_kcore_setup.txt`](./gbbs/scripts/approx_kcore_setup.txt) file within
[`gbbs/scripts/`](./gbbs/scripts).

The inputs within the file are as follows:

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

Reading Outputs
--------

You can read the outputs of the experiments by running within [`gbbs/scripts/`](./gbbs/scripts):


```
python3 read_approx_kcore_results.py
```

The result of the experiments is printed in your terminal as a comma-separated list of values.
A value is `N/A` if a setting is set to False. Note that this script reads the results of the
experiments configured by the `approx_kcore_setup.txt` file.

The order of the output is as follows:

```
Program, Dynamic graph file name, Delta, Lambda, Batch size, Number of hyper-threads, Levels divisor, Average runtime, Max runtime, Total runtime, Space (in Bytes)[or N/A], Average error [or N/A], Max error [or N/A]
```

Obtaining Input Graphs
--------

You can download the graphs we used in our experiments using the following link:

```
wget https://storage.googleapis.com/graph-files/<GRAPH_NAME>
```

The following graphs are available for download:

```
brain_batch
brain_deletion_edges
brain_insertion_edges
brain_initial
ctr_batch
ctr_deletion_edges
ctr_insertion_edges
ctr_initial
dblp_batch
dblp_deletion_edges
dblp_insertion_edges
dblp_initial
friendster_batch
friendster_deletion_edges
friendster_insertion_edges
friendster_initial
livejournal_batch
livejournal_deletion_edges
livejournal_insertion_edges
livejournal_initial
orkut_batch
orkut_deletion_edges
orkut_insertion_edges
orkut_initial
stackoverflow_batch
stackoverflow_deletion_edges
stackoverflow_insertion_edges
stackoverflow_initial
twitter_batch
twitter_deletion_edges
twitter_insertion_edges
twitter_initial
usa_batch
usa_deletion_edges
usa_insertion_edges
usa_initial
wiki_batch
wiki_deletion_edges
wiki_insertion_edges
wiki_initial
youtube_batch
youtube_deletion_edges
youtube_insertion_edges
youtube_initial
```

`_insertion_edges` contain edges formatted for insertions. `_deletion_edges`
contain edges formatted for deletions. `_initial` contain the initial graph for our
mixed batch experiments. `_batch` contain the batch of edge updates used
in our mixed batch (1/2 insertions and 1/2 deletions) experiments.


CLDS Algorithm and Experiments
--------

A description for how to run the experiments using the CLDS code is given in the
[ConcurrentLDS directory](https://github.com/qqliu/batch-dynamic-kcore-decomposition/tree/master/gbbs/benchmarks/EdgeOrientation/ConcurrentPLDS).
