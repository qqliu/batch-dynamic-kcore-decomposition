# Concurrent Batch-Dynamic K-Core Decomposition Algorithms
--------

Run

```
sh run_experiments.sh
```

from within the [CLDS scripts](../../../scripts/cplds_experiments/)
directory to see a sample of the experiments we ran in our paper.

Concurrent PLDS Algorithm and Experiments
--------

This repository contains code for our concurrent, parallel batch-dynamic k-core
decomposition algorithms. Our code for our concurrent, parallel batch-dynamic data
structures uses the framework of the [Graph-Based Benchmark Suite (GBBS)](https://github.com/ParAlg/gbbs).
The `gbbs/benchmarks/EdgeOrientation/` directory ([EdgeOrientation Directory Link](./gbbs/benchmarks/EdgeOrientation/README.md)) contains all relevant information
to our k-core decomposition algorithms, including how to run the experiments
from our paper, while the README within `gbbs` explains details about GBBS overall.

To run our script within [`gbbs/scripts/`](./gbbs/scripts):

```
python3 test_approx_kcore.py
```

The script outputs three different output files for each of the three different
programs: `CLDS_dblp_insertion_edges_0.4_3_1000000_60_.out`,
`Nonlin_dblp_insertion_edges_0.4_3_1000000_60_.out`,
`Delayed_dblp_insertion_edges_0.4_3_1000000_60_.out`.

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
are `CPLDS`, `NonSync`, and `SyncReads`.

`Numbers of workers:` Comma separated list of number of threads you want
to test, e.g., `8, 15`.

`Num Reader Threads:` Comma separated list of number of reader threads you want to
test, e.g. `8, 15`.

`Deltas:` Comma separated list of deltas you want to test, e.g., `0.4, 0.8,
1.6`. (Can only accomodate positive values.)

`Lambdas:` Comma separated list of lambdas you want to test, e.g., `3, 6, 12`.
(Can only accomodate positive values.)

`Batch sizes:` Comma separated list of batch sizes you want to test, e.g.,
`1000000, 10000000`.

`Output stats:` `True` if you want to output the error ratios for both reads and
writes; `False` if not.

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

`Percentile:` Comma separated list of percentiles for returning the latency, e.g.,
`0.95, 0.999` for returning the 95th and 99.9th percentile latencies.

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
wget https://storage.googleapis.com/k-core-decomposition/edge_orientation_exps/<GRAPH_NAME>
```

The following graphs are available for download:

```
brain_deletion_edges
brain_insertion_edges
ctr_deletion_edges
ctr_insertion_edges
dblp_deletion_edges
dblp_insertion_edges
livejournal_deletion_edges
livejournal_insertion_edges
orkut_deletion_edges
orkut_insertion_edges
stackoverflow_deletion_edges
stackoverflow_insertion_edges
twitter_deletion_edges
twitter_insertion_edges
usa_deletion_edges
usa_insertion_edges
wiki_deletion_edges
wiki_insertion_edges
youtube_deletion_edges
youtube_insertion_edges
```

`_insertion_edges` contain edges formatted for insertions. `_deletion_edges`
contain edges formatted for deletions.
