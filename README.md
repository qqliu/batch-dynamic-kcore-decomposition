#Dynamic K-Core Decomposition Algorithms
--------

PLDS Algorithm and Experiments
--------

This repository contains code for our parallel batch-dynamic k-core
decomposition algorithms. Our code for our parallel batch-dynamic data
structures uses the framework of the [Graph-Based Benchmark Suite (GBBS)](https://github.com/ParAlg/gbbs).
The `gbbs/benchmarks/EdgeOrientation/` directory ([EdgeOrientation Directory Link](./gbbs/benchmarks/EdgeOrientation)) contains all relevant information
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

Customizing Experiments
--------

You can change the parameters of the experiments by changing the
[`approx_kcore_setup.txt`](./gbbs/scripts/approx_kcore_setup.txt) file within
[`gbbs/scripts/`](./gbbs/scripts).

Obtaining Input Graphs
--------

You can download the graphs we used in our experiments using the following link:

```
wget https://storage.googleapis.com/k-core-decomposition/edge_orientation_exps/<GRAPH_NAME>
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
contain edges formatted for deletions. `_edges` contain a mixture of insertion
and deletion edges. We included the two graphs `dblp_edges` and
`livejournal_edges` for which we tested our mixed insertion/deletion
experiments.
