#Dynamic K-Core Decomposition Algorithms
--------

PLDS Algorithm and Experiments
--------

This repository contains code for our parallel batch-dynamic k-core
decomposition algorithms. Our code for our parallel batch-dynamic data
structures uses the framework of the [Graph-Based Benchmark Suite (GBBS)](https://github.com/ParAlg/gbbs).
The `gbbs-dynamic-k-core/benchmarks/EdgeOrientation/` directory [EdgeOrientation Directory Link](https://github.com/qqliu/batch-dynamic-kcore-decomposition/tree/master/gbbs-dynamic-k-core/benchmarks/EdgeOrientation) contains all relevant information
to our k-core decomposition algorithms, including how to run the experiments
from our paper, while the README within `gbbs-dynamic-k-core' explains details about GBBS overall, notably
the required graph input formats for static benchmarks.

To run our script within [`gbbs-dynamic-k-core/scripts/'](https://github.com/qqliu/batch-dynamic-kcore-decomposition/tree/master/gbbs-dynamic-k-core/scripts):

```
python3 test_approx_kcore.py
```

The script outputs three different output files for each of the three different
programs: `ParallelLDS_dblp_insertion_edges_0.4_3_1000000_60_.out`,
`LDS_dblp_insertion_edges_0.4_3_1000000_60_.out`,
`KCore_dblp_insertion_edges_0.4_3_1000000_60_.out`.

Hua et al. Code
--------

The first author of the work [``Faster Parallel Core Maintenance Algorithms in
Dynamic Graph''](https://ieeexplore.ieee.org/document/8935160) provided us with
the experimental code they used to run their experiments. We included their code
in the `hua-et-all-code' directory and also the small number of modifications we
made to their original code to fit our experimental environment.

Sun et al. Code
--------

The authors of [``Fully Dynamic Approximate k-Core Decomposition in
Hypergraphs''](https://dl.acm.org/doi/10.1145/3385416) has a public repository
containing their code. We included our fork of the directory (as well as
additional changes to fit our experimental environment) in the
`sun-et-all-sequential-code' directory.
