// Usage (dynamic):
// numactl -i all ./LDS -rounds 3 -s -eps 0.4 -delta 3 -i <dynamic graph> -b 1000 <static graph>
// flags:
//   required:
//     -s : indicates that the graph is symmetric
//   optional:
//     -m : indicate that the graph should be mmap'd
//     -c : indicate that the graph is compressed
//     -rounds : the number of times to run the algorithm
//     -i : file path to the dynamic graph, with a single dynamic edge on each
//          line, in the format "<+/-> u v", where + denotes insertion, -
//          denotes deletion, and (u, v) are the vertex ids denoting the edge
//          (e.g., "+ 0 1")
//     -b : batch size (output includes time per batch, from start_size to
//          end_size)
//     -eps: epsilon
//     -delta: lambda
//     -stats : indicates whether to output comparisons to exact coreness
//              values
//     -ins-opt : indicates whether to set lambda such that
//                (2 + 3 / lambda) = 1.1
//     -size : option for getting the size of the data structure in bytes
//     -opt : indicates whether to divide the number of levels by 50; faster
//     runtime, subtle decrease in accuracy
//     -init_graph_file: initial graph for which the additional edge updates
//     occur on; initial graph should follow the same format as the input to -i
//     (i.e. all edges should be prepended with '+')
// Note: The static graph is ignored if -i is specified.

#include "ETTree.h"

namespace gbbs {
template <class Graph>
double LDS_runner(Graph& G, commandLine P) {
  // Run LDS
  timer t; t.start();

  using K = std::pair<uintE, uintE>;
  using V = std::pair<bool, uintE>;
  using KV = std::pair<K, V>;

  KV empty =
       std::make_pair(std::make_pair(UINT_E_MAX, UINT_E_MAX), std::make_pair(false, UINT_E_MAX));

  auto hash_pair = [](const std::pair<uintE, uintE>& t) {
    size_t l = std::min(std::get<0>(t), std::get<1>(t));
    size_t r = std::max(std::get<0>(t), std::get<1>(t));
    size_t key = (l << 32) + r;
    return parlay::hash64_2(key);
  };

  auto edge_table =
    pbbslib::make_sparse_table<K, V>(16, empty, hash_pair);


  RunETTree(2, 2, 8, edge_table);
  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;
  return tt;
}
}  // namespace gbbs

generate_symmetric_main(gbbs::LDS_runner, false);
