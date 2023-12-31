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
// Note: The static graph is ignored if -i is specified.

#include "SkipList.h"

namespace gbbs {
template <class Graph>
double LDS_runner(Graph& G, commandLine P) {
  std::cout << "### Application: LDS" << std::endl;
  std::cout << "### Graph: " << P.getArgument(0) << std::endl;
  std::cout << "### Threads: " << num_workers() << std::endl;
  std::cout << "### n: " << G.n << std::endl;
  std::cout << "### m: " << G.m << std::endl;
  std::cout << "### Params: " <<  std::endl;
  std::cout << "### ------------------------------------" << std::endl;
  assert(P.getOption("-s"));

  // Load the dynamic graph
  timer t; t.start();
  uintE n = 1000;
  RunSkipList(n);
  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;

  return tt;
}
}  // namespace gbbs

generate_symmetric_main(gbbs::LDS_runner, false);
