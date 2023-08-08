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

#include "LDS.h"

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

  // Load dynamic graph options
  const std::string kInputFlag{"-i"};
  const char* const input_file{P.getOptionValue(kInputFlag)};
  long batch_size = P.getOptionLongValue("-b", 1);
  bool compare_exact = P.getOption("-stats");
  bool nonlinearizable = P.getOption("-nonlin");
  bool optimized_insertion = P.getOption("-ins-opt");
  bool get_size = P.getOption("-size");

  // Options for the approximation algorithm
  double eps = P.getOptionDoubleValue("-eps", 3);
  double delta = P.getOptionDoubleValue("-delta", 9);
  double percentile = P.getOptionDoubleValue("-percentile", 0.9);
  size_t optimized_all = P.getOptionIntValue("-opt", 1);
  size_t reader_threads = P.getOptionIntValue("-readers", 1);

  // Option for starting with a non-empty graph
  const char* const init_graph_file(P.getOptionValue("-init_graph_file"));

  // Option for computing approximation statistics using ground truth file
  const char* const core_ground_truth(P.getOptionValue("-ground_truth_file"));

  // Load the dynamic graph
  using W = typename Graph::weight_type;
  bool use_dynamic = (input_file && input_file[0]);
  BatchDynamicEdges<W> batch_edge_list = use_dynamic ?
    read_batch_dynamic_edge_list<W>(input_file) : BatchDynamicEdges<W>{};
  if (use_dynamic && batch_size == 0) batch_size = batch_edge_list.edges.size();

  // Prepend the initial graph if specified, and store the offset of the
  // initial graph in the resulting dynamic edges list
  size_t offset = 0;
  if (use_dynamic && init_graph_file) {
    BatchDynamicEdges<W> init_graph_list = read_batch_dynamic_edge_list<W>(init_graph_file);
    offset = prepend_dynamic_edge_list(batch_edge_list, init_graph_list);
  }

  // Get the ground truth numbers
  /*if (core_ground_truth) {
    BatchDynamicEdges<W> core_ground = read_batch_dynamic_edge_list<W>(core_ground_truth);
    auto core_ground_truth_vals = core_ground.edges;
    uintE node_cores[core_ground_truth_vals.size()/G.n][G.n];

    for (size_t batch_num = 0; batch_num < core_ground_truth_vals.size()/G.n; batch_num++) {
        for (size_t node = 0; node < G.n; node++) {
            node_cores[batch_num][node] = core_ground_truth_vals[batch_num * G.n + node].to;
        }
    }

    uintE* startRows[core_ground_truth_vals.size()/G.n];
    for (size_t b = 0; b < core_ground_truth_vals.size()/G.n; b++)
        startRows[b] = node_cores[b];

    RunLDS(G, batch_edge_list, batch_size, compare_exact, eps, delta,
          optimized_insertion, offset, get_size, optimized_all, reader_threads, startRows);
  }*/

  // Run LDS
  uintE node_cores[1][1];
  timer t; t.start();
  uintE* startRows[1] = {node_cores[0]};
  RunLDS(G, batch_edge_list, batch_size, compare_exact, eps, delta,
          optimized_insertion, offset, get_size, optimized_all, reader_threads, nonlinearizable,
          percentile);
  double tt = t.stop();

  std::cout << "### Running Time: " << tt << std::endl;

  return tt;
}
}  // namespace gbbs

generate_symmetric_main(gbbs::LDS_runner, false);
