#include <unordered_set>
#include <stack>
#include <limits>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
#include "gbbs/pbbslib/sparse_table.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "benchmarks/BatchDynamicConnectivity/EulerTourTree/ETTree.h"

namespace gbbs {

// A union-find node that stores its parent and rank
struct UFNode {
  uintE parent;
  uintE rank;
  UFNode() {}

  UFNode(int p, int r) : parent(p), rank(r) {}
};

// A union-find data structure that supports parallel find and union operations
class UnionFind {
  public:
  // An array of union-find nodes
  sequence<UFNode> nodes;

  UnionFind() {}
  // Constructor that initializes each node as a singleton set
  UnionFind(size_t n) {
      nodes = sequence<UFNode>(n);
      parallel_for(0, n, [&] (size_t i) {
        nodes[i] = UFNode(i, 0);
      });
  }

  // Find the root of the set containing x
  uintE find(uintE x) {
    // Use path compression to update the parent pointers along the path
    if (nodes[x].parent != x) {
      nodes[x].parent = find(nodes[x].parent);
    }
    return nodes[x].parent;
  }

  // Union the sets containing x and y
  void unionf(uintE x, uintE y) {
    // Find the roots of x and y
    int xRoot = find(x);
    int yRoot = find(y);
    // If they are already in the same set, do nothing
    if (xRoot == yRoot) return;
    // Use union by rank to merge the smaller set into the larger one
    if (nodes[xRoot].rank < nodes[yRoot].rank) {
      nodes[xRoot].parent = yRoot;
    } else if (nodes[xRoot].rank > nodes[yRoot].rank) {
      nodes[yRoot].parent = xRoot;
    } else {
      nodes[yRoot].parent = xRoot;
      nodes[xRoot].rank++;
    }
  }

  // Parallel find operation that returns an array of roots for each element
  parlay::sequence<std::pair<uintE, uintE>> parallel_find(sequence<std::pair<uintE, uintE>> nodes) {
    // Use parallel reduce to find the roots of each element in parallel
    sequence<std::pair<uintE, uintE>> parents = sequence<std::pair<uintE, uintE>>(nodes.size());
    parallel_for(0, nodes.size(),  [&](size_t i) {
      parents[i] = std::make_pair(find(nodes[i].first), find(nodes[i].second));
    });

    return parents;
  }
};

// A type alias for an edge represented by a pair of uintE
using edge_mst_type = std::pair<uintE, uintE>;

// A function that returns a sequence of edges in the spanning forest
sequence<edge_mst_type> parallel_spanning_forest(sequence<edge_mst_type> edges, size_t n) {
  // Initialize a union-find data structure
  UnionFind uf(n);
  // Initialize a sequence of edges in the spanning forest
  sequence<edge_mst_type> forest;
  // Repeat until all the edges are processed
  while (!edges.empty()) {
    // Find the roots of each endpoint in parallel
    auto not_selected = sequence<bool>(edges.size(), true);

    auto endpoints = sequence<std::pair<uintE, uintE>>(edges.size());
    parallel_for(0, edges.size(), [&](size_t i){
        endpoints[i] = std::make_pair(edges[i].first, edges[i].second);
    });

    auto roots = uf.parallel_find(endpoints);

    using K = std::pair<uintE, uintE>;
    using V = std::pair<std::pair<uintE, uintE>, uintE>;
    using KV = std::pair<K, V>;

    KV empty =
        std::make_pair(std::make_pair(UINT_E_MAX, UINT_E_MAX), std::make_pair(
                std::make_pair(UINT_E_MAX, UINT_E_MAX), UINT_E_MAX));

    auto hash_pair = [](const std::pair<uintE, uintE>& t) {
        size_t l = std::min(std::get<0>(t), std::get<1>(t));
        size_t r = std::max(std::get<0>(t), std::get<1>(t));
        size_t key = (l << 32) + r;
        return parlay::hash64_2(key);
    };

    auto root_edge_to_original =
         pbbslib::make_sparse_table<K, V>(roots.size(), empty, hash_pair);

    parallel_for(0, roots.size(), [&](size_t i){
        auto edge = roots[i];
        auto u = std::get<0>(edge);
        auto v = std::get<1>(edge);
        //std::cout << "root edges: " << u << ", " << v << std::endl;

        bool abort = false;
        root_edge_to_original.insert_check(std::make_pair(std::make_pair(u,
            v), std::make_pair(endpoints[i], i)), &abort);
    });

    // Filter out the edges that connect vertices in the same component

    auto bool_filter = sequence<bool>(roots.size());
    parallel_for(0, bool_filter.size(), [&](size_t i){
        bool_filter[i] = roots[i].first != roots[i].second;
        not_selected[i] = bool_filter[i];
    });
    auto filtered = parlay::pack(roots, bool_filter);

    auto compare_tup = [&] (const std::pair<uintE, uintE> l, const std::pair<uintE, uintE> r) {
        return l.first < r.first;
    };
    parlay::sort_inplace(parlay::make_slice(filtered), compare_tup);

    auto bool_seq = sequence<bool>(filtered.size() + 1);
    parallel_for(0, filtered.size() + 1, [&](size_t i) {
        bool_seq[i] = (i == 0) || (i == filtered.size()) ||
            (filtered[i-1].first != filtered[i].first);
    });
    auto starts = parlay::pack_index(bool_seq);

    auto selected = sequence<edge_mst_type>(starts.size()-1);

    parallel_for(0, starts.size()-1, [&](size_t i){
        auto original_edge = root_edge_to_original.find(filtered[starts[i]],
                std::make_pair(std::make_pair(UINT_E_MAX, UINT_E_MAX), UINT_E_MAX));
        selected[i] = original_edge.first;
        not_selected[original_edge.second] = false;
    });

    // Append the selected edges to the forest
    forest.append(selected);

    // Union the components connected by the selected edges
    parallel_for(0, selected.size(), [&](size_t i) {
      uf.unionf(selected[i].first, selected[i].second);
    });

    // Remove the selected edges from the original sequence
    edges = parlay::pack(edges, not_selected);
  }
  // Return the forest
  return forest;
}

struct Connectivity {
    size_t n;
    ETTree tree;

    Connectivity() {};

    Connectivity(size_t n_, int copies_, size_t m_, double pb_) {
        n = n_;
        tree = ETTree(n_, copies_, m_, pb_);
    }

    template <class Seq>
    void batch_insertion(const Seq& insertions) {
        auto non_empty_spanning_tree = true;
        auto first = true;
        sequence<std::pair<uintE, uintE>> edges_both_directions;
        sequence<size_t> starts;
        sequence<SkipList::SkipListElement*> representative_nodes;

        while(non_empty_spanning_tree) {

            //std::cout << "number of edges in insertions: " << insertions.size() << std::endl;
            if (first) {
                edges_both_directions =  sequence<std::pair<uintE, uintE>>(2 * insertions.size());

                parallel_for(0, insertions.size(), [&](size_t i) {
                    //std::cout << "original edge: " << insertions[i].first << ", " << insertions[i].second << std::endl;
                    auto [u, v] = insertions[i];
                    tree.create_edge_in_edge_table(u, v);
                    edges_both_directions[2 * i] = std::make_pair(u, v);
                    edges_both_directions[2 * i + 1] = std::make_pair(v, u);
                });

                auto compare_tup = [&] (const std::pair<uintE, uintE> l, const std::pair<uintE, uintE> r) {
                    return l.first < r.first;
                };
                parlay::sort_inplace(parlay::make_slice(edges_both_directions), compare_tup);

                auto bool_seq = sequence<bool>(edges_both_directions.size() + 1);
                parallel_for(0, edges_both_directions.size() + 1, [&](size_t i) {
                    bool_seq[i] = (i == 0) || (i == edges_both_directions.size()) ||
                       (edges_both_directions[i-1].first != edges_both_directions[i].first);
                });
                starts = parlay::pack_index(bool_seq);

                representative_nodes =
                    sequence<SkipList::SkipListElement*>(starts.size()-1);
                sequence<SkipList::SkipListElement*> nodes_to_update =
                    sequence<SkipList::SkipListElement*>(starts.size()-1);
                auto update_seq =
                    sequence<std::pair<SkipList::SkipListElement*,
                    sequence<sequence<std::pair<uintE, uintE>>>>>(starts.size()- 1);

            /*std::cout << "starts size: " << starts.size() << std::endl;
            std::cout << "edges both directions size: " << edges_both_directions.size() << std::endl;
            for (int i = 0; i < starts.size(); i++) {
                    std::cout << "starts: " << starts[i] << std::endl;
            }

            for (int i = 0; i < edges_both_directions.size(); i++) {
                    std::cout << "edge: " << edges_both_directions[i].first << ", " << edges_both_directions[i].second
                        << std::endl;
            }*/

                parallel_for(0, starts.size() - 1, [&](size_t i) {
                    for (size_t j = starts[i]; j < starts[i+1]; j++) {
                    /*std::cout << "starts for " << i <<
                        ": " << edges_both_directions[j].first << ", "
                        << edges_both_directions[j].second << std::endl;*/

                        tree.add_edge_to_cutsets(edges_both_directions[j]);
                    }
                    SkipList::SkipListElement* our_vertex = &tree.vertices[edges_both_directions[starts[i]].first];
                    representative_nodes[i] =
                        tree.skip_list.find_representative(our_vertex);
                    update_seq[i] = std::make_pair(our_vertex, our_vertex->values[0]);
                });
                tree.skip_list.batch_update(&update_seq);
            }
            first = false;
            if (!first) {
                parallel_for(0, starts.size() - 1, [&](size_t i) {
                    SkipList::SkipListElement* our_vertex = &tree.vertices[edges_both_directions[starts[i]].first];
                    representative_nodes[i] =
                        tree.skip_list.find_representative(our_vertex);
                });
            }

            //std::cout << "doing update" << std::endl;

            sequence<std::pair<uintE, uintE>> found_possible_edges =
                    sequence<std::pair<uintE, uintE>>(representative_nodes.size());
            sequence<bool> is_edge = sequence<bool>(representative_nodes.size());

            parallel_for(0, representative_nodes.size(), [&](size_t i) {
                auto id = representative_nodes[i]->id.first;
                auto sum = tree.get_tree_sum(id);
                bool is_present = false;
                bool is_real_edge = false;

                for (size_t ii = 0; ii < sum.size(); ii++) {
                if (!is_real_edge) {
                    for (size_t ij = 0; ij < sum[ii].size(); ij++) {
                        if (!is_real_edge) {
                            //std::cout << "rep values for " << id << ": " << sum[ii][ij].first << ", " << sum[ii][ij].second
                                //<< std::endl;
                            is_present = sum[ii][ij].first != UINT_E_MAX && sum[ii][ij].second != UINT_E_MAX;
                            auto u = sum[ii][ij].first;
                            auto v = sum[ii][ij].second;
                            if (is_present) {
                                if (u < n && v < n) {
                                    auto element = &tree.edge_table[u][v];
                                    /*std::cout << "element: " << element->id.first << ", " << element->id.second <<
                                        std::endl;*/

                                    if (element->id.first != UINT_E_MAX && element->id.second != UINT_E_MAX) {
                                        is_real_edge = true;
                                        found_possible_edges[i] = std::make_pair(u, v);
                                    }
                                }
                            }
                        }
                    }
                }
                }

                auto possible_edge =
                    tree.edge_table[found_possible_edges[i].first][found_possible_edges[i].second];
                is_edge[i] = is_real_edge;
            });

            auto real_edges = parlay::pack(found_possible_edges, is_edge);
            auto representative_edges = sequence<std::pair<uintE, uintE>>(real_edges.size());

            size_t max_index = 0;
            parallel_for(0, real_edges.size(), [&](size_t i) {
                auto u = tree.vertices[real_edges[i].first];
                auto v = tree.vertices[real_edges[i].second];
                auto node_a = tree.skip_list.find_representative(&u);
                auto node_b = tree.skip_list.find_representative(&v);
                auto node_a_id = std::min(node_a->id.first, node_a->id.second);
                auto node_b_id = std::min(node_b->id.first, node_b->id.second);
                if (node_a_id > max_index)
                    max_index = node_a_id;
                if (node_b_id > max_index)
                    max_index = node_b_id;

                representative_edges[i] = std::pair(node_a_id, node_b_id);
            });

            using K = std::pair<uintE, uintE>;
            using V = std::pair<uintE, uintE>;
            using KV = std::pair<K, V>;

            KV empty =
                std::make_pair(std::make_pair(UINT_E_MAX, UINT_E_MAX), std::make_pair(UINT_E_MAX, UINT_E_MAX));
            auto hash_pair = [](const std::pair<uintE, uintE>& t) {
                    size_t l = std::min(std::get<0>(t), std::get<1>(t));
                    size_t r = std::max(std::get<0>(t), std::get<1>(t));
                    size_t key = (l << 32) + r;
                return parlay::hash64_2(key);
            };
            auto representative_edge_to_original =
                pbbslib::make_sparse_table<K, V>(representative_edges.size(), empty, hash_pair);

            parallel_for(0, representative_edges.size(), [&](size_t i){
                auto edge = representative_edges[i];
                auto u = std::get<0>(edge);
                auto v = std::get<1>(edge);
                //std::cout << "rep edges: " << u << ", " << v << std::endl;

                bool abort = false;
                representative_edge_to_original.insert_check(std::make_pair(std::make_pair(u,
                        v), real_edges[i]), &abort);
            });

            auto spanning_forest = parallel_spanning_forest(representative_edges, n);
            if(spanning_forest.size() == 0)
                non_empty_spanning_tree = false;

            auto original_edges = sequence<std::pair<uintE, uintE>>(spanning_forest.size());
            parallel_for(0, spanning_forest.size(), [&](size_t i) {
                auto u = std::get<0>(spanning_forest[i]);
                auto v = std::get<1>(spanning_forest[i]);
                //std::cout << "spanning tree edges: " << u << ", " << v << std::endl;

                auto original_edge = representative_edge_to_original.find(std::make_pair(u, v), std::make_pair(UINT_E_MAX, UINT_E_MAX));
                //std::cout << "original edges: "<< original_edge.first << ", " << original_edge.second << std::endl;
                if (original_edge.first == UINT_E_MAX || original_edge.second == UINT_E_MAX)
                    std::cout << "FAILURE: THERE IS A BUG" << std::endl;
                original_edges[i] = std::make_pair(original_edge.first, original_edge.second);
            });
            tree.batch_link(original_edges);
         }
    }

    template <class Seq>
    void batch_deletion(const Seq& deletions) {
    }
};

void RunConnectivityTest() {
        std::cout << "Connectivity Test" << std::endl;
}

template <class W>
inline void RunConnectivity(BatchDynamicEdges<W>& batch_edge_list, long batch_size, bool compare_exact,
        size_t offset, size_t n, int copies, size_t m, double pb) {
        auto cutset = Connectivity(n, copies, m, pb);
        auto batch = batch_edge_list.edges;
        std::cout << "batch size: " << batch.size() << std::endl;

        if (offset != 0) {
            for (size_t i = 0; i < offset; i += 1000000) {
                auto end_size = std::min(i + 1000000, offset);
                auto insertions = parlay::filter(parlay::make_slice(batch.begin() + i,
                        batch.begin() + end_size), [&] (const DynamicEdge<W>& edge){
                    return edge.insert;
                });
                auto deletions = parlay::filter(parlay::make_slice(batch.begin() + i,
                        batch.begin() + end_size), [&] (const DynamicEdge<W>& edge){
                    return !edge.insert;
                });
                auto batch_insertions = parlay::delayed_seq<std::pair<uintE, uintE>>(insertions.size(),
                    [&] (size_t i) {
                    uintE vert1 = insertions[i].from;
                    uintE vert2 = insertions[i].to;
                    return std::make_pair(vert1, vert2);
                });
                auto batch_deletions = parlay::delayed_seq<std::pair<uintE, uintE>>(deletions.size(),
                [&] (size_t i) {
                    uintE vert1 = deletions[i].from;
                    uintE vert2 = deletions[i].to;
                    return std::make_pair(vert1, vert2);
                });
                cutset.batch_insertion(batch_insertions);
                cutset.batch_deletion(batch_deletions);
            }
        }

        for (size_t i = offset; i < batch.size(); i += batch_size) {
            std::cout << "batch: " << i << std::endl;
            timer t; t.start();
            auto end_size = std::min(i + batch_size, batch.size());
            auto insertions = parlay::filter(parlay::make_slice(batch.begin() + i,
                    batch.begin() + end_size), [&] (const DynamicEdge<W>& edge){
                return edge.insert;
            });

            auto deletions = parlay::filter(parlay::make_slice(batch.begin() + i,
                    batch.begin() + end_size), [&] (const DynamicEdge<W>& edge){
                return !edge.insert;
            });

            auto batch_insertions = parlay::delayed_seq<std::pair<uintE, uintE>>(insertions.size(),
                [&] (size_t i) {
                uintE vert1 = insertions[i].from;
                uintE vert2 = insertions[i].to;
                return std::make_pair(vert1, vert2);
            });


            auto batch_deletions = parlay::delayed_seq<std::pair<uintE, uintE>>(deletions.size(),
                [&] (size_t i) {
                uintE vert1 = deletions[i].from;
                uintE vert2 = deletions[i].to;
                return std::make_pair(vert1, vert2);
            });

            cutset.batch_insertion(batch_insertions);

        /*num_insertion_flips += layers.batch_insertion(batch_insertions);
        double insertion_time = t.stop();

        t.start();
        num_deletion_flips += layers.batch_deletion(batch_deletions);

        max_degree = layers.max_degree();

        double deletion_time = t.stop();
        double tt = insertion_time + deletion_time;

        std::cout << "### Batch Running Time: " << tt << std::endl;
        std::cout << "### Insertion Running Time: " << insertion_time << std::endl;
        std::cout << "### Deletion Running Time: " << deletion_time << std::endl;
        std::cout << "### Batch Num: " << end_size - offset << std::endl;
        std::cout << "### Coreness Estimate: " << layers.max_coreness() << std::endl;
        std::cout << "### Number Insertion Flips: " << num_insertion_flips << std::endl;
        std::cout << "### Number Deletion Flips: " << num_deletion_flips << std::endl;
        std::cout << "### Max Outdegree: " << max_degree << std::endl;
        if (get_size) {
            auto size = layers.get_size();
            std::cout << "### Size: " << size << std::endl;
        }
        if (compare_exact) {
            auto graph = dynamic_edge_list_to_symmetric_graph(batch_edge_list, std::min(batch.size(),
                        i + batch_size));

            // Run kcore on graph
            auto cores = KCore(graph, 16);

            auto max_core = parlay::reduce(cores, parlay::maxm<uintE>());
            std::cout << "### Coreness Exact: " << max_core << std::endl;

            // Compare cores[v] to layers.core(v)
            auto approximation_error = parlay::delayed_seq<float>(batch_edge_list.max_vertex,
                    [&] (size_t j) -> float {
                auto exact_core = j >= graph.n ? 0 : cores[j];
                auto approx_core = layers.core(j);
                if (exact_core == 0 || approx_core == 0) {
                    return 0;
                }
                return (exact_core > approx_core) ? (float) exact_core / (float) approx_core :
                       (float) approx_core / (float) exact_core;
            });

            double mult_appx = (2 + 2*layers.eps);
            float bad = parlay::reduce(parlay::delayed_seq<float>(batch_edge_list.max_vertex, [&](size_t j) -> float{
                auto true_core = j >= graph.n ? 0 : cores[j];
                auto appx_core = layers.core(j);
                return (appx_core > (mult_appx * true_core)) + (appx_core < (true_core/mult_appx));
            }), parlay::addm<float>());

            // Output min, max, and average error
            float sum_error = parlay::reduce(approximation_error, parlay::addm<float>());
            float max_error = parlay::reduce(approximation_error, parlay::maxm<float>());
            float min_error = parlay::reduce(approximation_error,
              parlay::make_monoid([](float l, float r){
                if (l == 0) return r;
                if (r == 0) return l;
                return std::min(r, l);
            }, (float) 0));
            float denominator = parlay::reduce(parlay::delayed_seq<float>(batch_edge_list.max_vertex,
                        [&] (size_t j) -> float{
                auto exact_core = j >= graph.n ? 0 : cores[j];
                auto approx_core = layers.core(j);
                return (exact_core != 0) && (approx_core != 0);
            }), parlay::addm<float>());
            auto avg_error = (denominator == 0) ? 0 : sum_error / denominator;
            std::cout << "### Num Bad: " << bad << std::endl;
            std::cout << "### Per Vertex Average Coreness Error: " << avg_error << std::endl; fflush(stdout);
            std::cout << "### Per Vertex Min Coreness Error: " << min_error << std::endl; fflush(stdout);
            std::cout << "### Per Vertex Max Coreness Error: " << max_error << std::endl; fflush(stdout);
        }*/
    }
}

}  // namespace gbbs
