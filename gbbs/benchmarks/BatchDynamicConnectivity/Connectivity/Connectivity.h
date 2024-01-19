#include <unordered_set>
#include <stack>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
#include "gbbs/pbbslib/sparse_table.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "benchmarks/EdgeOrientation/ParallelLDS/sparse_set.h"
#include "benchmarks/BatchDynamicConnectivity/EulerTourTree/ETTree.h"

namespace gbbs {

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
            sequence<std::pair<uintE, uintE>> edges_both_directions =
                sequence<std::pair<uintE, uintE>>(2 * insertions.size());

            parallel_for(0, insertions.size(), [&](size_t i) {
                auto [u, v] = insertions[i];
                edges_both_directions[2 * i] = std::make_pair(u, v);
                edges_both_directions[2 * i + 1] = std::make_pair(v, u);
            });

            auto compare_tup = [&] (const std::pair<uintE, uintE> l, const std::pair<uintE, uintE> r) {
                    return l.first < r.first;
            };
            parlay::sort_inplace(parlay::make_slice(edges_both_directions), compare_tup);

            auto bool_seq = sequence<bool>(edges_both_directions.size());
            parallel_for(0, edges_both_directions.size() + 1, [&](size_t i) {
                bool_seq[i] = (i == 0) || (i == edges_both_directions.size()) ||
                       (insertions[i-1].first != insertions[i].first);
            });
            auto starts = parlay::pack_index(bool_seq);

            sequence<SkipList::SkipListElement*> representative_nodes =
                sequence<SkipList::SkipListElement*>(starts.size()-1);
            sequence<SkipList::SkipListElement*> nodes_to_update =
                sequence<SkipList::SkipListElement*>(starts.size()-1);
            auto update_seq =
                sequence<std::pair<SkipList::SkipListElement*,
                sequence<sequence<std::pair<uintE, uintE>>>>>(starts.size()- 1);

            parallel_for(0, starts.size() - 1, [&](size_t i) {
                for (size_t j = starts[i]; j < starts[i+1]; j++) {
                    tree.add_edge_to_cutsets(insertions[j]);
                }
                SkipList::SkipListElement* our_vertex = tree.vertices[edges_both_directions[starts[i]].first];
                representative_nodes[i] =
                    tree.skip_list.find_representative(our_vertex);
                update_seq[i] = std::make_pair(our_vertex, our_vertex->values[0]);
            });

            tree.skip_list.batch_update(&update_seq);
    }

    template <class Seq>
    void batch_deletion(const Seq& insertions) {
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
}

}  // namespace gbbs
