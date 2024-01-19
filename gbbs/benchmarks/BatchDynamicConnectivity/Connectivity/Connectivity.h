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
            sequence<SkipList::SkipListElement*> representative_nodes =
                sequence<SkipList::SkipListElement*>(2 * insertions.size());

            parallel_for(0, insertions.size(), [&](size_t i) {
                tree.add_edge_to_cutsets(insertions[i]);
                representative_nodes[2 * i] =
                    tree.skip_list.find_representative(tree.vertices[insertions[i].first]);
                representative_nodes[2 * i + 1] =
                    tree.skip_list.find_representative(tree.vertices[insertions[i].second]);
            });
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
