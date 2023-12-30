#include <unordered_set>
#include <stack>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
#include "gbbs/pbbslib/sparse_table.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"

#include "sparse_set.h"

namespace gbbs {

struct SkipList {

    struct SkipListElement {
        size_t height;
        size_t lowest_needs_update = 0;
        uintE val;

        using pointers = std::pair<SkipListElement*, SkipListElement*>;
        using height_array = parlay::sequence<pointers>;

        height_array elements;

        SkipListElement(): height(0), lowest_needs_update(0) {}

        SkipListElement(size_t _h, SkipListElement* _r, SkipListElement* _l, uintE _val):
            height(_h), lowest_needs_update(_h) {
                elements.resize(_h);
                elements[0].first = _l;
                elements[0].second = _r;
                val = _val;
        }

        inline void set_left_pointer(size_t height, SkipListElement* left) {
            elements[height].first = left;
        }

        inline void set_right_pointer(size_t height, SkipListElement* right) {
            elements[height].second = right;
        }

        inline SkipListElement* get_left(size_t height) {
                return elements[height].first;
        }

        inline SkipListElement* get_right(size_t height) {
                return elements[height].second;
        }
    };

    size_t n;
    parlay::sequence<pbbslib::sparse_table> vertex_edge_elements;

    SkipList(): n(0) {}

    SkipList(size_t _n): n(_n) {
        vertex_edge_elements.resize(n);
        using KV = std::tuple<uintE, SkipListElement*>;

        auto hash_pair = [](const std::tuple<uintE, uintE>& t) {
            size_t l = std::min(std::get<0>(t), std::get<1>(t));
            size_t r = std::max(std::get<0>(t), std::get<1>(t));
            size_t key = (l << 32) + r;
            return parlay::hash64_2(key);
        };

        parallel_for(0, n, [&] (size_t i){
            KV empty = std::make_tuple(UINT_E_MAX, gbbs::empty());
            auto edge_table = pbbslib::make_sparse_table<uintE, SkipListElement*>(n, empty, hash_pair);
            vertex_edge_elements[i] = edge_table;
        });
    }

    /* The augmented skip list can take any arbitrary associative, commutative function */
};
template <class Graph, class W>
    inline void RunSkipList(Graph& G, BatchDynamicEdges<W> batch_edge_list) {
        size_t num_vertices = G.n;
        auto skip_list = SkipList(n);
}

}  // namespace gbbs
