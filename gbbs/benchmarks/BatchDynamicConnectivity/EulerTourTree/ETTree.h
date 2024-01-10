#include <unordered_set>
#include <stack>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
#include "gbbs/pbbslib/sparse_table.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "benchmarks/EdgeOrientation/ParallelLDS/sparse_set.h"
#include "benchmarks/BatchDynamicConnectivity/SkipList/SkipList.h"

namespace gbbs {

template <class H>
struct ETTree {
   pbbslib::sparse_table<std::tuple<uintE, uintE>, SkipList::SkipListElement*, H> edge_table;
   SkipList skip_list;
   sequence<SkipList::SkipListElement> vertices;

   ETTree() {}

   ETTree(size_t n) {
        skip_list = SkipList(n);
        std::tuple<std::tuple<uintE, uintE>, ETTreeBase*> empty =
            std::make_tuple(std::make_tuple(UINT_E_MAX, UINT_E_MAX), nullptr);
        auto hash_pair = [](const std::tuple<uintE, uintE>& t) {
            size_t l = std::get<0>(t);
            size_t r = std::get<1>(t);
            size_t key = (l << 32) + r;

            return parlay::hash64_2(key);
        };

        edge_table = pbbslib::make_sparse_table<std::tuple<uintE, uintE>,
                   ETTreeBase*>(n * n, empty, hash_pair);

        vertices = sequence<SkipList::SkipListElement>(n);
        parallel_for(0, n, [&] (size_t i) {
            vertices[i] = skip_list.create_node(i, nullptr, nullptr, std::make_pair(i, i));
            skip_list.join(&vertices[i], &vertices[i]);
        });
    }

    struct ETTreeBase {
        SkipList::SkipListElement* edge; // edge (u, v)
        SkipList::SkipListElement* twin; // edge (v, u)

        bool skip_mark;

        ETTreeBase() {
            edge = SkipList::SkipListElement();
            twin = nullptr;

            skip_mark = false;
        }

        ETTreeBase(SkipList::SkipListElement* orig, SkipList::SkipListElement* twin_) {
                edge = orig;
                skip_mark = false;
                twin = twin_;
        }
    };

    bool is_connected(int u, int v) {
        return skip_list.find_representative(&vertices[u]) == skip_list.find_representative(&vertices[v]);
    }

    void link(uintE u, uintE v) {
        auto uv = skip_list.create_node(u, nullptr, nullptr, std::make_pair(u, v));
        auto vu = skip_list.create_node(v, nullptr, nullptr, std::make_pair(v, u));

        auto uv_et = ETTreeBase(&uv, &vu);
        auto vu_et = ETTreeBase(&vu, &uv);

        edge_table.insert(std::make_tuple(u, v), &uv_et);
        edge_table.insert(std::make_tuple(v, u), &vu_et);

        auto u_left = &vertices[u];
        auto v_left = &vertices[v];

        auto u_right = skip_list.split(u_left);
        auto v_right = skip_list.split(v_left);

        skip_list.join(u_left, &uv);
        skip_list.join(&uv, v_right);
        skip_list.join(v_left, &vu);
        skip_list.join(&vu, u_right);
    }

    void cut(int u, int v){
            auto uv = edge_table.find(std::make_pair(u, v));
            auto vu = uv->twin;
            edge_table.remove(std::make_pair(u, v));
            edge_table.remove(std::make_pair(v, u));

    }
};

void RunETTree() {
        std::cout << "ET tree" << std::endl;
}

}  // namespace gbbs
