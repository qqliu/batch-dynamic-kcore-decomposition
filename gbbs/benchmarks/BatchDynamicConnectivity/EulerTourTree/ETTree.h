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
        std::tuple<std::tuple<uintE, uintE>, SkipList::SkipListElement*> empty =
            std::make_tuple(std::make_tuple(UINT_E_MAX, UINT_E_MAX), nullptr);
        auto hash_pair = [](const std::tuple<uintE, uintE>& t) {
            size_t l = std::get<0>(t);
            size_t r = std::get<1>(t);
            size_t key = (l << 32) + r;

            return parlay::hash64_2(key);
        };

        edge_table = pbbslib::make_sparse_table<std::tuple<uintE, uintE>,
                   SkipList::SkipListElement*>(n * n, empty, hash_pair);

        vertices = sequence<SkipList::SkipListElement>(n);
        parallel_for(0, n, [&] (size_t i) {
            vertices[i] = skip_list.create_node(i, nullptr, nullptr, std::make_pair(i, i));
            skip_list.join(&vertices[i], &vertices[i]);
        });
    }

    bool is_connected(int u, int v) {
        return skip_list.find_representative(&vertices[u]) == skip_list.find_representative(&vertices[v]);
    }

    void link(uintE u, uintE v) {
        auto uv = skip_list.create_node(u, nullptr, nullptr, std::make_pair(u, v));
        auto vu = skip_list.create_node(v, nullptr, nullptr, std::make_pair(v, u), &uv);
        uv.twin = &vu;

        edge_table.insert(std::make_tuple(u, v), &uv);
        edge_table.insert(std::make_tuple(v, u), &vu);

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

            auto u_left = uv->get_left(0);
            auto v_left = vu->get_left(0);
            auto v_right = skip_list.split(uv);
            auto u_right = skip_list.split(vu);

            skip_list.split(u_left);
            skip_list.split(v_left);
            uv->edge = nullptr;
            vu->edge = nullptr;
            uv->twin = nullptr;
            vu->twin = nullptr;

            skip_list.join(u_left, u_right);
            skip_list.join(v_left, v_right);
    }

    void batch_link_sequential(sequence<std::pair<uintE, uintE>>links) {
            for(size_t i = 0; i < links.size(); i++) {
                    link(links[i].first, links[i].second);
            }
    }

    void batch_link(sequence<std::pair<uintE, uintE>>links) {
        if (links.size() <= 75) {
                batch_link_sequential(links);
                return;
        }

        sequence<std::pair<uintE, uintE>> links_both_dirs = sequence<std::pair<uintE, uintE>>(2 * links.size());
        parallel_for(0, links.size(), [&] (size_t i) {
                links_both_dirs[2 * i] = links[i];
                links_both_dirs[2 * i + 1] = std::make_pair(links[i].second, links[i].first);
        });

        auto get_key = [&] (const std::pair<uintE, uintE>& elm) { return elm.first; };
        parlay::integer_sort_inplace(parlay::make_slice(links_both_dirs), get_key);

        auto split_successors = sequence<SkipList::SkipListElement*>(2 * links.size());

        parallel_for(0, 2 * links.size(), [&] (size_t i) {
            uintE u, v;
            u = links_both_dirs[i].first;
            v = links_both_dirs[i].second;

            if (i == 2 * links.size() - 1 || u != links_both_dirs[i+1].first) {
                    split_successors[i] = skip_list.split(&vertices[u]);
            }

            if (u < v) {
                    auto uv = skip_list.create_node(u, nullptr, nullptr, std::make_pair(u, v));
                    auto vu = skip_list.create_node(v, nullptr, nullptr, std::make_pair(v, u), &uv);
                    uv.twin = &vu;

                    edge_table.insert(std::make_tuple(u, v), &uv);
                    edge_table.insert(std::make_tuple(v, u), &vu);
            }
        });

        parallel_for(0, 2 * links.size(), [&] (size_t i) {
            uintE u, v;
            u = links_both_dirs[i].first;
            v = links_both_dirs[i].second;

            SkipList::SkipListElement* uv = edge_table.find(std::make_tuple(u, v));
            SkipList::SkipListElement* vu = uv -> twin;

            if (i == 0 || u != links_both_dirs[i-1].first) {
                skip_list.join(&vertices[u], uv);
            }

            if (i == 2 * links.size() - 1 || u != links_both_dirs[i+1].first) {
                skip_list.join(vu, split_successors[i]);
            } else {
                uintE u2, v2;
                u2 = links_both_dirs[i+1].first;
                v2 = links_both_dirs[i+1].second;

                skip_list.join(vu, edge_table.find(std::make_tuple(u2, v2)).orig);
            }
        });
    }

    void batch_cut_sequential(sequence<std::pair<uintE, uintE>> cuts) {
            for (size_t i = 0; i < cuts.size(); i++)
                cut(cuts[i].first, cuts[i].second);
    }

    void batch_cut_recurse(sequence<std::pair<uintE, uintE>> cuts) {
            sequence<SkipList::SkipListElement*> join_targets =
                sequence<SkipList::SkipListElement*>(4 * cuts.size(), nullptr);
            sequence<SkipList::SkipListElement*> edge_elements =
                sequence<SkipList::SkipListElement*>(cuts.size(), nullptr);

            parlay::random rng = parlay::random(time(0));

            if (cuts.size() <= 75) {
                    batch_cut_sequential(cuts);
                    return;
            }
            sequence<bool> ignored = sequence<bool>(cuts.size(), true);

            parallel_for(0, cuts.size(), [&] (size_t i) {
                rng.fork(i);
                bool rand_val = rng.rand() % 100 == 0;

                ignored[i] = rand_val;
                if (!ignored[i]) {
                    uintE u, v;
                    u = cuts[i].first;
                    v = cuts[i].second;

                    SkipList::SkipListElement* uv = edge_table.find(std::make_pair(u, v));
                    SkipList::SkipListElement* vu = uv->twin;

                    edge_elements[i] = uv;
                    uv -> split_mark = true;
                    vu -> split_mark = true;
                }

                rng = rng.next();

            });

            parallel_for(0, cuts.size(), [&] (size_t i) {
                    if (!ignored[i]) {
                        SkipList::SkipListElement* uv = edge_elements[i];
                        SkipList::SkipListElement* vu = uv->twin;

                        SkipList::SkipListElement* left_target = uv->get_left(0);
                        if (left_target -> split_mark) {
                            join_targets[4 * i] = nullptr;
                        } else {
                            SkipList::SkipListElement* right_target = vu->get_right(0);
                            while (right_target->split_mark) {
                                right_target = right_target->twin->get_right(0);
                                join_targets[4 * i] = left_target;
                                join_targets[4 * i + 1] = right_target;
                            }
                        }
                        left_target = vu->get_left(0);
                        if(left_target -> split_mark) {
                            join_targets[4 * i + 2] = nullptr;
                        } else {
                            SkipList::SkipListElement* right_target = uv -> get_right(0);
                            while (right_target -> split_mark) {
                                    right_target = right_target->twin->get_right(0);
                            }
                            join_targets[4 * i + 2] = left_target;
                            join_targets[4 * i + 3] = right_target;
                        }
                    }
            });

            parallel_for(0, cuts.size(), [&] (size_t i) {
                    if (!ignored[i]) {
                        SkipList::SkipListElement* uv = edge_elements[i];
                        SkipList::SkipListElement* vu = uv->twin;

                        skip_list.split(uv);
                        skip_list.split(vu);

                        SkipList::SkipListElement* predecessor = uv->get_left(0);
                        if (predecessor != nullptr) {
                            skip_list.split(predecessor);
                        }

                        predecessor = vu->get_left(0);
                        if (predecessor != nullptr) {
                            skip_list.split(predecessor);
                        }
                    }
           });

           parallel_for(0, cuts.size(), [&] (size_t i) {
                    if (!ignored[i]) {
                        uintE u, v;
                        u = cuts[i].first;
                        v = cuts[i].second;

                        edge_table.remove(std::make_pair(u, v));

                        if (join_targets[4 * i] != nullptr) {
                            skip_list.join(join_targets[4 * i], join_targets[4 * i + 1]);
                        }

                        if (join_targets[4 * i + 2] != nullptr) {
                            skip_list.join(join_targets[4 * i + 2], join_targets[4 * i + 3]);
                        }
                    }
            });

            auto element_indices = parlay::pack_index(ignored);
            auto next_cuts_seq = sequence<std::pair<uintE, uintE>>(element_indices.size());
            parallel_for(0, next_cuts_seq.size(), [&] (size_t i){
                next_cuts_seq = cuts[element_indices[i]];
            });
            batch_cut_recurse(next_cuts_seq);
    }

    void batch_cut(sequence<std::pair<uintE, uintE>> cuts) {
            if (cuts.size() <=  75) {
                    batch_cut_sequential(cuts);
                    return;
            }

            batch_cut_recurse(cuts);
    }
};

void RunETTree() {
        std::cout << "ET tree" << std::endl;
}

}  // namespace gbbs
