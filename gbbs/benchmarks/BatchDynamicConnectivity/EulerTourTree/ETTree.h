#include <unordered_set>
#include <stack>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
#include "gbbs/pbbslib/sparse_table.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "benchmarks/EdgeOrientation/ParallelLDS/sparse_set.h"
#include "benchmarks/BatchDynamicConnectivity/SkipList/SkipList.h"

namespace gbbs {

  uintE hash_function_parlay(std::tuple<uintE, uintE>& t) {
            size_t l = std::get<0>(t);
            size_t r = std::get<1>(t);
            size_t key = (l << 32) + r;

            return parlay::hash64_2(key);
  }

struct ETTree {
   std::tuple<std::tuple<uintE, uintE>, SkipList::SkipListElement*> empty =
          std::make_tuple(std::make_tuple(UINT_E_MAX, UINT_E_MAX), nullptr);
   using table_type = decltype(pbbslib::make_sparse_table<std::tuple<uintE, uintE>,
          SkipList::SkipListElement*>(10, empty, hash_function_parlay));
   table_type edge_table;
   //pbbslib::sparse_table<std::tuple<uintE, uintE>, SkipList::SkipListElement*, H> edge_table;
   SkipList skip_list;
   sequence<SkipList::SkipListElement> vertices;

   ETTree() {}

   ETTree(size_t n) {
        skip_list = SkipList(n);
        /*auto hash_pair = [](const std::tuple<uintE, uintE>& t) {
            size_t l = std::get<0>(t);
            size_t r = std::get<1>(t);
            size_t key = (l << 32) + r;

            return parlay::hash64_2(key);
        };*/

        edge_table = pbbslib::make_sparse_table<std::tuple<uintE, uintE>,
                   SkipList::SkipListElement*>(2 * n * n, empty, hash_function_parlay);

        vertices = sequence<SkipList::SkipListElement>(n);
        auto joins = sequence<std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>>(n);
        parallel_for(0, n, [&] (size_t i) {
            vertices[i] = skip_list.create_node(i, nullptr, nullptr, std::make_pair(i, i));
            joins[i] = std::make_pair(&vertices[i], &vertices[i]);
        });

        skip_list.batch_join(&joins);
    }

    bool is_connected(int u, int v) {
        return skip_list.find_representative(&vertices[u]) == skip_list.find_representative(&vertices[v]);
    }

    void link(uintE u, uintE v) {
        auto uv = skip_list.create_node(u, nullptr, nullptr, std::make_pair(u, v));
        auto vu = skip_list.create_node(v, nullptr, nullptr, std::make_pair(v, u), &uv);
        uv.twin = &vu;

        edge_table.insert(std::make_tuple(std::make_tuple(u, v), &uv));
        edge_table.insert(std::make_tuple(std::make_tuple(v, u), &vu));

        auto u_left = &vertices[u];
        auto v_left = &vertices[v];

        auto splits = sequence<SkipList::SkipListElement*>(2);
        splits[0] = u_left;
        splits[1] = v_left;
        auto results = skip_list.batch_split(&splits);

        auto u_right = results[0];
        auto v_right = results[1];

        auto joins = sequence<std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>>(4);

        joins[0] = std::make_pair(u_left, &uv);
        joins[1] = std::make_pair(&uv, v_right);
        joins[2] = std::make_pair(v_left, &vu);
        joins[3] = std::make_pair(&vu, u_right);

        skip_list.batch_join(&joins);
    }

    void cut(int u, int v){
            auto uv = edge_table.find(std::make_tuple(u, v), nullptr);
            auto vu = uv->twin;
            //edge_table.remove(std::make_pair(u, v));
            //edge_table.remove(std::make_pair(v, u));

            auto u_left = uv->get_left(0);
            auto v_left = vu->get_left(0);

            auto splits = sequence<SkipList::SkipListElement*>(2);
            splits[0] = uv;
            splits[1] = vu;
            auto results = skip_list.batch_split(&splits);

            auto v_right = splits[0];
            auto u_right = splits[1];

            splits = sequence<SkipList::SkipListElement*>(2);
            splits[0] = u_left;
            splits[1] = v_left;
            results = skip_list.batch_split(&splits);

            /*uv->orig = nullptr;
            vu->orig = nullptr;*/
            uv->twin = nullptr;
            vu->twin = nullptr;

            auto joins = sequence<std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>>(2);

            joins[0] = std::make_pair(u_left, u_right);
            joins[1] = std::make_pair(v_left, v_right);

            skip_list.batch_join(&joins);
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
        auto splits = sequence<SkipList::SkipListElement*>(2 * links.size(), nullptr);

        parallel_for(0, 2 * links.size(), [&] (size_t i) {
            uintE u, v;
            u = links_both_dirs[i].first;
            v = links_both_dirs[i].second;

            if (i == 2 * links.size() - 1 || u != links_both_dirs[i+1].first) {
                    splits[i] = &vertices[u];
            }

            if (u < v) {
                    auto uv = skip_list.create_node(u, nullptr, nullptr, std::make_pair(u, v));
                    auto vu = skip_list.create_node(v, nullptr, nullptr, std::make_pair(v, u), &uv);
                    uv.twin = &vu;

                    edge_table.insert(std::make_tuple(std::make_tuple(u, v), &uv));
                    edge_table.insert(std::make_tuple(std::make_tuple(v, u), &vu));
            }
        });

        auto bool_seq = parlay::delayed_seq<bool>(splits.size(), [&] (size_t i) {
                return (splits[i] != nullptr);
        });

        auto element_indices = parlay::pack_index(bool_seq);
        auto filtered_splits = sequence<SkipList::SkipListElement*>(element_indices.size());
        parallel_for(0, filtered_splits.size(), [&] (size_t i) {
            filtered_splits[i] = splits[element_indices[i]];
        });

        auto results = skip_list.batch_split(&filtered_splits);

        parallel_for(0, results.size(), [&] (size_t i) {
            auto split_index = element_indices[i];
            split_successors[split_index] = results[i];
        });

        auto joins = sequence<std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>>(2 * links.size(),
                std::make_pair(nullptr, nullptr));
        parallel_for(0, 2 * links.size(), [&] (size_t i) {
            uintE u, v;
            u = links_both_dirs[i].first;
            v = links_both_dirs[i].second;

            SkipList::SkipListElement* uv = edge_table.find(std::make_tuple(u, v), nullptr);
            SkipList::SkipListElement* vu = uv -> twin;

            if (i == 0 || u != links_both_dirs[i-1].first) {
                joins[2*i] = std::make_pair(&vertices[u], uv);
            }

            if (i == 2 * links.size() - 1 || u != links_both_dirs[i+1].first) {
                joins[2*i + 1] = std::make_pair(vu, split_successors[i]);
            } else {
                uintE u2, v2;
                u2 = links_both_dirs[i+1].first;
                v2 = links_both_dirs[i+1].second;

                auto found_element = edge_table.find(std::make_tuple(u2, v2), nullptr);
                joins[2*i + 1] = std::make_pair(vu, found_element);
            }
        });

        sequence<std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>> filtered =
            parlay::filter(joins, [&] (const std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>& e) {
                return e.first != nullptr || e.second != nullptr;
            });

        skip_list.batch_join(&filtered);
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

                    SkipList::SkipListElement* uv = edge_table.find(std::make_pair(u, v), nullptr);
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

            auto splits = sequence<SkipList::SkipListElement*>(4 * cuts.size());
            parallel_for(0, cuts.size(), [&] (size_t i) {
                    if (!ignored[i]) {
                        SkipList::SkipListElement* uv = edge_elements[i];
                        SkipList::SkipListElement* vu = uv->twin;

                        splits[4 * i] = uv;
                        splits[4 * i + 1] = vu;

                        SkipList::SkipListElement* predecessor = uv->get_left(0);
                        if (predecessor != nullptr) {
                            splits[4 * i + 2] = predecessor;
                        }

                        predecessor = vu->get_left(0);
                        if (predecessor != nullptr) {
                            splits[4 * i + 3] = predecessor;
                        }

                    }
           });

           sequence<SkipList::SkipListElement*> filtered =
                parlay::filter(splits, [&] (const SkipList::SkipListElement* e) {
                    return e != nullptr;
                });
           skip_list.batch_split(&filtered);

           auto joins = sequence<std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>>(
                join_targets.size() / 2,
                std::make_pair(nullptr, nullptr));
           parallel_for(0, cuts.size(), [&] (size_t i) {
                    if (!ignored[i]) {
                        /*uintE u, v;
                        u = cuts[i].first;
                        v = cuts[i].second;*/

                        //edge_table.remove(std::make_pair(u, v));

                        if (join_targets[4 * i] != nullptr) {
                            joins[2 * i] = std::make_pair(join_targets[4 * i], join_targets[4 * i + 1]);
                        }

                        if (join_targets[4 * i + 2] != nullptr) {
                            joins[2 * i + 1] = std::make_pair(join_targets[4 * i + 2], join_targets[4 * i + 3]);
                        }
                    }
            });
            sequence<std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>> filtered_joins =
                parlay::filter(joins, [&] (const std::pair<SkipList::SkipListElement*,
                            SkipList::SkipListElement*>& e) {
                    return e.first != nullptr || e.second != nullptr;
            });

            skip_list.batch_join(&filtered_joins);

            auto element_indices = parlay::pack_index(ignored);
            auto next_cuts_seq = sequence<std::pair<uintE, uintE>>(element_indices.size());
            parallel_for(0, next_cuts_seq.size(), [&] (size_t i){
                next_cuts_seq[i] = cuts[element_indices[i]];
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

    uintE get_subsequence_sum(uintE v, uintE parent) {
        auto edge = edge_table.find(std::make_tuple(parent, v), nullptr);
        auto twin = edge->twin;

        auto result =  skip_list.get_subsequence_sum(edge, twin);
        auto final_result = result.first ^ result.second;
        return final_result;
    }

    bool is_connected(uintE u, uintE v) {
            auto uu = vertices[u];
            auto vv = vertices[v];

            return skip_list.find_representative(&uu) == skip_list.find_representative(&vv);
    }
};

void RunETTree() {
        std::cout << "ET tree" << std::endl;
        auto tree = ETTree((size_t) 10);

        sequence<std::pair<uintE, uintE>> links = sequence<std::pair<uintE, uintE>>(5);
        links[0] = std::make_pair(2, 3);
        links[1] = std::make_pair(3, 4);
        links[2] = std::make_pair(0, 1);
        links[3] = std::make_pair(0, 8);
        links[4] = std::make_pair(7, 9);

        tree.batch_link(links);
        std::cout << "Connected 2, 3: " << tree.is_connected(2, 3) << std::endl;

}

}  // namespace gbbs
