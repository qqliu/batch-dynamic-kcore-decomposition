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

sequence<sequence<std::pair<uintE, uintE>>> make_values(uintE a, uintE b, double pb, int copies, size_t m) {
        return sequence<sequence<std::pair<uintE, uintE>>>(copies, sequence<std::pair<uintE, uintE>>(
                    ceil(log(m)/log(pb)), std::make_pair(a, b)));
}

struct ETTree {
   sequence<SkipList::SkipListElement> edge_table;
   SkipList skip_list;
   sequence<SkipList::SkipListElement> vertices;
   int copies;
   size_t m;
   parlay::random cutset_rng = parlay::random(time(0));
   double pb;

   ETTree() {}

   ETTree(size_t n, int copies_, size_t m_, double pb_) {
        skip_list = SkipList(n);
        edge_table = sequence<SkipList::SkipListElement>(2 * m_);
        pb = pb_;
        m = m_;
        copies = copies_;

        vertices = sequence<SkipList::SkipListElement>(n);
        auto joins = sequence<std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>>(n);
        parallel_for(0, n, [&] (size_t i) {
            vertices[i] = skip_list.create_node(i, nullptr, nullptr, make_values(0, 0, pb, copies, m), nullptr, true,
                    std::make_pair(i, i), pb, copies, m);
            joins[i] = std::make_pair(&vertices[i], &vertices[i]);
        });

        skip_list.batch_join(&joins);
    }

    void print_value(std::string label, SkipList::SkipListElement* v) {
            std::cout << label << ", id: " << v->id.first << ", " << v->id.second <<
            ", values: " << v->values[0][0][0].first << ", " << v->values[0][0][0].second << std::endl;
    }

    /*void create_edge_in_edge_table(uintE u, uintE v) {
        if (edge_table[u][v].id.first == UINT_E_MAX) {
            edge_table[u][v] = skip_list.create_node(u, nullptr, nullptr, make_values(0, 0, pb, copies, m), nullptr, false,
                std::make_pair(u, v), pb, copies, m);
            auto uv = &edge_table[u][v];
            edge_table[v][u] = skip_list.create_node(v, nullptr, nullptr, make_values(0, 0, pb, copies, m), uv, false,
                std::make_pair(v, u), pb, copies, m);
            auto vu = &edge_table[v][u];
        }
    }*/

    template <class KY, class VL, class HH>
    void link(uintE u, uintE v, pbbslib::sparse_table<KY, VL, HH> edge_index_table) {
        auto index_uv = edge_index_table.find(std::make_pair(u, v), UINT_E_MAX); //, UINT_E_MAX);
        auto index_vu = edge_index_table.find(std::make_pair(v, u), UINT_E_MAX); //, UINT_E_MAX);

        if (index_uv == UINT_E_MAX || index_vu == UINT_E_MAX)
            std::cout << "THERE IS AN ERROR in edge_index_table" << std::endl;

        edge_table[index_uv] = skip_list.create_node(u, nullptr, nullptr, make_values(0, 0, pb, copies, m), nullptr, false,
                std::make_pair(u, v), pb, copies, m);
        auto uv = &edge_table[index_uv];
        edge_table[index_vu] = skip_list.create_node(v, nullptr, nullptr, make_values(0, 0, pb, copies, m), uv, false,
                std::make_pair(v, u), pb, copies, m);
        auto vu = &edge_table[index_vu];
        uv->twin = vu;
        vu->twin = uv;
        /*print_value("uv value: ", uv);
        print_value("vu value: ", vu);
        print_value("uv twin value: ", uv->twin);*/

        auto u_left = &vertices[u];
        auto v_left = &vertices[v];

       /*print_value("u_left: ", u_left);
       print_value("v_left: ", v_left);*/

        auto splits = sequence<SkipList::SkipListElement*>(2);
        splits[0] = u_left;
        splits[1] = v_left;
        //std::cout << "starting split" << std::endl;
        auto results = skip_list.batch_split(&splits);

        auto u_right = results[0];
        auto v_right = results[1];
        /*print_value("u_right ", u_right);
        print_value("v_right ", v_right);*/

        auto joins = sequence<std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>>(4);

        joins[0] = std::make_pair(u_left, uv);
        joins[1] = std::make_pair(uv, v_right);
        joins[2] = std::make_pair(v_left, vu);
        joins[3] = std::make_pair(vu, u_right);

        skip_list.batch_join(&joins);
        auto new_u_right = u_left->get_right(0);
        //print_value("u_left right: ", new_u_right);
        auto new_rr = new_u_right->get_right(0);
        //print_value("new right right: ", new_rr);
        auto rrr = new_rr->get_right(0);
        auto rrrr = rrr->get_right(0);
        /*print_value("rrr: ", rrr);
        print_value("rrrr: ", rrrr);*/
        auto rrrrr = rrrr->get_right(0);
        //print_value("rrrrr: ", rrrrr);
        auto rrrrrr = rrrrr->get_right(0);
        //print_value("rrrrrr: ", rrrrrr);
        auto rrrrrrr = rrrrrr->get_right(0);
        //print_value("rrrrrrr: ", rrrrrrr);
    }

    template <class KY, class VL, class HH>
    void cut(int u, int v, pbbslib::sparse_table<KY, VL, HH> edge_index_table){
            auto index_uv = edge_index_table.find(std::make_pair(u, v), UINT_E_MAX); //, UINT_E_MAX);
            if (index_uv == UINT_E_MAX)
                std::cout << "There is an error in edge_index_table" << std::endl;

            auto uv = &edge_table[index_uv];
            auto vu = uv->twin;
            /*print_value("cut uv: ", uv);
            print_value("cut vu: ", vu);*/
            //edge_table.remove(std::make_pair(u, v));
            //edge_table.remove(std::make_pair(v, u));

            //std::cout << "get lefts" << std::endl;
            auto u_left = uv->get_left(0);
            auto v_left = vu->get_left(0);
            /*print_value("u_left", u_left);
            print_value("v_left", v_left);*/

            /*std::cout << "correct rights" << std::endl;
            print_value("v_right: ", uv->get_right(0));
            print_value("u_right: ", vu->get_right(0));*/


            auto splits = sequence<SkipList::SkipListElement*>(2);
            splits[0] = uv;
            splits[1] = vu;
            auto results = skip_list.batch_split(&splits);

            //std::cout << "get rights" << std::endl;
            auto v_right = results[0];
            auto u_right = results[1];
            /*print_value("v_right", v_right);
            print_value("u_right", u_right);*/

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

            //std::cout << "inserted first" << std::endl;
            edge_table[index_uv] = SkipList::SkipListElement();
            auto index_vu = edge_index_table.find(std::make_pair(v, u), UINT_E_MAX); //, UINT_E_MAX);
            //std::cout << "inserted second" << std::endl;
            edge_table[index_vu] = SkipList::SkipListElement();
    }

    template <class KY, class VL, class HH>
    void batch_link_sequential(sequence<std::pair<uintE, uintE>>links,
            pbbslib::sparse_table<KY, VL, HH> edge_index_table) {
            for(size_t i = 0; i < links.size(); i++) {
                    std::cout << "u: " << links[i].first << ", "
                        << "v: " << links[i].second << std::endl;

                    link(links[i].first, links[i].second, edge_index_table);
                    std::cout << "end current link" << std::endl;

                    /*auto u_right = vertices[links[i].first].get_right(0);
                    auto v_right = vertices[links[i].second].get_right(0);
                    auto uv_right = u_right->get_right(0);
                    auto vu_right = v_right->get_right(0);
                    if (vu_right == nullptr)
                        std::cout << "NOOOO" << std::endl;

                    std::cout << "u right: " << u_right->values[0].first << ", " << u_right->values[0].second << std::endl;
                    std::cout << "vu right: " << vu_right->values[0].first << ", " << vu_right->values[0].second << std::endl;
                    std::cout << "v right: " << v_right->values[0].first << ", " << v_right->values[0].second << std::endl;
                    std::cout << "uv right: " << uv_right->values[0].first << ", " << uv_right->values[0].second << std::endl;*/

            }
    }

    template <class KY, class VL, class HH>
    void batch_link(sequence<std::pair<uintE, uintE>>links, pbbslib::sparse_table<KY, VL, HH>& edge_index_table) {
        if (links.size() <= 75) {
                batch_link_sequential(links, edge_index_table);
                return;
        }

        //std::cout << "start batch links" << std::endl;
        sequence<std::pair<uintE, uintE>> links_both_dirs = sequence<std::pair<uintE, uintE>>(2 * links.size());
        parallel_for(0, links.size(), [&] (size_t i) {
                links_both_dirs[2 * i] = links[i];
                links_both_dirs[2 * i + 1] = std::make_pair(links[i].second, links[i].first);
        });
        // std::cout << "get bidirectional links" << std::endl;

        auto get_key = [&] (const std::pair<uintE, uintE>& elm) { return elm.first; };
        parlay::integer_sort_inplace(parlay::make_slice(links_both_dirs), get_key);

        auto split_successors = sequence<SkipList::SkipListElement*>(2 * links.size());
        auto splits = sequence<SkipList::SkipListElement*>(2 * links.size(), nullptr);
        //std::cout << "initialize successors" << std::endl;

        parallel_for(0, 2 * links.size(), [&] (size_t i) {
            uintE u, v;
            u = links_both_dirs[i].first;
            v = links_both_dirs[i].second;

            if (i == 2 * links.size() - 1 || u != links_both_dirs[i+1].first) {
                    splits[i] = &vertices[u];
            }

            if (u < v) {
                    auto index_uv = edge_index_table.find(std::make_pair(u, v), UINT_E_MAX); //, UINT_E_MAX);
                    if (index_uv == UINT_E_MAX)
                        std::cout << "This is an error in edge_index_table in recursive insertions" << std::endl;

                    edge_table[index_uv] = skip_list.create_node(u, nullptr, nullptr, make_values(0, 0, pb, copies, m),
                            nullptr, false, std::make_pair(u, v), pb, copies, m);
                    auto uv = &edge_table[index_uv];

                    auto index_vu = edge_index_table.find(std::make_pair(v, u), UINT_E_MAX); //, UINT_E_MAX);
                    if (index_vu == UINT_E_MAX)
                        std::cout << "This is an error in edge_index_table in recursive insertions" << std::endl;

                    edge_table[index_vu] = skip_list.create_node(v, nullptr, nullptr, make_values(0, 0, pb, copies, m),
                            uv, false, std::make_pair(v, u), pb, copies, m);
                    auto vu = &edge_table[index_vu];
                    vu->twin = uv;
                    uv->twin = vu;
            }
        });
        //std::cout << "successfully inserted splits" << std::endl;

        auto bool_seq = parlay::delayed_seq<bool>(splits.size(), [&] (size_t i) {
                return (splits[i] != nullptr);
        });

        auto element_indices = parlay::pack_index(bool_seq);
        auto filtered_splits = sequence<SkipList::SkipListElement*>(element_indices.size());
        parallel_for(0, filtered_splits.size(), [&] (size_t i) {
            filtered_splits[i] = splits[element_indices[i]];
        });
        //std::cout << "succesfully obtain filtered splits" << std::endl;

        auto results = skip_list.batch_split(&filtered_splits);

        parallel_for(0, filtered_splits.size(), [&] (size_t i) {
            auto split_index = element_indices[i];
            split_successors[split_index] = results[i];
        });
        //std::cout << "succesfully obtained successors" << std::endl;

        auto joins = sequence<std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>>(4 * links.size(),
                std::make_pair(nullptr, nullptr));
        parallel_for(0, 2 * links.size(), [&] (size_t i) {
            uintE u, v;
            u = links_both_dirs[i].first;
            v = links_both_dirs[i].second;

            auto index_uv = edge_index_table.find(std::make_pair(u, v), UINT_E_MAX); //, UINT_E_MAX);
            if (index_uv == UINT_E_MAX)
                std::cout << "This is an error in edge_index_table in recursive insertions" << std::endl;

            SkipList::SkipListElement* uv = &edge_table[index_uv];
            //print_value("uv value: ", uv);
            SkipList::SkipListElement* vu = uv -> twin;

            if (i == 0 || u != links_both_dirs[i-1].first) {
                /*if (&vertices[u] == nullptr)
                    std::cout << "One has error" << std::endl;
                if (uv == nullptr)
                    std::cout << "Two has error" << std::endl;*/
                joins[2*i] = std::make_pair(&vertices[u], uv);
            }

            if (i == 2 * links.size() - 1 || u != links_both_dirs[i+1].first) {
                /*if (vu == nullptr)
                    std::cout << "Three has error" << std::endl;
                if (split_successors[i] == nullptr)
                    std::cout << "Four has error" << std::endl;*/
                joins[2*i + 1] = std::make_pair(vu, split_successors[i]);
            } else {
                uintE u2, v2;
                u2 = links_both_dirs[i+1].first;
                v2 = links_both_dirs[i+1].second;
                auto index_uv2 = edge_index_table.find(std::make_pair(u2, v2), UINT_E_MAX); //, UINT_E_MAX);
                if (index_uv2 == UINT_E_MAX)
                    std::cout << "This is an error in edge_index_table in recursive insertions uv2" << std::endl;


                auto found_element = &edge_table[index_uv2];
                /*if (vu == nullptr)
                    std::cout << "Five has error" << std::endl;
                if (found_element == nullptr)
                    std::cout << "Six has error" << std::endl;*/
                joins[2*i + 1] = std::make_pair(vu, found_element);
            }
        });

        //std::cout << "successfully compute joins" << std::endl;

        sequence<std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>> filtered =
            parlay::filter(joins, [&] (const std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>& e) {
                    /*if (e.first != nullptr && e.second != nullptr) {
                        std::cout << "works" << std::endl;
                        print_value("first " , e.first);
                        print_value("second ", e.second);
                     } else if ((e.first == nullptr || e.second == nullptr)
                             && (e.first != nullptr || e.second != nullptr)) {
                        std::cout << "IMPOSSIBLE" << std::endl;
                     }*/

                return e.first != nullptr && e.second != nullptr;
            });

        skip_list.batch_join(&filtered);
    }

    template <class KY, class VL, class HH>
    void batch_cut_sequential(sequence<std::pair<uintE, uintE>> cuts,
            pbbslib::sparse_table<KY, VL, HH> edge_index_table) {
            //std::cout << "cutting sequentially" << std::endl;
            for (size_t i = 0; i < cuts.size(); i++)
                cut(cuts[i].first, cuts[i].second, edge_index_table);
    }

    template <class KY, class VL, class HH>
    void batch_cut_recurse(sequence<std::pair<uintE, uintE>> cuts,
        pbbslib::sparse_table<KY, VL, HH> edge_index_table) {
            sequence<SkipList::SkipListElement*> join_targets =
                sequence<SkipList::SkipListElement*>(4 * cuts.size(), nullptr);
            sequence<SkipList::SkipListElement*> edge_elements =
                sequence<SkipList::SkipListElement*>(cuts.size(), nullptr);

            parlay::random rng = parlay::random(time(0));
            //std::cout << "cuts size: !!!!!!!!!!!" << cuts.size() << std::endl;

            if (cuts.size() <= 75) {
                //std::cout << "PASSED IF 2" << std::endl;
                    batch_cut_sequential(cuts, edge_index_table);
                    return;
            }
            sequence<bool> ignored = sequence<bool>(cuts.size(), true);

            parallel_for(0, cuts.size(), [&] (size_t i) {
                rng.fork(i);
                rng = rng.next();
                bool rand_val = rng.rand() % 100 == 0;

                ignored[i] = rand_val;
                if (!ignored[i]) {
                    uintE u, v;
                    u = cuts[i].first;
                    v = cuts[i].second;

                    auto index_uv = edge_index_table.find(std::make_pair(u, v), UINT_E_MAX); //, UINT_E_MAX);
                    if (index_uv == UINT_E_MAX)
                        std::cout << "This is an error in edge_index_table in batch deletions" << std::endl;

                    SkipList::SkipListElement* uv = &edge_table[index_uv];
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
                            }
                            join_targets[4 * i] = left_target;
                            join_targets[4 * i + 1] = right_target;
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
                        uintE u, v;
                        u = cuts[i].first;
                        v = cuts[i].second;

                        auto index_uv = edge_index_table.find(std::make_pair(u, v), UINT_E_MAX); //, UINT_E_MAX);
                        if (index_uv == UINT_E_MAX)
                            std::cout << "This is an error in edge_index_table in batch deletions" << std::endl;

                        //edge_table.remove(std::make_pair(u, v));
                        edge_table[index_uv] = SkipList::SkipListElement();
                        auto index_vu = edge_index_table.find(std::make_pair(v, u), UINT_E_MAX); //, UINT_E_MAX);
                        if (index_vu == UINT_E_MAX)
                            std::cout << "This is an error in edge_index_table in batch deletions" << std::endl;

                        edge_table[index_vu] = SkipList::SkipListElement();

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
                        /*if ((e.first != nullptr || e.second != nullptr) && (e.first == nullptr || e.second == nullptr))
                            std::cout << "something is WRONNNNNNNGGGGGG" << std::endl;*/
                    return e.first != nullptr || e.second != nullptr;
            });

            skip_list.batch_join(&filtered_joins);

            auto element_indices = parlay::pack_index(ignored);
            auto next_cuts_seq = sequence<std::pair<uintE, uintE>>(element_indices.size());
            parallel_for(0, next_cuts_seq.size(), [&] (size_t i){
                next_cuts_seq[i] = cuts[element_indices[i]];
            });
            batch_cut_recurse(next_cuts_seq, edge_index_table);
    }

    template <class KY, class VL, class HH>
    void batch_cut(sequence<std::pair<uintE, uintE>> cuts, pbbslib::sparse_table<KY, VL, HH> edge_index_table) {
            if (cuts.size() <=  75) {
                //std::cout << "PASSED IF 1" << std::endl;
                    batch_cut_sequential(cuts, edge_index_table);
                    return;
            }

            batch_cut_recurse(cuts, edge_index_table);
    }

    void add_edge_to_cutsets(std::pair<uintE, uintE> edge) {
        auto min_edge = std::min(edge.first, edge.second);
        auto max_edge = std::max(edge.first, edge.second);
        auto u = &vertices[edge.first];

        parallel_for(0, u->values[0].size(), [&](size_t ii) {
            parallel_for(0, u->values[0][ii].size(), [&](size_t ij) {
                    cutset_rng.fork(edge.first + ii * ij);
                    cutset_rng = cutset_rng.next();
                    auto base = uintE(floor(pow(pb, ij)));
                    if (base == 1)
                        base = pb;
                    auto rand_val_u = cutset_rng.rand() % base;
                    //std::cout << "rand_val " << rand_val_u << ", base: " << base << std::endl;
                    if (rand_val_u <= 0)
                        u->values[0][ii][ij] = std::make_pair(u->values[0][ii][ij].first ^ min_edge,
                                u->values[0][ii][ij].second ^ max_edge);

                    cutset_rng = cutset_rng.next();
                    /*std::cout << "new values: " << u->values[0][ii][ij].first << ", " << u->values[0][ii][ij].second
                        << std::endl;*/
            });
        });
    }

    template <class KY, class VL, class HH>
    sequence<sequence<std::pair<uintE, uintE>>> get_subtree_sum(uintE v, uintE parent,
            pbbslib::sparse_table<KY, VL, HH> edge_index_table) {
        auto index_uv = edge_index_table.find(std::make_pair(parent, v), UINT_E_MAX); //, UINT_E_MAX);
        if (index_uv == UINT_E_MAX)
            std::cout << "This is an error in edge_index_table in get_subtree_sum" << std::endl;

        auto edge = &edge_table[index_uv];
        auto twin = edge->twin;

        /*print_value("v: ", edge);
        print_value("parent: ", twin);
        auto right = edge->get_right(0);
        while (right != twin) {
                print_value("right: ", right);
                right = right->get_right(0);
        }*/

        auto result =  skip_list.get_subsequence_sum(edge, twin);
        /*auto left_v = edge->get_left(0);
        while(left_v->is_vertex && left_v != edge) {
                left_v = left_v -> get_left(0);
        }*/
        /*result.first = result.first ^ left_v->id.first ^ twin->id.first;
        result.second = result.second ^ left_v->id.second ^ twin->id.second;*/
        /*print_value("left_v ", left_v);
        print_value("twin ", twin);*/
        /*print_value("edge: ", edge);
        print_value("twin: ", twin);
        auto right = edge->get_right(0);
        while (right != twin) {
                print_value("right: ", right);
                right = right->get_right(0);
        }*/
        return result;
        //auto final_result = result.first ^ result.second;
        //return final_result;
    }

    sequence<sequence<std::pair<uintE, uintE>>> get_tree_sum(uintE v) {
        auto edge = &vertices[v];
        auto edge_sum = skip_list.get_sum(edge);
        return edge_sum;
        //return edge_sum.first ^ edge_sum.second;
    }

    bool is_connected(uintE u, uintE v) {
            auto uu = vertices[u];
            auto vv = vertices[v];
            /*std::cout << "finding representative: " << skip_list.find_representative(&uu)->values[0].first
                << ", " << skip_list.find_representative(&uu)->values[0].second << std::endl;
            std::cout << "finding representative: " << skip_list.find_representative(&vv)->values[0].first
                << ", " << skip_list.find_representative(&vv)->values[0].second << std::endl;*/

            return skip_list.find_representative(&uu) == skip_list.find_representative(&vv);
    }
};

void RunETTree(double pb, int copies, size_t m) {
        std::cout << "ET tree" << std::endl;
        using K = std::pair<uintE, uintE>;
        using V = uintE;
        using KV = std::pair<K, V>;

        KV empty =
            std::make_pair(std::make_pair(UINT_E_MAX, UINT_E_MAX), UINT_E_MAX);

        auto hash_pair = [](const std::pair<uintE, uintE>& t) {
            size_t l = std::min(std::get<0>(t), std::get<1>(t));
            size_t r = std::max(std::get<0>(t), std::get<1>(t));
            size_t key = (l << 32) + r;
            return parlay::hash64_2(key);
        };

        auto edge_index_table =
            pbbslib::make_sparse_table<K, V>(2 * 8, empty, hash_pair);
        bool abort = false;

        auto tree = ETTree((size_t) 10, copies, m, pb);

        sequence<std::pair<uintE, uintE>> links = sequence<std::pair<uintE, uintE>>(8);
        links[0] = std::make_pair(2, 3);
        edge_index_table.insert_check(std::make_pair(std::make_pair(2, 3), 0), &abort);
        edge_index_table.insert_check(std::make_pair(std::make_pair(3, 2), 1), &abort);
        links[1] = std::make_pair(3, 4);
        edge_index_table.insert_check(std::make_pair(std::make_pair(3, 4), 2), &abort);
        edge_index_table.insert_check(std::make_pair(std::make_pair(4, 3), 3), &abort);
        links[2] = std::make_pair(0, 1);
        edge_index_table.insert_check(std::make_pair(std::make_pair(0, 1), 4), &abort);
        edge_index_table.insert_check(std::make_pair(std::make_pair(1, 0), 5), &abort);
        links[3] = std::make_pair(0, 8);
        edge_index_table.insert_check(std::make_pair(std::make_pair(0, 8), 6), &abort);
        edge_index_table.insert_check(std::make_pair(std::make_pair(8, 0), 7), &abort);
        links[4] = std::make_pair(7, 2);
        edge_index_table.insert_check(std::make_pair(std::make_pair(7, 2), 8), &abort);
        edge_index_table.insert_check(std::make_pair(std::make_pair(2, 7), 9), &abort);
        links[5] = std::make_pair(1, 2);
        edge_index_table.insert_check(std::make_pair(std::make_pair(1, 2), 10), &abort);
        edge_index_table.insert_check(std::make_pair(std::make_pair(2, 1), 11), &abort);
        links[6] = std::make_pair(6, 5);
        edge_index_table.insert_check(std::make_pair(std::make_pair(6, 5), 12), &abort);
        edge_index_table.insert_check(std::make_pair(std::make_pair(5, 6), 13), &abort);
        links[7] = std::make_pair(6, 9);
        edge_index_table.insert_check(std::make_pair(std::make_pair(6, 9), 14), &abort);
        edge_index_table.insert_check(std::make_pair(std::make_pair(9, 6), 15), &abort);

        tree.batch_link(links, edge_index_table);
        std::cout << "Connected 2, 3: " << tree.is_connected(2, 3) << std::endl;
        std::cout << "Connected 2, 0: " << tree.is_connected(2, 0) << std::endl;
        std::cout << "Connected 2, 4: " << tree.is_connected(2, 4) << std::endl;
        std::cout << "Connected 3, 7: " << tree.is_connected(3, 7) << std::endl;
        std::cout << "Connected 3, 8: " << tree.is_connected(3, 8) << std::endl;
        std::cout << "Connected 1, 8: " << tree.is_connected(1, 8) << std::endl;
        std::cout << "Connected 1, 5: " << tree.is_connected(1, 5) << std::endl;
        std::cout << "Connected 5, 9: " << tree.is_connected(5, 9) << std::endl;

        std::cout << "Getting subsequence sums: " << std::endl;
        auto a = tree.get_subtree_sum(2, 3, edge_index_table);
        auto b = tree.get_subtree_sum(6, 5, edge_index_table);
        auto c = tree.get_tree_sum(2);
        auto d = tree.get_tree_sum(5);

        std::cout << "v: 2, parent: 3: " << a[0][0].first << ", " << a[0][0].second << std::endl;
        std::cout << "v: 6, parent: 5: " << b[0][0].first << ", " << b[0][0].second << std::endl;
        std::cout << "v: 2, parent: 2: " << c[0][0].first << ", " << c[0][0].second << std::endl;
        std::cout << "v: 5, parent: 5: " << d[0][0].first << ", " << d[0][0].second << std::endl;

        sequence<std::pair<uintE, uintE>> cuts = sequence<std::pair<uintE, uintE>>(4);
        cuts[0] = std::make_pair(0, 8);
        cuts[1] = std::make_pair(2, 3);
        cuts[2] = std::make_pair(6, 9);
        cuts[3] = std::make_pair(0, 1);

        tree.batch_cut(cuts, edge_index_table);
        std::cout << "After Batch Cut" << std::endl;
        std::cout << "Connected 2, 3: " << tree.is_connected(2, 3) << std::endl;
        std::cout << "Connected 2, 0: " << tree.is_connected(2, 0) << std::endl;
        std::cout << "Connected 2, 4: " << tree.is_connected(2, 4) << std::endl;
        std::cout << "Connected 3, 7: " << tree.is_connected(3, 7) << std::endl;
        std::cout << "Connected 3, 4: " << tree.is_connected(3, 4) << std::endl;
        std::cout << "Connected 1, 8: " << tree.is_connected(1, 8) << std::endl;
        std::cout << "Connected 1, 5: " << tree.is_connected(1, 5) << std::endl;
        std::cout << "Connected 5, 9: " << tree.is_connected(5, 9) << std::endl;
        std::cout << "Connected 5, 6: " << tree.is_connected(5, 6) << std::endl;

        std::cout << "Getting subsequence sums: " << std::endl;
        a = tree.get_subtree_sum(4, 3, edge_index_table);
        b = tree.get_subtree_sum(6, 5, edge_index_table);

        std::cout << "v: 4, parent: 3: " << a[0][0].first << ", " << a[0][0].second << std::endl;
        std::cout << "v: 6, parent: 5: " << b[0][0].first << ", " << b[0][0].second << std::endl;
}

}  // namespace gbbs
