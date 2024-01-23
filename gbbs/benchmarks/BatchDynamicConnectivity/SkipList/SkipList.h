#include <unordered_set>
#include <stack>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
#include "gbbs/pbbslib/sparse_table.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "benchmarks/EdgeOrientation/ParallelLDS/sparse_set.h"


namespace gbbs {

struct SkipList {

    struct SkipListElement {
        size_t height;
        size_t lowest_needs_update = 0;

        using pointers = std::pair<SkipListElement*, SkipListElement*>;
        using height_array = sequence<pointers>;
        using values_array = sequence<sequence<sequence<std::pair<uintE, uintE>>>>;

        height_array elements;
        values_array values;
        uintE update_level;
        bool split_mark;
        SkipListElement* twin;
        bool is_vertex;
        std::pair<uintE, uintE> id;
        double probability_base;
        int num_duplicates;

        SkipListElement(): height(0), lowest_needs_update(0) { update_level = UINT_E_MAX; split_mark = false;
            twin = nullptr;
            is_vertex = false;
            id = std::make_pair(UINT_E_MAX, UINT_E_MAX);
            probability_base = 2;
            num_duplicates = 2;
        }

        SkipListElement(size_t _h, SkipListElement* _r, SkipListElement* _l,
                sequence<sequence<std::pair<uintE, uintE>>> _vals,
                SkipListElement* twin_ = nullptr, bool is_vertex_ = false, std::pair<uintE, uintE>id_ =
                std::make_pair(UINT_E_MAX, UINT_E_MAX), double pb = 2, int num_dup = 2, size_t m = 10):
            height(_h), lowest_needs_update(_h) {
                elements.resize(_h);
                update_level = UINT_E_MAX;
                elements[0].first = _l;
                elements[0].second = _r;
                values.resize(_h);
                values[0] = _vals;
                split_mark = false;
                twin = twin_;
                is_vertex = is_vertex_;
                id = id_;
                probability_base = pb;
                num_duplicates = num_dup;

                parallel_for(1, _h, [&](size_t i){
                    values[i] = sequence<sequence<std::pair<uintE, uintE>>>(num_duplicates,
                            sequence<std::pair<uintE, uintE>>(ceil(log(m)/log(pb)), std::make_pair(0, 0)));
                });
        }

        inline void set_left_pointer(size_t height, SkipListElement* left) {
            elements[height].first = left;
        }

        inline void set_right_pointer(size_t height, SkipListElement* right) {
            elements[height].second = right;
        }

        inline bool CASleft(size_t level,
                SkipListElement* old_left, SkipListElement* new_left) {
               return pbbslib::atomic_compare_and_swap(&elements[level].first, old_left, new_left);
        }

        inline bool CASright(size_t level,
                SkipListElement* old_right, SkipListElement* new_right) {
               return pbbslib::atomic_compare_and_swap(&elements[level].second, old_right, new_right);
        }

        /*inline bool CASvalue(std::pair<uintE, uintE> old_value, std::pair<uintE, uintE> new_value) {
               return pbbslib::atomic_compare_and_swap(&values[0], old_value, new_value);
               values[0] = new_value;
               return true;
        }*/

        inline SkipListElement* get_left(size_t height) {
                return elements[height].first;
        }

        inline SkipListElement* get_right(size_t height) {
                return elements[height].second;
        }
    };

    size_t n;
    parlay::random rng = parlay::random(time(0));

    SkipList(): n(0) {}

    SkipList(size_t _n): n(_n) {}

    SkipListElement create_node(size_t index, SkipListElement* left, SkipListElement* right,
            sequence<sequence<std::pair<uintE, uintE>>> vals,
            SkipListElement* twin = nullptr, bool is_vertex = false,
            std::pair<uintE, uintE> id = std::make_pair(UINT_E_MAX, UINT_E_MAX),
            double pb = 2, int num_dup = 2, size_t m = 10) {
        rng.fork(index);
        rng = rng.next();
        auto rand_val = rng.rand() % UINT_E_MAX;
        rng = rng.next();

        size_t cur_height = 1;
        while (rand_val & 1) {
                rand_val >>= 1;
                cur_height++;
        }

        auto height = std::min(cur_height, (size_t) 32);
        //std::cout << "initializing node" << std::endl;

        auto node = SkipListElement(height, left, right, vals, twin, is_vertex, id);
        return node;
    }

    SkipListElement* find_left_parent(size_t level, SkipListElement* this_element) {
            SkipListElement* cur_element = this_element;
            SkipListElement* start_element = this_element;
            if (cur_element->height > level + 1)
                return cur_element;
            cur_element = cur_element->elements[level].first;

            while (cur_element != nullptr && cur_element != start_element) {
                if (cur_element->height > level + 1)
                    return cur_element;
                cur_element = cur_element->elements[level].first;
            }
            return nullptr;
    }

    SkipListElement* find_right_parent(size_t level, SkipListElement* this_element) {
            SkipListElement* cur_element = this_element;
            SkipListElement* start_element = this_element;

            if (cur_element->height > level+1)
                return cur_element;
            cur_element = cur_element->elements[level].second;

            while (cur_element != nullptr && cur_element != start_element) {
                if (cur_element->height > level + 1)
                    return cur_element;
                cur_element = cur_element->elements[level].second;
            }
            return nullptr;
    }

    SkipListElement* find_representative(SkipListElement* this_element) {
        ///std::cout << "started find representative" << std::endl;
        SkipListElement* cur_element = this_element;
        SkipListElement* seen_element = nullptr;

        size_t current_level = cur_element->height - 1;

        while(cur_element->elements[current_level].second != nullptr &&
                seen_element != cur_element) {
                //std::cout << "started while" << std::endl;
                if (seen_element == nullptr || cur_element < seen_element)
                    seen_element = cur_element;
                /*std::cout << "start middle while" << std::endl;
                if (cur_element->elements[current_level].second == nullptr)
                    std::cout << "found null ptr in while" << std::endl;
                std::cout << "height: " << cur_element->elements[current_level].second->id.first << std::endl;*/
                cur_element = cur_element->elements[current_level].second;
                //std::cout << "end middle while" << std::endl;
                auto top_level = cur_element->height-1;
                //std::cout << "top_level: " << top_level << std::endl;
                if (current_level < top_level) {
                        current_level = top_level;
                        seen_element = nullptr;
                }
                //std::cout << "end while" << std::endl;
        }
        //std::cout << "finished find representative while" << std::endl;

        if (seen_element == cur_element) {
        //std::cout << "finished find representative" << std::endl;
                return seen_element;
        } else {
                while (cur_element->elements[current_level].first != nullptr) {
                        cur_element = cur_element->elements[current_level].first;
                        current_level = cur_element->height - 1;
                }
        //std::cout << "finished find representative" << std::endl;
                return cur_element;
        }
        //std::cout << "finished find representative" << std::endl;
    }

    /*void update_edge_values(sequence<SkipListElement*> nodes) {
            parallel_for(0, nodes.size(), [&](size_t i) {
                if (!nodes[i]->is_vertex && nodes[i] != nullptr) {
                    auto left = nodes[i]->get_left(0);
                    while (left->is_vertex && left != nodes[i]) {
                        left = left->get_left(0);
                    }
                    if (left != nodes[i]) {
                        nodes[i] -> CASvalue(nodes[i]->values[0],
                           std::make_pair(left->id.first ^ nodes[i]->id.first,
                                left->id.second ^ nodes[i]->id.second));
                   }
                }
            });
    }*/

    void join(SkipListElement* left, SkipListElement* right) {
            size_t level = 0;
            while(left != nullptr && right != nullptr) {
                    if (left->elements[level].second == nullptr &&
                                left->CASright(level, nullptr, right)) {
                            right->CASleft(level, nullptr, left);
                            left = find_left_parent(level, left);
                            right = find_right_parent(level, right);
                            level++;
                    } else {
                            return;
                    }
            }
    }

    SkipListElement* split(SkipListElement* this_element) {
            SkipListElement* successor = nullptr; //this_element->get_right(0);
            SkipListElement* cur_element = this_element;

            size_t level = 0;

            //SkipListElement* next;
            while(cur_element != nullptr) { // && (next = cur_element->elements[level].second) != nullptr) {
                SkipListElement* next = cur_element->elements[level].second;
                if (next != nullptr && //(next = cur_element->elements[level].second) != nullptr) {
                        cur_element->CASright(level, next, nullptr)) {
                        if (level == 0) {
                           successor = next;
                           /*if (cur_element->is_vertex)
                                cur_element->values[0] = std::make_pair(0, 0);*/
                        }
                        //cur_element->elements[level].second = nullptr;
                        next->elements[level].first = nullptr;
                        cur_element = find_left_parent(level, cur_element);
                        level++;
               }
               else {
                        break;
               }
            }
            return successor;
    }

    /* The augmented skip list can take any arbitrary associative, commutative function. Here, it is
     * implemented using the XOR function. */
    void update_top_down(size_t level, SkipListElement* this_element) {
            if (level == 0) {
                    if (this_element->height == 1) {
                        this_element->update_level = UINT_E_MAX;
                    }
                    return;
            }

            if (this_element->update_level < level) {
                    update_top_down(level - 1, this_element);
            }

            sequence<sequence<std::pair<uintE, uintE>>> xor_total = this_element->values[level-1];
            SkipListElement* curr = this_element->elements[level-1].second;
            while (curr != nullptr && curr->height < level + 1) {
                    if (curr->update_level != UINT_E_MAX && curr->update_level < level) {
                            update_top_down(level-1, curr);
                    }

                    /*auto first_level_values = curr->values[0];
                    if (level-1 != 0 || (level-1 == 0 && first_level_values.first != first_level_values.second)) {*/
                    parallel_for(0, xor_total.size(), [&] (size_t ii) {
                            parallel_for(0, xor_total[ii].size(), [&] (size_t ij) {
                                xor_total[ii][ij].first ^= curr->values[level-1][ii][ij].first;
                                xor_total[ii][ij].second ^= curr->values[level-1][ii][ij].second;
                            });
                    });
                    //}
                    curr = curr->elements[level-1].second;
            }
            this_element->values[level] = xor_total;

            if(this_element->height == level+1) {
                    this_element->update_level = UINT_E_MAX;
            }
    }

    void batch_update(sequence<std::pair<SkipListElement*,
            sequence<sequence<std::pair<uintE, uintE>>>>>* new_values) {
        auto top_nodes = sequence<SkipListElement*>(new_values->size(), nullptr);
        sequence<std::pair<SkipListElement*, sequence<sequence<std::pair<uintE, uintE>>>>>&
            new_values_ref = *new_values;
        if (new_values != nullptr) {
            parallel_for(0, new_values->size(), [&](size_t i) {
                    SkipListElement* this_element = new_values_ref[i].first;
                    sequence<sequence<std::pair<uintE, uintE>>>
                        this_element_values = new_values_ref[i].second;

                    parallel_for(0, this_element_values.size(), [&](size_t ii) {
                        parallel_for(0, this_element_values[ii].size(), [&](size_t ij) {
                            this_element->values[0][ii][ij] = this_element_values[ii][ij];
                        });
                    });

                    size_t level = 0;
                    SkipListElement* curr = this_element;
                    while(true) {
                        uintE curr_update_level = curr->update_level;
                        if (curr_update_level == UINT_E_MAX && pbbslib::atomic_compare_and_swap(&curr->update_level,
                                UINT_E_MAX, (uintE) level)) {
                            level = curr->height-1;
                            SkipListElement* parent = find_left_parent(level, curr);

                            if (parent == nullptr) {
                                top_nodes[i] = curr;
                                break;
                            } else {
                                curr = parent;
                                level++;
                            }
                        } else {
                            // Some other execution claimed this ancestor
                            if (curr_update_level > level) {
                                uintE c = curr->update_level;
                                while(c > level && !pbbslib::atomic_compare_and_swap(&curr->update_level, (uintE)c,
                                            (uintE)level))
                                    c = curr->update_level;
                            }
                            top_nodes[i] = nullptr;
                            break;
                        }
                    }
            });

            parallel_for(0, new_values->size(), [&](size_t i){
                    if (top_nodes[i] != nullptr) {
                        update_top_down(top_nodes[i]->height-1, top_nodes[i]);
                    }
            });
        }
    }

    void batch_join(sequence<std::pair<SkipListElement*, SkipListElement*>>* joins) {
            sequence<std::pair<SkipListElement*, SkipListElement*>>& joins_ref = *joins;
            //auto update_nodes = sequence<SkipListElement*>(2 * joins->size());
            auto join_lefts = sequence<std::pair<SkipListElement*,
                 sequence<sequence<std::pair<uintE, uintE>>>>>(joins->size());
            parallel_for(0, joins->size(), [&] (size_t i) {
                /*auto left = joins_ref[i].first;
                auto right = joins_ref[i].second;*/
                join(joins_ref[i].first, joins_ref[i].second);
                //update_nodes[2 * i] = left;
                //update_nodes[2 * i + 1] = right;
                join_lefts[i] = std::make_pair(joins_ref[i].first, joins_ref[i].first->values[0]);
            });

            //update_edge_values(update_nodes);

            /*parallel_for(0, joins->size(), [&] (size_t i) {
                join_lefts[i] = std::make_pair(joins_ref[i].first, joins_ref[i].first->values[0]);
            });*/

            batch_update(&join_lefts);
    }

    sequence<SkipListElement*> batch_split(sequence<SkipListElement*>* splits) {
            sequence<SkipListElement*>& splits_ref = *splits;
            sequence<SkipListElement*> results = sequence<SkipListElement*>(splits->size());
            parallel_for(0, splits->size(), [&](size_t i){
                results[i] = split(splits_ref[i]);
            });

            // Perform updates but only if some other thread hasn't already performed the update
            parallel_for(0, splits->size(), [&](size_t i){
                SkipListElement* curr = splits_ref[i];
                bool can_proceed = curr->update_level == UINT_E_MAX
                    && pbbslib::atomic_compare_and_swap(&curr->update_level, UINT_E_MAX, (uintE)0);

                if (can_proceed) {
                    sequence<sequence<std::pair<uintE, uintE>>> xor_sums = curr->values[0];
                    size_t level = 0;
                    while(true) {
                        if(level < curr->height - 1) {
                            level++;
                            curr->values[level] = xor_sums;
                        } else {
                            curr = curr->elements[level].first;
                            if (curr == nullptr) {
                                break;
                            } else {
                                /*auto first_level_values = curr->values[0];
                                if (level != 0 || (level == 0 &&
                                            first_level_values.first != first_level_values.second)) {*/
                                parallel_for(0, xor_sums.size(), [&](size_t ii) {
                                        parallel_for(0, xor_sums[ii].size(), [&](size_t ij) {
                                            xor_sums[ii][ij].first ^= curr->values[level][ii][ij].first;
                                            xor_sums[ii][ij].second ^= curr->values[level][ii][ij].second;
                                        });
                                });
                                //}
                            }
                        }
                    }
                }
            });

            parallel_for(0, splits->size(), [&](size_t i) {
                splits_ref[i]->update_level = UINT_E_MAX;
            });
            return results;
    }

    sequence<sequence<std::pair<uintE, uintE>>> get_subsequence_sum(SkipListElement* left, SkipListElement* right) {
            size_t level = 0;
            sequence<sequence<std::pair<uintE, uintE>>> xor_sums = right->values[level];

            while(left != right) {
                    level = std::min(left->height, right->height) - 1;
                    if (level == left->height-1) {
                        /*auto first_level_values = left->values[0];
                        if (level != 0 || (level == 0 && first_level_values.first != first_level_values.second)) {*/
                        parallel_for(0, xor_sums.size(), [&](size_t ii) {
                                parallel_for(0, xor_sums[ii].size(), [&](size_t ij) {
                                    xor_sums[ii][ij].first ^= left->values[level][ii][ij].first;
                                    xor_sums[ii][ij].second ^= left->values[level][ii][ij].second;
                            });
                        });
                        // }
                        left = left->elements[level].second;
                    } else {
                        right = right->elements[level].first;
                        /*auto first_level_values = right->values[0];
                        if (level != 0 || (level == 0 && first_level_values.first != first_level_values.second)) {*/
                        parallel_for(0, xor_sums.size(), [&](size_t ii) {
                                parallel_for(0, xor_sums[ii].size(), [&](size_t ij) {
                                    xor_sums[ii][ij].first ^= right->values[level][ii][ij].first;
                                    xor_sums[ii][ij].second ^= right->values[level][ii][ij].second;
                            });
                        });
                        //}
                    }
            }
            return xor_sums;
    }

    // Get the sum of the entire sequence
    sequence<sequence<std::pair<uintE, uintE>>> get_sum(SkipListElement* this_element) {
            SkipListElement* root = find_representative(this_element);
            size_t level = root->height-1;
            sequence<sequence<std::pair<uintE, uintE>>> xor_sums = root->values[level];

            SkipListElement* curr = root->elements[level].second;
            while (curr != nullptr && curr != root) {
                /*auto first_level_values = curr->values[0];
                if (level != 0 || (level == 0 && first_level_values.first != first_level_values.second)) {*/
                parallel_for(0, xor_sums.size(), [&](size_t ii) {
                        parallel_for(0, xor_sums[ii].size(), [&](size_t ij) {
                            xor_sums[ii][ij].first ^= curr->values[level][ii][ij].first;
                            xor_sums[ii][ij].second ^= curr->values[level][ii][ij].second;
                        });
                });
                //}
                curr = curr->elements[level].second;
            }

            if (curr == nullptr) { // the list is not circular
                    curr = root;
                    while(true) {
                            while(level > 0 && curr->elements[level].first == nullptr) {
                                    level--;
                            }

                            if(level == 0 && curr->elements[level].first == nullptr)
                                break;

                            while(curr->elements[level].first != nullptr) {
                                curr = curr->elements[level].first;
                                /*auto first_level_values = curr->values[0];
                                if (level != 0 || (level == 0
                                            && first_level_values.first != first_level_values.second)) {*/
                                parallel_for(0, xor_sums.size(), [&](size_t ii) {
                                        parallel_for(0, xor_sums[ii].size(), [&](size_t ij) {
                                            xor_sums[ii][ij].first ^= curr->values[level][ii][ij].first;
                                            xor_sums[ii][ij].second ^= curr->values[level][ii][ij].second;
                                        });
                                });
                                //}
                            }
                    }
            }
            return xor_sums;
    }
};

sequence<sequence<std::pair<uintE, uintE>>> default_values(uintE a, uintE b) {
        auto values_seq = sequence<sequence<std::pair<uintE, uintE>>>(2, sequence<std::pair<uintE, uintE>>(
                    ceil(log(10)/log(2)),
                    std::make_pair(a, b)));
        return values_seq;
}

inline void RunSkipList(uintE n) {
    std::cout << "Creating skip list" << std::endl;
    auto skip_list = SkipList(6);
    sequence<SkipList::SkipListElement> skip_list_elements = sequence<SkipList::SkipListElement>(10);

    std::cout << "creating nodes" << std::endl;
    auto curr_node = skip_list.create_node(2, nullptr, nullptr, default_values(2, 2));
    std::cout << "created first node" << std::endl;

    skip_list_elements[1] = curr_node;
    auto curr_node2 = skip_list.create_node(3, nullptr, nullptr, default_values(3, 3));
    skip_list_elements[2] = curr_node2;
    skip_list_elements[0] = skip_list.create_node(1, nullptr, nullptr, default_values(1, 1));
    skip_list_elements[3] = skip_list.create_node(4, nullptr, nullptr, default_values(4, 4));
    skip_list_elements[4] = skip_list.create_node(5, nullptr, nullptr, default_values(5, 5));
    skip_list_elements[5] = skip_list.create_node(6, nullptr, nullptr, default_values(6, 6));

    std::cout << "created nodes" << std::endl;

    std::cout << "joining nodes" << std::endl;
    sequence<std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>> join_updates
        = sequence<std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>>(5);
    join_updates[0] = std::make_pair(&skip_list_elements[1], &skip_list_elements[2]);
    join_updates[1] = std::make_pair(&skip_list_elements[0], &skip_list_elements[1]);
    join_updates[2] = std::make_pair(&skip_list_elements[2], &skip_list_elements[0]);
    join_updates[3] = std::make_pair(&skip_list_elements[3], &skip_list_elements[4]);
    join_updates[4] = std::make_pair(&skip_list_elements[4], &skip_list_elements[3]);
    skip_list.batch_join(&join_updates);

    std::cout << "printing answers" << std::endl;

    std::cout << "node 1 height: " << skip_list_elements[0].height << std::endl;
    std::cout << "node 1 value: " << skip_list_elements[0].values[0][0][0].first << ", "
        << skip_list_elements[0].values[0][0][0].second << std::endl;

    std::cout << "node 2 height: " << skip_list_elements[1].height << std::endl;
    std::cout << "node 2 value: " << skip_list_elements[1].values[0][0][0].first << ", " <<
       skip_list_elements[1].values[0][0][0].second << std::endl;

    std::cout << "node 3 height: " << skip_list_elements[2].height << std::endl;
    std::cout << "node 3 value: " << skip_list_elements[2].values[0][0][0].first << ", " <<
       skip_list_elements[2].values[0][0][0].second << std::endl;

    if (skip_list_elements[1].elements[0].first != nullptr)
        std::cout << "node 2 left: " <<
            skip_list_elements[1].elements[0].first -> values[0][0][0].first <<
            skip_list_elements[1].elements[0].first -> values[0][0][0].second << std::endl;
    if (skip_list_elements[1].elements[0].second != nullptr)
        std::cout << "node 2 right: " <<
            skip_list_elements[1].elements[0].second -> values[0][0][0].first <<
            skip_list_elements[1].elements[0].second -> values[0][0][0].second << std::endl;
    if (skip_list_elements[2].elements[0].first != nullptr)
        std::cout << "node 3 left: " <<
            skip_list_elements[2].elements[0].first -> values[0][0][0].first <<
            skip_list_elements[2].elements[0].first -> values[0][0][0].second << std::endl;
    if (skip_list_elements[2].elements[0].second != nullptr)
        std::cout << "node 3 right: " <<
            skip_list_elements[2].elements[0].second -> values[0][0][0].first <<
            skip_list_elements[2].elements[0].second -> values[0][0][0].second << std::endl;

    std::cout << "height values" << std::endl;

    auto node_1_height = skip_list_elements[0].height - 1;
    auto node_2_height = skip_list_elements[1].height - 1;
    auto node_3_height = skip_list_elements[2].height - 1;

    std::cout << "printing height values" << std::endl;

    if (skip_list_elements[0].elements[node_1_height].first != nullptr) {
       std::cout << "node 1 left pointer" << std::endl;
       auto left_neighbor = skip_list_elements[0].elements[node_1_height].first;
       std::cout << "node 1 height left: " <<
           skip_list_elements[0].elements[node_1_height].first -> values[left_neighbor->height - 1][0][0].first <<
           skip_list_elements[0].elements[node_1_height].first -> values[left_neighbor->height - 1][0][0].second
           << std::endl;
    }

    if (skip_list_elements[0].elements[node_1_height].second != nullptr) {
       auto right_neighbor = skip_list_elements[0].elements[node_1_height].second;
       std::cout << "node 1 height right: " <<
           skip_list_elements[0].elements[node_1_height].second -> values[right_neighbor->height - 1][0][0].first <<
           skip_list_elements[0].elements[node_1_height].second -> values[right_neighbor->height - 1][0][0].second
           << std::endl;
    }

    std::cout << "node 1 done" << std::endl;

    if (skip_list_elements[1].elements[node_2_height].first != nullptr) {
       auto left_neighbor = skip_list_elements[1].elements[node_2_height].first;
       std::cout << "node 2 height left: " <<
           skip_list_elements[1].elements[node_2_height].first -> values[left_neighbor->height - 1][0][0].first <<
           skip_list_elements[1].elements[node_2_height].first -> values[left_neighbor->height - 1][0][0].second
           << std::endl;
    }

    if (skip_list_elements[1].elements[node_2_height].second != nullptr) {
       auto right_neighbor = skip_list_elements[1].elements[node_2_height].second;
       std::cout << "node 2 height right: " <<
           skip_list_elements[1].elements[node_2_height].second -> values[right_neighbor->height - 1][0][0].first <<
           skip_list_elements[1].elements[node_2_height].second -> values[right_neighbor->height - 1][0][0].second
           << std::endl;
    }

    if (skip_list_elements[2].elements[node_3_height].first != nullptr) {
       auto left_neighbor = skip_list_elements[2].elements[node_3_height].first;
       std::cout << "node 3 height left: " <<
           skip_list_elements[2].elements[node_3_height].first -> values[left_neighbor->height - 1][0][0].first <<
           skip_list_elements[2].elements[node_3_height].first -> values[left_neighbor->height - 1][0][0].second
           << std::endl;
    }

    if (skip_list_elements[2].elements[node_3_height].second != nullptr) {
       auto right_neighbor = skip_list_elements[2].elements[node_3_height].second;
       std::cout << "node 3 height right: " <<
           skip_list_elements[2].elements[node_3_height].second -> values[right_neighbor->height - 1][0][0].first <<
           skip_list_elements[2].elements[node_3_height].second -> values[right_neighbor->height - 1][0][0].second
           << std::endl;
    }

    auto one_left_parent = skip_list.find_left_parent(0, &skip_list_elements[0]);
    if (one_left_parent != nullptr)
        std::cout << "node 1 left parent: " <<
            one_left_parent->values[0][0][0].first <<
            one_left_parent->values[0][0][0].second << std::endl;

    auto one_right_parent = skip_list.find_right_parent(0, &skip_list_elements[0]);
    if (one_right_parent != nullptr)
        std::cout << "node 1 right parent: " <<
            one_right_parent->values[0][0][0].first <<
            one_right_parent->values[0][0][0].second << std::endl;

    auto two_left_parent = skip_list.find_left_parent(0, &skip_list_elements[1]);
    if (two_left_parent != nullptr)
        std::cout << "node 2 left parent: " <<
            two_left_parent->values[0][0][0].first <<
            two_left_parent->values[0][0][0].second << std::endl;

    auto two_right_parent = skip_list.find_right_parent(0, &skip_list_elements[1]);
    if (two_right_parent != nullptr)
        std::cout << "node 2 right parent: " <<
            two_right_parent->values[0][0][0].first <<
            two_right_parent->values[0][0][0].second << std::endl;

    auto three_right_parent = skip_list.find_right_parent(0, &skip_list_elements[2]);
    if (three_right_parent != nullptr)
        std::cout << "node 3 right parent: " <<
            three_right_parent->values[0][0][0].first <<
            three_right_parent->values[0][0][0].second << std::endl;

    auto three_left_parent = skip_list.find_left_parent(0, &skip_list_elements[2]);
    if (three_left_parent != nullptr)
        std::cout << "node 3 left parent: " <<
            three_left_parent->values[0][0][0].first <<
            three_left_parent->values[0][0][0].second << std::endl;

    std::cout << "total sum 2, 3: " <<
        skip_list.get_subsequence_sum(&skip_list_elements[1],
            &skip_list_elements[2])[0][0].first << ", "
        << skip_list.get_subsequence_sum(&skip_list_elements[1],
            &skip_list_elements[2])[0][0].second << ", "
            << std::endl;
    std::cout << "total sum 1, 2: "
        << skip_list.get_subsequence_sum(&skip_list_elements[0],
            &skip_list_elements[1])[0][0].first
        << skip_list.get_subsequence_sum(&skip_list_elements[0],
            &skip_list_elements[1])[0][0].second << std::endl;
    std::cout << "total sum 1, 3: "
        << skip_list.get_subsequence_sum(&skip_list_elements[0],
            &skip_list_elements[2])[0][0].first
        << skip_list.get_subsequence_sum(&skip_list_elements[0],
            &skip_list_elements[2])[0][0].second << std::endl;
    std::cout << "total sum shouldn't work: "
        << skip_list.get_subsequence_sum(&skip_list_elements[2],
            &skip_list_elements[0])[0][0].first
        << skip_list.get_subsequence_sum(&skip_list_elements[2],
            &skip_list_elements[0])[0][0].second << std::endl;
    std::cout << "total sum shouldn't work 4, 5: "
        << skip_list.get_subsequence_sum(&skip_list_elements[3],
            &skip_list_elements[4])[0][0].first
        << skip_list.get_subsequence_sum(&skip_list_elements[3],
            &skip_list_elements[4])[0][0].second << std::endl;

    std::cout << "representative node 1: " << skip_list.find_representative(&skip_list_elements[0])->values[0][0][0].first
        << std::endl;
    std::cout << "representative node 2: " << skip_list.find_representative(&skip_list_elements[1])->values[0][0][0].first
        << std::endl;
    std::cout << "representative node 3: " << skip_list.find_representative(&skip_list_elements[2])->values[0][0][0].first
        << std::endl;
    std::cout << "representative node 4: " << skip_list.find_representative(&skip_list_elements[3])->values[0][0][0].first
        << std::endl;
    std::cout << "representative node 5: " << skip_list.find_representative(&skip_list_elements[4])->values[0][0][0].first
        << std::endl;
    std::cout << "representative node 6: " << skip_list.find_representative(&skip_list_elements[5])->values[0][0][0].first
        << std::endl;

    std::cout << "total sum subtree 1: "
        << skip_list.get_sum(&skip_list_elements[0])[0][0].first
        << skip_list.get_sum(&skip_list_elements[0])[0][0].second
        << "; "
        << skip_list.get_sum(&skip_list_elements[1])[0][0].first << ", "
        << skip_list.get_sum(&skip_list_elements[1])[0][0].second
        << "; "
        << skip_list.get_sum(&skip_list_elements[2])[0][0].first << ", "
        << skip_list.get_sum(&skip_list_elements[2])[0][0].second
        << std::endl;

    std::cout << "total sum subtree 2: "
        << skip_list.get_sum(&skip_list_elements[3])[0][0].first << ", "
        << skip_list.get_sum(&skip_list_elements[3])[0][0].second << "; "
        << skip_list.get_sum(&skip_list_elements[4])[0][0].first << ", "
        << skip_list.get_sum(&skip_list_elements[4])[0][0].second
        << std::endl;

    std::cout << "total sum subtree 3: "
        << skip_list.get_sum(&skip_list_elements[5])[0][0].first << ", "
        << skip_list.get_sum(&skip_list_elements[5])[0][0].second
        << std::endl;

    std::cout << "splitting some nodes" << std::endl;
    sequence<SkipList::SkipListElement*> splits
        = sequence<SkipList::SkipListElement*>(4);
    splits[0] = &skip_list_elements[2];
    splits[1] = &skip_list_elements[4];
    splits[2] = &skip_list_elements[5];
    splits[3] = &skip_list_elements[1];
    skip_list.batch_split(&splits);

    std::cout << "representative node 1: " << skip_list.find_representative(&skip_list_elements[0])->values[0][0][0].second
        << std::endl;
    std::cout << "representative node 2: " << skip_list.find_representative(&skip_list_elements[1])->values[0][0][0].second
        << std::endl;
    std::cout << "representative node 3: " << skip_list.find_representative(&skip_list_elements[2])->values[0][0][0].second
        << std::endl;
    std::cout << "representative node 4: " << skip_list.find_representative(&skip_list_elements[3])->values[0][0][0].second
        << std::endl;
    std::cout << "representative node 5: " << skip_list.find_representative(&skip_list_elements[4])->values[0][0][0].second
        << std::endl;
    std::cout << "representative node 6: " << skip_list.find_representative(&skip_list_elements[5])->values[0][0][0].second
        << std::endl;

    std::cout << "total sum subtree 1: " << skip_list.get_sum(&skip_list_elements[0])[0][0].first << ", "
        << skip_list.get_sum(&skip_list_elements[1])[0][0].first << std::endl;

    std::cout << "total sum subtree 2: "  << skip_list.get_sum(&skip_list_elements[2])[0][0].first
        << std::endl;

    std::cout << "total sum subtree 3: " << skip_list.get_sum(&skip_list_elements[3])[0][0].second << ", "
        << skip_list.get_sum(&skip_list_elements[4])[0][0].second
        << std::endl;

    std::cout << "total sum subtree 6: " << skip_list.get_sum(&skip_list_elements[5])[0][0].first << std::endl;

    std::cout << "total sum 1, 2: " << skip_list.get_subsequence_sum(&skip_list_elements[0],
            &skip_list_elements[1])[0][0].first << std::endl;
    std::cout << "total sum 4, 5: " << skip_list.get_subsequence_sum(&skip_list_elements[3],
            &skip_list_elements[4])[0][0].first << std::endl;

    std::cout << "joining more nodes" << std::endl;
    sequence<std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>> join_updates_1
        = sequence<std::pair<SkipList::SkipListElement*, SkipList::SkipListElement*>>(3);
    join_updates_1[0] = std::make_pair(&skip_list_elements[1], &skip_list_elements[2]);
    join_updates_1[1] = std::make_pair(&skip_list_elements[5], &skip_list_elements[0]);
    join_updates_1[2] = std::make_pair(&skip_list_elements[2], &skip_list_elements[5]);
    skip_list.batch_join(&join_updates_1);

    std::cout << "representative node 1: " << skip_list.find_representative(&skip_list_elements[0])->values[0][0][0].first
        << std::endl;
    std::cout << "representative node 2: " << skip_list.find_representative(&skip_list_elements[1])->values[0][0][0].first
        << std::endl;
    std::cout << "representative node 3: " << skip_list.find_representative(&skip_list_elements[2])->values[0][0][0].first
        << std::endl;
    std::cout << "representative node 4: " << skip_list.find_representative(&skip_list_elements[3])->values[0][0][0].first
        << std::endl;
    std::cout << "representative node 5: " << skip_list.find_representative(&skip_list_elements[4])->values[0][0][0].first
        << std::endl;
    std::cout << "representative node 6: " << skip_list.find_representative(&skip_list_elements[5])->values[0][0][0].first
        << std::endl;

    std::cout << "total sum subtree 1: " << skip_list.get_sum(&skip_list_elements[0])[0][0].first << ", "
        << skip_list.get_sum(&skip_list_elements[1])[0][0].first
        << ", " << skip_list.get_sum(&skip_list_elements[2])[0][0].first
        << ", " << skip_list.get_sum(&skip_list_elements[5])[0][0].first << std::endl;

    std::cout << "total sum subtree 2: " << skip_list.get_sum(&skip_list_elements[3])[0][0].first << ", "
        << skip_list.get_sum(&skip_list_elements[4])[0][0].first
        << std::endl;

}

}  // namespace gbbs
