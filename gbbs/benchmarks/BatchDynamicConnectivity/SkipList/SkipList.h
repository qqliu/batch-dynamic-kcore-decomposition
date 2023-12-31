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
        uintE val;

        using pointers = std::pair<SkipListElement*, SkipListElement*>;
        using height_array = parlay::sequence<pointers>;
        using values_array = parlay::sequence<uintE>;

        height_array elements;
        values_array values;
        uintE update_level;

        SkipListElement(): height(0), lowest_needs_update(0) { update_level = UINT_E_MAX; }

        SkipListElement(size_t _h, SkipListElement* _r, SkipListElement* _l, uintE _val):
            height(_h), lowest_needs_update(_h) {
                elements.resize(_h);
                values.resize(_h);
                update_level = UINT_E_MAX;
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

        inline bool CASleft(size_t level,
                SkipListElement* old_left, SkipListElement* new_left) {
               return pbbslib::atomic_compare_and_swap(&elements[level].first, old_left, new_left);
        }

        inline bool CASright(size_t level,
                SkipListElement* old_right, SkipListElement* new_right) {
               return pbbslib::atomic_compare_and_swap(&elements[level].second, old_right, new_right);
        }

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

    SkipListElement create_node(size_t index, SkipListElement* left, SkipListElement* right, uintE val) {
        rng.fork(index);
        auto rand_val = rng.rand() % UINT_E_MAX;

        size_t cur_height = 1;
        while (rand_val & 1) {
                rand_val >>= 1;
                cur_height++;
        }

        auto height = std::min(cur_height, (size_t) 32);

        auto node = SkipListElement(height, left, right, val);
        return node;
    }

    SkipListElement* find_left_parent(size_t level, SkipListElement* this_element) {
            SkipListElement* cur_element = this_element;
            SkipListElement* start_element = this_element;
            while (cur_element != nullptr && cur_element != start_element) {
                if (cur_element->height > level + 1)
                    return cur_element;
                cur_element->elements[level].first;
            }
            return nullptr;
    }

    SkipListElement* find_right_parent(size_t level, SkipListElement* this_element) {
            SkipListElement* cur_element = this_element;
            SkipListElement* start_element = this_element;
            while (cur_element != nullptr && cur_element != start_element) {
                if (cur_element->height > level + 1)
                    return cur_element;
                cur_element->elements[level].second;
            }
            return nullptr;
    }

    SkipListElement* find_representative(SkipListElement* this_element) {
        SkipListElement* cur_element = this_element;
        SkipListElement* seen_element = nullptr;

        size_t current_level = cur_element->height - 1;

        while(cur_element->elements[current_level].second != nullptr &&
                seen_element != cur_element) {
                if (seen_element == nullptr || cur_element < seen_element)
                    seen_element = cur_element;
                cur_element = cur_element->elements[current_level].second;
                auto top_level = cur_element->height-1;
                if (current_level < top_level) {
                        current_level = top_level;
                        seen_element = nullptr;
                }
        }

        if (seen_element == cur_element) {
                return seen_element;
        } else {
                while (cur_element->elements[current_level].first != nullptr) {
                        cur_element = cur_element->elements[current_level].first;
                        current_level = cur_element->height - 1;
                }
                return cur_element;
        }
    }

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
            SkipListElement* successor = nullptr;
            SkipListElement* cur_element = this_element;

            size_t level = 0;

            while(cur_element != nullptr) {
                SkipListElement* next = cur_element->elements[level].second;
                if (next != nullptr && cur_element->CASright(level, next, nullptr)) {
                        if (level == 0)
                            successor = next;
                        next->elements[level].first = nullptr;
                        cur_element = find_left_parent(level, cur_element);
                        level++;
                } else {
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

            uintE xor_total = this_element->values[level-1];
            SkipListElement* curr = this_element->elements[level-1].second;
            while (curr != nullptr && curr->height < level + 1) {
                    if (curr->update_level != UINT_E_MAX && curr->update_level < level) {
                            update_top_down(level-1, curr);
                    }
                    xor_total ^= curr->values[level-1];
                    curr = curr->elements[level-1].second;
            }
            this_element->values[level] = xor_total;

            if(this_element->height == level+1) {
                    this_element->update_level = UINT_E_MAX;
            }
    }

    void batch_update(sequence<std::pair<SkipListElement*, uintE>>* new_values) {
        auto top_nodes = sequence<SkipListElement*>(new_values->size(), nullptr);
        sequence<std::pair<SkipListElement*, uintE>>& new_values_ref = *new_values;
        if (new_values != nullptr) {
            parallel_for(0, new_values->size(), [&](size_t i){
                SkipListElement* this_element = new_values_ref[i].first;
                uintE this_element_value = new_values_ref[i].second;

                this_element->values[0] = this_element_value;

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
            auto join_lefts = sequence<std::pair<SkipListElement*, uintE>>(joins->size());
            parallel_for(0, joins->size(), [&] (size_t i){
                join(joins_ref[i].first, joins_ref[i].second);
                join_lefts[i] = std::make_pair(joins_ref[i].first, joins_ref[i].first->values[0]);
            });

            batch_update(&join_lefts);
    }

    void batch_split(sequence<SkipListElement*>* splits) {
            sequence<SkipListElement*>& splits_ref = *splits;
            parallel_for(0, splits->size(), [&](size_t i){
                split(splits_ref[i]);
            });

            // Perform updates but only if some other thread hasn't already performed the update
            parallel_for(0, splits->size(), [&](size_t i){
                SkipListElement* curr = splits_ref[i];
                bool can_proceed = curr->update_level == UINT_E_MAX
                    && pbbslib::atomic_compare_and_swap(&curr->update_level, UINT_E_MAX, (uintE)0);

                if (can_proceed) {
                    size_t xor_sum = curr->values[0];
                    size_t level = 0;
                    while(true) {
                        if(level < curr->height - 1) {
                            level++;
                            curr->values[level] = xor_sum;
                        } else {
                            curr = curr->elements[level].first;
                            if (curr == nullptr) {
                                break;
                            } else {
                                xor_sum ^= curr->values[level];
                            }
                        }
                    }
                }
            });

            parallel_for(0, splits->size(), [&](size_t i) {
                splits_ref[i]->update_level = UINT_E_MAX;
            });
    }

    uintE get_subsequence_sum(SkipListElement* left, SkipListElement* right) {
            size_t level = 0;
            size_t xor_sum = right->values[level];

            while(left != right) {
                    level = std::min(left->height, right->height) - 1;
                    if (level == left->height-1) {
                            xor_sum ^= left->values[level];
                            left = left->elements[level].second;
                    } else {
                            right = right->elements[level].first;
                            xor_sum ^= right->values[level];
                    }
            }
            return xor_sum;
    }

    // Get the sum of the entire sequence
    uintE get_sum(SkipListElement* this_element) {
            SkipListElement* root = find_representative(this_element);
            size_t level = root->height-1;
            uintE xor_sum = root->values[level];

            SkipListElement* curr = root->elements[level].second;
            while (curr != nullptr && curr != root) {
                    xor_sum ^= curr->values[level];
                    curr = curr->elements[level].second;
            }

            if (curr == nullptr) {
                // the list is not circular
                    curr = root;
                    while(true) {
                            while(level > 0 && curr->elements[level].first == nullptr) {
                                    level--;
                            }

                            if(level == 0 && curr->elements[level].first == nullptr)
                                break;

                            while(curr->elements[level].first != nullptr) {
                                curr = curr->elements[level].first;
                                xor_sum ^= curr->values[level];
                            }
                    }
            }
            return xor_sum;
    }
};

inline void RunSkipList(uintE n) {
    std::cout << "Creating skip list" << std::endl;
    auto skip_list = SkipList(n);
    std::cout << "Creating node" << std::endl;
    auto curr_node = skip_list.create_node(10, nullptr, nullptr, (uintE) 10);
    std::cout << "Node height" << std::endl;
    std::cout << curr_node.height << std::endl;
    std::cout << "Printing node value" << std::endl;
    std::cout << curr_node.values[0] << std::endl;
}

}  // namespace gbbs
