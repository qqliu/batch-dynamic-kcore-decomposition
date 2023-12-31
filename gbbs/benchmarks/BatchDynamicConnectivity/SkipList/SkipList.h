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

    SkipListElement* find_left_parent(int level, SkipListElement* this_element) {
            SkipListElement* cur_element = this_element;
            SkipListElement* start_element = this_element;
            while (cur_element != nullptr && cur_element != start_element) {
                if (cur_element->height > level + 1)
                    return cur_element;
                cur_element->elements[level].first;
            }
            return nullptr;
    }

    SkipListElement* find_right_parent(int level, SkipListElement* this_element) {
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
    void update_top_down_sequential(size_t level, SkipListElement* this_element) {
            if (level == 0) {
                    if (this_element->height == 1) {
                        this_element->update_level = UINT_E_MAX;
                    }
                    return;
            }

            if (this_element->update_level < level) {
                    update_top_down_sequential(level - 1, this_element);
            }

            uintE xor_total = this_element->values[level-1];
            SkipListElement* curr = this_element->elements[level-1].second;
            while (curr != nullptr && curr->height < level + 1) {
                    if (curr->update_level != UINT_E_MAX && curr->update_level < level) {
                            update_top_down_sequential(level-1, curr);
                    }
                    xor_total ^= curr->values[level-1];
                    curr = curr->elements[level-1].second;
            }
            this_element->values[level] = xor_total;

            if(this_element->height == level+1) {
                    this_element->update_level = UINT_E_MAX;
            }
    }

    void update_top_down(size_t level, SkipListElement* this_element) {
            if (level <= 6) {
                    update_top_down_sequential(level, this_element);
                    return;
            }

            SkipListElement* curr = this_element;
            while(curr != nullptr && curr->height < level + 1) {
                    if (curr->update_level != UINT_E_MAX && curr->update_level < level) {
                            // can do this in parallel somehow; STOP HERE
                    }
            }
    }

};

template <class Graph, class W>
    inline void RunSkipList(Graph& G, BatchDynamicEdges<W> batch_edge_list) {
        size_t num_vertices = G.n;
        auto skip_list = SkipList(num_vertices);
}

}  // namespace gbbs
