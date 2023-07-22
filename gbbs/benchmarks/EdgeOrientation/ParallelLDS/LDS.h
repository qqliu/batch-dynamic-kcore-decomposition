#include <unordered_set>
#include <stack>
#include <pthread.h>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"

#include "sparse_set.h"

namespace gbbs {

struct ground_truth_struct{
    sequence<sequence<uintE>> ground_truth;

    ground_truth_struct() {}

    ground_truth_struct(size_t batches, uintE n) {
        ground_truth = sequence<sequence<uintE>>(batches, sequence<uintE>(n, 0));
    }
};

struct LDS {

  double delta = 48.0;
  double UpperConstant = 2 + ((double) 3 / delta);
  double eps = 1.6;
  double OnePlusEps = 1 + eps;
  bool optimized_insertion = false;

  volatile uintE batch_num = 0;

  static constexpr uintE kUpLevel = UINT_E_MAX;
  // Used to indicate that this vertex is not moving.
  static constexpr uintE kNotMoving = UINT_E_MAX;

  using levelset = gbbs::sparse_set<uintE>;
  using down_neighbors = parlay::sequence<levelset>;
  using up_neighbors = levelset;
  using edge_type = std::pair<uintE, uintE>;

  // The descriptor contains the parent for our dependency DAG and the old level
  struct descriptor {
    uintE root;
    uintE old_level;

    descriptor(): root(UINT_E_MAX), old_level(UINT_E_MAX) {}

    descriptor(uintE r, uintE ol): root(r), old_level(ol) {}
  };

  struct LDSVertex {
    uintE level;         // The level of this vertex.
    uintE desire_level;  // The desire level of this vertex (set only when moving this vertex).
    down_neighbors down; // The neighbors in levels < level, bucketed by their level.
    up_neighbors up;     // The neighbors in levels >= level.

    LDSVertex() : level(0), desire_level(kNotMoving) {}

    // Used when Invariant 1 (upper invariant) is violated.
    template <class Levels>
    inline uintE get_desire_level_upwards(uintE vtx_id, Levels& L, const size_t levels_per_group,
            double upper_constant, double eps) const {
      using LV = std::pair<uintE, uintE>;
      auto LV_seq = parlay::delayed_seq<LV>(up.size(), [&] (size_t i) {
        uintE v = up.table[i];
        uintE l_v = UINT_E_MAX;
        if (levelset::valid(v)) l_v = L[v].level;
        return std::make_pair(l_v, v);
      });

      auto LVs = parlay::filter(LV_seq, [&] (const LV& lv) {
        return lv.first != UINT_E_MAX;
      });

      parlay::sort_inplace(parlay::make_slice(LVs));

      uintE new_level = UINT_E_MAX;
      uintE prev_level = level;
      for (size_t i=0; i<LVs.size(); i++) {
        // Start of a new level among our up-neighbors
        if ((i == 0) || (LVs[i].first != LVs[i-1].first)) {
          uintE cur_level = LVs[i].first;
          if (cur_level == level) continue;  // we are surely not moving to the same level

          size_t degree_at_level = LVs.size() - i;
          bool done = false;
          uintE prev_level_group = (prev_level+1)/levels_per_group;
          uintE cur_level_group = (cur_level)/levels_per_group;

          for (uintE j = prev_level_group; j <= cur_level_group; j++) {
                uintE up_degree = upper_constant * group_degree(j, eps);

                if (degree_at_level <= up_degree) {
                        new_level = std::max(uintE{j * levels_per_group}, uintE{prev_level + 1});
                        done = true;
                        break;
                }
          }
          if (done) break;

          prev_level = cur_level;
        }
      }

      if (new_level == UINT_E_MAX)
        new_level = prev_level + 1;


      return new_level;
    }

    // Used when Invariant 2 (up* degree invariant) is violated.
    template <class Levels>
    inline uintE get_desire_level_downwards(uintE vtx_id, Levels& L, const size_t levels_per_group,
            double upper_constant, double eps) const {
        assert(!lower_invariant(levels_per_group, eps));
        assert(level > 0);

        // This is the current level of the node.
        uintE cur_level = desire_level;
        if (desire_level == kNotMoving)
            cur_level = level - 1;
        if (cur_level == 0) return cur_level;
        uintE num_up_neighbors = num_neighbors_higher_than_level((uintE) cur_level);
        uintE num_up_star_neighbors = num_up_neighbors;
        while (cur_level > 0) {
            num_up_neighbors = num_up_star_neighbors;
            num_up_star_neighbors += down[cur_level-1].num_elms();
            if (upstar_degree_satisfies_invariant(cur_level, levels_per_group, num_up_star_neighbors,
                        num_up_neighbors, upper_constant, eps)) {
                break;
            }
            cur_level -= 1;
        }

        return cur_level;
    }

    inline uintE num_up_neighbors() const {
      return up.num_elms();
    }

    inline uintE num_up_star_neighbors() const {
      if (level == 0) return up.num_elms();
      return up.num_elms() + down[level - 1].num_elms();
    }

    // Get the sum of the number of neighbors at or higher than start_level
    inline uintE num_neighbors_higher_than_level(uintE start_level) const {
        uintE num_flipped_neighbors = 0;
        while (start_level < level) {
            num_flipped_neighbors += down[start_level].num_elms();
            start_level++;
        }
        return num_flipped_neighbors + up.num_elms();
    }

    template <class OutputSeq>
    inline uintE emit_up_neighbors(OutputSeq output, uintE our_id) const {
      size_t off = 0;
      for (size_t i=0; i<up.table_seq.size(); i++) {
        uintE v = up.table[i];
        if (levelset::valid(v))
          output[off++] = std::make_pair(v, our_id);
      }
      return off;
    }

    template <class OutputSeq>
    inline uintE emit_neighbors_higher_than_level(OutputSeq output, uintE our_id, uintE start_level) const {
      size_t off = 0;
      for (size_t j = start_level; j < level; j++) {
        for (size_t i=0; i < down[j].table_seq.size(); i++) {
            uintE v = down[j].table[i];
            if (levelset::valid(v))
                output[off++] = std::make_pair(v, our_id);
        }
      }

      for (size_t i = 0; i < up.table_seq.size(); i++) {
        uintE v = up.table[i];
        if (levelset::valid(v))
            output[off++] = std::make_pair(v, our_id);
      }
      return off;
    }


    template <class Levels>
    inline void filter_up_neighbors(uintE vtx_id, Levels& L) {
      uintE removed = 0;
      auto all_up = up.entries();
      assert(desire_level != kNotMoving);

      auto resize_sizes_seq = parlay::sequence<uintE>(desire_level - level, (uintE)0);
      auto resize_sizes = resize_sizes_seq.begin();

      for (size_t i=0; i<up.table_seq.size(); i++) {
        uintE v = up.table[i];
        if (levelset::valid(v)) {
          if (L[v].level == level && L[v].desire_level == kNotMoving) {
            up.table[i] = levelset::kTombstone;
            uintE norm_v = L[v].level - level;
            resize_sizes[norm_v]++;
            removed++;
          } else if (L[v].level > level && L[v].level < desire_level) {
            // Another case: L[v].level < desire_level
            up.table[i] = levelset::kTombstone;
            uintE norm_v = L[v].level - level;
            resize_sizes[norm_v]++;
            removed++;
          }
        }
      }

      // Perform resizing.
      up.resize_down(removed);
      assert(down.size() > level);
      for (size_t i=0; i<resize_sizes_seq.size(); i++) {
        auto l = level + i;
        down[l].resize(resize_sizes[i]);
      }

      uintE inserted = 0;
      for (size_t i=0; i<all_up.size(); i++) {
        uintE v = all_up[i];
        if (L[v].level == level && L[v].desire_level == kNotMoving) {
          down[level].insert(v);
          inserted++;
        } else if (L[v].level > level && L[v].level < desire_level) {
          down[L[v].level].insert(v);
          inserted++;
        }
      }
      assert(removed == inserted);
    }

    inline double group_degree(size_t group, double eps) const {
      return pow((1 + eps), group);
    }

    inline bool upper_invariant(const size_t levels_per_group, double upper_constant,
             double eps, bool optimized_insertion) const {
      uintE group = level / levels_per_group;
      uintE up_degree = 0;
      if (!optimized_insertion)
        up_degree = upper_constant * group_degree(group, eps);
      else
        up_degree = 1.1 * group_degree(group, eps);
      return up.num_elms() <= up_degree;
    }

    inline bool lower_invariant(const size_t levels_per_group, double eps) const {
      if (level == 0) return true;
      uintE lower_group = (level - 1) / levels_per_group;
      size_t our_group_degree = static_cast<size_t>(group_degree(lower_group, eps));
      return num_up_star_neighbors() >= our_group_degree;
    }

    // Method for computing whether the updegree from cur_level is enough to
    // satisfy Invariant 2.
    inline bool upstar_degree_satisfies_invariant(uintE cur_level, const size_t
            levels_per_group, uintE num_up_star_neighbors, uintE num_up_neighbors,
            double upper_constant, double eps) const {
        if (cur_level == 0) return true;
        uintE group = (cur_level - 1) / levels_per_group;
        size_t our_group_degree = static_cast<size_t>(group_degree(group, eps));
        uintE upper_group = cur_level / levels_per_group;
        size_t upper_group_degree = static_cast<size_t>(group_degree(upper_group, eps));

        return num_up_star_neighbors >= our_group_degree && num_up_neighbors
            < upper_constant * upper_group_degree;
    }

    inline bool is_dirty(const size_t levels_per_group, double upper_constant, double eps,
            bool optimized_insertion) const {
      bool upper = upper_invariant(levels_per_group, upper_constant, eps, optimized_insertion);
      bool lower = lower_invariant(levels_per_group, eps);
      return !(upper && lower);
    }

  };

  size_t n;  // number of vertices
  size_t levels_per_group;  // number of inner-levels per group,  O(\log n) many.
  parlay::sequence<LDSVertex> L_seq;
  LDSVertex* L;
  parlay::sequence<descriptor> descriptor_array; // keep a descriptor for every node in the graph

  LDS(size_t _n, bool _optimized_insertion, size_t _optimize_all) : n(_n),
    optimized_insertion(_optimized_insertion) {
    if (optimized_insertion)
        UpperConstant = 1.1;
    if (_optimize_all > 1)
        levels_per_group = ceil(log(n) / size_t(_optimize_all) * log(OnePlusEps));
    else
        levels_per_group = ceil(log(n) / log(OnePlusEps));
    L_seq = parlay::sequence<LDSVertex>(n);
    L = L_seq.begin();
    descriptor_array = parlay::sequence<descriptor>(_n);
  }

  LDS(size_t _n, double _eps, double _delta, bool _optimized_insertion,
          size_t _optimize_all) : n(_n), eps(_eps),
        delta(_delta), optimized_insertion(_optimized_insertion){
    OnePlusEps = (1 + _eps);
    if (optimized_insertion)
        UpperConstant = 1.1;
    else
        UpperConstant = 2 + ((double) 3 / _delta);
    if (_optimize_all > 1)
        levels_per_group = ceil(log(n) / size_t(_optimize_all) * log(OnePlusEps));
    else
        levels_per_group = ceil(log(n) / log(OnePlusEps));
    L_seq = parlay::sequence<LDSVertex>(_n);
    L = L_seq.begin();
    descriptor_array = parlay::sequence<descriptor>(_n);
  }

  uintE get_level(uintE ngh) {
    return L[ngh].level;
  }

  bool edge_exists(edge_type e) {
    auto[u, v] = e;
    auto l_u = L[u].level;
    auto l_v = L[v].level;
    if (l_u < l_v) {  // look in up(u)
      if (L[u].up.contains(v)) {
        assert(L[v].down[l_u].contains(u));
      }
      return L[u].up.contains(v);
    } else {  // look in up(v)
      if (L[v].up.contains(u)) {
          assert(L[u].up.contains(v) || L[u].down[l_v].contains(v));
      }
      return L[v].up.contains(u);
    }
  }

  // Find the root of the dependency DAG for a node
  inline uintE find_compress(uintE cur_node, descriptor* parents) {
    uintE j = cur_node;

    // Stop if the parent node is the same as the current node or if it is unmarked
    if (parents[j].root == j || parents[j].root == UINT_E_MAX)
        return parents[j].root;

    // Otherwise, continue looking for the root
    while (j != UINT_E_MAX && parents[j].root != j) {
        j = parents[j].root;
    }

    // Do path compression up to the root
    uintE tmp = parents[cur_node].root;
    while (cur_node != UINT_E_MAX && tmp != UINT_E_MAX && (tmp > j || j == UINT_E_MAX)) {
            parents[cur_node].root = j;
            cur_node = tmp;
            tmp = parents[cur_node].root;
    }

    if (tmp == UINT_E_MAX)
        j = UINT_E_MAX;
    return j;
  }

  // Merging two dependency trees
  inline bool unite_impl(uintE u_orig, uintE v_orig, descriptor* parents) {
        auto u = u_orig;
        auto v = v_orig;

        u = find_compress(u, parents);
        v = find_compress(v, parents);

        while (u != v && u != UINT_E_MAX && v != UINT_E_MAX) {
                u = find_compress(u, parents);
                v = find_compress(v, parents);

                if (u == UINT_E_MAX || v == UINT_E_MAX)
                    return false;

                if (u > v && parents[u].root == u
                    && pbbslib::atomic_compare_and_swap(&parents[u].root, u, v)) {
                        return true;
                } else if (v > u && parents[v].root == v
                    && pbbslib::atomic_compare_and_swap(&parents[v].root, v, u)) {
                        return true;
                }
        }

        return false;
  }

  // Invariant checking for an edge e that we expect to exist
  bool check_both_directions(edge_type e) {
    auto [u, v] = e;
    auto l_u = L[u].level;
    auto l_v = L[v].level;
    bool ok = true;
    if (l_u < l_v) {  // look in up(u)
      ok &= L[u].up.contains(v); assert(ok);
      ok &= L[v].down[l_u].contains(u); assert(ok);
    } else if (l_v < l_u) {  // look in up(v)
      ok &= L[v].up.contains(u); assert(ok);
      ok &= L[u].down[l_v].contains(v); assert(ok);
    } else {  // (l_v == l_u)
      ok &= L[v].up.contains(u); assert(ok);
      ok &= L[u].up.contains(v); assert(ok);
    }
    return ok;
  }

  // Input: dirty vertices
  // Output: Levels structure containing vertices at their desire-level
  template <class Seq, class Levels>
  void update_desire_levels(Seq&& possibly_dirty, Levels& levels,
          bool initialize_roots = false) {
    using level_and_vtx = std::pair<uintE, uintE>;

    // Compute dirty vertices that violate lower threshold
    auto level_and_vtx_seq = parlay::delayed_seq<level_and_vtx>(possibly_dirty.size(), [&] (size_t i){
        uintE v = possibly_dirty[i];
        uintE desire_level = UINT_E_MAX;
        assert(L[v].upper_invariant(levels_per_group, UpperConstant, eps, optimized_insertion));
        if (!L[v].lower_invariant(levels_per_group, eps)) {
            auto our_desire_level = L[v].get_desire_level_downwards(v, L, levels_per_group, UpperConstant,
                    eps);
            // Only add it if it's not in the level structure
            if (our_desire_level >= levels.size() || !levels[our_desire_level].contains(v)) {
                desire_level = our_desire_level;
            }
        }
        return std::make_pair(desire_level, v);
    });

    auto dirty = parlay::filter(level_and_vtx_seq, [&] (const level_and_vtx& lv) {
      return lv.first != UINT_E_MAX;
    });

    if (dirty.size() == 0) return;

    if  (initialize_roots) {
        parallel_for(0, dirty.size(), [&](size_t i) {
            if (descriptor_array[dirty[i].second].root == UINT_E_MAX)
                descriptor_array[dirty[i].second].root = dirty[i].second;
            descriptor_array[dirty[i].second].old_level = L[dirty[i].second].level;
        });
    }

    parallel_for(0, dirty.size(), [&] (size_t i){
        auto v = dirty[i].second;
        descriptor_array[v].old_level = L[v].level;
        if (descriptor_array[v].root == UINT_E_MAX)
            descriptor_array[v].root = v;
        else {
            for (size_t j = 0; j < L[v].down.size(); j++) {
                auto cur_neighbors = L[v].down[j].entries();
                for (size_t k = 0; k < cur_neighbors.size(); k++) {
                    if (cur_neighbors[k] != UINT_E_MAX) {
                        unite_impl((uintE) v, (uintE) cur_neighbors[k],
                                descriptor_array.begin());
                    }
                }
            }
        }
    });

    // Compute the removals from previous desire_levels
    auto previous_desire_levels = parlay::delayed_seq<level_and_vtx>(dirty.size(), [&] (size_t i){
        auto lv = dirty[i];
        uintE v = lv.second;
        uintE removal_level = UINT_E_MAX;
        uintE prev_desire_level = L[v].desire_level;
        if (prev_desire_level != lv.first && !(prev_desire_level >= levels.size())
                && levels[prev_desire_level].contains(v)) {
           removal_level = prev_desire_level;
           // Remove from the previous level
           levels[prev_desire_level].remove(v);
        }
        return std::make_pair(removal_level, v);
    });

    auto non_empty_levels = parlay::filter(previous_desire_levels, [&] (const level_and_vtx& old_level) {
        return old_level.first != UINT_E_MAX;
    });

    // Sort by level.
    auto get_key = [&] (const level_and_vtx& elm) { return elm.first; };
    parlay::integer_sort_inplace(parlay::make_slice(dirty), get_key);

    // Get unique levels.
    auto bool_seq = parlay::delayed_seq<bool>(dirty.size() + 1, [&] (size_t i) {
     if (i < dirty.size())
        assert(dirty[i].first != UINT_E_MAX);
      return (i == 0) || (i == dirty.size()) || (dirty[i].first != dirty[i-1].first);
    });
    auto level_starts = parlay::pack_index(bool_seq);
    assert(level_starts[level_starts.size() - 1] == dirty.size());
    assert(level_starts.size() >= 2);

    parlay::integer_sort_inplace(parlay::make_slice(non_empty_levels), get_key);
    auto bool_seq_2 = parlay::delayed_seq<bool>(non_empty_levels.size() + 1, [&] (size_t i) {
      return (i == 0) || (i == non_empty_levels.size()) ||
           (non_empty_levels[i].first != non_empty_levels[i-1].first);
    });
    auto prev_level_starts = parlay::pack_index(bool_seq_2);

    uintE max_current_level = dirty[level_starts[level_starts.size() - 2]].first + 1;
    if (levels.size() < max_current_level) {
      levels.resize(max_current_level);
    }

    // Resize the previous levels
    parallel_for(0, prev_level_starts.size() - 1, [&] (size_t i){
        uintE idx = prev_level_starts[i];
        uintE level = non_empty_levels[idx].first;
        uintE num_in_level = level_starts[i+1] - idx;

        if (level < levels.size())
            levels[level].resize_down(num_in_level);
    });

    parallel_for(0, level_starts.size() - 1, [&] (size_t i) {
      uintE idx = level_starts[i];
      uintE level = dirty[idx].first;
      uintE num_in_level = level_starts[i+1] - idx;

      auto stuff_in_level = parlay::delayed_seq<uintE>(num_in_level, [&] (size_t j) {
        L[dirty[idx+j].second].desire_level = dirty[idx + j].first;
        assert(L[dirty[idx+j].second].desire_level != UINT_E_MAX);
        return dirty[idx + j].second;
      });

      levels[level].append(stuff_in_level);
    });
  }

  // Input: sequence of vertex_ids
  template <class Seq, class Levels>
  void update_levels(Seq&& possibly_dirty, Levels& levels, sparse_set<uintE>& roots,
          bool initialize_roots = false) {
    using level_and_vtx = std::pair<uintE, uintE>;

    // Compute dirty vertices, which have either lower / upper threshold
    // violated. Output this vertex as dirty only if it is not already in the
    // bucketing structure.
    auto level_and_vtx_seq = parlay::delayed_seq<level_and_vtx>(possibly_dirty.size(), [&] (size_t i) {
      uintE v = possibly_dirty[i];
      uintE level = UINT_E_MAX;
      if (L[v].is_dirty(levels_per_group, UpperConstant, eps, optimized_insertion)) {
        auto our_level = L[v].level;
        // Return only if the vertex is not in the bucketing structure.
        if (our_level >= levels.size() || !levels[our_level].contains(v)) {
          level = our_level;
        }
      }
      return std::make_pair(level, v);
    });

    auto dirty = parlay::filter(level_and_vtx_seq, [&] (const level_and_vtx& lv) {
      return lv.first != UINT_E_MAX;
    });

    if  (initialize_roots) {
        parallel_for(0, dirty.size(), [&](size_t i) {
            if (descriptor_array[dirty[i].second].root == UINT_E_MAX)
                descriptor_array[dirty[i].second].root = dirty[i].second;
            descriptor_array[dirty[i].second].old_level = L[dirty[i].second].level;
        });
    }

    parallel_for(0, dirty.size(), [&] (size_t i) {
        auto v = dirty[i];
        descriptor_array[v].old_level = L[v].level;
        auto my_up_neighbors = L[v].up.entries();

        // do merging the dependency trees using
        // compare and swap for all up neighbors until it succeeds for every pair
        parallel_for(j = 0, my_up_neighbors.size(), [&] (size_t j) {
            unite_impl((uintE) v, (uintE) my_up_neighbors[j], descriptor_array.begin());
        });
    });

    if (dirty.size() == 0) return;

    // Sort by level.
    auto get_key = [&] (const level_and_vtx& elm) { return elm.first; };
    parlay::integer_sort_inplace(parlay::make_slice(dirty), get_key);

    // Get unique levels.
    auto bool_seq = parlay::delayed_seq<bool>(dirty.size() + 1, [&] (size_t i) {
      return (i == 0) || (i == dirty.size()) || (dirty[i].first != dirty[i-1].first);
    });
    auto level_starts = parlay::pack_index(bool_seq);
    assert(level_starts[level_starts.size() - 1] == dirty.size());
    assert(level_starts.size() >= 2);

    assert(level_starts.size() >= 2);
    uintE max_current_level = dirty[level_starts[level_starts.size() - 2]].first + 1;
    if (levels.size() < max_current_level) {
      levels.resize(max_current_level);
    }

    parallel_for(0, level_starts.size() - 1, [&] (size_t i) {
      uintE idx = level_starts[i];
      uintE level = dirty[idx].first;
      uintE num_in_level = level_starts[i+1] - idx;

      auto stuff_in_level = parlay::delayed_seq<uintE>(num_in_level, [&] (size_t j) {
        return dirty[idx + j].second;
      });

      levels[level].append(stuff_in_level);
    });
  }

  // Neighbors is a slice of sorted (level, neighbor_id) pairs.
  template <class Neighbors>
  void insert_neighbors(uintE vtx, Neighbors neighbors) {
    // Resize level containers. Can do this in parallel in many ways (e.g.,
    // pack), but a simple way that avoids memory allocation is to map over
    // everyone, and have the first index for each level search for the level.
    // If we use a linear search the algorithm is work eff. and has depth
    // proportional to the max incoming size of a level. We can also improve
    // to log(n) depth by using a doubling search.
    auto bool_seq = parlay::delayed_seq<bool>(neighbors.size() + 1, [&] (size_t i) {
      return (i == 0) || (i == neighbors.size()) || (std::get<0>(neighbors[i-1]) != std::get<0>(neighbors[i]));
    });
    auto starts = parlay::pack_index(bool_seq);

    parallel_for(0, starts.size() - 1, [&] (size_t i){
        auto index = starts[i];
        auto next_index = starts[i+1];

        uintE level_id = neighbors[index].first;
        uintE num_in_level = next_index - index;

        if (level_id != kUpLevel) {
          L[vtx].down[level_id].resize(num_in_level);
        } else {
          L[vtx].up.resize(num_in_level);
        }
    });

    parallel_for(0, neighbors.size(), [&] (size_t i) {
      auto [level_id, v] = neighbors[i];
      if (level_id != kUpLevel) {
        bool inserted = L[vtx].down[level_id].insert(v);
        assert(inserted);
      } else {
        bool inserted = L[vtx].up.insert(v);
        assert(inserted);
      }
    });
  }

  // Resize the tables
  // Can probably be optimized more
  template <class Neighbors>
  void delete_neighbors(uintE vtx, Neighbors neighbors) {
       // Delete neighbors from the adjacency list structures in parallel
       parallel_for(0, neighbors.size(), [&] (size_t i){
          auto [level_id, v] = neighbors[i];
          if (level_id != kUpLevel) {
              bool deleted = L[vtx].down[level_id].remove(v);
              assert(deleted);
          } else {
              bool deleted = L[vtx].up.remove(v);
              assert(deleted);
          }
      });

       auto bool_seq = parlay::delayed_seq<bool>(neighbors.size() + 1, [&] (size_t i) {
           return (i == 0) || (i == neighbors.size()) || (std::get<0>(neighbors[i-1]) != std::get<0>(neighbors[i]));
       });
       auto starts = parlay::pack_index(bool_seq);

       parallel_for(0, starts.size() - 1, [&] (size_t i){
           auto index = starts[i];
           auto next_index = starts[i+1];

           uintE level_id = neighbors[index].first;
           uintE num_in_level = next_index - index;

           if (level_id != kUpLevel) {
               L[vtx].down[level_id].resize_down(num_in_level);
           } else {
               L[vtx].up.resize_down(num_in_level);
           }
       });
  }

  // returns the total number of moved vertices
  template <class Levels>
  size_t rebalance_insertions(Levels&& levels, size_t actual_cur_level_id) {
    size_t total_moved = 0;
    size_t cur_level_id = actual_cur_level_id;
    if (cur_level_id >= levels.size())
      return total_moved;

    while (cur_level_id < levels.size()) {
        auto& cur_level = levels[cur_level_id];
        if (cur_level.num_elms() == 0) {
            cur_level_id++;
            continue;
        }

    // Figure out the desire_level for each vertex in cur_level.
    // Sets the desire level for each vertex with a valid desire_level in L.
    auto moved_vertices = parlay::sequence<size_t>(cur_level.size(), 0);
    parallel_for(0, cur_level.size(), [&] (size_t i) {
      uintE v = cur_level.table[i];
      uintE desire_level = UINT_E_MAX;
      if (levelset::valid(v) && L[v].is_dirty(levels_per_group, UpperConstant, eps, optimized_insertion)) {
        assert(L[v].lower_invariant(levels_per_group, eps));
        assert(!L[v].upper_invariant(levels_per_group, UpperConstant, eps, optimized_insertion));

        desire_level = L[v].get_desire_level_upwards(v, L, levels_per_group, UpperConstant, eps);
        L[v].desire_level = desire_level;
        moved_vertices[i] = 1;

        assert(L[v].level == cur_level_id);
        assert(desire_level > cur_level_id);
      }
    });
    total_moved += parlay::scan_inplace(parlay::make_slice(moved_vertices));

    size_t outer_level_sizes = 0;
    for (size_t i = 0; i < cur_level.size(); i++) {
      uintE u = cur_level.table[i];
      if (levelset::valid(u)) {
        if (L[u].desire_level == kNotMoving) {
          cur_level.remove(u);
          outer_level_sizes++;
        }
      }
    }
    cur_level.resize_down(outer_level_sizes);

    auto vertex_seq = parlay::delayed_seq<uintE>(cur_level.size(), [&] (size_t i) {
      uintE u = cur_level.table[i];
      if (levelset::valid(u)) {
        if (L[u].desire_level == kNotMoving) {
          u = UINT_E_MAX;  // filter out, otherwise leave in
        }
      }
      return u;
    });
    // Dirty now contains (dl(v), v) pairs for all v \in the current level
    // moving to dl(v) (their desire level).
    auto dirty_seq = parlay::filter(vertex_seq, [&] (const uintE& u) {
      return u != UINT_E_MAX;
    });
    auto dirty = dirty_seq.begin();

    parallel_for(0, dirty_seq.size(), [&] (size_t i) {
      uintE v = dirty[i];
      uintE desire_level = L[v].desire_level;
      L[v].down.resize(desire_level);
    });

    // Consider the neighbors N(u) of a vertex u moving from cur_level to dl(u) > cur_level.
    // - {v \in N(u) | dl(u) < l(v)}: move u from cur_level -> dl(u) in v's Down
    // - {v \in N(u) | cur_level < l(v) <= dl(u)}:
    //     remove u from v's Down and insert into v's Up
    // - {v \in N(u) | l(v) == cur_level}:
    //     if dl(v) <= dl(u): unaffected (u is already in v's Up)
    //     else (dl(u) < dl(v)): remove u from v's Up and insert into v's Down

    // Compute the number of neighbors we affect.
    auto degrees = parlay::map(parlay::make_slice(dirty_seq), [&] (auto v) {
      return L[v].num_up_neighbors();
    });
    size_t sum_degrees = parlay::scan_inplace(parlay::make_slice(degrees));

    // Write the affected (flipped) edges into an array of
    //   (affected_neighbor, moved_id)
    // edge pairs.
    auto flipped = sequence<edge_type>::uninitialized(sum_degrees);
    parallel_for(0, dirty_seq.size(), [&] (size_t i) {
      uintE v = dirty[i];
      size_t offset = degrees[i];
      size_t end_offset = (i == dirty_seq.size()-1) ? sum_degrees : degrees[i+1];

      auto output = flipped.cut(offset, end_offset);
      size_t written = L[v].emit_up_neighbors(output, v);
      assert(written == (end_offset - offset));

      // Remove neighbors in up with level < desire_level that are not moving.
      L[v].filter_up_neighbors(v, L);
    });

    // Sort based on the affected_neighbor. Note that there are no dup edges.
    auto compare_tup = [&] (const edge_type& l, const edge_type& r) { return l < r; };
    parlay::sort_inplace(parlay::make_slice(flipped), compare_tup);

    // Compute the starts of each (modified) vertex's new edges.
    auto bool_seq = parlay::delayed_seq<bool>(flipped.size() + 1, [&] (size_t i) {
      return (i == 0) || (i == flipped.size()) || (std::get<0>(flipped[i-1]) != std::get<0>(flipped[i]));
    });
    auto starts = parlay::pack_index(bool_seq);

    // Save the vertex ids (we will use this in update_levels).
    auto affected = sequence<uintE>::from_function(starts.size() - 1, [&] (size_t i) {
      size_t idx = starts[i];
      uintE vtx_id = flipped[idx].first;
      return vtx_id;
    });

    // Map over the vertices being modified. Overwrite their incident edges with
    // (level_id, neighbor_id) pairs, sort by level_id, resize each level to the
    // correct size, and then insert in parallel.
    auto num_flips = parlay::sequence<size_t>(starts.size(), (size_t) 0);
    parallel_for(0, starts.size() - 1, [&] (size_t i) {
      size_t idx = starts[i];
      uintE u = std::get<0>(flipped[idx]);
      uintE incoming_degree = starts[i+1] - starts[i];
      auto neighbors = parlay::make_slice(flipped.begin() + idx,
          flipped.begin() + idx + incoming_degree);

      uintE l_u = L[u].level;

      // (1) vtx (u) is a vertex moving from the current level.
      if (l_u == cur_level_id && L[u].desire_level != UINT_E_MAX) {
        // Get up neighbors.
        auto my_up_neighbors = L[u].up.entries();

        // parallel for loop for the merge using compare and swap with all dependency DAGs of
        // pairs of up neighbors
        parallel_for (0, my_up_neighbors.size(), [&] (size_t i) {
            unite_impl(u, my_up_neighbors[i], descriptor_array.begin());
        });

        descriptor_array[u].old_level = L[u].level;

        uintE dl_u = L[u].desire_level;
        assert(dl_u != UINT_E_MAX);
        assert(dl_u > l_u);

        // Map the incident edges to (level, neighbor_id).
        parallel_for(0, incoming_degree, [&] (size_t off) {
          auto [u_dup, v] = neighbors[off];
          assert(u == u_dup);
          uintE dl_v = L[v].desire_level;  // using desire_level not level.
          if (dl_v >= dl_u) { dl_v = kUpLevel; }  // stay in up
          neighbors[off] = {dl_v, v};  // send to dl_v
        });

        // Sort neighbors by level.
        auto compare_neighbors_tup = [&] (const edge_type& l, const edge_type& r) { return l < r; };
        parlay::internal::quicksort_serial(neighbors.begin(), neighbors.size(), compare_neighbors_tup);

        auto level_seq = parlay::delayed_seq<uintE>(neighbors.size(), [&] (size_t i) {
          return neighbors[i].first;
        });
        // first key greater than or equal to kUpLevel
        size_t upstart = parlay::internal::binary_search(level_seq, kUpLevel, std::less<uintE>());

        if (upstart > 0) {  // stuff to delete from L[u].up
          parallel_for(0, upstart, [&] (size_t j) {
            uintE v = neighbors[j].second;
            bool removed = L[u].up.remove(v);
            assert(removed);
          });
          L[u].up.resize_down(upstart);

          // No need to update stuff in [upstart, end), since they are already
          // in u's Up.
          auto lower_neighbors = neighbors.cut(0, upstart);

          // Insert the neighbors to their new locations.
          insert_neighbors(u, lower_neighbors);
        }
        num_flips[i] = upstart;
        //std::cout << "finished all of 2" << std::endl;
      } else if (l_u > cur_level_id) {

        // Map the incident edges to (level, neighbor_id).
        parallel_for(0, incoming_degree, [&] (size_t off) {
          auto [u_dup, v] = neighbors[off];
          assert(u == u_dup);
          uintE dl_v = L[v].desire_level;  // using desire_level not level.
          if (dl_v >= l_u) { dl_v = kUpLevel; }  // move to up
          neighbors[off] = {dl_v, v};  // send to dl_v
        });

        // Sort neighbors by level.
        auto compare_neighbors_tup = [&] (const edge_type& l, const edge_type& r) { return l < r; };
        parlay::internal::quicksort_serial(neighbors.begin(), neighbors.size(), compare_neighbors_tup);

        parallel_for(0, neighbors.size(), [&] (size_t i) {
          uintE v = neighbors[i].second;
          bool removed = L[u].down[cur_level_id].remove(v);
          assert(removed);
        });
        // Every updated edge is removed from cur_level.
        L[u].down[cur_level_id].resize_down(incoming_degree);

        // Insert the neighbors to their new locations.
        insert_neighbors(u, neighbors);
        num_flips[i] = neighbors.size();
      } else {
        assert(l_u == cur_level_id && L[u].desire_level == UINT_E_MAX);
        // Don't have to do anything for these vertices. They are staying put at
        // l_u, but their neighbors (already in u's Up) are moving to higher
        // levels (and thus staying in u's Up).
      }
    });

    //total_moved += parlay::scan_inplace(parlay::make_slice(num_flips));
    // Update current level for the dirty vertices, and reset
    // the desire_level.
    parallel_for(0, dirty_seq.size(), [&] (size_t i) {
      uintE v = dirty[i];
      L[v].level = L[v].desire_level;
      L[v].desire_level = UINT_E_MAX;
    });

    sparse_set<uintE> empty_roots = sparse_set<uintE>();
    update_levels(std::move(affected), levels, empty_roots);

    cur_level_id++;
    }

    return total_moved;
  }

  // Used to balance the level data structure for deletions
  // Returns the total number of moved vertices
  template <class Levels>
  size_t rebalance_deletions(Levels&& levels, size_t actual_cur_level_id, size_t total_moved = 0) {
    size_t cur_level_id = actual_cur_level_id;
    if (cur_level_id >= levels.size()) {
      return total_moved;
    }

    while (cur_level_id < levels.size()) {
      // Get vertices that have desire level equal to the current desire level

      // For deletion, we need to figure out all vertices that want to move to
      // the current level. This is maintained in the levels data structure.
      if (levels[cur_level_id].size() == 0) {
        cur_level_id++;
        continue;
      }
      auto nodes_to_move = levels[cur_level_id];

      // Turn nodes_to_move into a sequence
      auto nodes_to_move_seq = nodes_to_move.entries();

      // Compute the number of neighbors we affect
      auto degrees = parlay::map(parlay::make_slice(nodes_to_move_seq), [&] (auto v) {
          auto desire_level = L[v].desire_level;
          assert(desire_level < L[v].level);
          return L[v].num_neighbors_higher_than_level(desire_level);
      });
      size_t sum_degrees = parlay::scan_inplace(parlay::make_slice(degrees));

      // Write the affected neighbors into an array of (affected_neighbor,
      // moved_id) edge pairs.
      auto flipped = sequence<edge_type>::uninitialized(sum_degrees);
      parallel_for(0, nodes_to_move_seq.size(), [&] (size_t i){
        uintE v = nodes_to_move_seq[i];
        size_t offset = degrees[i];
        size_t end_offset = ((i == nodes_to_move_seq.size() - 1) ? sum_degrees : degrees[i + 1]);

        auto output = flipped.cut(offset, end_offset);
        size_t written = L[v].emit_neighbors_higher_than_level(output, v, L[v].desire_level);
        assert(written == (end_offset - offset));
      });

      // Sort by neighbor vertex index
      auto compare_tup = [&] (const edge_type& l, const edge_type& r) {return l < r;};
      parlay::sort_inplace(parlay::make_slice(flipped), compare_tup);

      // Compute the starts of the neighbor vertex indices
      auto bool_seq = parlay::delayed_seq<bool>(flipped.size() + 1, [&] (size_t i) {
        return (i == 0) || (i == flipped.size()) ||
            (std::get<0>(flipped[i-1]) != std::get<0>(flipped[i]));
      });
      auto starts = parlay::pack_index(bool_seq);

      // Save the vertex ids of the vertices which did not move
      auto affected = sequence<uintE>::from_function(starts.size() - 1, [&] (size_t i) {
        size_t idx = starts[i];
        uintE vtx_id = flipped[idx].first;
        return vtx_id;
      });

      // First, update the down-levels of vertices that do not move (not in
      // nodes_to_move). Then, all vertices which moved to the same level should
      // be both in each other's up adjacency lists.
      parallel_for(0, starts.size() - 1, [&] (size_t i) {
        size_t idx = starts[i];
        uintE u = std::get<0>(flipped[idx]);
        if (!nodes_to_move.contains(u) && L[u].level > cur_level_id) {
            uintE incoming_degree = starts[i+1] - starts[i];
            auto neighbors = parlay::make_slice(flipped.begin() + idx,
                    flipped.begin() + idx + incoming_degree);

            // All vertices in neighbors are moving to level cur_level
            // Update down[cur_level] to contain these vertices
            assert(cur_level_id <= L[u].level);
            L[u].down[cur_level_id].resize(neighbors.size());
            parallel_for(0, neighbors.size(), [&] (size_t q) {
                L[u].down[cur_level_id].insert(neighbors[q].second);
            });

            // Remove the vertices that moved from their previous levels in
            // L[u].down.
            //
            auto my_level = L[u].level;
            auto neighbor_levels = sequence<size_t>::uninitialized(neighbors.size());
            for (size_t j = 0; j < neighbors.size(); j++) {
                auto neighbor_id = neighbors[j].second;
                auto neighbor = L[neighbor_id];
                auto neighbor_level = neighbor.level;
                assert(neighbor_level >= cur_level_id);

                if (neighbor_level < my_level) {
                    assert(L[u].down[neighbor_level].contains(neighbor_id));
                    L[u].down[neighbor_level].remove(neighbor_id);
                    neighbor_levels[j] = neighbor_level;
                } else {
                    assert(L[u].up.contains(neighbor_id));
                    L[u].up.remove(neighbor_id);
                    neighbor_levels[j] = my_level;
                }
            }

            // Get the num deleted by sorting the levels of all neighbors
            auto compare_tup = [&] (const size_t l, const size_t r) { return l < r; };
            parlay::integer_sort_inplace(parlay::make_slice(neighbor_levels));
            auto new_bool_seq = parlay::delayed_seq<bool>(neighbor_levels.size() + 1, [&] (size_t i) {
                      return (i == 0) || (i == neighbor_levels.size()) ||
                             (neighbor_levels[i-1] != neighbor_levels[i]);
                 });
            auto new_starts = parlay::pack_index(new_bool_seq);

            for (size_t j = 0; j < new_starts.size() - 1; j++) {
                size_t idx = new_starts[j];
                size_t n_level = neighbor_levels[idx];
                size_t num_deleted = new_starts[j+1] - idx;
                if (n_level == my_level)
                    L[u].up.resize_down(num_deleted);
                else
                    L[u].down[n_level].resize_down(num_deleted);
            }
        }
      });

      // Update descriptors of every node in nodes_to_move
      parallel_for(0, nodes_to_move.size(), [&] (size_t cur_node_to_move){
        descriptor_array[cur_node_to_move].old_level = L[cur_node_to_move].level;
        if (descriptor_array[cur_node_to_move].root == UINT_E_MAX)
            descriptor_array[cur_node_to_move].root = cur_node_to_move;
        else {
            for (size_t down_level = 0; down_level < L[cur_node_to_move].down.size(); down_level++) {
                auto cur_neighbors = L[cur_node_to_move].down[down_level].entries();
                for (size_t neighbor_index = 0; neighbor_index < cur_neighbors.size(); neighbor_index++) {
                    if (cur_neighbors[neighbor_index] != UINT_E_MAX) {
                        unite_impl((uintE) cur_node_to_move,
                                (uintE) cur_neighbors[neighbor_index], descriptor_array.begin());
                    }
                }
            }
        }
      });

      // Move vertices in nodes_to_move to cur_level. Update the data structures
      // of each moved vertex and neighbors in flipped.
      //
      // Re-sort the flipped edges by endpoint which moved
      auto flipped_reverse = sequence<edge_type>::uninitialized(sum_degrees);
      parallel_for(0, flipped_reverse.size(), [&] (size_t i){
        flipped_reverse[i] = std::make_pair(flipped[i].second, flipped[i].first);
      });
      auto compare_flipped = [&] (const edge_type& l, const edge_type& r) { return l < r; };
      parlay::sort_inplace(parlay::make_slice(flipped_reverse), compare_flipped);

      // Compute the starts of the reverse flipped edges
      auto bool_seq_reverse = parlay::delayed_seq<bool>(flipped_reverse.size() + 1, [&] (size_t i){
        return (i == 0) || (i == flipped_reverse.size()) || (std::get<0>(flipped_reverse[i-1])
                != std::get<0>(flipped_reverse[i]));
      });
      auto reverse_starts = parlay::pack_index(bool_seq_reverse);

      // Update the data structures (vertices kept at each level) of each vertex
      // that moved.
      parallel_for(0, reverse_starts.size() - 1, [&] (size_t i) {
        size_t idx = reverse_starts[i];
        uintE moved_vertex_v = std::get<0>(flipped_reverse[idx]);
        size_t idx_plus = reverse_starts[i+1];
        uintE next_v = (idx_plus == flipped_reverse.size()) ? UINT_E_MAX : std::get<0>(flipped_reverse[idx_plus]);
        assert(moved_vertex_v < next_v);
        uintE num_neighbors = reverse_starts[i+1] - reverse_starts[i];

        auto neighbors = parlay::make_slice(flipped_reverse.begin() + idx,
                flipped_reverse.begin() + idx + num_neighbors);

        auto my_level = L[moved_vertex_v].level;
        size_t indegree_sum = 0;
        for (size_t j = 0; j < neighbors.size(); j++) {
            assert(moved_vertex_v == neighbors[j].first);

            if (L[neighbors[j].second].level < my_level) {
                indegree_sum++;
            }
        }

        L[moved_vertex_v].up.resize(indegree_sum);
        for (size_t k = 0; k < neighbors.size(); k++) {
            if (L[neighbors[k].second].level < my_level) {
                L[moved_vertex_v].up.insert(neighbors[k].second);
            }
        }
      });

      // Update current level of the moved vertices and reset the desired level
      parallel_for(0, nodes_to_move_seq.size(), [&] (size_t i){
        uintE v = nodes_to_move_seq[i];
        L[v].level = cur_level_id;
        L[v].desire_level = kNotMoving;
        L[v].down.resize(cur_level_id);
        assert(L[v].upper_invariant(levels_per_group, UpperConstant, eps, optimized_insertion));
        assert(L[v].lower_invariant(levels_per_group, eps));
      });

      size_t size_down = levels[cur_level_id].num_elms();
      levels[cur_level_id].clear();
      levels[cur_level_id].resize_down(size_down);
      // update the levels with neighbors
      update_desire_levels(std::move(affected), levels);

      total_moved += nodes_to_move_seq.size();
      cur_level_id++;
    }
    return total_moved;
  }

  size_t get_size() {
    size_t size = 0;
    size += sizeof(delta) + sizeof(UpperConstant) + sizeof(eps) + sizeof(OnePlusEps) +
        sizeof(optimized_insertion) + sizeof(n) + sizeof(levels_per_group);
    for (size_t i = 0; i < n; i++) {
        auto vertex = L[i];
        size += sizeof(vertex.level);
        size += sizeof(vertex.desire_level);
        for (size_t j = 0; j < vertex.down.size(); j++) {
            for (size_t k = 0; k < vertex.down[j].size(); k++) {
                size += sizeof(vertex.down[j].table[k]);
            }
        }

        for (size_t j = 0; j < vertex.up.size(); j++) {
            size += sizeof(vertex.up.table[j]);
        }
    }
    return size;
  }

  template <class Seq>
  size_t batch_insertion(const Seq& insertions_unfiltered) {
    // Remove edges that already exist from the input.
    auto insertions_filtered = parlay::filter(parlay::make_slice(insertions_unfiltered),
        [&] (const edge_type& e) { return !edge_exists(e); });

    // Duplicate the edges in both directions and sort.
    auto insertions_dup = sequence<edge_type>::uninitialized(2*insertions_filtered.size());
    parallel_for(0, insertions_filtered.size(), [&] (size_t i) {
      auto [u, v] = insertions_filtered[i];
      insertions_dup[2*i] = {u, v};
      insertions_dup[2*i + 1] = {v, u};
    });
    auto compare_tup = [&] (const edge_type& l, const edge_type& r) { return l < r; };
    parlay::sort_inplace(parlay::make_slice(insertions_dup), compare_tup);

    // Remove duplicate edges to get the insertions.
    auto not_dup_seq = parlay::delayed_seq<bool>(insertions_dup.size(), [&] (size_t i) {
      auto [u, v] = insertions_dup[i];
      bool not_self_loop = u != v;
      return not_self_loop && ((i == 0) || (insertions_dup[i] != insertions_dup[i-1]));
    });
    auto insertions = parlay::pack(parlay::make_slice(insertions_dup), not_dup_seq);


    // Compute the starts of each (modified) vertex's new edges.
    auto bool_seq = parlay::delayed_seq<bool>(insertions.size() + 1, [&] (size_t i) {
      return (i == 0) || (i == insertions.size()) || (std::get<0>(insertions[i-1]) != std::get<0>(insertions[i]));
    });
    auto starts = parlay::pack_index(bool_seq);

    // Save the vertex ids. The next step will overwrite the edge pairs to store
    // (neighbor, current_level).  (saving is not necessary if we modify + sort
    // in a single parallel loop)
    auto affected = sequence<uintE>::from_function(starts.size() - 1, [&] (size_t i) {
      size_t idx = starts[i];
      uintE vtx_id = std::get<0>(insertions[idx]);
      return vtx_id;
    });

    // Map over the vertices being modified. Overwrite their incident edges with
    // (level_id, neighbor_id) pairs, sort by level_id, resize each level to the
    // correct size, and then insert in parallel.
    parallel_for(0, starts.size() - 1, [&] (size_t i) {
      size_t idx = starts[i];
      uintE vtx = std::get<0>(insertions[idx]);
      uintE our_level = L[vtx].level;
      uintE incoming_degree = starts[i+1] - starts[i];
      auto neighbors = parlay::make_slice(insertions.begin() + idx,
          insertions.begin() + idx + incoming_degree);

      // Map the incident edges to (level, neighbor_id).
      for (size_t off = 0; off < incoming_degree; off++) {
        auto [u, v] = neighbors[off];
        assert(vtx == u);
        uintE neighbor_level = L[v].level;
        if (neighbor_level >= our_level) { neighbor_level = kUpLevel; }
        neighbors[off] = {neighbor_level, v};
      }

      // Sort neighbors by level.
      auto compare_neighbors_tup = [&] (const edge_type& l, const edge_type& r) { return l < r; };
      parlay::internal::quicksort_serial(neighbors.begin(), neighbors.size(), compare_neighbors_tup);

      // Insert moving neighbors.
      insert_neighbors(vtx, neighbors);

    }, 1);

    // New edges are done being deleted. Update the level structure.
    // Interface: supply vertex seq -> process will settle everything.
    using dirty_elts = sparse_set<uintE>;
    sequence<dirty_elts> levels;
    dirty_elts root_set = sparse_set<uintE>();

    // Place the affected vertices into levels based on their current level.
    update_levels(std::move(affected), levels, root_set, 1);

    // Update the level structure (basically a sparse bucketing structure).
    size_t total_moved = rebalance_insertions(std::move(levels), 0);

    parallel_for(0, n, [&] (size_t i) {
        if (descriptor_array[i].root == i) {
            descriptor_array[i].root = UINT_E_MAX;
            descriptor_array[i].old_level = UINT_E_MAX;
        }
    });

    parallel_for(0, n, [&] (size_t i) {
        descriptor_array[i].root = UINT_E_MAX;
        descriptor_array[i].old_level = UINT_E_MAX;
    });
    return total_moved;
  }

  template <class Seq>
  size_t batch_deletion(const Seq& deletions_unfiltered) {
    // Remove edges that do not exist in the graph.
    auto deletions_filtered = parlay::filter(parlay::make_slice(deletions_unfiltered),
            [&] (const edge_type& e) { return edge_exists(e); });


    // Duplicate the edges in both directions and sort.
    auto deletions_dup = sequence<edge_type>::uninitialized(2*deletions_filtered.size());
    parallel_for(0, deletions_filtered.size(), [&] (size_t i) {
        auto [u, v] = deletions_filtered[i];
        deletions_dup[2*i] = {u, v};
        deletions_dup[2*i + 1] = {v, u};
    });
    auto compare_tup = [&] (const edge_type& l, const edge_type& r) {return l < r;};
    parlay::sort_inplace(parlay::make_slice(deletions_dup), compare_tup);

    // Remove duplicate deletions.
    auto not_dup_seq = parlay::delayed_seq<bool>(deletions_dup.size(), [&](size_t i){
        auto [u, v] = deletions_dup[i];
        bool not_self_loop = u != v;
        return not_self_loop && ((i == 0) || (deletions_dup[i] != deletions_dup[i-1]));
    });
    auto deletions = parlay::pack(parlay::make_slice(deletions_dup), not_dup_seq);

    // Compute the starts of each (modified) vertex's deleted edges.
    auto bool_seq = parlay::delayed_seq<bool>(deletions.size() + 1, [&] (size_t i) {
      return (i == 0) || (i == deletions.size()) ||
        (std::get<0>(deletions[i-1]) != std::get<0>(deletions[i]));
    });
    auto starts = parlay::pack_index(bool_seq);

    // Save the vertex ids. The next step will overwrite the edge pairs to store
    // (neighbor, current_level).  (saving is not necessary if we modify + sort
    // in a single parallel loop)
    auto affected = sequence<uintE>::from_function(starts.size() - 1, [&] (size_t i) {
      size_t idx = starts[i];
      uintE vtx_id = std::get<0>(deletions[idx]);
      return vtx_id;
    });

    // Map over the vertices being modified. Overwrite their incident edges with
    // (level_id, neighbor_id) pairs, sort by level_id, resize each level to the
    // correct size, and then insert in parallel.
    parallel_for(0, starts.size() - 1, [&] (size_t i) {
      size_t idx = starts[i];
      uintE vtx = std::get<0>(deletions[idx]);
      uintE our_level = L[vtx].level;

      // Number of edges deleted that are adjacent to vtx
      uintE outgoing_degree = starts[i+1] - starts[i];

      // Neighbors incident to deleted edges
      auto neighbors = parlay::make_slice(deletions.begin() + idx,
          deletions.begin() + idx + outgoing_degree);

      // Map the incident edges to (level, neighbor_id).
      //parallel_for(0, outgoing_degree, [&] (size_t off) {
      for (size_t off = 0; off < outgoing_degree; off++) {
        auto [u, v] = neighbors[off];
        assert(vtx == u);
        uintE neighbor_level = L[v].level;
        if (neighbor_level >= our_level) { neighbor_level = kUpLevel; }
        neighbors[off] = {neighbor_level, v};
      }

      // Sort neighbors by level.
      parlay::sort_inplace(neighbors);

      // Delete moving neighbors.
      delete_neighbors(vtx, neighbors);

    }, 1);

    // New edges are done being inserted. Update the level structure.
    // Interface: supply vertex seq -> process will settle everything.

    using dirty_elts = sparse_set<uintE>;

    // Maintains the dirty vertices by their level
    sequence<dirty_elts> levels;

    // Place the affected vertices into levels based on their current level.
    // For deletion, need to move only vertices that are going to the lowest
    // level and then re-update the dirty vertices and repeat.
    // First, start by finding all the dirty vertices via the below method that
    // are adjacent to an edge deletion.
    update_desire_levels(std::move(affected), levels, true);

    // Update the level structure (basically a sparse bucketing structure).
    size_t total_moved = rebalance_deletions(std::move(levels), 0);

    return total_moved;
  }

  void check_invariants() {
    bool invs_ok = true;
    for (size_t i=0; i<n; i++) {
      bool upper_ok = L[i].upper_invariant(levels_per_group, UpperConstant, eps, optimized_insertion);
      bool lower_ok = L[i].lower_invariant(levels_per_group, eps);
      assert(upper_ok);
      assert(lower_ok);
      invs_ok &= upper_ok;
      invs_ok &= lower_ok;
    }
    assert(invs_ok);
  }

  inline uintE group_for_level(uintE level) const {
    return level / levels_per_group;
  }

  uintE core(uintE v) const {
    auto l = L[v].level;
    uintE group = group_for_level(l);
    if (l % levels_per_group != levels_per_group - 1 && group != 0) group--;
    return ceil(L[v].group_degree(group, eps));
  }

  uintE get_core_from_level(uintE l) const {
    uintE group = group_for_level(l);
    if (l % levels_per_group != levels_per_group - 1 && group != 0) group--;
    return ceil(pow((1 + eps), group));
  }

  uintE max_degree() const {
    auto outdegrees = parlay::delayed_seq<uintE>(n, [&] (size_t i) {
        return L[i].up.size();
    });
    uintE max_degree = pbbslib::reduce_max(outdegrees);
    return max_degree;
  }

  uintE max_coreness() const {
    auto levels = parlay::delayed_seq<uintE>(n, [&] (size_t i) {
        return L[i].level;
    });
    uintE max_level = pbbslib::reduce_max(levels);
    uintE max_group = group_for_level(max_level);
    return ceil(L[0].group_degree(max_group, eps));
  }
};

struct barrier_t {
  pthread_cond_t complete;
  pthread_mutex_t mutex;
  int count;
  int crossing;
};

static inline void barrier_init(barrier_t *b, int n) {
    pthread_cond_init(&b->complete, NULL);
    pthread_mutex_init(&b->mutex, NULL);
    b->count = n;
    b->crossing = 0;
}

static inline void barrier_cross(barrier_t *b) {
    pthread_mutex_lock(&b->mutex);
    /* One more thread through */
    b->crossing++;
    /* If not all here, wait */
    if (b->crossing < b->count) {
        pthread_cond_wait(&b->complete, &b->mutex);
    } else {
        pthread_cond_broadcast(&b->complete);
        /* Reset for next time */
        b->crossing = 0;
    }
    pthread_mutex_unlock(&b->mutex);
}


template <class W>
inline void RunLDS (BatchDynamicEdges<W>& batch_edge_list, long batch_size, bool compare_exact, LDS& layers, bool optimized_insertion, size_t offset, bool get_size,
        size_t num_reader_threads, bool nonlinearizable = false) {
    auto batch = batch_edge_list.edges;
    size_t num_insertion_flips = 0;
    size_t num_deletion_flips = 0;
    size_t max_degree = 0;
    auto ground_truth_container = ground_truth_struct(
            (size_t) ceil(batch.size()/(batch_size * 1.0)), layers.n);

    if (compare_exact) {
        for (size_t batch_i = 0; batch_i < batch.size(); batch_i += batch_size) {
            auto graph = dynamic_edge_list_to_symmetric_graph(batch_edge_list,
                    std::min(batch.size(),
                    batch_i + batch_size));

            // Run kcore on graph
            auto cores = KCore(graph, 16);

            parallel_for (0, layers.n, [&] (size_t cur_node_index) {
                auto b = (size_t) ceil(batch_i/(batch_size * 1.0));
                ground_truth_container.ground_truth[b][cur_node_index] = cores[cur_node_index];
            });
        }
    }

    parlay::random rng = parlay::random();

    std::cout << "Barrier is initializing" << std::endl;
    volatile struct {
        volatile uint64_t flag = 0;
        char padding[64 - sizeof(uint64_t)];
    } stop;

    std::atomic<int> counter_seq[num_reader_threads];
    std::atomic<double> error_seq[num_reader_threads];
    std::atomic<double> latency_seq[num_reader_threads];
    std::thread read_thread_seq[num_reader_threads];

    size_t main_counter = 0;
    barrier_t* sync_point = new barrier_t();
    barrier_t* sync_point_end = new barrier_t();
    barrier_init(sync_point, num_reader_threads + 1);
    barrier_init(sync_point_end, num_reader_threads + 1);

    for(size_t thread_i = 0; thread_i < num_reader_threads; thread_i++) {
        read_thread_seq[thread_i] = std::thread([sync_point, sync_point_end,
                &layers, &rng, &stop, &counter_seq, &ground_truth_container, &compare_exact,
                &error_seq, thread_i] {
            std::cout << "Thread " << thread_i << " is waiting" << std::endl;
            barrier_cross(sync_point);
            std::cout << "Thread " << thread_i << " is running" << std::endl;
            size_t my_counter = 0;

            float approx_error_total = 0.0;
            std::vector<double> all_latency;

            rng.fork(thread_i);
            while (stop.flag == 0) {
                timer read_timer; read_timer.start();
                auto random_vertex = rng.rand() % layers.n;
                rng = rng.next();

                auto retry = true;
                if (nonlinearizable) {
                    retry = false;
                    auto approx_core =
                        layers.get_core_from_level(layers.L[random_vertex].level);
                    if (compare_exact) {
                        uintE exact_core =
                            ground_truth_container.ground_truth[b1][random_vertex];

                        if (exact_core > 0 && approx_core > 0) {
                            approx_error_total += (exact_core > approx_core) ?
                            (float) exact_core / (float) approx_core :
                            (float) approx_core / (float) exact_core;
                        }
                    }
                }

                while(retry && stop.flag == 0) {
                    auto b1 = layers.batch_num;
                    auto l1 = layers.L[random_vertex].level;

                    auto root = layers.find_compress(random_vertex, layers.descriptor_array.begin());
                    auto l2 = layers.L[random_vertex].level;
                    auto b2 = layers.batch_num;

                    if (b1 != b2)
                        continue;

                    if (root != UINT_E_MAX) {
                        auto approx_core =
                            layers.get_core_from_level(
                                    layers.descriptor_array[random_vertex].old_level);

                        if (compare_exact) {
                                uintE exact_core = ground_truth_container.ground_truth[b1][random_vertex];
                                if (exact_core > 0 && approx_core > 0) {
                                    approx_error_total += (exact_core > approx_core) ?
                                        (float) exact_core / (float) approx_core :
                                        (float) approx_core / (float) exact_core;
                                }
                        }
                        retry = false;
                    }

                    else {
                        if (l1 == l2) {
                            auto approx_core = layers.get_core_from_level(l1);
                            if (compare_exact) {
                                auto exact_core = ground_truth_container.ground_truth[b1][random_vertex];
                                if (exact_core > 0 && approx_core > 0) {
                                    approx_error_total += (exact_core > approx_core) ?
                                        (float) exact_core / (float) approx_core :
                                        (float) approx_core / (float) exact_core;
                                }
                            }

                            retry = false;
                        } else
                            continue;
                    }
                }
                double read_latency = read_timer.end();
                all_latency.push_back(read_latency);
                my_counter++;
            }

            if (my_counter > 0) {
                double average_error = ((double)approx_error_total)/((double)my_counter);
                std::cout << average_error << std::endl;
                error_seq[thread_i] = (double) average_error; //, std::memory_order_release);
            } else {
                error_seq[thread_i] = (double) 1000000; //, std::memory_order_release);
            }
            counter_seq[thread_i].store(my_counter, std::memory_order_release);

            barrier_cross(sync_point_end);
        });
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(thread_i, &cpuset);
        pthread_setaffinity_np(read_thread_seq[thread_i].native_handle(), sizeof(cpu_set_t), &cpuset);
    }

    // First, insert / delete everything up to offset
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
            num_insertion_flips += layers.batch_insertion(batch_insertions);
            num_deletion_flips += layers.batch_deletion(batch_deletions);
        }
    }

    stop.flag = 0;
    std::cout << "Main thread batch is waiting" << std::endl;
    timer overall_timer; overall_timer.start();
    barrier_cross(sync_point);
    std::cout << "Main thread batch started" << std::endl;

    for (size_t i = offset; i < batch.size(); i += batch_size) {
        timer t; t.start();
        auto end_size = std::min(i + batch_size, batch.size());
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

        num_insertion_flips += 0;

        layers.batch_insertion(batch_insertions);
        double insertion_time = t.stop();

        t.start();
        num_deletion_flips += layers.batch_deletion(batch_deletions);

        max_degree = layers.max_degree();

        double deletion_time = t.stop();
        double tt = insertion_time + deletion_time;

        std::cout << "### Batch Running Time: " << tt << std::endl;
        std::cout << "### Insertion Running Time: " << insertion_time << std::endl;
        std::cout << "### Deletion Running Time: " << deletion_time << std::endl;
        std::cout << "### Batch Num: " << end_size - offset << std::endl;
        std::cout << "### Coreness Estimate: " << layers.max_coreness() << std::endl;
        std::cout << "### Number Insertion Flips: " << num_insertion_flips << std::endl;
        std::cout << "### Number Deletion Flips: " << num_deletion_flips << std::endl;
        std::cout << "### Max Outdegree: " << max_degree << std::endl;

        if (layers.batch_num * batch_size <= batch.size() - batch_size - 2)
            layers.batch_num++;

        if (get_size) {
            auto size = layers.get_size();
            std::cout << "### Size: " << size << std::endl;
        }
    }

    stop.flag = 1;
    double overall_time = overall_timer.stop();
    barrier_cross(sync_point_end);
    for (size_t t_i = 0; t_i < num_reader_threads; t_i++) {
        read_thread_seq[t_i].join();
    }

    std::cout << "### Total Time: " << overall_time << std::endl;
    double total_error = 0;
    for (size_t t_i = 0; t_i < num_reader_threads; t_i++) {
        main_counter += counter_seq[t_i].load(std::memory_order_acquire);
        total_error += error_seq[t_i].load(std::memory_order_acquire);
    }

    if (num_reader_threads > 0) {
        std::cout << "### Average Error: " << total_error/num_reader_threads << std::endl;
        std::cout << "### Main Counter: " << main_counter << std::endl;
        std::cout << "### Num Reader Threads: " << num_reader_threads << std::endl;
        std::cout << "### Throughput: " << main_counter/overall_time << std::endl;
    }
}

template <class Graph, class W>
inline void RunLDS(Graph& G, BatchDynamicEdges<W> batch_edge_list, long batch_size,
        bool compare_exact, double eps, double delta, bool optimized_insertion,
        size_t offset, bool get_size, size_t optimized_all, size_t num_reader_threads) {
        //uintE** core_ground_truth) {
    uintE max_vertex = std::max(uintE{G.n}, batch_edge_list.max_vertex);
    auto layers = LDS(max_vertex, eps, delta, optimized_insertion, optimized_all);
    if (batch_edge_list.max_vertex > 0)
        RunLDS(batch_edge_list, batch_size, compare_exact, layers, optimized_insertion,
                offset, get_size, num_reader_threads);
}

}  // namespace gbbs
