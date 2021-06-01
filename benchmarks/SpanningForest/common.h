#pragma once

#include "gbbs/graph.h"

namespace gbbs {
using edge = std::pair<uintE, uintE>;
using parent = uintE;

constexpr edge empty_edge = std::make_pair(UINT_E_MAX, UINT_E_MAX);

namespace spanning_forest {
  uintE largest_comp = UINT_E_MAX;

  auto parents_to_edges(sequence<parent>& parents) -> sequence<edge> {
    auto all_edges = pbbslib::make_delayed<edge>(parents.size(), [&] (uintE i) {
      return std::make_pair(i,parents[i]);
    });
    return pbbslib::filter(all_edges, [&] (const edge& e) {
      return (e.first != e.second) && (e.second != UINT_E_MAX);
    });
  }
}
}  // namespace gbbs
