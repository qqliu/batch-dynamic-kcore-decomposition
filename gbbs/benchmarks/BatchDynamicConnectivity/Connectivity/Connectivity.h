#include <unordered_set>
#include <stack>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
#include "gbbs/pbbslib/sparse_table.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "benchmarks/EdgeOrientation/ParallelLDS/sparse_set.h"
#include "benchmarks/BatchDynamicConnectivity/EulerTourTree/ETTree.h"

namespace gbbs {

struct Connectivity {
    Connectivity() {};
};

void RunConnectivityTest() {
        std::cout << "Connectivity Test" << std::endl;
}

}  // namespace gbbs
