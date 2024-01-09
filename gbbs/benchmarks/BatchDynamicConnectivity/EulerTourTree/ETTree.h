#include <unordered_set>
#include <stack>

#include "gbbs/gbbs.h"
#include "gbbs/dynamic_graph_io.h"
#include "benchmarks/KCore/JulienneDBS17/KCore.h"
#include "benchmarks/EdgeOrientation/ParallelLDS/sparse_set.h"
#include "benchmarks/BatchDynamicConnectivity/SkipList/SkipList.h"

namespace gbbs {

struct ETTree {

    struct ETTreeBase {
        SkipList::SkipListElement edge; // edge (u, v)
        SkipList::SkipListElement* twin; // edge (v, u)

        bool skip_mark;

        ETTreeBase() {
            edge = SkipList::SkipListElement();
            twin = nullptr;

            skip_mark = false;
        }

        ETTreeBase(uintE u, uintE v, SkipList::SkipListElement* twin_) {
                // TODO: make the value take two values: u and v
                edge = SkipList::SkipListElement(u, nullptr, nullptr, v);
                skip_mark = false;
                twin = twin_;
        }
    };

    void cut(uintE u, uintE v) {
            SkipList::SkipListElement
    }
};

void RunETTree() {
        std::cout << "ET tree" << std::endl;
}

}  // namespace gbbs
