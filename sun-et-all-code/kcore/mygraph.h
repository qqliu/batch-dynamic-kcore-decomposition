#ifndef GRAPH_H
#define GRAPH_H

#include <vector>

struct vertex_descriptor
{
    int core;
    int mcd;
    int pcd;
    int cd;

    vertex_descriptor()
    {
        core = mcd = pcd = cd = 0;
    }
};

class Graph
{
public:
    int vertexNum;
    int edgeNum;
    std::vector<int> allcores;

    std::vector<vertex> vertices;
    std::vector<std::vector<int>> edges;
};

#endif // GRAPH_H




