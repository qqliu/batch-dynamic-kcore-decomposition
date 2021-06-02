#ifndef KCORE_GRAPH_H
#define KCORE_GRAPH_H

#include <utility>
#include <vector>
#include <unordered_map>
#include <map>

typedef std::pair<int,int> Edge;
typedef std::vector<Edge> EdgeSeq;

class Graph
{
public:
    Graph();
    bool read(const char* filename);

    void kcore_decomposition();

    void insert_edges(EdgeSeq& edges);
    bool insert_edges(const char* edges_filename);

    void saveCores(const char* filename)const;
    std::unordered_map<int,int> cores();
private:
    void remove_edge(Edge edge);
    std::vector<int> centerPoints();
    std::vector<int> twoHopNeighbors(std::map<int,std::vector<int>>::iterator);
    std::map<int,std::vector<int>>::iterator maxIncreaseDegree();
    std::vector<Edge> computeInsertEdgeSet();
    int precore(int node);
    std::vector<int> neighborsWithInserted(int node);
private:
    std::unordered_map<int,std::unordered_map<int,bool>> matrix;
    std::map<int,std::vector<int>> inserted_graph;
    std::vector<int> inserted_size;
    int maxNodeId;
    std::vector<int> nodeInfos;

    std::unordered_map<int,int> kcores;
};

#endif //KCORE_GRAPH_H
