#ifndef _UNDIRECTEDGRAPH_HPP_
#define _UNDIRECTEDGRAPH_HPP_


#include <vector>
#include <set>
#include <utility>

typedef std::pair<int,int> Edge;
typedef std::vector<Edge> EdgeSeq;


class UndirectedGraph
{
public:
    UndirectedGraph();
    ~UndirectedGraph();

    // read only operation
    bool empty()const {return edgeNum == 0;}
    const std::vector<int>& neighborsOf(int node)const;
    int degreeOf(int node)const;
    int size()const {return (int)edges.size();}
    // bool isConnected(int from,int to) const;
    int getEdgeNum()const{return edgeNum;}

    bool addEdge(int node1,int node2);

    template <class InputIt>
        void addEdges(InputIt first,InputIt last)
    {
        while(first != last) {
            addEdge(first->first,first->second);
            ++first;
        }
    }

    bool removeEdge(int node1,int node2);

    template <class InputIt>
        void removeEdges(InputIt first,InputIt last)
    {
        while(first != last) {
            removeEdge(first->first,first->second);
            ++first;
        }
    }

private:
    class EdgeIterator {
    public:
        explicit EdgeIterator(const UndirectedGraph *g) : graph(g), i(-1), j(-1) {}

        explicit EdgeIterator(const UndirectedGraph *g, int ii, int jj) : graph(g), i(ii), j(jj) {}

        bool operator==(const EdgeIterator &other) {
            if (i < 0 || j < 0 || other.i < 0 || other.j < 0) {
                return true;
            }
            return graph->edges[i][j] == other.graph->edges[i][j];
        }

        bool operator!=(const EdgeIterator &other) {
            return i != other.i && j != other.j;

        }

        EdgeIterator operator++(int) {
            ++j;
            while (i < graph->edges.size()) {
                for (; j < graph->edges[i].size(); j++) {
                    if (i < graph->edges[i][j]) {
                        return EdgeIterator(graph, i, j);
                    }
                }

                ++i;
                j = 0;
            }
            return EdgeIterator(graph);
        }

        EdgeIterator &operator++() {
            ++j;
            while (i < graph->edges.size()) {
                for (; j < graph->edges[i].size(); j++) {
                    if (i < graph->edges[i][j]) {
                        return *this;
                    }
                }

                ++i;
                j = 0;
            }
            i = -1;
            j = -1;
            return *this;
        }

        Edge operator*() {
            return Edge(i, graph->edges[i][j]);
        }

    private:
        const UndirectedGraph *graph;
        int i, j;
    };

public:
    typedef EdgeIterator iterator;

    iterator begin() const {
        int i = 0, j = 0;
        while (i < edges.size()) {
            while (j < edges[i].size()) {
                if (i < edges[i][j]) {
                    return iterator(this, i, j);
                }
            }
            ++i;
            j = 0;
        }

        return iterator(this);
    }

    iterator end() const { return iterator(this); }
public:
    int vertexNum;
    int edgeNum;
    std::vector<std::vector<int>> edges;
    std::vector<int> emptySet;

    std::vector<std::vector<int>> twoHopEdges;
};

#endif // _UNDIRECTEDGRAPH_HPP_