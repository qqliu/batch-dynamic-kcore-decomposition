#ifndef _KCORE_HPP_
#define _KCORE_HPP_

#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <set>
#include <fstream>
#include <random>
#include "UndirectedGraph.h"
#include "timer.h"
#include "bitset.h"
#include "thread_pool.h"
#include "queue.h"

using namespace std;

class Decomposition {
public:
    Decomposition();

    Decomposition(int nThreads);

    bool Load(const char *graphfile, const char *edgefile);

    float LoadBatch(std::vector<string>& vstrings, int start, int batch_size);

    float LoadBatchDeletes(std::vector<string>& vstrings, int start, int batch_size);

    long long GetSize();

    bool LoadStatic(const char *staticfile);

    void Insert(float load_time);

    void Delete(float delete_time);

    void dumpCores() const;

    bool calcCores(const char *filename);
private:
    void init();

    void computeCores();

    void computeInsertEdgeSet();

    void computeDeleteEdgeSet();

    std::vector<int> insertCenterPoints();

    std::vector<int> deleteCenterPoints();


    int insertPrecore(int node);

    int deletePrecore(int node);

    void superiorInsertMCD(std::unordered_map<int, std::vector<Edge> >::const_iterator);

    void superiorDelete(int k,std::vector<int>& root);

    void deleteRemove(int *cd, char *removed, int node, int kcore);

    void insertRemoveMCD(int *cd, char *removed, int node, int kcore);

    int computeMCD(int node, int kcore);

    void insertPreprocess();

    void deletePreprocess();

    void compare();

    void addDeleteRoot(int node) {
      if (kcores[node] < deleteRoots.size()) {
        deleteRoots[kcores[node]].push_back(node);
      } else {
        deleteRoots.resize(kcores[node] * 2);
        deleteRoots[kcores[node]].push_back(node);
      }
    }
private:
    std::vector<int> dummy;
    std::vector<int> mcd;
    std::vector<int> pcd;
    std::vector<char> visited;
    std::vector<char> removed;
    std::vector<int> cd;

    UndirectedGraph graph;
    UndirectedGraph insertedGraph;
    UndirectedGraph deletedGraph;

    int maxNodeId;

    std::vector<int> nodeDegree;
    std::vector<int> kcores;
    std::vector<int> oldcores;
    std::vector<int> cps; // center points
    std::unordered_map<int, std::vector<Edge> > core2Edges;

    std::vector<Bitset> palette;
    std::vector<int> colors;
    std::vector<std::vector<int> > centerPoints;
    std::vector<int> delete_vertices;
    std::vector<int> original_cores;

    std::vector<int> insertedNodes;
    std::vector<int> isSelected;
    std::vector<int> deletedNodes;
    std::vector<std::vector<int>> deleteRoots;
    int searchIndex;

    Timer timer;
    std::mt19937 random;
    ThreadPool pool;
    lqueue<int> queue;
};


#endif // _KCORE_H
