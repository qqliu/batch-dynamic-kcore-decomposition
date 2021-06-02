#include "UndirectedGraph.h"
#include "timer.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <stack>
#include <unordered_map>
#include <vector>

using namespace std;

UndirectedGraph::UndirectedGraph() { vertexNum = edgeNum = 0; }

UndirectedGraph::~UndirectedGraph() {
  edges.clear();
  edgeNum = 0;
  vertexNum = 0;
}

const std::vector<int> &UndirectedGraph::neighborsOf(int node) const {
  if (node >= edges.size()) {
    return emptySet;
  }
  return edges[node];
}

int UndirectedGraph::degreeOf(int node) const {
  return (int)edges[node].size();
}

bool UndirectedGraph::addEdge(int node1, int node2) {
//  edges[node1].push_back(node2);
//  edges[node2].push_back(node1);
//  ++edgeNum;
      if(node1 == node2) {
          return false;
      }

      int maxNode = max(node1,node2);
      while(maxNode + 1 > vertexNum) {
          ++vertexNum;
          edges.emplace_back(vector<int>());
      }

      bool ret = false;
      auto it = lower_bound(edges[node1].begin(),edges[node1].end(),node2);
      if(it == edges[node1].end() || *it != node2) {
          edges[node1].insert(it,node2);
          ret = true;
      }
      it = lower_bound(edges[node2].begin(),edges[node2].end(),node1);
      if(it == edges[node2].end() || *it != node1) {
          edges[node2].insert(it,node1);
          ret = true;
      }
      if(ret) {
          edgeNum++;
      }
      return ret;
}

bool UndirectedGraph::removeEdge(int node1, int node2) {
  if (node1 >= vertexNum || node2 >= vertexNum) {
    return false;
  }

  bool found = false;
  auto it = lower_bound(edges[node1].begin(), edges[node1].end(), node2);
  if (it != edges[node1].end() && *it == node2) {
    edges[node1].erase(it);
    found = true;
  }

  it = lower_bound(edges[node2].begin(), edges[node2].end(), node1);
  if (it != edges[node2].end() && *it == node1) {
    edges[node2].erase(it);
    found = true;
  }

  if (found) {
    --edgeNum;
  }
  return found;
}

