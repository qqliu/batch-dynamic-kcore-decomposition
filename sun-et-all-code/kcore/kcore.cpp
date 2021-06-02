#include <algorithm>
#include <cassert>
#include <fstream>
#include <future>
#include <iostream>
#include <iterator>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "kcore.h"

using namespace std;

Decomposition::Decomposition()
    : pool(16),
      random(std::random_device{}()) {
  maxNodeId = -1;
  searchIndex = 0;
  deleteRoots.resize(64);
}

Decomposition::Decomposition(int nThreads)
    : pool(nThreads),
      random(std::random_device{}()) {
  maxNodeId = -1;
  searchIndex = 0;
  deleteRoots.resize(64);
}

long long Decomposition::GetSize() {
    long long size = 0;
    for (int i = 0; i < dummy.size(); i++) {
        size += sizeof(dummy[i]);
    }
    for (int i = 0; i < mcd.size(); i++) {
        size += sizeof(mcd[i]);
    }
    for (int i = 0; i < pcd.size(); i++) {
        size += sizeof(pcd[i]);
    }
    for (int i = 0; i < visited.size(); i++) {
        size += sizeof(visited[i]);
    }
    for (int i = 0; i < removed.size(); i++) {
        size += sizeof(removed[i]);
    }
    for (int i = 0; i < cd.size(); i++) {
        size += sizeof(cd[i]);
    }

    size += sizeof(maxNodeId);

    for (int i = 0; i < nodeDegree.size(); i++) {
        size += sizeof(nodeDegree[i]);
    }
    for (int i = 0; i < kcores.size(); i++) {
        size += sizeof(kcores[i]);
    }
    for (int i = 0; i < oldcores.size(); i++) {
        size += sizeof(oldcores[i]);
    }
    for (int i = 0; i < cps.size(); i++) {
        size += sizeof(cps[i]);
    }
    for (int i = 0; i < palette.size(); i++) {
        size += sizeof(palette[i]);
    }
    for (int i = 0; i < colors.size(); i++) {
        size += sizeof(colors[i]);
    }
    for (int i = 0; i < centerPoints.size(); i++) {
        for (int j = 0; j < centerPoints[i].size(); j++)
            size += sizeof(centerPoints[i][j]);
    }
    for (int i = 0; i < delete_vertices.size(); i++) {
        size += sizeof(delete_vertices[i]);
    }
    for (int i = 0; i < original_cores.size(); i++) {
        size += sizeof(original_cores[i]);
    }
    for (int i = 0; i < insertedNodes.size(); i++) {
        size += sizeof(insertedNodes[i]);
    }
    for (int i = 0; i < isSelected.size(); i++) {
        size += sizeof(isSelected[i]);
    }
    for (int i = 0; i < deletedNodes.size(); i++) {
        size += sizeof(deletedNodes[i]);
    }
    for (int i = 0; i < deleteRoots.size(); i++) {
        for (int j = 0; j < deleteRoots[i].size(); j++)
            size += sizeof(deleteRoots[i][j]);
    }

    size += sizeof(searchIndex);
    return size;

    /*for (int i = 0; i < queue.size(); i++) {
        size += sizeof(queue[i]);
    }*/
}

float Decomposition::LoadBatch(vector<string>& vstrings, int start, int batch_size) {
    auto insert_edge_time = 0;
    {
      set<int> delete_set;
      assert(!(vstrings.size() & 1));
      int node1, node2;
      Timer insert_edge_timer;
      for (int i = start; i <= min((int)vstrings.size() - 2, 2 * batch_size + start - 2); i += 2) {
        node1 = stoi(vstrings[i]);
        node2 = stoi(vstrings[i + 1]);
        delete_set.insert(node1);
        delete_set.insert(node2);
        insert_edge_timer.reset();
        insertedGraph.addEdge(node1, node2);
        insertedGraph.addEdge(node2, node1);
        insert_edge_time += insert_edge_timer.elapsed();
        //deletedGraph.addEdge(node1, node2);
        //deletedGraph.addEdge(node2, node1);
        maxNodeId = max(max(node1, node2), maxNodeId);
      }
      //copy(delete_set.begin(),delete_set.end(),back_inserter(delete_vertices));
    }

  dummy.resize(maxNodeId + 1);
  mcd.resize(maxNodeId + 1);
  pcd.resize(maxNodeId + 1);
  removed.resize(maxNodeId + 1);
  visited.resize(maxNodeId + 1);
  cd.resize(maxNodeId + 1);
  palette.resize(maxNodeId + 1);
  colors.resize(maxNodeId + 1);
  nodeDegree.resize(maxNodeId + 1);
  isSelected.resize(maxNodeId + 1);
  fill(colors.begin(), colors.end(), -1);

  kcores.resize(maxNodeId + 1);
  oldcores.resize(maxNodeId + 1);

  //computeCores();
  original_cores = kcores;

  /*for (int i = 0; i < (int)mcd.size(); i++) {
    mcd[i] = 0;
    for (auto u : graph.neighborsOf(i)) {
      if (kcores[u] >= kcores[i]) {
        ++mcd[i];
      }
    }
  }*/
  return insert_edge_time;
}

float Decomposition::LoadBatchDeletes(vector<string>& vstrings, int start, int batch_size) {
    auto delete_edge_time = 0;
    {
      set<int> delete_set;
      assert(!(vstrings.size() & 1));
      int node1, node2;
      for (int i = start; i <= min((int)vstrings.size() - 2, 2 * batch_size + start - 2); i += 2) {
        node1 = stoi(vstrings[i]);
        node2 = stoi(vstrings[i + 1]);
        Timer delete_edge_timer;
        delete_edge_timer.reset();
        delete_set.insert(node1);
        delete_set.insert(node2);
        delete_edge_time += delete_edge_timer.elapsed();
        //insertedGraph.addEdge(node1, node2);
        //insertedGraph.addEdge(node2, node1);
        deletedGraph.addEdge(node1, node2);
        deletedGraph.addEdge(node2, node1);
        maxNodeId = max(max(node1, node2), maxNodeId);
      }
      copy(delete_set.begin(),delete_set.end(),back_inserter(delete_vertices));
    }

  /*
  dummy.resize(maxNodeId + 1);
  mcd.resize(maxNodeId + 1);
  pcd.resize(maxNodeId + 1);
  removed.resize(maxNodeId + 1);
  visited.resize(maxNodeId + 1);
  cd.resize(maxNodeId + 1);
  palette.resize(maxNodeId + 1);
  colors.resize(maxNodeId + 1);
  nodeDegree.resize(maxNodeId + 1);
  isSelected.resize(maxNodeId + 1);
  fill(colors.begin(), colors.end(), -1);

  kcores.resize(maxNodeId + 1);
  oldcores.resize(maxNodeId + 1);*/

  //computeCores();
  //original_cores = kcores;

  /*for (int i = 0; i < (int)mcd.size(); i++) {
    mcd[i] = 0;
    for (auto u : graph.neighborsOf(i)) {
      if (kcores[u] >= kcores[i]) {
        ++mcd[i];
      }
    }
  }*/
  return delete_edge_time;
}

bool Decomposition::LoadStatic(const char *staticfile) {
    // construct graph
    ifstream fin;
    fin.open(staticfile, ios_base::in);
    if (fin.fail()) {
      return false;
    }
    string content;
    content.assign((istreambuf_iterator<char>(fin)),
                   (istreambuf_iterator<char>()));
    {
      stringstream ss(content);
      istream_iterator<string> begin(ss);
      istream_iterator<string> end;
      vector<string> vstrings(begin, end);
      assert(!(vstrings.size() & 1));
      int node1, node2;
      for (int i = 0; i <= vstrings.size() - 2; i += 2) {
        node1 = stoi(vstrings[i]);
        node2 = stoi(vstrings[i + 1]);

        graph.addEdge(node1, node2);
        graph.addEdge(node2, node1);
        maxNodeId = max(max(node1, node2), maxNodeId);
      }
    }
    fin.close();

    dummy.resize(maxNodeId + 1);
    mcd.resize(maxNodeId + 1);
    pcd.resize(maxNodeId + 1);
    removed.resize(maxNodeId + 1);
    visited.resize(maxNodeId + 1);
    cd.resize(maxNodeId + 1);
    palette.resize(maxNodeId + 1);
    colors.resize(maxNodeId + 1);
    nodeDegree.resize(maxNodeId + 1);
    isSelected.resize(maxNodeId + 1);
    fill(colors.begin(), colors.end(), -1);

    kcores.resize(maxNodeId + 1);
    oldcores.resize(maxNodeId + 1);

    computeCores();
    original_cores = kcores;

    for (int i = 0; i < (int)mcd.size(); i++) {
        mcd[i] = 0;
        for (auto u : graph.neighborsOf(i)) {
            if (kcores[u] >= kcores[i]) {
                ++mcd[i];
            }
        }
    }

    return true;
}

bool Decomposition::Load(const char *graphfile, const char *edgefile) {
  {
    // construct graph
    ifstream fin;
    fin.open(graphfile, ios_base::in);
    if (fin.fail()) {
      return false;
    }
    ifstream fedges;
    fedges.open(edgefile);
    if (fedges.fail()) {
      fin.close();
      return false;
    }

    string content;
    content.assign((istreambuf_iterator<char>(fin)),
                   (istreambuf_iterator<char>()));
    {
      stringstream ss(content);
      istream_iterator<string> begin(ss);
      istream_iterator<string> end;
      vector<string> vstrings(begin, end);
      assert(!(vstrings.size() & 1));
      int node1, node2;
      for (int i = 0; i <= vstrings.size() - 2; i += 2) {
        node1 = stoi(vstrings[i]);
        node2 = stoi(vstrings[i + 1]);

        graph.addEdge(node1, node2);
        graph.addEdge(node2, node1);
        maxNodeId = max(max(node1, node2), maxNodeId);
      }

    }

    content.clear();
    content.assign((istreambuf_iterator<char>(fedges)),
                   (istreambuf_iterator<char>()));
    {
      set<int> delete_set;
      stringstream ss(content);
      istream_iterator<string> begin(ss);
      istream_iterator<string> end;
      vector<string> vstrings(begin, end);
      assert(!(vstrings.size() & 1));
      int node1, node2;
      for (int i = 0; i <= vstrings.size() - 2; i += 2) {
        node1 = stoi(vstrings[i]);
        node2 = stoi(vstrings[i + 1]);
        delete_set.insert(node1);
        delete_set.insert(node2);
        insertedGraph.addEdge(node1, node2);
        insertedGraph.addEdge(node2, node1);
        deletedGraph.addEdge(node1, node2);
        deletedGraph.addEdge(node2, node1);
        maxNodeId = max(max(node1, node2), maxNodeId);
      }
      copy(delete_set.begin(),delete_set.end(),back_inserter(delete_vertices));
    }
    fin.close();
    fedges.close();
  }


  dummy.resize(maxNodeId + 1);
  mcd.resize(maxNodeId + 1);
  pcd.resize(maxNodeId + 1);
  removed.resize(maxNodeId + 1);
  visited.resize(maxNodeId + 1);
  cd.resize(maxNodeId + 1);
  palette.resize(maxNodeId + 1);
  colors.resize(maxNodeId + 1);
  nodeDegree.resize(maxNodeId + 1);
  isSelected.resize(maxNodeId + 1);
  fill(colors.begin(), colors.end(), -1);

  kcores.resize(maxNodeId + 1);
  oldcores.resize(maxNodeId + 1);

  computeCores();
  original_cores = kcores;

  for (int i = 0; i < (int)mcd.size(); i++) {
    mcd[i] = 0;
    for (auto u : graph.neighborsOf(i)) {
      if (kcores[u] >= kcores[i]) {
        ++mcd[i];
      }
    }
  }


  // cout << edgefile << '\t' << insertedGraph.getEdgeNum() << '\t';
  return true;
}

void Decomposition::init() {
  fill(mcd.begin(), mcd.end(), 0);
  fill(cd.begin(), cd.end(), 0);
  fill(visited.begin(), visited.end(), 0);
  fill(removed.begin(), removed.end(), 0);
}

void Decomposition::insertPreprocess() {

  insertedNodes.clear();
  fill(isSelected.begin(), isSelected.end(), 0);
  fill(colors.begin(), colors.end(), -1);
  fill(palette.begin(), palette.end(), Bitset());

  int maxSize = insertedGraph.size();
  unordered_map<int, int> node2inc;

  for (int i = 0; i < maxSize; i++) {
    if (insertedGraph.degreeOf(i) > 1) {
      int c = 0;
      for (auto n : insertedGraph.neighborsOf(i)) {
        if (kcores[n] >= kcores[i]) {
          c++;
        }
      }

      if (c > 9) {
        node2inc[i] = c;
        insertedNodes.push_back(i);
      }
    }
  }

  sort(insertedNodes.begin(), insertedNodes.end(),
       [this, &node2inc](int lhs, int rhs) {
        return node2inc[lhs] > node2inc[rhs];
       });

  centerPoints.clear();
  int maxColor = -1;

  timer.reset();

  for (auto n : insertedNodes) {
    if (maxColor == -1) {
      maxColor = 0;
      colors[n] = 0;
      if (insertedNodes.size() == 1)
        break;
      vector<int> neighbors(graph.neighborsOf(n).begin(),
                            graph.neighborsOf(n).end());
      copy(insertedGraph.neighborsOf(n).begin(),
           insertedGraph.neighborsOf(n).end(), back_inserter(neighbors));

      for (auto v : neighbors) {
        palette[v].set(0);
      }
      continue;
    } else {

      Bitset bt;
      for (auto v : graph.neighborsOf(n)) {
        bt.mix(palette[v]);
      }
      for (auto v : insertedGraph.neighborsOf(n)) {
        bt.mix(palette[v]);
      }
      int color = bt.nextAvailable();
      colors[n] = color;
      maxColor = max(maxColor, color);

      for (auto v : graph.neighborsOf(n)) {
        palette[v].set(color);
      }
      for (auto v : insertedGraph.neighborsOf(n)) {
        palette[v].set(color);
      }
    }
  }

  centerPoints.resize(maxColor + 1);
  for (auto n : insertedNodes) {
    centerPoints[colors[n]].push_back(n);
  }
}

void Decomposition::Insert(float load_time) {

  vector<future<void>> futures;
  int round = 0;
  long findTime = 0;

  // init + compute k-Joint Edge Set
  Timer correct_time;
  correct_time.reset();
  searchIndex = 0;

  insertPreprocess();
  auto preprocess = correct_time.elapsed();

  //Timer t;
  timer.reset();

  while (!insertedGraph.empty()) {
    ++round;

    fill(mcd.begin(),mcd.end(),0);
    init();

    //t.reset();
    computeInsertEdgeSet();
    //findTime += t.elapsed();

    futures.clear();

    for (auto it = core2Edges.begin(); it != core2Edges.end(); ++it) {
      // superiorInsertMCD(it);
      futures.emplace_back(
          pool.enqueue(&Decomposition::superiorInsertMCD, this, it));
    }

    for (auto &f : futures) {
      f.get();
    }

    for (int i = 0; i < maxNodeId + 1; i++) {
      if (visited[i] && !removed[i]) {
        kcores[i]++;
      }
    }

    for (auto u : cps) {
      kcores[u] = insertPrecore(u);
    }

  }

  auto time_elapsed = timer.elapsed();
  cout << "Insert\t" << time_elapsed << '\t' << round << '\n';
  cout << "Insert Correct\t" << time_elapsed + preprocess << '\t' << round << '\n';
  cout << "Insert Batch\t" << time_elapsed + preprocess + load_time << '\t' << round << '\n';
  //compare();
}

void Decomposition::deletePreprocess() {

  deletedNodes.clear();
  fill(isSelected.begin(), isSelected.end(), 0);
  fill(colors.begin(), colors.end(), -1);
  fill(palette.begin(), palette.end(), Bitset());

  unordered_map<int, int> node2dec;

  for (auto i = 0; i < deletedGraph.size(); i++) {
    if (deletedGraph.degreeOf(i) > 1) {
      int c = 0;
      for (auto n : deletedGraph.neighborsOf(i)) {
        if (kcores[n] >= kcores[i]) {
          c++;
        }
      }
      if (c > 9) {
        node2dec[i] = c;
        deletedNodes.push_back(i);
      }
    }
  }

  sort(deletedNodes.begin(), deletedNodes.end(),
       [this, &node2dec](int lhs, int rhs) {
         return node2dec[lhs] > node2dec[rhs];
       });

  centerPoints.clear();
  int maxColor = -1;
  timer.reset();

  for (auto n : deletedNodes) {
    if (maxColor == -1) {
      maxColor = 0;
      colors[n] = 0;
      if (deletedNodes.size() == 1)
        break;
      const vector<int> &neighbors = graph.neighborsOf(n);
      for (auto v : neighbors) {
        palette[v].set(0);
      }
      continue;
    } else {
      const vector<int> &neighbors = graph.neighborsOf(n);
      Bitset bt;
      for (auto v : neighbors) {
        bt.mix(palette[v]);
      }
      int color = bt.nextAvailable();
      colors[n] = color;
      maxColor = max(maxColor, color);

      for (auto v : neighbors) {
        palette[v].set(color);
      }
    }
  }

  centerPoints.resize(maxColor + 1);
  for (auto n : deletedNodes) {
    centerPoints[colors[n]].push_back(n);
  }
}

void Decomposition::Delete(float delete_time) {

  vector<future<void>> futures;
  long findTime = 0;
  int round = 0;

  fill(dummy.begin(),dummy.end(),-1);
  // init
  Timer correct_time;
  correct_time.reset();
  searchIndex = 0;
  deletePreprocess();
  auto preprocess = correct_time.elapsed();

  Timer t;
  timer.reset();

  while (!deletedGraph.empty()) {
    ++round;
    init();

    t.reset();
    computeDeleteEdgeSet();
    findTime += t.elapsed();

    futures.clear();

    /*
    for (auto it = core2Edges.begin(); it != core2Edges.end(); ++it) {
      futures.emplace_back(
          pool.enqueue(&Decomposition::superiorDelete, this, it));
    }
     */
    vector<int> affected;
    for (int i = 0; i < (int)deleteRoots.size(); i++) {
      if (!deleteRoots[i].empty()) {
        // superiorDelete(i,deleteRoots[i]);
        affected.push_back(i);
        futures.emplace_back(
            pool.enqueue(&Decomposition::superiorDelete,this,i,deleteRoots[i])
        );
      }
    }



    for (auto &f : futures) {
      f.get();
    }

    for (auto k : affected) {
      deleteRoots[k].clear();
    }

    /*
    for (int i = 0; i < maxNodeId + 1; i++) {
      if (visited[i] && removed[i]) {
        kcores[i]--;
      }
    }
     */
    int v;
    while (queue.pop(v)) {

      if (visited[v] && removed[v]) {
        kcores[v]--;
      }
    }

    for (auto u : cps) {
      kcores[u] = deletePrecore(u);
    }
  }

  auto time_elapsed = timer.elapsed();

  cout << "Delete\t" << time_elapsed << '\t' << round << '\n';
  cout << "Delete Correct\t" << preprocess + time_elapsed << '\t' << round << '\n';
  cout << "Delete Batch\t" << preprocess + time_elapsed + delete_time << '\t' << round << '\n';
}

void Decomposition::superiorDelete(int k,std::vector<int>& root) {

  for (auto r : root) {
    if (!visited[r]) {
      visited[r] = 1;
      if (!mcd[r]) {
        mcd[r] = computeMCD(r, kcores[r]);
      }
      cd[r] = mcd[r];
    }
    if (!removed[r] && cd[r] < k) {
      deleteRemove(&(*cd.begin()), &(*removed.begin()), r, k);
    }
  }
  root.clear();
}

void Decomposition::superiorInsertMCD(
    std::unordered_map<int, std::vector<Edge>>::const_iterator input) {
  int k = input->first;

  stack<int> stk;

  for (const auto &item : input->second) {

    int r = item.first;
    if (!visited[r] && !removed[r]) {
      if (mcd[r] == 0) {
        mcd[r] = computeMCD(r, kcores[r]);
      }
      if (cd[r] >= 0) {
        cd[r] = mcd[r];
      } else {
        cd[r] += mcd[r];
      }

      stk.push(r);
      visited[r] = 1;

      while (!stk.empty()) {
        int v = stk.top();
        stk.pop();
        if (cd[v] > k) {
          for (auto w : graph.neighborsOf(v)) {
            if (kcores[w] == k && !visited[w]) {

              if (mcd[w] == 0) {
                mcd[w] = computeMCD(w, kcores[w]);
              }

              if (mcd[w] > k) {
                stk.push(w);
                visited[w] = 1;
                cd[w] += mcd[w];
              } else {
                if (!removed[w]) {
                  insertRemoveMCD(&(*cd.begin()), &(*removed.begin()), w, k);
                }
              }
            }
          }
        } else {
          if (!removed[v]) {
            insertRemoveMCD(&(*cd.begin()), &(*removed.begin()), v, k);
          }
        }
      }
    }
  }
}


void Decomposition::deleteRemove(int *cd, char *removed, int node, int kcore) {
  stack<int> stk;
  removed[node] = 1;
  queue.push(node);

  stk.push(node);

  while (!stk.empty()) {
    int v = stk.top();
    stk.pop();
    const vector<int> &neighbors = graph.neighborsOf(v);
    for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
      int w = *it;

      if (kcores[w] == kcore) {
        if (!visited[w]) {
          visited[w] = 1;
          if (!mcd[w]) {
            mcd[w] = computeMCD(w, kcores[w]);
          }
          cd[w] += mcd[w];
        }

        cd[w]--;
        if (cd[w] < kcore && !removed[w]) {
          stk.push(w);
          removed[w] = 1;
          queue.push(w);
        }
      }
    }
  }
}

void Decomposition::insertRemoveMCD(int *cd, char *removed, int node,
                                    int kcore) {
  stack<int> stk;
  *(removed + node) = 1;
  stk.push(node);

  while (!stk.empty()) {
    int v = stk.top();
    stk.pop();
    for (auto w : graph.neighborsOf(v)) {
      if (kcores[w] == kcore) {
        cd[w]--;
        if (cd[w] > 0 && cd[w] <= kcore && !removed[w]) {
          stk.push(w);
          removed[w] = 1;
        }
      }
    }
  }
}

void Decomposition::computeDeleteEdgeSet() {

  static vector<Edge> Edel;
  static vector<int> garbage;

  // init
  garbage.clear();
  Edel.clear();


  cps = deleteCenterPoints();

	// std::cout << "delete cps.size() = " << cps.size() << '\t';

  for (auto v : cps) {
    dummy[v] = -2;
    garbage.push_back(v);
    oldcores[v] = kcores[v];
    kcores[v] = deletePrecore(v);
  }

  // remove center vertices' neighbors
  for (int i = 0; i < (int)cps.size(); i++) {
    int v = cps[i];

    for (auto u : graph.neighborsOf(v)) {
      if (kcores[u] >= kcores[v] && kcores[u] <= oldcores[v]) {
        dummy[u] = 0;
        garbage.push_back(u);
        addDeleteRoot(u);
      }
    }

    for (auto u : deletedGraph.neighborsOf(v)) {
      if (kcores[u] <= kcores[v]) {
        dummy[u] = 0;
        garbage.push_back(u);
        addDeleteRoot(u);
      }

      if (v < u) {
        Edel.emplace_back(Edge(v, u));
      } else {
        Edel.emplace_back(Edge(u, v));
      }
    }
  }

  graph.removeEdges(Edel.begin(), Edel.end());
  deletedGraph.removeEdges(Edel.begin(), Edel.end());
  Edel.clear();

  // superior
  vector<int> indices;

  for (int i = 0; i < (int)delete_vertices.size(); i++) {
    int w = delete_vertices[i];

    if (deletedGraph.degreeOf(w) > 0) {
      for (auto u : deletedGraph.neighborsOf(w)) {
        if (u > w) {
          if (kcores[u] < kcores[w] && dummy[u] == -1) {
            dummy[u] = 0;
            garbage.push_back(u);
            Edel.emplace_back(make_pair(u,w));
            addDeleteRoot(u);
          } else if (kcores[w] < kcores[u] && dummy[w] == -1) {
            dummy[w] = 0;
            garbage.push_back(w);
            Edel.emplace_back(make_pair(w,u));
            addDeleteRoot(w);
          } else if(kcores[w] == kcores[u] && dummy[w] == -1 && dummy[u] == -1) {
            dummy[w] = dummy[u] = 0;
            garbage.push_back(u);
            garbage.push_back(w);
            addDeleteRoot(u);
            addDeleteRoot(w);
            Edel.emplace_back(make_pair(u,w));
          }
        }
      }
    } else {
      indices.push_back(i);
    }
  }

  if (!indices.empty()) {
    for (int j = (int)indices.size() - 1; j >= 0; j--) {
      swap(delete_vertices.back(), delete_vertices[indices[j]]);
      delete_vertices.pop_back();
    }
  }

  for (auto v : garbage) {
    dummy[v] = -1;
  }

  graph.removeEdges(Edel.begin(), Edel.end());
  deletedGraph.removeEdges(Edel.begin(), Edel.end());

  // std::cout << Edel.size() << std::endl;
}

void Decomposition::computeInsertEdgeSet() {

  vector<Edge> Eins;
  // init
  core2Edges.clear();

  cps = insertCenterPoints();

	// std::cout << "insert cps.size() = " << cps.size() << '\t';

  fill(dummy.begin(), dummy.end(), -1);

  for (auto v : cps) {
    dummy[v] = -2;
    oldcores[v] = kcores[v];
    kcores[v] = insertPrecore(v);
  }

  for (auto v : cps) {
    // original neighbors
    for (auto n : graph.neighborsOf(v)) {
      if (kcores[n] >= oldcores[v] && kcores[n] <= kcores[v]) {
        dummy[n] = 0;
        core2Edges[kcores[n]].emplace_back(Edge(n, v));
      }
    }

    // new added neighbors
    for (auto u : insertedGraph.neighborsOf(v)) {

      if (kcores[u] <= kcores[v]) {
        core2Edges[kcores[u]].emplace_back(Edge(u, v));
        dummy[u] = 0;
      }

      if (v < u) {
        Eins.emplace_back(make_pair(v, u));
      } else {
        Eins.emplace_back(make_pair(u, v));
      }
    }
  }


  graph.addEdges(Eins.begin(), Eins.end());
  insertedGraph.removeEdges(Eins.begin(), Eins.end());
  Eins.clear();

  // add superior edge
  for (auto it = insertedGraph.begin(); it != insertedGraph.end(); ++it) {
    Edge edge = *it;
    int small = edge.first, large = edge.second;
    int core1 = kcores[small], core2 = kcores[large];

    if (core1 > core2) {
      small = edge.second;
      large = edge.first;
      int tmp = core1;
      core1 = core2;
      core2 = tmp;
    }

    if (core1 < core2 && dummy[small] == -1) {
      dummy[small] = 0;
      Eins.emplace_back(Edge(small, large));
      core2Edges[core1].emplace_back(Edge(small, large));
    } else if (core1 == core2 && dummy[small] == -1 && dummy[large] == -1) {
      dummy[small] = 0;
      dummy[large] = 0;
      Eins.emplace_back(Edge(small, large));
      core2Edges[core1].emplace_back(Edge(small, large));
    }
  }


  graph.addEdges(Eins.begin(), Eins.end());
  insertedGraph.removeEdges(Eins.begin(), Eins.end());
	// std::cout << Eins.size() << std::endl;
}

int Decomposition::computeMCD(int node, int kcore) {
  assert(kcores[node] == kcore);
  int c = 0;

  const std::vector<int> &neighbors = graph.neighborsOf(node);
  for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
    if (kcores[*it] >= kcore) {
      ++c;
    }
  }

  return c;
}

vector<int> Decomposition::deleteCenterPoints() {
  vector<int> points;

  for (; searchIndex < (int)centerPoints.size(); searchIndex++) {
    bool found = false;
    for (auto n : centerPoints[searchIndex]) {
      if (deletedGraph.degreeOf(n) > 1) {
        found = true;
        break;
      }
    }
    if (found) {
      points = centerPoints[searchIndex++];
      break;
    }
  }

  return points;
}

vector<int> Decomposition::insertCenterPoints() {

  vector<int> points;

  for (; searchIndex < (int)centerPoints.size(); ++searchIndex) {
    bool found = false;
    for (auto n : centerPoints[searchIndex]) {
      if (insertedGraph.degreeOf(n) > 1) {
        found = true;
        break;
      }
    }
    if (found) {
      points = centerPoints[searchIndex++];
      break;
    }
  }

  return points;
}

int Decomposition::deletePrecore(int node) {
  const vector<int> &neighbors = graph.neighborsOf(node);
  vector<int> core2count(64,0);
  int maxCore = -1;

  for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
    int core = kcores[*it];
    if (core > maxCore) {
      maxCore = core;
    }
    if (core < core2count.size()) {
      ++core2count[core];
    } else {
      core2count.resize(core * 2);
      ++core2count[core];
    }
  }

  const vector<int> &deletedNeighbors = deletedGraph.neighborsOf(node);
  for (auto it = deletedNeighbors.begin(); it != deletedNeighbors.end(); ++it) {
    int core = kcores[*it];
    if (core > maxCore) {
      maxCore = core;
    }
    if (core < core2count.size()) {
      --core2count[core];
    } else {
      core2count.resize(core * 2);
      --core2count[core];
    }
  }

  for (int i = maxCore; i >= 0; i--) {
    if (core2count[i] >= i) {
      return i;
    } else if (i > 0) {
      core2count[i - 1] += core2count[i];
    }
  }

  return 0;
}

int Decomposition::insertPrecore(int node) {

  vector<int> core2count(64, 0);
  int maxCore = -1;

  const vector<int> &neighbors = graph.neighborsOf(node);
  for (auto n : neighbors) {
    int core = kcores[n];

    if (core > maxCore) {
      maxCore = core;
      if (maxCore >= (int)core2count.size()) {
        core2count.resize(maxCore + 1);
      }
    }
    core2count[core]++;
  }
  const vector<int> &insertedNeighbors = insertedGraph.neighborsOf(node);
  for (auto n : insertedNeighbors) {
    int core = kcores[n];

    if (core > maxCore) {
      maxCore = core;
      if (maxCore >= (int)core2count.size()) {
        core2count.resize(maxCore + 1);
      }
    }
    core2count[core]++;
  }

  int acc = 0;
  for (int i = maxCore; i >= 1; i--) {
    acc += core2count[i];
    if (acc >= i) {
      return i;
    }
  }
  return 0;
}

void Decomposition::computeCores() {
  int d, md, start, num;
  int v, u, w, du, pu, pw;

  vector<int> bin;
  unordered_map<int, int> pos;
  map<int, int> vert;

  md = 0;
  for (int i = 0; i < graph.size(); i++) {
    d = graph.degreeOf(i);
    kcores[i] = d;
    if (d > md) {
      md = d;
    }
  }

  bin.resize(md + 1);
  fill(bin.begin(), bin.end(), 0);
  for (int i = 0; i < graph.size(); i++) {
    bin[kcores[i]]++;
  }

  start = 0;
  for (d = 0; d <= md; d++) {
    num = bin[d];
    bin[d] = start;
    start += num;
  }
  for (int i = 0; i < graph.size(); i++) {
    pos[i] = bin[kcores[i]];
    vert[pos[i]] = i;
    bin[kcores[i]]++;
  }
  for (d = md; d >= 1; d--) {
    bin[d] = bin[d - 1];
  }
  bin[0] = 0;
  for (auto it = vert.begin(); it != vert.end(); ++it) {
    v = it->second;
    const vector<int> &neighbors = graph.neighborsOf(v);
    for (auto it2 = neighbors.begin(); it2 != neighbors.end(); ++it2) {
      u = *it2;
      if (kcores[u] > kcores[v]) {
        du = kcores[u];
        pu = pos[u];
        pw = bin[du];
        w = vert[pw];
        if (u != w) {
          pos[u] = pw;
          vert[pu] = w;
          pos[w] = pu;
          vert[pw] = u;
        }

        bin[du]++;
        kcores[u]--;
      }
    }
  }
}

void Decomposition::dumpCores() const {
  for (int i = 0; i < (int)kcores.size(); i++) {
    cout << i << " " << kcores[i] << endl;
  }
}

bool Decomposition::calcCores(const char *filename) {
  ifstream fin;
  fin.open(filename, ios_base::in);
  if (fin.fail()) {
    return false;
  }
  ifstream fedges;

  int node1, node2;
  while (fin >> node1 >> node2) {
    graph.addEdge(node1, node2);
    graph.addEdge(node2, node1);
    int node = max(node1, node2);
    if (node > maxNodeId) {
      maxNodeId = node;
    }
  }

  kcores.resize(maxNodeId + 1);

  computeCores();
  dumpCores();
  return true;
}

void Decomposition::compare() {
  for (int i = 0; i < (int)kcores.size(); i++) {
    if (kcores[i] != original_cores[i]) {
      cerr << "[" << i << "] " <<  kcores[i] << " " << original_cores[i] << '\n';
    }
  }
}
