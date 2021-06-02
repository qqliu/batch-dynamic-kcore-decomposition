#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <functional>
#include <utility>
#include <set>
#include "graph.h"

using namespace std;

Graph::Graph()
{
    maxNodeId = 0;
}

bool Graph::read(const char *filename)
{
    ifstream fin(filename);
    if(fin.fail()){
        return false;
    }
    int node1,node2;
    while(fin >> node1 >> node2) {
        matrix[node1][node2] = true;
        matrix[node2][node1] = true;
        maxNodeId = max(max(maxNodeId,node1),node2);
    }
    nodeInfos.resize(maxNodeId + 1);
    fill(nodeInfos.begin(),nodeInfos.end(),0);
    fin.close();
}

void Graph::kcore_decomposition()
{
    int iter = 0;
    int total = 0;
    set<Edge> allEdges;

    while(!inserted_graph.empty()) {
        vector<Edge> insertEdges = computeInsertEdgeSet();

        // insert Eins into G
        for(int i = 0;i < (int)insertEdges.size();i++) {


            matrix[insertEdges[i].first][insertEdges[i].second] = true;
            matrix[insertEdges[i].second][insertEdges[i].first] = true;
        }

        // delete E ins from E'
        for(int i = 0;i < (int)insertEdges.size();i++) {
            int node1 = insertEdges[i].first,node2 = insertEdges[i].second;


            if(node1 < node2) {
                allEdges.insert(Edge(node1,node2));
            } else {
                allEdges.insert(Edge(node2,node1));
            }




            auto it = find(inserted_graph[node1].begin(),inserted_graph[node1].end(),node2);
            if(it != inserted_graph[node1].end()){
                inserted_graph[node1].erase(it);
            }

            it = remove(inserted_graph[node2].begin(),inserted_graph[node2].end(),node1);
            if(it != inserted_graph[node2].end()) {
                inserted_graph[node2].erase(it);
            }

            if(inserted_graph[node1].empty()) {
                inserted_graph.erase(node1);
            }
            if(inserted_graph[node2].empty()) {
                inserted_graph.erase(node2);
            }
        }
        cout << "第" << ++iter << " 次迭代中插入了" << insertEdges.size() << "条边" << endl;
        total += insertEdges.size();
    }
    cout << "total = " << total << endl;
    cout << "allEdges.size() = " << allEdges.size() << endl;
}

vector<Edge> Graph::computeInsertEdgeSet()
{
    // nodeInfos:
    // -2
    // -1 没有连向superior edge
    // 0  已经连向一条superior edge
    // >0 pre-core
    vector<Edge> inserted;
    vector<int> cps = centerPoints();
    fill(nodeInfos.begin(),nodeInfos.end(),-1);
    for(int i = 0;i < (int)cps.size();i++) {
        nodeInfos[cps[i]] = -2;
    }

    for(int i = 0;i < (int)cps.size();i++) {
        int v = cps[i];
        nodeInfos[v] = precore(v);
        vector<int> newNeighbors = inserted_graph[v];

        for(auto u : newNeighbors) {
            inserted.push_back(make_pair(v,u));

            if(kcores[u] - nodeInfos[v] <= 1) {
                nodeInfos[u] = 0;
            }
            vector<int> insertEdges = inserted_graph[u];
            for(auto w : insertEdges) {
                if(nodeInfos[w] == -2)continue;

                if(kcores[w] < kcores[u] && nodeInfos[w] == -1) {
                    inserted.push_back(make_pair(u,w));
                    nodeInfos[w] = 0;
                }
                if(nodeInfos[u] == -1 && kcores[u] <= kcores[w]) {
                    inserted.push_back(make_pair(u,w));
                    nodeInfos[u] = 0;
                }
            }
        }
    }

    return inserted;
}

int Graph::precore(int node)
{
    vector<int> neighbors = neighborsWithInserted(node);
    unordered_map<int,int> core2count;
    int pc = 0;
    for(int i = 0;i < (int)neighbors.size();i++) {
        int k = kcores[neighbors[i]];
        core2count[k]++;
        if(core2count[k] >= k && k > pc) {
            pc = k;
        }
    }
    return pc;
}

vector<int> Graph::neighborsWithInserted(int node)
{
    vector<int> neighbors;
    auto nit = matrix.find(node);
    if(nit != matrix.end()) {
        for(auto it = nit->second.begin();it != nit->second.end();++it) {
            neighbors.push_back(it->first);
        }
    }
    auto iit = inserted_graph.find(node);
    if(iit != inserted_graph.end()) {
        for(auto it = iit->second.begin();it != iit->second.end();++it) {
            neighbors.push_back(*it);
        }
    }
    return neighbors;
}

std::map<int,std::vector<int>>::iterator Graph::maxIncreaseDegree()
{
    if(inserted_graph.empty()) return inserted_graph.end();

    auto maxit = inserted_graph.end();
    int maxSize = -1;

    for(auto it = inserted_graph.begin();it != inserted_graph.end();++it) {
        if(nodeInfos[it->first] == 1 &&
           int(it->second.size()) > maxSize) {
            maxit = it;
            maxSize = it->second.size();
        }
    }

    return maxit;
};

vector<int> Graph::centerPoints()
{
    fill(nodeInfos.begin(),nodeInfos.end(),1);
    vector<int> points;
    int total = 0;

    do {
        auto it = maxIncreaseDegree();

        if(it == inserted_graph.end())break;
        points.push_back(it->first);
        nodeInfos[it->first] = 0;
        ++total;
        vector<int>twoHopNodes = twoHopNeighbors(it);
        total += twoHopNodes.size();
    }while(total < maxNodeId + 1);

    return points;
}

vector<int> Graph::twoHopNeighbors(std::map<int,std::vector<int>>::iterator nodeit)
{
    vector<int> neighbors;

    // get one-hop neighbors of node `nodeit->first` in new graph
    vector<int> oneHopNeighbors;
    auto neighborsIterator = matrix.find(nodeit->first);
    if(neighborsIterator == matrix.end()) {
        return neighbors;
    }
    for(auto it = neighborsIterator->second.begin();
            it != neighborsIterator->second.end();++it) {
        if(nodeInfos[it->first]) {
            oneHopNeighbors.push_back(it->first);
            nodeInfos[it->first] = 0;
        }
    }
    for(auto node : nodeit->second) {
        if(nodeInfos[node]) {
            oneHopNeighbors.push_back(node);
            nodeInfos[node] = 0;
        }
    }

    // calculate two-hop neighbors
    for(int i = 0;i < (int)oneHopNeighbors.size();i++) {
        neighborsIterator = matrix.find(oneHopNeighbors[i]);
        for(auto it = neighborsIterator->second.begin();
                it != neighborsIterator->second.end();++it) {
            if(nodeInfos[it->first]) {
                neighbors.push_back(it->first);
                nodeInfos[it->first] = 0;
            }
        }

        auto insertedNeighbor = inserted_graph.find(oneHopNeighbors[i]);
        if(insertedNeighbor != inserted_graph.end()) {
            for(auto it = insertedNeighbor->second.begin();
                    it != insertedNeighbor->second.end();++it) {
                if(nodeInfos[*it]) {
                    neighbors.push_back(*it);
                    nodeInfos[*it] = 0;
                }
            }
        }
        neighbors.push_back(oneHopNeighbors[i]);
    }

    return neighbors;
}

void Graph::insert_edges(EdgeSeq &edges)
{
    fill(inserted_size.begin(),inserted_size.end(),0);

    for(int i = 0;i < edges.size();i++) {
        int node1 = edges[i].first,node2 = edges[i].second;
        inserted_graph[node1].push_back(node2);
        inserted_graph[node2].push_back(node1);
        //inserted_size[node1]++;
        //inserted_size[node2]++;

        remove_edge(edges[i]);
    }
    kcores = cores();
    //sort(inserted_size.begin(),inserted_size.end(),greater<int>());
}

bool Graph::insert_edges(const char *edges_filename)
{
    ifstream fin(edges_filename);
    if(fin.fail()) {
        return false;
    }

    vector<Edge> edges;
    int node1,node2;
    while(fin >> node1 >> node2) {
        Edge edge = make_pair(node1,node2);
        edges.push_back(edge);
        //remove_edge(edge);
    }
    cerr << "一共插入" << edges.size() << "条新边" << endl;
    insert_edges(edges);
}

void Graph::remove_edge(Edge edge)
{
    // <node1,node2>
    int node1 = edge.first,node2 = edge.second;
    auto it1 = matrix.find(node1);
    if(it1 != matrix.end()) {
        auto it2 = it1->second.find(node2);
        if(it2 != it1->second.end()) {
            it1->second.erase(it2);
        }
    }


    it1 = matrix.find(node2);
    if(it1 != matrix.end()) {
        auto it2 = it1->second.find(node1);
        if(it2 != it1->second.end()) {
            it1->second.erase(it2);
        }
    }
}

unordered_map<int,int> Graph::cores()
{
    int d, md, i, start, num;
    int v, u, w, du, pu, pw;

    unordered_map<int,unordered_map<int,bool> >& g = matrix;
    vector<int> bin;
    unordered_map<int,int> pos;
    unordered_map<int,int> deg;
    map<int,int> vert;

    md = 0; // md 是图中最大的度
    // 计算每个节点初始时刻的度
    for (unordered_map<int,unordered_map<int,bool> >::iterator it = g.begin();it != g.end();++it){
        v = it->first;
        d = it->second.size();
        deg[v] = d;
        if (d > md){
            md = d;
        }
    }


    bin.resize(md + 1);
    fill(bin.begin(),bin.end(),0);
    // bin[i] 中存储度为 i 的节点个数
    for (unordered_map<int,unordered_map<int,bool> >::iterator it = g.begin();it != g.end();++it){
        v = it->first;
        bin[deg[v]]++;
    }

    start = 0;
    for (d = 0; d <= md; d++){
        num = bin[d];
        bin[d] = start;
        start += num;
    }
    for (unordered_map<int,unordered_map<int,bool> >::iterator it = g.begin();it != g.end();++it){
        v = it->first;
        pos[v] = bin[deg[v]];
        vert[pos[v]] = v;
        bin[deg[v]]++;
    }
    for(d = md;d >= 1;d--){
        bin[d] = bin[d - 1];
    }
    bin[0] = 0;
    for(map<int,int>::iterator it = vert.begin();it != vert.end();++it){
        v = it->second;
        unordered_map<int,unordered_map<int,bool> >::iterator inducedit = g.find(v);
        if(inducedit != g.end()){
            unordered_map<int,bool>& neighbors = inducedit->second;
            for(unordered_map<int,bool>::iterator neighbor = neighbors.begin();neighbor != neighbors.end();++neighbor){
                u = neighbor->first;
                if(deg[u] > deg[v]){
                    du = deg[u];pu = pos[u];
                    pw = bin[du];w = vert[pw];
                    if (u != w){
                        pos[u] = pw;vert[pu] = w;
                        pos[w] = pu;vert[pw] = u;
                    }
                    bin[du]++;deg[u]--;
                }
            }
        }
    }

    return deg;
}

void Graph::saveCores(const char *filename) const
{

}