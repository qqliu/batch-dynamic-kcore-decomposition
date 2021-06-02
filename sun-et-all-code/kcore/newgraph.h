
#ifndef NEWGRAPH_H_INCLUDED
#define NEWGRAPH_H_INCLUDED
#include<vector>
#include<map>
#include<iostream>
#include<stack>
#include<algorithm>
#include<pthread.h>
#include<time.h>
using namespace std;

typedef vector<int> VECINT;

struct params{
    int k;
    int i;
    params(){
        k = i = 0;
    }
    params(int K,int index){
        k = K;
        i = index;
    }
};

struct newvertex{
    int id;
    int core;

    //adjacent list of this vertex, store the neighbor`s index
    vector<int> adjacent;
    newvertex(){
        id = core = 0;
        adjacent.resize(0);
    }

};

class newGraph
{

/** Number of vertexes and edges */
public: int vertexNum;
        int edgeNum;
        map<int,int> index_map;
        vector<newvertex> vertices;
		
//--------------------------------------------------------------------------
//Constructor
//--------------------------------------------------------------------------
public: newGraph()
{
	vertexNum = edgeNum = 0;
}

public: void init(int vNum, int eNum){
    vertexNum=vNum;
    edgeNum=eNum;
    for(int i=0;i<vertexNum;i++){
        newvertex temp=newvertex();
        vertices.push_back(temp);
	}
}

public: void clear(){
    vertices.clear();
    index_map.clear();
    //vector<newvertex>(v).swap(v);
}
// --------------------------------------------------------------------------
// Methods
// --------------------------------------------------------------------------
/**
*
*construct the index map
**/
public: void Map_index(vector<pair<int,int> > allNewEdges){
    int from, to, index;
    from = to = index = 0;
    for(size_t i=0;i<allNewEdges.size();i++){
        from = allNewEdges[i].first;
        to = allNewEdges[i].second;
        if(index_map.find(from) == index_map.end()){
            index_map[from] = index;
            //fout<<from<<" "<<index<<endl;
            newvertex temp=newvertex();
            temp.id = from;
            vertices.push_back(temp);
            index++;
        }
        if(index_map.find(to) == index_map.end()){
            index_map[to] = index;
            //fout<<to<<" "<<index<<endl;
            newvertex temp=newvertex();
            temp.id = to;
            vertices.push_back(temp);
            index++;
        }
    }
    vertexNum = index_map.size();
	for(int i=0;i<allNewEdges.size();i++){
		int a = allNewEdges[i].first;
		int b = allNewEdges[i].second;
		addEdge(a,b);
		addEdge(b,a);
	}
}
/**
 * Adds an edge to the graph.
 *
 * @param from
 *          the starting node of an edge
 * @param to
 *          the final node of an edge
 */

public: bool addEdge(int from, int to)
{
	if(from == to )
        return false;
    //find the index of vertex from
	int index_f = index_map[from];
	int index_t = index_map[to];
	// Search if the edge is already present; if so, stop here and return
	// false.
	int len = vertices[index_f].adjacent.size();
	for (int i = 0; i < len; i++) {
		if (vertices[index_f].adjacent[i] == index_t)
			return false;
	}

	// Finally store the edge
	vertices[index_f].adjacent.push_back(index_t);
	edgeNum++;
	//cout<<"("<<from<<","<<to<<") added!"<<endl;
	return true;
}

public: bool deleteEdge(int from, int to){
    bool flag = true;
    int index_f = index_map[from];
	int index_t = index_map[to];
    if(index_f > vertexNum){
        flag = false;
    }
    //delete <from, to>
    vector <int>::iterator Iter = std::find(vertices[index_f].adjacent.begin(),vertices[index_f].adjacent.end(), index_t);
    if(Iter != vertices[index_f].adjacent.end()){
        vertices[index_f].adjacent.erase(Iter);
        edgeNum--;
        //cout<<"delete "<<from<<","<<to;
    }
    else{
        flag = false;
    }
    //delete <to, from>
    VECINT::iterator it1 = std::find(vertices[index_t].adjacent.begin(),vertices[index_t].adjacent.end(), index_f);
    if(it1 != vertices[index_t].adjacent.end()){
        vertices[index_t].adjacent.erase(it1);
        edgeNum--;
        flag = true;
        //cout<<"delete "<<to<<","<<from;
    }
    else{
        flag = false;
    }
    return flag;
}

/**
 * Returns the number of nodes in the graph.
 */
public: int GetvertexNum()
{
	return vertexNum;
}

/**
 * Returns the number of edges in the graph.
 */
public: int GetedgeNum()
{
	return edgeNum;
}

//set cores for new graph, if its a new vertex, core is default 0
public: void SetCores(vector<int> allcores){
	for(int u=0;u<vertexNum;u++){
		int idu = vertices[u].id;
		if(idu+1 > allcores.size()){
			vertices[u].core = 0;
		}else{
			vertices[u].core = allcores[idu];
		}
	}
}

//get root core numbers and 
public: vector<params> GetParams(){
	VECINT cores;
	vector<params> param;
	int index_p = 0;
	for(int u=0;u<vertexNum;u++){
		int coreu = vertices[u].core;
		int corer = coreu;
		int uSize = vertices[u].adjacent.size();
		//vertex u has a neighbor with larger core number
		bool flag = false;
		for(int i=0;i<uSize;i++){
			int v = vertices[u].adjacent[i];
			int corev = vertices[v].core;
			if(corev >= coreu){
				flag = true;
				break;
			}
		}
		if(flag){
			VECINT::iterator iter=std::find(cores.begin(),cores.end(),corer);
			// if the core is not in rootcores, add it
			if(iter == cores.end()){
				cores.push_back(corer);
				param.push_back(params(corer,index_p));
				index_p++;
			}
		}
	}
	return param;
}
				
public: int GetVertexId(int u)
{
	return vertices[u].id;
}

int GetVertexCore(int u)
{
	return vertices[u].core;
}

vector<int> GetVertexAdj(int u)
{
	return vertices[u].adjacent;
}


};

#endif // GRAPH_H_INCLUDED
