#ifndef GRAPH_H_INCLUDED
#define GRAPH_H_INCLUDED
#include<vector>
#include<map>
#include<stack>
#include<pthread.h>
#include<time.h>
#include<algorithm>
#include<fstream>

using namespace std;


//sort the corepair according to core number by ascent
bool cmp_by_core(const pair<int,int>& v1, const pair<int,int>& v2) {
    if(v1.second == v2.second)
        return v1.first < v2.first;
    else
        return v1.second > v2.second;
}

struct vertex{
    bool visited;
    bool removed;
    int core;
    int mcd;
    int pcd;
    int cd;
    vertex(){
        visited=removed=false;
        core=mcd=pcd=cd=0;
    }
};

class Graph
{

//--------------------------------------------------------------------------
//Variables
//--------------------------------------------------------------------------

/** Number of vertexes and edges */
public: int vertexNum;
        int edgeNum;
		vector<int> allcores;

vector<vertex> vertices;
/** Edges associated with the vertexes stored at this node
edges[u][i] means the i-th neighbor of vertex u*/
vector<vector<int> > edges;

//--------------------------------------------------------------------------
//Constructor
//--------------------------------------------------------------------------
Graph()
{
	vertexNum = 0;
	edgeNum = 0;
}

void init(int vNum, int eNum){
    vertexNum=vNum;
    edgeNum=eNum;
    for(int i=0;i<vertexNum;i++){
        vertex temp=vertex();
        vertices.push_back(temp);
        edges.push_back(vector<int>(0));
    }
}

void clear(){
    vertices.clear();
    edges.clear();
    vertexNum = 0;
	edgeNum = 0;
}

/**
 * Adds an edge to the graph.
 *
 * @param from
 *          the starting node of an edge
 * @param to
 *          the final node of an edge
 */

 bool addEdge(int from, int to)
{
    if(from == to)
        return false;
	// Update size,add another  vertex
	while(from+1 > vertexNum) {
		vertexNum ++;
		vertex temp=vertex();
        vertices.push_back(temp);
        edges.push_back(vector<int>(0));
    }
        // Search if the edge is already present; if so, stop here and return
	for (int i = 0; i < edges[from].size(); i++) {
		if (edges[from][i] == to)
			return false;
	}
	// Finally store the edge
	edges[from].push_back(to);
	edgeNum++;
	return true;
}


 bool deleteEdge(int from, int to){
    if(from > vertexNum){
        return false;
    }
    vector <int>::iterator Iter;
    for(Iter=edges[from].begin();Iter!=edges[from].end();Iter++){
      if (*Iter== to){
       edges[from].erase(Iter);
       edgeNum--;
       Iter = edges[from].begin();
       return true;
      }
    }
    //if there is no edge <from, to>
    if(Iter == edges[from].end()){
        return false;
    }
}

/**
 * Returns the number of nodes in the graph.
 */
 int GetvertexNum()
{
	return vertexNum;
}

 int GetedgeNum()
{
	return edgeNum;
}
/**
*compute core number for the graph, static method
*/
void ComputeCores(){
    int a,b;
    vector<int> bin,pos,vert,deg;
    a = b = 0;
    int n = vertexNum;
    int maxdegree = 0;
    for(int u = 0;u < n;u++){
        pos.push_back(0);
        vert.push_back(0);
        int usize = edges[u].size();
        deg.push_back(usize);
        if(usize > maxdegree){
            maxdegree = usize;
        }
    }
    //cout<<"maxdegree: "<<maxdegree<<endl;
    for(int k = 0;k <= maxdegree;k ++){
        bin.push_back(0);
    }
    for(int u = 0;u < n;u ++){
        bin[deg[u]]++;
    }
    int start  = 0;
    for(int k = 0;k <= maxdegree;k ++){
        int num = bin[k];
        bin[k] = start;
        start += num;
    }
    for(int u = 0;u < n ;u ++){
        pos[u] = bin[deg[u]];
        vert[pos[u]] = u;
        bin[deg[u]]++;
    }
    for(int d = maxdegree;d>0;d--){
        bin[d]=bin[d-1];
    }
    bin[0] = 0;
    int du,pu,pw,w;
    for(int i = 0;i < n;i++){
        int v = vert[i];
        for(int j=0;j<edges[v].size();j++){
            int u = edges[v][j];
            if(deg[u] > deg[v]){
                du = deg[u];
                pu = pos[u];
                pw = bin[du];
                w = vert[pw];
                if(u != w){
                    pos[u] = pw;
                    vert[pu] = w;
                    pos[w] = pu;
                    vert[pw] = u;
                }
                bin[du]++;
                deg[u]--;
            }
        }
    }
    vector<pair<int,int> > corepair;
    for(int i=0;i < n;i++){
        vertices[i].core = deg[i];
    }
	allcores = deg;
}

vector<int> GetAllcores(){
	return allcores;
}

// --------------------------------------------------------------------------
// main Methods for centralized algorithm
// --------------------------------------------------------------------------

/**
*delete all selected edges and conduct the centralized algorithm
*/
void Deletion(vector<pair<int,int> > allNewEdges){//delete edges	
	for(int k=0;k < allNewEdges.size();k ++){
		int a = allNewEdges[k].first;
		int b = allNewEdges[k].second;
		if(deleteEdge(a,b) && deleteEdge(b,a)){
			computeMcd();
			TravelDelete(a,b);
			resetVertex();
		}
	}		
}

/**
*insert all selected edges and conduct the centralized algorithm
*/
void Insertion(vector<pair<int,int> > allNewEdges){
	for(int k=0;k < allNewEdges.size();k ++){
		int a = allNewEdges[k].first;
		int b = allNewEdges[k].second;
		if(addEdge(a,b) && addEdge(b,a)){
			computeMcd();
			computePcd();
			TravelInsert(a,b);
			for(int i=0;i<vertexNum;i++){						
				if((vertices[i].visited)&&(!vertices[i].removed)){
					vertices[i].core++;
				}
			}
			resetVertex();
		}
	}
}

//deal with one edge insertion 
void  TravelInsert(int u1, int u2)
{
    int r=u1;
    int coreu1=vertices[u1].core;
    int coreu2=vertices[u2].core;
    if(coreu1 > coreu2){
        r=u2;
    }
    stack<int> s;
    s.push(r);
    int K=vertices[r].core;
    vertices[r].visited=true;
    vertices[r].cd=vertices[r].pcd;
    int cdr = vertices[r].cd;
    while(!s.empty()){
        int v=s.top();
        s.pop();
        int cdv = vertices[v].cd;
        if(cdv > K){
            for(int j=0;j<edges[v].size();j++){
                int w=edges[v][j];
                int corew = vertices[w].core;
                int mcdw = vertices[w].mcd;
                if(corew == K && mcdw > K &&(!vertices[w].visited)){
                    s.push(w);
                    vertices[w].visited = true;
                    vertices[w].cd += vertices[w].pcd;
                }
            }
        }
        else{
            if(!  vertices[v].removed){
                InsertRemovement(K,v);
            }
        }
    }
}

void  InsertRemovement(int k, int v )
{
    vertices[v].removed=true;
    for(int i=0;i<edges[v].size();i++){
        int w=edges[v][i];
        int corew = vertices[w].core;
        if(corew == k){
            vertices[w].cd--;
            int cdw = vertices[w].cd;
            if(cdw == k && !vertices[w].removed){
                InsertRemovement(k,w);
            }
		}
    }
}

//deal with one edge deletion
void TravelDelete(int u1, int u2)
{
    int r = u1;
    int coreu1 = vertices[u1].core;
    int coreu2 = vertices[u2].core;
    int k = coreu1;
    if(coreu1 > coreu2){ 
		r=u2;
        k=coreu2;
    }
    if(coreu1 != coreu2){
        vertices[r].visited = true;
        vertices[r].cd = vertices[r].mcd;
        int cdr = vertices[r].cd;
        if(cdr < k){
            DeleteRemove(k,r);
        }
    }
    else{
       vertices[u1].visited = true;
       vertices[u1].cd = vertices[u1].mcd;
       int cdu1 = vertices[u1].cd;
       if(cdu1 < k){
            DeleteRemove(k,u1);
       }
       vertices[u2].visited = true;
       vertices[u2].cd = vertices[u2].mcd;
       int cdu2 = vertices[u2].cd;
       if(!vertices[u2].removed && cdu2 < k ){
            DeleteRemove(k,u2);
       }
    }
}

void  DeleteRemove(int k, int v)
{
    vertices[v].removed = true;
    vertices[v].core--;
    for(int i=0;i< edges[v].size();i++){
        int w= edges[v][i];
        int corew =   vertices[w].core;
        if(corew == k){
            if(!  vertices[w].visited){
                  vertices[w].cd += vertices[w].mcd;
                  vertices[w].visited = true;
            }
              vertices[w].cd--;
            int cdw = vertices[w].cd;
            if(cdw < k && ! vertices[w].removed){
                DeleteRemove(k,w);
            }
        }
    }
}

/**
*compute MCD for all vertices
*/
void computeMcd(){
    for(int v=0;v<vertexNum;v++){
        for(int j=0;j<edges[v].size();j++)
        {
            int w=edges[v][j];
            int corev=vertices[v].core;
            int corew=vertices[w].core;
            if(corew >= corev){
                vertices[v].mcd++;
            }
        }
    }
}

/**
*compute PCD for all vertices
*/
void computePcd(){
    for(int v = 0;v < vertexNum;v++){
        for(int j=0;j<edges[v].size();j++)
        {
            int w=edges[v][j];
            int corev=vertices[v].core;
            int corew=vertices[w].core;
            int mcdw=vertices[w].mcd;
            if(corew > corev ||
               (corew == corev && mcdw> corev)){
                    vertices[v].pcd++;
            }
        }
    }
}


// --------------------------------------------------------------------------
// Methods for parallel algorithm
// --------------------------------------------------------------------------

/**
*decrease core numbers for vertices that are visited and removed
*/
void delCores(){
	for(int i=0;i<vertexNum;i++){
		if((vertices[i].visited)&&(vertices[i].removed)){
			vertices[i].core--;
			allcores[i]--;
		}
	}
}
/**
*increase core numbers for vertices that are visited but not removed
*/
void insCores(){
	for(int i=0;i<vertexNum;i++){
		if((vertices[i].visited)&&(!vertices[i].removed)){
			vertices[i].core++;
			allcores[i]--;
		}
	}
}

/**
*given a superior edge set, compute the MCD for vertices in the EXPTree
*/
 void computeInsertMcd(vector<pair<int,int> > superioredges){
	int edgenums = superioredges.size();
    map<int,bool> visited;	
    for(int i = 0;i < edgenums;i++){
        stack<int> s;
        int a = superioredges[i].first;
        int b = superioredges[i].second;
        int sizea = edges[a].size();
        int sizeb = edges[b].size();
        int corea = vertices[a].core;
        int coreb = vertices[b].core;
        int r = a;
        int corer = corea;
        int sizer = sizea;
        if(corea > coreb){
            r = b;
            corer = coreb;
            sizer = sizeb;
        }
		if(!visited[r]){
			s.push(r);
			//vertices[r].visited = true;
			visited[r] = true;
			for(int j=0;j<sizer;j++)
			{
				int w = edges[r][j];
				int corew = vertices[w].core;
				if(corew >= corer){
					vertices[r].mcd++;
				}
			}
		   while(!s.empty()){
				int v=s.top();
				s.pop();
				int sizev = edges[v].size();
				int corev = vertices[v].core;
				for(int j=0;j<sizev;j++){
					int p = edges[v][j];
					if(vertices[p].core == corev && !visited[p]){
						s.push(p);
						//vertices[p].visited = true;
						visited[p] = true;
						if(!vertices[p].mcd){
							int sizep = edges[p].size();
							int corep = vertices[p].core;
							for(int k = 0;k < sizep;k ++){
								int w = edges[p][k];
								int corew = vertices[w].core;
								if(corew >= corep){
									vertices[p].mcd++;
								}
							}
						}
					}
				}
			}
		}
	}

}

/**
*compute MCD for a vertex v
*/
 void computeMcd(int v){
    int corev = vertices[v].core;
    int sizev = edges[v].size();
    for(int j=0;j<sizev;j++)
    {
        int w = edges[v][j];
        int corew = vertices[w].core;
        if(corew >= corev){
            vertices[v].mcd++;
        }
    }
}

/**
*compute PCD for a vertex v, all MCDs are already known
*/
 void computePcd(int v){
    int corev = vertices[v].core;
    int sizev = edges[v].size();
    for(int j=0;j< sizev;j++){
        int w = edges[v][j];
        int corew = vertices[w].core;
        int mcdw = vertices[w].mcd;
        if(corew > corev || (corew == corev && mcdw > corev)){
            vertices[v].pcd++;
        }
    }
}


 void resetVertex(){
    for(int v = 0;v < vertexNum;v++){
        vertices[v].visited=false;
        vertices[v].removed=false;
        vertices[v].mcd=0;
        vertices[v].pcd=0;
        vertices[v].cd=0;
    }
}

void WriteCores(string corefile){
	ofstream fcore(corefile.data());
	vector<pair<int,int> > corepair;
	for(int i=0;i<vertexNum;i++){
		corepair.push_back(make_pair(i,vertices[i].core));
	}
	sort(corepair.begin(),corepair.end(),cmp_by_core);
	for(int i=0;i<corepair.size();i++){
		fcore<<corepair[i].first<<","<<corepair[i].second<<endl;
	}
	fcore.close();
}

};


#endif // GRAPH_H_INCLUDED
