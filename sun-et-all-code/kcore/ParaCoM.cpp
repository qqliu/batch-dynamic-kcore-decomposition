#include<iostream>
#include<stdio.h>
#include <stdlib.h>
#include <cstring>
#include <dirent.h>
#include <sys/time.h>
#include"ParaCoM.h"
#include"newgraph.h"

using namespace std;

#define CLOCK_PER_MS 1000

void FindThreads(vector<params> param);

void DeleteThread();

void InsertThread();

void *GetKSuperiorEdges(void *param);

void *TravelInsert(void *sameCoreEdges);

void InsertRemove(int k, int r);

void *TravelDelete(void *sameCoreEdges);

void DeleteRemove(int k, int r);

Graph graph;
newGraph newgraph;
vector<vector<pair<int, int> > > superiorEdges;

//create finding threads to find superior edges
void FindThreads(vector<params> param) {
    int num_threads = param.size();
    pthread_t tids[num_threads];
    for (int i = 0; i < num_threads; i++) {
        int ret = pthread_create(&tids[i], NULL, GetKSuperiorEdges, &param[i]);
        if (ret != 0) {
            cout << "pthread_create error:error_code=" << ret << endl;
        } else {
            //cout <<tids[i]<<" thread_create"<<endl;
        }
    }
    int edges = 0;
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(tids[i], NULL);//(void**)&superiorEdges[i]);

        edges += superiorEdges[i].size();
    }
}

//delete superior edges from new graph and insert them into original graph
void InsertEdgesIntoGraph() {
    int superiorEdgesSize = superiorEdges.size();

    int total = 0;
    for (int i = 0; i < superiorEdgesSize; i++) {
        //int ksize = superiorEdges[i].size();
        vector<pair<int, int> > kedges = superiorEdges[i];
        int ksize = kedges.size();
        total += ksize;
        for (int k = 0; k < ksize; k++) {
            int u = kedges[k].first;
            int v = kedges[k].second;
            if (!newgraph.deleteEdge(u, v)) {
                cout << "delete " << u << "," << v << endl;
                cout << "wrong delete edges in new graph!\n";
                return;
            }
            if (!graph.addEdge(u, v) || !graph.addEdge(v, u)) {
                cout << "insert " << u << "," << v << endl;
                cout << "wrong inserting edges in original graph!\n";
                return;
            }
        }
    }
}

//delete superior edges from new graph and original graph
void DeleteEdgesFromGraph() {
    int total = 0;

    for (int i = 0; i < superiorEdges.size(); i++) {
        vector<pair<int, int> > kedges = superiorEdges[i];
        int ksize = kedges.size();
        total += ksize;

        for (int k = 0; k < ksize; k++) {
            int u = kedges[k].first;
            int v = kedges[k].second;
            if (!newgraph.deleteEdge(u, v)) {
                cout << "delete " << u << "," << v << endl;
                cout << "wrong delete edges in new graph!\n";
                return;
            }
            if (!graph.deleteEdge(u, v) || !graph.deleteEdge(v, u)) {
                cout << "delete " << u << "," << v << endl;
                cout << "wrong delete edges in original graph!\n";
                return;
            }
        }
    }
    std::cout << "size: " << total << std::endl;
}

//start the deletion thread
void DeleteThread() {
    int num_threads = superiorEdges.size();
    pthread_t threads[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        int ret = pthread_create(&threads[i], NULL, TravelDelete, &superiorEdges[i]);
        if (ret != 0) {
            cout << "pthread_create error:error_code=" << ret << endl;
        } else {
            //cout <<tids[i]<<" thread_create"<<endl;
        }
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], NULL);
    }
}

//start the insertion thread
void InsertThread() {
    int num_threads = superiorEdges.size();
    pthread_t threads[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        int ret = pthread_create(&threads[i], NULL, TravelInsert, &superiorEdges[i]);
        if (ret != 0) {
            cout << "pthread_create error:error_code=" << ret << endl;
        } else {
            //cout <<tids[i]<<" thread_create"<<endl;
        }
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(threads[i], NULL);
    }

}

/**
* given a core number k, find the k-superiorEdges in newgraph(the inserted/deleted edges) and add to superiorEdges
*/
void *GetKSuperiorEdges(void *param) {
    params par = *(params *) param;
    int k = par.k;
    int index = par.i;
    VECINT rootvertex;
    vector<bool> IsPicked;
    int verNumber = 0;
    for (int i = 0; i < newgraph.vertexNum; i++) {
        if (newgraph.vertices[i].core == k) {
            rootvertex.push_back(newgraph.vertices[i].id);
        }
        IsPicked.push_back(false);
    }
    verNumber = rootvertex.size();
    for (int i = 0; i < verNumber; i++) {
        int idu = rootvertex[i];
        int index_u = newgraph.index_map[idu];
        if (!IsPicked[index_u]) {
            vector<int> adjacent = newgraph.vertices[index_u].adjacent;
            int uSize = adjacent.size();
            for (int j = 0; j < uSize; j++) {
                int index_v = adjacent[j];
                int idv = newgraph.vertices[index_v].id;
                int corev = newgraph.vertices[index_v].core;
                if (corev >= k && !IsPicked[index_v]) {
                    IsPicked[index_u] = true;
                    //if corev == k, set v picked
                    if (corev == k) {
                        IsPicked[index_v] = true;
                    }
                    superiorEdges[index].push_back(make_pair(idu, idv));
                    break;
                }
            }
        }
    }
}

/**
*after inserting a superior edge set, search vertices whose core numbers are same as the roots and
*find vertices that will change cores
*executed by a thread
*@param sameCoreEdges:the roots of new inserted edges, for each root conducts the InsertOneEdge algorithm
*
**/
void *TravelInsert(void *sameCoreEdges) {
    vector<pair<int, int> > rootEdges = *(vector<pair<int, int> > *) sameCoreEdges;
    stack<int> s;
    int k = 0;
    int r = 0;
    int rootsize = rootEdges.size();
    graph.computeInsertMcd(rootEdges);
    for (int i = 0; i < rootsize; i++) {
        int u = rootEdges[i].first;
        int v = rootEdges[i].second;
        int coreu = graph.vertices[u].core;
        int corev = graph.vertices[v].core;
        r = u;
        k = coreu;
        if (coreu > corev) {
            r = v;
            k = corev;
        }
        //only vertex that was not updated need to be travel again
        if (!(graph.vertices[r].visited && !graph.vertices[r].removed)) {
            graph.vertices[r].visited = true;
            if (!graph.vertices[r].pcd) {
                graph.computePcd(r);
            }
            if (graph.vertices[r].cd >= 0)
                graph.vertices[r].cd = graph.vertices[r].pcd;//pcdpair[r];//
                //if r has been attacked before, its cd should be updated not initialized
            else
                graph.vertices[r].cd += graph.vertices[r].pcd;//pcdpair[r];//
            s.push(r);
            int cdr = graph.vertices[r].cd;
            while (!s.empty()) {
                int v = s.top();
                s.pop();
                int cdv = graph.vertices[v].cd;
                if (cdv > k) {
                    int sizev = graph.edges[v].size();
                    for (int j = 0; j < sizev; j++) {
                        int w = graph.edges[v][j];
                        int corew = graph.vertices[w].core;
                        /*if(!graph.vertices[w].mcd){
                            graph.computeMcd(w);
                        }*/
                        int mcdw = graph.vertices[w].mcd;//mcdpair[w];//
                        if (corew == k && mcdw > k && (!graph.vertices[w].visited)) {
                            s.push(w);
                            graph.vertices[w].visited = true;
                            if (!graph.vertices[w].pcd) {
                                graph.computePcd(w);
                            }
                            graph.vertices[w].cd += graph.vertices[w].pcd;//pcdpair[w];//
                        }
                    }
                } else if (!graph.vertices[v].removed) {
                    InsertRemove(k, v);
                }
            }
        }
    }
}

void InsertRemove(int k, int r) {
    stack<int> s;
    s.push(r);
    graph.vertices[r].removed = true;
    while (!s.empty()) {
        int v = s.top();
        s.pop();
        for (int i = 0; i < graph.edges[v].size(); i++) {
            int w = graph.edges[v][i];
            int corew = graph.vertices[w].core;
            if (corew == k) {
                graph.vertices[w].cd--;
                int cdw = graph.vertices[w].cd;
                if (cdw == k && !graph.vertices[w].removed) {
                    graph.vertices[w].removed = true;
                    s.push(w);
                }
            }
        }
    }
}

/**
*after deleting a superior edge set, search vertices whose core numbers are same as the roots and
*find vertices that will change cores
*executed by a thread
*@param sameCoreEdges:the roots of deleted edges, for each root conducts the DeleteOneEdge algorithm
*
**/
void *TravelDelete(void *sameCoreEdges) {
    vector<pair<int, int> > rootEdges = *(vector<pair<int, int> > *) sameCoreEdges;
    int rootsize = rootEdges.size();
    for (int i = 0; i < rootsize; i++) {
        int u1 = rootEdges[i].first;
        int u2 = rootEdges[i].second;
        int coreu1 = graph.vertices[u1].core;
        int coreu2 = graph.vertices[u2].core;
        int r = u1;
        int k = coreu1;
        if (coreu1 > coreu2) {
            r = u2;
            k = coreu2;
        }
        if (coreu1 != coreu2) {
            if (!graph.vertices[r].visited) {
                graph.vertices[r].visited = true;
                if (!graph.vertices[r].mcd) {
                    graph.computeMcd(r);
                }
                graph.vertices[r].cd = graph.vertices[r].mcd;
            }
            if (!graph.vertices[r].removed) {
                int cdr = graph.vertices[r].cd;
                if (cdr < k) {
                    DeleteRemove(k, r);
                }
            }
        } else {
            if (!graph.vertices[u1].visited) {
                graph.vertices[u1].visited = true;
                if (!graph.vertices[u1].mcd) {
                    graph.computeMcd(u1);
                }
                graph.vertices[u1].cd = graph.vertices[u1].mcd;
            }
            if (!graph.vertices[u1].removed) {
                int cdu1 = graph.vertices[u1].cd;
                if (cdu1 < k) {
                    DeleteRemove(k, u1);
                }
            }
            if (!graph.vertices[u2].visited) {
                graph.vertices[u2].visited = true;
                if (!graph.vertices[u2].mcd) {
                    graph.computeMcd(u2);
                }
                graph.vertices[u2].cd = graph.vertices[u2].mcd;
            }
            if (!graph.vertices[u2].removed) {
                int cdu2 = graph.vertices[u2].cd;
                if (cdu2 < k) {
                    DeleteRemove(k, u2);
                }
            }
        }
    }

}

void DeleteRemove(int k, int r) {
    stack<int> s;
    s.push(r);
    graph.vertices[r].removed = true;
    while (!s.empty()) {
        int v = s.top();
        s.pop();
        for (int i = 0; i < graph.edges[v].size(); i++) {
            int w = graph.edges[v][i];
            int corew = graph.vertices[w].core;
            if (corew == k) {
                if (!graph.vertices[w].visited) {
                    if (!graph.vertices[w].mcd) {
                        graph.computeMcd(w);
                    }
                    int mcdw = graph.vertices[w].mcd;
                    graph.vertices[w].cd += mcdw;
                    graph.vertices[w].visited = true;
                }
                graph.vertices[w].cd--;
                int cdw = graph.vertices[w].cd;
                if (cdw < k && !graph.vertices[w].removed) {
                    graph.vertices[w].removed = true;
                    s.push(w);
                }
            }
        }
    }
}


int main(int argc, char *argv[]) {
    if (argc != 4) {
        cout << "Usage: -p(parallel)/-c(centralized) graph_filename edge_filename" << endl;
        return 0;
    }
    string fname = argv[2];
    string edge_file = argv[3];
    string graphfile = fname + ".txt";
    string edgefile = edge_file + ".txt";
    string delcorefile = fname + "_core_del.txt";
    string inscorefile = fname + "_core_ins.txt";
    ifstream fingraph(graphfile.data());
    ifstream finedge(edgefile.data());
    vector<pair<int, int> > allNewEdges;
    vector<int> allcores;
    //clock_t startall,endall,start,end;
    struct timeval t_start, t_end, f_start, f_end;
    long dur;
    {//open files, read graph and edge file and compute core
        int a, b;
        if (!fingraph) {
            cout << "Error opening " << graphfile << " for input" << endl;
            return 0;
        }
        while (fingraph >> a >> b) {
            graph.addEdge(a, b);
            graph.addEdge(b, a);
        }
        graph.ComputeCores();
        allcores = graph.GetAllcores();
        //get new graph and set cores
        if (!finedge) {
            cout << "Error opening " << edgefile << " for input" << endl;
            return 0;
        }
        while (finedge >> a >> b) {
            allNewEdges.push_back(make_pair(a, b));
        }
        newgraph.Map_index(allNewEdges);
        newgraph.SetCores(allcores);
    }

    //parallel algorithm
    if (strcmp(argv[1], "-p") == 0) { //delete the edges first and then insert them back

//        string output = fname + "_time_p.txt";
//        ofstream fout(output.data(), ios::app);
//        fout << fname << "\t" << edge_file << "\t";
        std::cout << fname << '\t' << edge_file << '\t';


        int newEdgeNum = newgraph.GetedgeNum();
        //cout<<newEdgeNum<<endl;
        {//delete the edges
            gettimeofday(&t_start, NULL);
            long startall_t = ((long) t_start.tv_sec) * 1000000 + (long) t_start.tv_usec;
            long findTime = 0.0;
            int round = 0;
            while (newEdgeNum) {
                round++;
                vector<params> param;
                gettimeofday(&f_start, NULL);
                long start_t = ((long) f_start.tv_sec) * 1000000 + (long) f_start.tv_usec;
                param = newgraph.GetParams();
                superiorEdges.resize(param.size());
                FindThreads(param);
                gettimeofday(&f_end, NULL);
                long end_t = ((long) f_end.tv_sec) * 1000000 + (long) f_end.tv_usec;
                dur = (end_t - start_t);
                findTime += dur;
                DeleteEdgesFromGraph();

                {//deleting threads and update cores for new graph
                    DeleteThread();
                    graph.delCores();
                    graph.resetVertex();
                    allcores = graph.GetAllcores();
                    newgraph.SetCores(allcores);
                    newEdgeNum = newgraph.GetedgeNum();
                    superiorEdges.clear();
                }
            }

            gettimeofday(&t_end, NULL);
            long endall_t = ((long) t_end.tv_sec) * 1000000 + (long) t_end.tv_usec;
            dur = endall_t - startall_t;

            //fout << dur << "\t" << findTime << "\t" << round << "\t";
            std::cout << dur << "\t" << findTime << "\t" << round << "\t";

            //write to new core file
            graph.WriteCores(delcorefile);

        }

        {//reconstruct new graph and set cores
            graph.ComputeCores();
            allcores = graph.GetAllcores();
            newgraph.clear();
            newgraph.Map_index(allNewEdges);
            newgraph.SetCores(allcores);
        }

        {//insertion
            gettimeofday(&t_start, NULL);
            long startall = ((long) t_start.tv_sec) * 1000000 + (long) t_start.tv_usec;
            newEdgeNum = newgraph.GetedgeNum();
            long findTime = 0;
            int round = 0;
            while (newEdgeNum) {
                round++;
                vector<params> param;
                gettimeofday(&f_start, NULL);
                long start = ((long) f_start.tv_sec) * 1000000 + (long) f_start.tv_usec;
                param = newgraph.GetParams();
                superiorEdges.resize(param.size());
                FindThreads(param);
                gettimeofday(&f_end, NULL);
                long end = ((long) f_end.tv_sec) * 1000000 + (long) f_end.tv_usec;
                dur = end - start;
                findTime += dur;
                InsertEdgesIntoGraph();

                {//insertion threads and update cores for new graph
                    InsertThread();
                    graph.insCores();
                    graph.resetVertex();
                    allcores = graph.GetAllcores();
                    newgraph.SetCores(allcores);
                    newEdgeNum = newgraph.GetedgeNum();
                    superiorEdges.clear();
                }
            }
            gettimeofday(&t_end, NULL);
            long endall = ((long) t_end.tv_sec) * 1000000 + (long) t_end.tv_usec;
            dur = endall - startall;
            //fout << dur << "\t" << findTime << "\t" << round << "\n";
            std::cout << dur << "\t" << findTime << "\t" << round << "\n";
            //write to new core file
            graph.WriteCores(inscorefile);
        }
        //fout.close();
    } else if (strcmp(argv[1], "-c") == 0) {
        string output = fname + "_time_c.txt";
        ofstream fout(output.data(), ios::app);
        fout << fname << "\t" << edge_file << "\t";
        {//delete the selected edges
            timeval t_start, t_end;
            gettimeofday(&t_start,NULL);
            graph.Deletion(allNewEdges);
            gettimeofday(&t_end,NULL);
            long endall = ((long) t_end.tv_sec) * 1000000 + (long) t_end.tv_usec;
            long startall = ((long)t_start.tv_sec) * 1000000 + (long) t_start.tv_usec;
            long dur = endall - startall;
            fout << dur << "\t";
            //write to core file
            graph.WriteCores(delcorefile);
        }
        {//insert edges back
            timeval t_start, t_end;
            gettimeofday(&t_start,NULL);
            graph.Insertion(allNewEdges);
            gettimeofday(&t_end,NULL);
            long endall = ((long) t_end.tv_sec) * 1000000 + (long) t_end.tv_usec;
            long startall = ((long)t_start.tv_sec) * 1000000 + (long) t_start.tv_usec;
            long dur = endall - startall;
            fout << dur << "\n";
            //write to new core file
            graph.WriteCores(inscorefile);
        }
        fout.close();
    }

    {//close the file
        fingraph.close();
        finedge.close();

    }
    cout << fname << " finished!" << endl;
    return 0;
}
