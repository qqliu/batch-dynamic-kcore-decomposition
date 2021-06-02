#include <ctime>
#include <iostream>
#include <vector>
#include "kcore.h"

using namespace std;

int main(int argc,char* argv[])
{
    if(argc != 6) {
        fprintf(stderr,"Usage: %s <graph-file> <select-core-lower>  <select-core-upper> <node-prob> <edge-prob>\n",argv[0]);
        return -1;
    }

    srand(static_cast<unsigned>(time(0)));

    vector<int> select_cores;
    Decomposition dec;

    double node_prob = stod(argv[argc - 2]);
    double edge_prob = stod(argv[argc - 1]);

    if(!dec.select(argv[1],stoi(argv[2]),stoi(argv[3]),node_prob,edge_prob)) {
        fprintf(stderr,"Load file failed\n");
        return -1;
    }
    return 0;
}