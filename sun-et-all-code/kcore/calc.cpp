#include <ctime>
#include <iostream>
#include "kcore.h"

using namespace std;

int main(int argc,char *argv[])
{
    if(argc != 2) {
        fprintf(stderr,"Usage: %s <graph-file>\n",argv[0]);
        return -1;
    }
    if(!Decomposition().calcCores(argv[1])) {
        cerr << "Open " << argv[1] << " failed\n";
    } 
    
    return 0;
}