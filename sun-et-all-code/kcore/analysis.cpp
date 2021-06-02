#include <iostream>
#include "kcore.h"

int main(int argc,char *argv[])
{
    if(argc != 2) {
        fprintf(stderr,"Usage: %s <graph-file>\n",argv[0]);
        return -1;
    }

    Decomposition().coreAnalysis(argv[1]);
    return 0;
}