#include <ctime>
#include <iostream>
#include "kcore.h"

using namespace std;

int main(int argc,char* argv[])
{

    if(argc != 3) {
        fprintf(stderr,"Usage: %s <graph-file> <edge-file>\n",argv[0]);
        return -1;
    }
    Decomposition dec;
    if(!dec.Load(argv[1],argv[2])) {
        fprintf(stderr,"Load file failed\n");
        return -1;
    }

    dec.SerialDelete();
    dec.SerialInsert();
}