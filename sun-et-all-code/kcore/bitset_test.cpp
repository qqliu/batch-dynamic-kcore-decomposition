#include <cassert>
#include "bitset.h"

int main(int argc,char *argv[])
{
    Bitset bitset;
    for(int i = 0;i < 100;i++) {
        if(i & 1) {
            bitset.set(i);
        }
    }
    
    for(int i = 0;i < 100;i++) {
        if(i & 1) {
            assert(bitset.get(i));
        } else {
            assert(!bitset.get(i));
        }
    }
    return 0;
}