#!/bin/bash

# let "bs  = 1000"
# let "s = $1 * $bs"
# let "e = $s + $bs"

# -blocksize 100
# ./Triangle -s -nb 10 -n 2048 -bo 0 ../../../inputs/empty ../../../inputs/rMatGraph_J_2K_10K

# ./Triangle -s -static -nb 10 -n 2048 -bo 0 ../../../inputs/empty ../../../inputs/rMatGraph_J_2K_10K
# 
# ./Triangle -s -w 2 -nb 10 -n 2048 -bo 0 ../../../inputs/empty ../../../inputs/rMatGraph_J_2K_10K
# 


# ./Triangle -s -w -sg -trict 3 -nb  1 -n 7 -bo 0 ../../../inputs/small/triangles4.txt ../../../inputs/small/update4.txt

# ./Triangle -s -w -shuffle -nb  1 -n 7 -bo 0 ../../../inputs/small/triangles4.txt ../../../inputs/small/empty

# ./Triangle -s -blocksize 5 -n 128 -bo $1 ../../../inputs/empty ../../../inputs/rMatEdgesShuffle

# ./Triangle -s -w -static -shuffle -nb  1 -n 7 -bo 0 ../../../inputs/small/triangles4.txt ../../../inputs/small/empty

# ./Triangle -s -w 1 -nb 10 -n 128 -bo 0 ../../../inputs/empty ../../../inputs/rMatEdgesShuffle

# ./Triangle -s -shuffle -static -nb 10 -n 1138499 -bo $1 ../../../inputs/youtube_sym.adj ../../../inputs/small/empty


./Triangle -s -shuffle -batchsize 1000 -n 1138499 -blocksize 16 -bo 0 ../../../inputs/youtube_sym.adj ../../../inputs/small/empty >> tmp

# ./Triangle -s -shuffle -nb 10 -n 1138499  -blocksize 16 -bo $1 -bp 10 ../../../inputs/youtube_sym.adj ../../../inputs/small/empty

# ./Triangle -s -w 2 -batchsize 500 -nb 10 -n 128 -bo 10 ../../../inputs/empty ../../../inputs/rMatEdgesShuffle

# ./Triangle -s -w 1 -blocksize 16 -batchsize 1000 -n 2048 -partial 5000 ../../../inputs/empty ../../../inputs/rMatGraph_J_2K_10K


# ./exp.sh 0 70
# ./exp.sh 71 140
