# Build  
  `make` in kcore fold.

# Run
````
    ./kcore LJ.txt LJ_10W.txt 16 
````
  LJ.txt: graph file  
  LJ_10W.txt: updated edges file  
  16: numbers of thread  

# Result
````
  Delete  193.976 3  
  Insert  240.639 3  
````
  Delete/Insert: deletion/insertion Results  
  192.973/240.639: total run time (unit: ms)  
  3/3: iteration times
