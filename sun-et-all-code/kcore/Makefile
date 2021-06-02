# git root directory
CXX = g++
CXXFLAGS = -std=c++0x -O3

all:kcore  paracom

kcore:timer.cpp UndirectedGraph.cpp kcore.cpp main.cpp
	$(CXX) $(CXXFLAGS) -o $@ timer.cpp UndirectedGraph.cpp kcore.cpp main.cpp -lpthread -latomic
calc:timer.cpp UndirectedGraph.cpp kcore.cpp calc.cpp
	$(CXX) $(CXXFLAGS) -o $@ timer.cpp UndirectedGraph.cpp kcore.cpp calc.cpp -lpthread -latomic
paracom:ParaCoM.cpp
	$(CXX) $(CXXFLAGS) -o $@ ParaCoM.cpp -lpthread
