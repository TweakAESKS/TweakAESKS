CXX = g++
CXXFLAGS = -O3 -march=native -Wall -Wextra -std=c++11

IFLAGS = 
LFLAGS =

search : search.o CPsolve.o allPathsNoEquations.o DicStateKey.o ShortestPath.o Subset.o
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o search search.o CPsolve.o allPathsNoEquations.o DicStateKey.o ShortestPath.o Subset.o $(LFLAGS)

search.o : search.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o search.o -c search.cpp $(LFLAGS)

CPsolve.o : CPsolve.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o CPsolve.o -c CPsolve.cpp $(LFLAGS)

allPathsNoEquations.o : allPathsNoEquations.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o allPathsNoEquations.o -c allPathsNoEquations.cpp $(LFLAGS)

DicStateKey.o : DicStateKey.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o DicStateKey.o -c DicStateKey.cpp $(LFLAGS)

ShortestPath.o : ShortestPath.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o ShortestPath.o -c ShortestPath.cpp $(LFLAGS)

Subset.o : Subset.cpp
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o Subset.o -c Subset.cpp $(LFLAGS)




clean :
	rm -rf *.o