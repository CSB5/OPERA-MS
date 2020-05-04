CC:=g++
CFLAGS:=-c -Wall -g -std=c++11
LDFLAGS:=-pthread -std=c++11
INC:=-I ./
SOURCES:=AbundanceLoader.cpp EManager.cpp EucDist.cpp fastaReader.cpp \
	KendallDist.cpp kmerMap.cpp logger.cpp main.cpp ManhattanDist.cpp \
	MinManhattanDist.cpp Profiler.cpp quickSort.cpp SpearmanDist.cpp \
	AbstractDist.cpp NormalDistribution.cpp
OBJECTS:=$(SOURCES:.cpp=.o)
EXECUTABLE:=MaxBin

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $(EXECUTABLE) $(OBJECTS) $(INC) $(LDFLAGS)

$(OBJECTS): %.o: %.cpp
	$(CC) $(CFLAGS) $< -o $@ $(INC)

clean:
	rm $(EXECUTABLE) $(OBJECTS)

