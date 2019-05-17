PROGRAM   := OPERA-LG

OBJS := $(wildcard *.o)
SOURCE := $(wildcard *.cpp)
CXX := g++

install : $(PROGRAM)

$(PROGRAM) : $(OBJS)     
#<<<<<<< HEAD
	$(CXX) -c -O3 -std=c++11 -Wall -Wno-deprecated $(SOURCE) $<
#	$(CXX) -c -g -std=c++11 -fno-inline -w $(SOURCE) $<
#	$(CXX) -c -g -fno-inline -w $(CXX) $< -pg
#	$(CXX) -o ../bin/$(PROGRAM) *.o -pg
	$(CXX) -o ../bin/$(PROGRAM) *.o
#=======
#	$(CXX) -c -O3 -std=c++98 -Wall -Wno-deprecated $(SOURCE) $<
#	$(CXX) -c -g -fno-inline -w $(CXX) $<
#	$(CXX) -c -g -fno-inline -w $(CXX) $< -pg
#	$(CXX) -o ../bin/$(PROGRAM) *.o -pg
#	$(CXX) -o ../bin/$(PROGRAM) *.o
#>>>>>>> a0221e115ddf9f91d3667cf2fd6a7fecdf15a8c7
#	rm *.o

clean:
	rm *.o
