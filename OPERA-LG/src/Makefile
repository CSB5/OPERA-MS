PROGRAM   := OPERA-LG

OBJS := $(wildcard *.o)
SOURCE := $(wildcard *.cpp)
#CXX := g++

install : $(PROGRAM)

$(PROGRAM) : $(OBJS)     
	$(CXX) -c -O3 -std=c++98 -Wall -Wno-deprecated $(SOURCE) $<
#	$(CXX) -c -g -fno-inline -w $(CXX) $<
#	$(CXX) -c -g -fno-inline -w $(CXX) $< -pg
#	$(CXX) -o ../bin/$(PROGRAM) *.o -pg
	$(CXX) -o ../bin/$(PROGRAM) *.o
#	rm *.o

clean:
	rm *.o
