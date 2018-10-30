PROGRAM   := OPERA-LG

OBJS := $(wildcard *.o)
CPP := $(wildcard *.cpp)

install : $(PROGRAM)

$(PROGRAM) : $(OBJS)     
	g++ -c -O3 -std=c++11 -Wall -Wno-deprecated $(CPP) $<
#	g++ -c -g -std=c++11 -fno-inline -w $(CPP) $<
#	g++ -c -g -fno-inline -w $(CPP) $< -pg
#	g++ -o ../bin/$(PROGRAM) *.o -pg
	g++ -o ../bin/$(PROGRAM) *.o
#	rm *.o

clean:
	rm *.o
