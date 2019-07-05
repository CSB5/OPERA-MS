CC = g++
short_read_analysis: CFLAGS = -std=c++11 -O3 -Wall -Wno-deprecated
debug: CFLAGS = -std=c++0x -g -O0 -Wall -Wconversion -fno-inline
#debug: CFLAGS = -std=c++11 -g -O0 -Wall -Wconversion -fno-inline

all: short_read_analysis

debug: Array.o ARTreeNode.o cluster_graph.o cluster.o CommonFunction.o Configure.o configureReader.o ContigConverter.o contig.o Contig.o contig_reader.o edge.o edge_reader.o FinalScaffold.o GapCorrecter.o Graph.o LibInfo.o MapConverter.o mapping_reader.o opera.o PartialScaffold.o PetLibrary.o PET.o probability_distribution.o QuadProg.o ScaffoldResult.o sigma.o SinglePet.o StartPoint.o TreeNode.o Tree.o
	$(CC) $(CFLAGS) Array.o ARTreeNode.o cluster_graph.o cluster.o CommonFunction.o Configure.o configureReader.o ContigConverter.o contig.o Contig.o contig_reader.o edge.o edge_reader.o FinalScaffold.o GapCorrecter.o Graph.o LibInfo.o MapConverter.o mapping_reader.o opera.o PartialScaffold.o PetLibrary.o PET.o probability_distribution.o QuadProg.o ScaffoldResult.o sigma.o SinglePet.o StartPoint.o TreeNode.o Tree.o -o debug

short_read_analysis: Array.o ARTreeNode.o cluster_graph.o cluster.o CommonFunction.o Configure.o configureReader.o ContigConverter.o contig.o Contig.o contig_reader.o edge.o edge_reader.o FinalScaffold.o GapCorrecter.o Graph.o LibInfo.o MapConverter.o mapping_reader.o opera.o PartialScaffold.o PetLibrary.o PET.o probability_distribution.o QuadProg.o ScaffoldResult.o sigma.o SinglePet.o StartPoint.o TreeNode.o Tree.o
	$(CC) $(CFLAGS) Array.o ARTreeNode.o cluster_graph.o cluster.o CommonFunction.o Configure.o configureReader.o ContigConverter.o contig.o Contig.o contig_reader.o edge.o edge_reader.o FinalScaffold.o GapCorrecter.o Graph.o LibInfo.o MapConverter.o mapping_reader.o opera.o PartialScaffold.o PetLibrary.o PET.o probability_distribution.o QuadProg.o ScaffoldResult.o sigma.o SinglePet.o StartPoint.o TreeNode.o Tree.o  -o ../../bin/short_read_analysis

Array.o: Array.cpp Array.h
	$(CC) $(CFLAGS) -c Array.cpp

ARTreeNode.o: TreeNode.h TreeNode.cpp ARTreeNode.cpp ARTreeNode.h
	$(CC) $(CFLAGS) -c ARTreeNode.cpp

cluster_graph.o: cluster.h cluster.cpp edge.h edge.cpp contig.h contig.cpp probability_distribution.h probability_distribution.cpp cluster_graph.h cluster_graph.cpp
	$(CC) $(CFLAGS) -c cluster_graph.cpp

cluster.o: contig.h contig.cpp cluster.h cluster.cpp
	$(CC) $(CFLAGS) -c cluster.cpp

CommonFunction.o: Configure.h Configure.cpp CommonFunction.h CommonFunction.cpp
	$(CC) $(CFLAGS) -c CommonFunction.cpp

Configure.o: LibInfo.h LibInfo.cpp Configure.h Configure.cpp
	$(CC) $(CFLAGS) -c Configure.cpp

configureReader.o: Configure.h Configure.cpp LibInfo.h LibInfo.cpp CommonFunction.h CommonFunction.cpp configureReader.h configureReader.cpp
	$(CC) $(CFLAGS) -c configureReader.cpp

ContigConverter.o: Configure.h Configure.cpp Contig.h Contig.cpp Graph.h Graph.cpp CommonFunction.h CommonFunction.cpp ContigConverter.h ContigConverter.cpp
	$(CC) $(CFLAGS) -c ContigConverter.cpp

contig.o: contig.h contig.cpp
	$(CC) $(CFLAGS) -c contig.cpp

Contig.o: PET.h PET.cpp CommonFunction.h CommonFunction.cpp Contig.h Contig.cpp
	$(CC) $(CFLAGS) -c Contig.cpp

contig_reader.o: contig.h contig.cpp contig_reader.h contig_reader.cpp
	$(CC) $(CFLAGS) -c contig_reader.cpp

edge.o: contig.h contig.cpp edge.h edge.cpp
	$(CC) $(CFLAGS) -c edge.cpp

edge_reader.o: contig.h contig.cpp edge.h edge.cpp edge_reader.h edge_reader.cpp
	$(CC) $(CFLAGS) -c edge_reader.cpp

FinalScaffold.o: FinalScaffold.h FinalScaffold.cpp
	$(CC) $(CFLAGS) -c FinalScaffold.cpp

GapCorrecter.o: CommonFunction.h CommonFunction.cpp GapCorrecter.h GapCorrecter.cpp
	$(CC) $(CFLAGS) -c GapCorrecter.cpp

Graph.o: CommonFunction.h CommonFunction.cpp ScaffoldResult.h ScaffoldResult.cpp Graph.h Graph.cpp
	$(CC) $(CFLAGS) -c Graph.cpp

LibInfo.o: Configure.h Configure.cpp LibInfo.h LibInfo.cpp
	$(CC) $(CFLAGS) -c LibInfo.cpp

MapConverter.o: Contig.h Contig.cpp Configure.h Configure.cpp Graph.h Graph.cpp SinglePet.h SinglePet.cpp PET.h PET.cpp CommonFunction.h CommonFunction.cpp PetLibrary.h PetLibrary.cpp GapCorrecter.h GapCorrecter.cpp MapConverter.h MapConverter.cpp
	$(CC) $(CFLAGS) -c MapConverter.cpp

mapping_reader.o: contig.h contig.cpp mapping_reader.h mapping_reader.cpp
	$(CC) $(CFLAGS) -c mapping_reader.cpp

opera.o: configureReader.h configureReader.cpp ContigConverter.h ContigConverter.cpp MapConverter.h MapConverter.cpp Configure.h Configure.cpp Graph.h Graph.cpp StartPoint.h StartPoint.cpp CommonFunction.h CommonFunction.cpp PartialScaffold.h PartialScaffold.cpp QuadProg.h QuadProg.cpp Tree.h Tree.cpp FinalScaffold.h FinalScaffold.cpp PetLibrary.h PetLibrary.cpp GapCorrecter.h GapCorrecter.cpp opera.h opera.cpp
	$(CC) $(CFLAGS) -c opera.cpp

PartialScaffold.o: Contig.h Contig.cpp PET.h PET.cpp PartialScaffold.h PartialScaffold.cpp
	$(CC) $(CFLAGS) -c PartialScaffold.cpp

PetLibrary.o: PET.h PET.cpp SinglePet.h SinglePet.cpp PetLibrary.h PetLibrary.cpp
	$(CC) $(CFLAGS) -c PetLibrary.cpp

PET.o: Contig.h Contig.cpp Configure.h Configure.cpp CommonFunction.h CommonFunction.cpp PetLibrary.h PetLibrary.cpp PET.h PET.cpp
	$(CC) $(CFLAGS) -c PET.cpp

probability_distribution.o: probability_distribution.h probability_distribution.cpp
	$(CC) $(CFLAGS) -c probability_distribution.cpp

QuadProg.o: Array.h Array.cpp QuadProg.h QuadProg.cpp
	$(CC) $(CFLAGS) -c QuadProg.cpp

ScaffoldResult.o: ScaffoldResult.h ScaffoldResult.cpp
	$(CC) $(CFLAGS) -c ScaffoldResult.cpp

sigma.o: sigma.h sigma.cpp contig_reader.h mapping_reader.h edge_reader.h contig.h cluster.h cluster_graph.h probability_distribution.h
	$(CC) $(CFLAGS) -c sigma.cpp

SinglePet.o: Configure.h Configure.cpp SinglePet.h SinglePet.cpp
	$(CC) $(CFLAGS) -c SinglePet.cpp

StartPoint.o: Contig.h Contig.cpp Configure.h Configure.cpp StartPoint.h StartPoint.cpp
	$(CC) $(CFLAGS) -c StartPoint.cpp

TreeNode.o: TreeNode.h TreeNode.cpp
	$(CC) $(CFLAGS) -c TreeNode.cpp

Tree.o: TreeNode.h TreeNode.cpp ARTreeNode.h ARTreeNode.cpp CommonFunction.h CommonFunction.cpp Tree.h Tree.cpp
	$(CC) $(CFLAGS) -c Tree.cpp

clean:
	-rm *.o short_read_analysis
