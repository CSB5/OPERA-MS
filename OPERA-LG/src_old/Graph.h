#pragma once

#include "Contig.h"
#include <list>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include "CommonFunction.h"
#include "ScaffoldResult.h"
#include <ext/hash_map>
//#include <hash_map>

using namespace std;
using namespace __gnu_cxx; 

struct eqName
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1,  s2) == 0;
  }
};

class Graph
{
public:
	Graph(void);
	~Graph(void);

	// Methods
public:
	// Add repeat contigs
	void AddRepeat( Contig *c );
	// Add small contigs
	void AddSmallContigs( Contig *c );
	// set total number of contigs
	void SetContigNum( int n );
	// Add non-repeat contigs
	void AddContig( Contig *c );
	// output the contigs into a new file
	int OutputContigs( string fileName );
	// get contig index according to contig name
	int GetContigIndex( string contigName );
	// get contig in specific index
	Contig * GetContig( int pos );
	// add edge to graph
	void AddEdge( PET *p );
	// generate border contig and singletons
	void GenerateBorderAndScaffold();
	// update the graph using results of subgraph
	void UpdateGraph( list<Contig*> *subgraph, ScaffoldResult **scaffold, FILE* discardEdgeFile );
	// find subgraph
	list<Contig*> * FindSubgraph( int &numOfContig, int &numOfBorderContig, int &minClusterSize );
	int GetNumOfContigs();
	// check if current graph still has edges
	bool HasEdges();
	// output final scaffolds to files
	int OutputScaffolds( string fileName );

	// generate the id map of all contigs
	void GenerateIDMap();
	// get the contig using velvet id
	Contig* GetContigUsingID( int id );
	// get the contig in position pos
	Contig* GetContigUsingPos( int pos );

	// multiple libraries related
	void AddEdgeMultiLib( PET *p );		// add edges from all libraries to graph
	// clear the edges number
	void ClearEdge();	
	// initialize contig list for a new library
	void InitializeContig(  int libNum );

	// dynamically change cluster threshold related functions
	// remove edge p from graph
	void RemoveOneEdge();
	// set a contig as singleton
	void SetAsSingleton( Contig *c );

	// contig graph related
	// get all contigs
	list<Contig*>* GetContigs();

	// save contig name
	void SaveContigName( string name );
	// check if a contig name exist
	bool IfExistInContig( string name );

	// count all the subgraph
	void CountSubgraph();
	void ClearSubgraph( list<Contig*> *subgraph );

	// insert repeat contig name
	void InsertRepeatContigName( Contig *c );
	
	// check if a certain contig is a repeat
	bool IsRepeat( string name );

	// check the number of edges of certain contig
	void CheckNumberOfEdges( string name );

	void CheckEdgesForAdjacentRepeat( PET *edge );	

	

private:
	void DeleteEdges( Contig *c, list<PET*> *edges, bool isFirst, ScaffoldResult *scaffold, bool ifSame );			// release the memory of edges
	// check if current scaffold has the same orientation in scaffold and edge
	bool IfSameOriInScaffoldAndEdge( Contig *c, int ori, ScaffoldResult *s );
	// traverse and find subgraph starting from contig c 
	list<Contig*> * TraverseSubgraph( Contig *c, int direction, int &numOfContig, int &numOfBorderContig, int &minClusterSize );
	// add list of edges into subgraph
	void AddEdgesToSubgraph( Contig *c, list<PET*> *edges, list<Contig*> *contigList, int &minClusterSize );
	inline void DeleteContigs( list<Contig*> *contigs );
	// check if contig c has the same orientation as in scaffold
	bool SameOri( Contig *c, string scaffold );

	// multiple libraries related
	void DeleteEdgesMultiLib( Contig *c, list<PET*> *edges, bool isFirst, ScaffoldResult *scaffold, bool ifSame,
				  FILE* discardEdgeFile);			// release the memory of edges

	// traverse and find subgraph starting from contig c 
	list<Contig*> * TraverseSubgraphToCount( Contig *c, int direction, int &numOfContig, int &numOfBorderContig, int &minClusterSize, int &numberOfRepeats );
	
	// check if two contigs contain the same repeats near each other, if yes, modify the distance
	void CheckEdgesForAdjacentRepeat( list<PET*> *edges );	

	


	// Attributes
private:
	int m_numberOfContigs;				// the number of contigs
	int m_numberOfEdges;				// the number of edges
	Contig **m_contigsArray;			// non repeat contigs array
	list<Contig*> *m_contigsList;		// non repeat contig list for contigs having edges;
	list<Contig*> *m_scaffoldsList;		// the list of all scaffolds( singletons )
	list<Contig*> *m_borderContigsList;	// the border contig list(with edges)
	//map<string, int> *m_contigNameMap;	// map saving contig names 
	hash_map<const char*, int, hash<const char*>, eqName> *m_contigNameHashMap;
	//hash_map<const char*, int, eqName> *m_contigNameHashMap;
	bool ifInitializeContig;

	map<int, Contig*> *m_contigID;		// record the relationship of contig ID of velvet and conigs
	set<string> *m_contigNames;             // all the contig names in contig file

	void AddEdgeToSubgraph( list<Contig*> *subgraph );          // check if there are missing edges in this subgraph

	hash_map<const char*, bool, hash<const char*>, eqName> *m_repeatContigSet;    // a set saving all the repetitive contigs' names
};

// comparison, not case sensitive.
bool compare_scaffold( Contig *first, Contig *second);
