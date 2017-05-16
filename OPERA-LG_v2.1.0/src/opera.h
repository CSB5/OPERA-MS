#pragma once

#include <list>
#include <set>
#include <map>
#include <time.h>
#include <algorithm>
#include <sys/time.h>
#include <unistd.h>
#include <stdio.h>

#include "configureReader.h"
#include "ContigConverter.h"
#include "MapConverter.h"
#include "Configure.h"
#include "Graph.h"
#include "StartPoint.h"
#include "CommonFunction.h"
#include "PartialScaffold.h"
#include "QuadProg.h"
#include "Tree.h"
#include "FinalScaffold.h"
#include "PetLibrary.h"
#include "GapCorrecter.h"

using namespace std;
using namespace QuadProgPP;

struct great_length {
  bool operator() (const FinalScaffold *s1, const FinalScaffold *s2) const
  {return (s1->GetLength() > s2->GetLength());}
};

// the structure for saving mapping information of contigs
struct reference_mapping {
	string m_contigName;                     // the name of this contig
	string m_referenceGenomeName;            // the reference genome name of this contig
	int m_referenceStartPos;                 // the 
	int m_referenceEndPos;                   // the start position on reference genome
	int m_referenceOri;                      // the orientation of this contig on reference genome
	int m_referenceIndex;                    // the index of this contig on reference genome

	int m_subgraphReferenceIndex;            // the index of this contig within one subgraph on reference genome
};

class opera
{
public:
	opera(void);
	~opera(void);

	// Attributes
public:
	Graph *m_graph;
	list<PetLibrary*> *m_libraries;		// all the libraries in the data
	list<PetLibrary*> *m_superLibrary;      // the library saving super clusters
	list<string> *m_conflictingEdges;        // all the conflicting edges
	FILE *m_discardEdgeFile;				// record all discarded edges due to increased threshold
	bool m_isSpecialGraph;
	int m_totalPartialScaffold;
        int m_numOfTopUnassignedNodes;
	Contig *m_repeatTraceBack;               // record the repeats need to trace back to
	
	// properties about increasing threshold
	long m_totalNumberOfSubgraph;                    // the total number of subgraphs
	long m_totalNumberOfIncreasedSubgraph;           // the total number of subgraphs which increase threshold
	long m_totalNumberOfNodes;                       // total number of contigs in the whole graph
	set<string> m_totalIncreasedNodes;                  // the total number of nodes which increase threshold

private:
	list<Contig*> *m_activeRegion;			// the global active region
	list<PET*> *m_happyDanglingEdges;		// the global happy dangling edges
	list<PET*> *m_unhappyDanglingEdges;		// the global unhappy dangling edges
	list< pair<Contig*, int> > *m_unassignedNodes;		// the global unassigned contigs
	list<string> *m_unhappyEdgeString;			// the string of all unhappy edges
	// initialize the tree saving all visited scaffolds
	Tree *m_visitedTree;
	FILE *logFile;						// the log file to record the output of program
	string m_stats;						// the statistics of assembly

	int m_graphIDIncreaseThreshold;
	map<string, reference_mapping*> *m_realMapping;         // the real mapping of all unique contigs
	int m_subgraphSize;
	
	

	// Methods
public:
	// general pipeline of scaffolding
	int StartScaffold();
	// scaffolding
	int Scaffolding( PartialScaffold *&currentScaffold, int maxOfUnhappyEdges, int numberOfContigs );
	// output scaffolds
	int OutputScaffold( string fileName );
	// output unhappy edges
	int OutputUnhappyEdges( string fileName );
	// sort all scaffold according to length
	int SortScaffold();
	// check if the names have confliction
	bool CheckNameConfliction();
	// check contig file format: velvet, soap or normal fasta file
	int CheckContigFormat();
	void OpenLogFile();
	// Get statistics of assembly
	string GetStats();
	// check if all files exist
	bool CheckFileExist();

	// multiple libraries related methods
	// move the edges in first library into graph
	void MoveEdgesToGraph( PetLibrary *lib );
	// initialize for another run
	void Initialize();

	// analyze assembly file to get kmer size
	bool GetKmer( string contigFileName, int type );
	
	// print the parameter into a file
	int PrintParameter( string fileName );

	// output conflicting edges
	int OutputConflictingEdges( string fileName );

	// read real position file
	int ReadRealPositionFile( string fileName );

private:
	// select the starting contig
	void FindStartingContig( list<Contig*> *subgraph, int numOfContig, int numOfBorderContig,
		Contig *&startContig, int &startContigOri );
	// create a scaffold containing a start contig
	PartialScaffold* CreatScaffold( Contig *c, int ori );
	// create the initial scaffold containing all nodes as unassigned nodes
	PartialScaffold* CreateInitialScaffold( list<Contig*> *subgraph, int numOfContig );

	// add happy dangling edges to first scaffolds
	void AddHappyDEToFirstScaffold( list<PET*> *edges, PartialScaffold *p );
	// add unhappy dangling edges to first scaffolds
	void AddUnhappyDEToFirstScaffold( list<PET*> *edges, PartialScaffold *p, Contig *c );
	// generate active region string
	void GenARString( PartialScaffold *p );
	// generate unhappy dangling edges string
	void GenUDEString( PartialScaffold *p );
	// generate unassigned contig set
	void GenUnassigndNodes( PartialScaffold *p, bool isStart, bool ifRecalculate );
	// generate unassigned contigs for initial scaffold
	void GenUnassigndNodesForInitialScaffold( StartPoint **allStart, int numOfStarts );
	// traverse edge and find unassigned node
	void TraverseEdges( Contig *c, list<PET*> *edges, list< pair<Contig*, int> > *&possibleContig,
			    set<Contig*> *&unassignedNodes, set<PET*> *&visitedEdges, int direction, bool isStart, bool ifRecalculate );
	// add one contig to partial scaffold
	bool AddContigToScaffold( PartialScaffold *&p, int maxUnhappyEdge, Contig *newContig, int newContigOri );
	// check edges of newly added contig
	void CheckEdgesOfNewContig( PartialScaffold *p, Contig *c, list<PET*> *edges, int ori );
	// check incident edge happiness
	bool CheckDanglingEdgeHappiness( PET *edge, Contig *newContig, int ori );
	// check if it is a left edge of current contig with certain orientation
	inline bool IsLeftEdge( PET *edge, Contig *newContig, int ori );
	// check connectivity and update scaffold active region & dangling edges
	bool CheckScaffold( PartialScaffold *p );
	// check if contig c could connect to active region
	bool CheckActiveRegionConnectivity( Contig *c );
	// trace back to parent partial scaffold
	void TraceBack( PartialScaffold *&p );
	// find start contig of new scaffold
	StartPoint* FindStartOfNewScaffold();
	// add connected contigs to corresponding list
	inline void AddToList( list<StartPoint*> *borderContigs, list<StartPoint*> *allContigs, Contig *c, int ori );
	// check edges to find start point
	inline void CheckEdgesForStartPoint( list<StartPoint*> *borderContigs, list<StartPoint*> *allContigs, 
		Contig *c, list<PET*> *edges, set<Contig*> *visitedContigs, list<Contig*> *possibleContigs );
	// Break scaffold
	void BreakScaffold( PartialScaffold *&p, int maxUnhappyEdges );
	// generate results of scaffold
	void GenerateResults( PartialScaffold *p, ScaffoldResult **&results, int num );
	// print scaffold result
	void PrintScaffoldResult( ScaffoldResult **results, int number );
	// clear scaffold variables
	inline void Clear();
	// delete all scaffold
	inline void DeleteAllScaffolds( PartialScaffold *p, PartialScaffold *head, bool clearRepeat );
	// calculate gap sizes
	void CalculateGap( vector<Contig*> *contigs, int numberOfGaps );

		// read result contig files
	int ReadContigFile( string fileName, multiset<FinalScaffold*, great_length> *scaffolds, double &totalLength );
	// read scaffolds file
	int ReadScaffoldFile( string fileName, multiset<FinalScaffold*, great_length> *scaffolds, double &totalLength, int &nonSingleton );
	// read contig fasta file
	int ReadContigFasta( string fileName, map<string, string> *contigs );

	        // sorting unassigned nodes
	// insert a contig into unassigned nodes set according to distance
	void InsertContigIntoUnassignedNodeSet( pair<Contig*, int> pa );

	        // methods about analyzing scaffolds (without checking connectivity)
	// separate the scaffolds into independent ones
	void SeparateScaffolds( list<Contig*> *subgraph, int &numberOfScaf, vector<vector<Contig*>*> *scafs, PartialScaffold *currentScaffold,
				list<Contig*> *allRepeats);
	// generate results of scaffold
	void GenerateResults( vector<vector<Contig*>*> *scafs, ScaffoldResult **&results );

	// recover unassigned nodes
	void RecoverUnassignedNodes( list< pair<Contig*, int> > *nodes );
	
	// generate the graph file for certain subgraph
	void GenerateDotFile( list<Contig*> *subgraph, int id );
	// print edges to file
	int PrintEdges( FILE *file, Contig *c, list<PET*> *edges, int type, list<Contig*> *subgraph );

	// label the order of contigs in the subgraph according to the real mapping
	void LabelContigsUsingMapping( list<Contig*> *subgraph );

	// print AR string
	string PrintARStringOfName();

	// create a new scaffold for all disconnected nodes
	void CreateNewScaffold();

	// check edges to find all unused contigs
	void CheckEdgesForUnusedContigs( list<Contig*> *subgraph, Contig *c, list<PET*> *edges, 
					 set<Contig*> *visitedContigs, list<Contig*> *possibleContigs );

	// check if a repeat can be removed from the active region without introducing any discordant edges
	bool BringDiscordantEdgeByRemoving( Contig *c, int ori, list<Contig*> erasedContigs, PartialScaffold *p, int sizeOfActiveRegion );
	
	// Given a repeat, collect the edges related with a set of contigs
	void CollectRelatedEdges( list<PET*> *collectedEdges, list<PET*> *possibleEdges, Contig *repeat, list<SUBCONTIG*> *region );
	void CollectRelatedLeftEdges( list<PET*> *collectedEdges, list<PET*> *possibleEdges, Contig *repeat, list<SUBCONTIG*> *region );
	void CollectRelatedRightEdges( list<PET*> *collectedEdges, list<PET*> *possibleEdges, Contig *repeat, list<SUBCONTIG*> *region );

	// check if all the edges can be satisfied
	bool CheckEdgeSatisfaction( list<PET*> *collectedEdges, list<SUBCONTIG*> *region, Contig *repeat );

	// check repeated related edge of a new contig
	bool CheckRepeatEdgeOfNewContig( Contig *repeatContig, int repeatOri, Contig *addedContig, PET *edge );

	// add the repetitive edges to repeats
	void AddRepetitiveEdges( list<PET*> leftEdges, list<PET*> rightEdges, Contig *repeat );
	void AddRepetitiveEdges( PET *edge, Contig *repeat );

	// remove the labeled repeats
	void RemoveRepeats( vector<Contig*> *scafs );

	// assign repetitive edges
	void AssignRepetitiveEdges( list<PET*> *edges, Contig *repeat );

	// create the unique edge of an edge connecting with a repeat
	void CreateUniqueEdge( Contig *startContig, int startOri, Contig *endContig, int endOri, double disDelta, double disDeltaWithGap, PET *oriEdge );

	// add scaffold string
	// addtach the new string before old string
	string AddScaffoldString( string newString, string oldString );

      
};

bool SortStartPoints( const StartPoint *s1, const StartPoint *s2 );
bool SortStartPointsUsingDis( const StartPoint *s1, const StartPoint *s2 );

// comparison, not case sensitive.
bool compare_real_position( Contig *first, Contig *second);
