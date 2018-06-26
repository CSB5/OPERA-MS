#pragma once

#include <vector>
#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <list>
#include <set>
#include <algorithm>
#include <ext/hash_map>
#include <map>
#include <time.h>
#include <sys/time.h>
#include "Contig.h"
#include "Configure.h"
#include "Graph.h"
#include "SinglePet.h"
#include "PET.h"
#include "CommonFunction.h"
#include "PetLibrary.h"
#include "GapCorrecter.h"

using namespace std;
using namespace __gnu_cxx; 

struct less_distance {
  bool operator() (const SinglePet* p1, const SinglePet* p2) const
  {return (p1->GetDistance() < p2->GetDistance());}
};

struct less_std {
  bool operator() (const PET* p1, const PET* p2) const
  {return (p1->GetStd() < p2->GetStd());}
};


class MapConverter
{
public:
	MapConverter(void);
	MapConverter( Graph *graph );
	~MapConverter(void);

	// Methods
public:
	// analyze mapping file or opera edge file
	int Analyze( string fileName );
	// start bundle
	int Bundle();

	// for multiple libraries
	// analyze mapping file or opera edge file
	int AnalyzeMultiLib( string fileName, list<PetLibrary*> *libs  );
	// start bundle
	int BundleMultiLib( list<PetLibrary*> *libs, FILE *discardedEdgeFile, list<PetLibrary*> *superCluster, list<string> *conflictingEdges );

private:
			// mapping file converter related functions
	// analyze bowtie file
	int AnalyzeBowtie( string fileName );
	// analyze opera edge file
	int AnalyzeOpera( string fileName );
	// calculate library mean and std
	int CalculateLibDis( string fileName );
	// check if two lines represent a pair of reads
	bool IsPair( string firstRead, string secondRead );
	// calculate the distance of two reads in the same contig
	int CalculateReadDisOnSameContig( vector<string> *firstColumn, vector<string> *secondColumn );
	// convert bowtie format into opera format
	int ConvertBowtieFile( string fileName );
	// convert bowtie format into opera format
	int ConvertBowtieFile();
	// convert bowtie format
	string ConvertBowtieFormat( vector<string> *firstColumn, vector<string> *secondColumn, int id );
	
			// bundle related functions
	// bundle a certain cluster
	string BundleCluster( multiset<SinglePet*, less_distance> *group );

	// for multiple libraries
	// analyze bowtie file
	int AnalyzeBowtieMultiLib( list<PetLibrary*> *libs );
	// calculate library mean and std
	int CalculateLibDisMultiLib( list<PetLibrary*> *libs );
	// convert bowtie format into opera format
	int ConvertBowtieFileMultiLib( PetLibrary *lib );
	// convert bowtie format
	string ConvertBowtieFormatMultiLib( mappingInfo *firstRead, mappingInfo *secondRead, PetLibrary *currentLib, int id ); 
	// bundle a certain cluster
	string BundleClusterMultiLib( multiset<SinglePet*, lessDistance> *group, PetLibrary *lib, ofstream *clusterInfoFile, FILE *discardedEdgeFile,
				      PetLibrary *currentLibrary );
	// analyze opera edge file
	int AnalyzeOperaMultiLib( list<PetLibrary*> *libs );

	// for contig graph
	// analyze contig graph
	int ConvertSoapContigGraph( string fileName, list<PetLibrary*> *libs );
	// get two kmers from the line
	void GetTwoKmer( string line, string &firstKmer, string &lastKmer );
	// insert the kmer into map
	void InsertKmerIntoMap( string kmer, map<string, list<int>* >* kmerList, int contigID, int ori );
	// find all overlap information of SOAPdenovo
	void FindOverlapOfSOAP();

	// debugging purpose
	string CalculateAnchorRegion( list<SinglePet*> *pets );

	// calculate mean and std using Interquartile Range
	int IQR( vector<double> *dis, double &mean, double &std, string libName, double &percentageOfOutliers, double &numberOfUsedPairs, 
		 PetLibrary *currentLibrary );

	// bundle multiple libraries to create super clusters
	void BundleSuperCluster( list<PetLibrary*> *libs, list<PetLibrary*> *superCluster, list<string> *conflictingEdges, int libType );
	// bundle a certain super cluster
	string BundleOneSuperCluster( multiset<PET*, less_std> *group, PetLibrary *lib, list<string> *conflictingEdges ); //, ofstream *clusterInfoFile, FILE *discardedEdgeFile );

	// correct the distance of clusters
	//void CorrectDistance( multiset<SinglePet*, lessDistance> *group, PetLibrary *lib );
	
	// create tables for gap correction
	void CreateGapTables( PetLibrary *lib );
	
	// refine the mean of a cluster
	void RefineMean( PET *e, Contig *c1, Contig *c2 );

	// Attributes
private:
	// graph
	Graph *m_graph;
	// number of single PET
	int m_numOfPet;
	// the array of single PET
	list<SinglePet*> *m_singlePets;
	// the map of single PET, first element is the contig pairs, second element is the set of all edges
	map<pair<int, int>, multiset<SinglePet*, less_distance>*> *m_singlePetsMap;
	// list of all paired reads mapped on different scaffold
	list<string> *m_pairMap;
	double m_totalTime;
	double m_pairTime;

	// for multiple libraries
	string m_libString;

	// for contig graph
	hash_map<const char*, list<int>*, hash<const char*>, eqName> *m_headHashMap;
	hash_map<const char*, list<int>*, hash<const char*>, eqName> *m_tailHashMap;
	map<string, list<int>* > *m_headList;			// the map to save all the first kmer of each contig, contigID orientation
	map<string, list<int>* > *m_tailList;			// the map to save all the last kmer of each contig

	// the map of all clusters, first element is the contig pairs, second element is the set of all clusters
	map<pair<string, string>, multiset<PET*, less_std>*> *m_ClustersMap;
	
	// the set of table to correct the edge distance
	vector<GapCorrecter*> m_gapTables;
};

