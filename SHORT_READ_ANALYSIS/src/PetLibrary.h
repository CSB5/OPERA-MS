#pragma once

#include <string>
#include <list>
#include <map>
#include <set>
#include <vector>
class PET;
#include "PET.h"
#include "SinglePet.h"

using namespace std;

struct lessDistance {
  bool operator() (const SinglePet* p1, const SinglePet* p2) const
  {return (p1->GetDistance() < p2->GetDistance());}
};

struct mappingInfo{
	Contig *contig;                // the mapped contig
	int position;             // the mapping position
	string orientation;          // the mapping orientation
	int readLength;           // the read length
};

class PetLibrary
{
public:
	PetLibrary( string name );
	virtual ~PetLibrary(void);

	// Attributes
private:
	int m_mean;		// the mean of this library
	int m_std;		// the standard deviation of this library
	int m_ori;		// orientation of paired reads
	string m_fileName;			// the name of mapping file
	string m_fileNameWithoutPath;		// the pure file name without path and extension
	// the map of single PET, first element is the contig pairs, second element is the set of all edges
	map<pair<int, int>, multiset<SinglePet*, lessDistance>*> *m_singlePetsMap;
	list<PET*> *m_edges;			// the edges of this library
	//int m_threshold;

	// variables related with gap correction
	vector<double> *m_possibilityOfDistance;    // the possibility of each distance for a paired-end reads
	int m_minDistance;
	int m_maxDistance;
	
	int m_readLength;                       // the average read length (both reads)

	list<mappingInfo*> *m_mapInfo;          // the mapping information of paired-end reads mapped on different contigs

	// Methods
public:	
	void SetMean( int mean );
	void SetStd( int std );
	void SetOri( int ori );
	int GetOri();	
	int GetMean();
	int GetStd();
	//void SetThreshold( int thre );
	//int GetThreshold();
	void SetNameWithoutPath( string name );
	string GetNameWithoutPath();
	// insert a contig pair into map
	void InsertPair( int c1, int c2, SinglePet *pet );
	// get the pair map
	map<pair<int, int>, multiset<SinglePet*, lessDistance>*>* GetSinglePetsMap();
	// add a cluster
	void AddCluster( PET *pet );
	// get the pet clusers
	list<PET*>* GetEdges();
	// get file name
	string GetFileName();

	// initialize the possibility array of all possible distance
	void InitPossibility( int min, int max );
	// add one occurrance of a certain distance
	void AddOccuOfDistance( int dis );
	// Calculate possibility
	void CalculatePossibility();
	// get the vector of possibility
	vector<double>* GetPossibilities();
	int GetMinDis();
	int GetMaxDis();

	// set read length
	void SetReadLength( int l );
	int GetReadLength();

	///// methods about mapping information of paired-end reads mapped on different contigs
	// add one mapping information 
	void AddMapInfo( Contig *c, int pos, string ori, int readLength );   
	// get the mapping information
	list<mappingInfo*>* GetMappingInfo();

	// remove an edge
	void RemoveEdge( list<PET*>::iterator iter );
};

