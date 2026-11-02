#pragma once

// the paired-end reads cluster, used as edge in graph

#include <string>
class Contig;
#include "Contig.h"
#include "Configure.h"
#include "CommonFunction.h"
class PetLibrary;
#include "PetLibrary.h"

using namespace std;

class PET
{
public:
	PET(void);
	PET( Contig *c1, int o1, Contig *c2, int o2, int dis, int std, int size );
	PET( Contig *c1, int o1, Contig *c2, int o2, int dis, int std, int size, PET *edge );
	~PET(void);

	// Method
public:
	Contig * GetStartContig();
	Contig * GetEndContig();
	int GetSize();
	int GetStartContigOri();
	int GetEndContigOri();
	friend ostream& operator<<( ostream& out, const PET& p );
	string ToString();
	bool IfUnhappy();
	bool IfInSubgraph();
	bool IfVisited();
	void VisitEdge();
	int GetPositionOfContig( Contig *c );		// get the positon of contig c: start or end
	int GetOrientationOfContig( Contig *c );
	void ReplaceContig( int pos, Contig *c );		// replace the contig in pos using c
	void SetOri( int pos, int ori );
	void SetDis( int dis );
	int GetDis();
	void SetDisWithGap( int dis );
	int GetDisWithGap();
	void SetInSubgraph();
	void SetNotInSubgraph();
	bool isInSubgraph();
	Contig* GetOtherContig( Contig *c );			// get the other contig rather than c
	void SetUnhappy( Contig *c );			// label the edge as unhappy
	void SetHappy();			// label the edge as happy
	int GetID();
	void SetID( int id ); 
	void SetDE( bool isDE );		// set if it is a dangling edge
	bool IsDE();			// check if it is a dangling edge
	int GetStd() const;
	void SetHappyDEIter( list<PET*>::iterator iter );
	list<PET*>::iterator GetHappyDEIter();
	void SetUnhappyDEIter( list<PET*>::iterator iter );
	list<PET*>::iterator GetUnhappyDEIter();
	// find which contig is in active region
	Contig* FindARContig();
	// set this edge as unhappy
	void SetUnhappy();
	// get unused contig
	Contig* GetUnusedContig();
	void SetOriString( string s );
	string GetOriString();
	// generate original cluster string
	void GenClusterString();

	// multiple libraries related
	// set the edge iterater, edges in all libraries of contig
	void SetEdgeIterMultiLib( list<PET*>::iterator iter, int type, int contigPos );
	list<PET*>::iterator GetEdgeIterMultiLib( int contigPos );
	int GetEdgeTypeMultiLib( int contigPos );
	// add the corresponding library and iterator of edges
	void AddIterInLib( PetLibrary *lib, list<PET*>::iterator iter );
	// delete this edge from library list
	void DeleteFromLib();

	// initialize all properties
	void Initialize();

	// calculate anchor region
	void AddAnchorString( string s );

	// add one cluster to super cluster
	void AddOneClusterToSuperCluster( PET *e );
	// set all variables about a cluster
	void SetProperties( Contig *c1, int o1, Contig *c2, int o2, int dis, int std, int size );
	// get all clusters forming this super cluster
	list<PET*>* GetClustersOfSuperCluster();

	void SetRepetitiveEdge( PET *edge );
	PET* GetRepetitiveEdge();
	void SetOriginalEdge( PET *edge );
	PET* GetOriginalEdge();
	void SetUniqueEdge( PET *edge );
	PET* GetUniqueEdge();

	// all original information
	void SetOriginalStartContig( Contig *c );
	Contig* GetOriginalStartContig();
	void SetOriginalEndContig( Contig *c );
	Contig* GetOriginalEndContig();
	void SetOriginalStartContigOri( int ori );
	int GetOriginalStartContigOri();
	void SetOriginalEndContigOri( int ori );
	int GetOriginalEndContigOri();
	void SetOriginalDistance( int dis );
	int GetOriginalDistance();
	Contig* GetOtherContigOriginal( Contig *c );	
	int GetOrientationOfContigOriginal( Contig *c );
	int GetPositionOfContigOriginal( Contig *c );

	void RemoveFromMultiLib();

	void SetRepetitiveEdgeOnly( PET *edge );
	void SetUniqueEdgeOnly( PET *edge );
	

	// Attributes
private:
	int m_id;
	Contig *m_startContig;
	int m_startContigOri;
	Contig *m_endContig;
	int m_endContigOri;
	int m_distance;
	int m_distanceWithGap;
	int m_std;
	int m_size;
	string m_oriCluster;		// the original cluster content
	bool m_isUnhappy;			// record if this edge is unhappy
	bool m_isInSubgraph;		// the flag to record if this edge is in subgraph
	bool m_isVisited;			// record if this edge is already been visited in deletion mode
	Contig *m_unusedContig;		// record which contig is not used if it is unhappy dangling edges
	bool m_isDE;				// record if it is a dangling edge
	list<PET*>::iterator m_happyDEIter;
	list<PET*>::iterator m_unhappyDEIter;
	string m_oriString;			// original cluster string

	// multiple libraries related
	list<PET*>::iterator m_allLibraryIter[2];		// save the iterator of the list saving all libraries information in contig
	int m_allLibraryEdgeType[2];					// record the type of the all library edge, right or left. Help for later search
	PetLibrary *m_lib;							// record the library this edge belonging to
	list<PET*>::iterator m_edgeIterMultiLib;	// record the iterator of this edge in library

	// super cluster
	list<PET*> *m_originalClusters;                 // the original clusters of this super cluster

	// corresponding repetitive edge
	PET *m_repetitiveEdge;
	// corresponding original edge
	PET *m_originalEdge;
	// corresponding unique edge (the repeat contig is combined with other nodes to form a unique node, update this repetitive contig in the edge)
	PET *m_uniqueEdge;

	// original information
	Contig *m_startContigOriginal;
	int m_startContigOriOriginal;
	Contig *m_endContigOriginal;
	int m_endContigOriOriginal;
	int m_distanceOriginal;
};
