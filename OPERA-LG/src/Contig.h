#pragma once

#include <string>
#include <list>
#include <vector>
#include <set>

class PET;
#include "PET.h"
#include "CommonFunction.h"

using namespace std;

struct SUBCONTIG{
	string m_contigName;
	int m_ori;
	double m_length;
	double m_startPos;
	int m_gapSize;
	bool m_isRepeat;
};

class Contig
{
public:
	Contig(void);
	Contig( string name, double length, double cov );
	~Contig(void);

	// Attributes
private:
	string m_name;		// contig name
	double m_length;		// contig length
	double m_cov;			// contig coverage
	int m_id;			// contig ID
	int m_type;                     // the type of the contig: unique, small or repeat
	list<PET*> *m_leftEdges;		// the left edges, suppose ori is +
	list<PET*> *m_rightEdges;		// the right edges, suppose ori is +
	list<Contig*>::iterator m_listPos;		// the list iterator of current contig
	list<Contig*>::iterator m_borderListPos;		// the border list iterator of current contig
	bool m_ifInSubgraph;		// the flag to record if this contig is in current subgraph
	double m_leftDistance;			// the left distance to the scaffold, suppose ori is +
	double m_rightDistance;		// the right distance to the scaffold, suppose ori is +
	string m_scaffold;			// the scaffold content of this contig
	int m_scaffoldID;
	bool m_isBorderContig;		// flag to record if it is a border contig
	int m_extentionDirection;		// the direction during extension

	bool m_isRepeat;                // flag if this contig is repeat or not
	
	// include gaps
	double m_lengthWithGap;          // the length of this contig including gap
	double m_leftDistanceWithGap;
	double m_rightDistanceWithGap;
	double m_startPositionWithGap;		// the start position in scaffold with gap

	// scaffold related
	int m_ori;			// orientation in scaffold
	double m_startPosition;			// the startPosition in scaffold without gap
	list< pair<Contig*, int> >::iterator m_unassignedNodeIter[2];		// iterator in unassigned node set, both plus and minus
	bool m_isUnassignedNode[2];			// record if it is unassigned node
	bool m_isInAR;				// recod if it is in active region
	list< pair<Contig*, int> >::iterator m_endOfList;	
	//*****unassigned nodes*****
	int m_unassignedOri;			// the untried orientation in unassigned nodes
	list< pair<Contig*, int> > *m_emptyList;			// the empty list, used for get end() iterator
	//*****end of unassigned nodes******
	//***********gap sizes related
	int m_scaffoldPos;			// the order in scaffolds
	int m_gapSize;				// the gap size after this scaffold
	//***********

	// heuristic parameters
	int m_step;			// the steps to find this contig

	// multiple libraries related
	list<PET*> *m_leftEdgesMultiLib;		// the left edges of all libraries, suppose ori is +
	list<PET*> *m_rightEdgesMultiLib;		// the right edges of all libraries, suppose ori is +

	double m_numberOfBases;                   // the number of bases mapped to this contig

	// calculate the distance from the first contig in the scaffold
	double m_distance;                       // the distance of this contig to the first contig in the scaffold
	int m_visitedTime;                         // record how many time this contig has been visited
	double m_weightedVisitedTime;            // record the weighted visited time
	int m_oriInExtension;                    // the orientation during extension
	bool m_ifCalculated;                     // if the distance of this contig has been calculated before
	
	// separate the scaffolds into independent components
	int m_scafID;

	// properties about the real positons on reference genome
	string m_referenceGenomeName;            // the reference genome name of this contig
	int m_referenceStartPos;                 // 
	int m_referenceEndPos;                   // the start position on reference genome
	int m_referenceOri;                      // the orientation of this contig on reference genome
	int m_referenceIndex;                    // the index of this contig on reference genome

	int m_subgraphReferenceIndex;            // the index of this contig within one subgraph on reference genome

	// super contig related properties
	list<SUBCONTIG*> m_subcontigs;           // subcontig information
	double m_removedDis;                     // record if it need to remove the first contig or not
	Contig *m_originalRepeat;                // the original repeat contig of this copy
	set<int> m_scaffoldIDSet;
	set<int> m_scafIDSet;

	list<PET*> *m_repeatLeftEdges;		// the left edges of this occurrence, suppose ori is +
	list<PET*> *m_repeatRightEdges;		// the right edges of this occurrence, suppose ori is +

	bool m_needToRemove;                    // record if this repeat needs to be removed from the final scaffold

	// for traversing graph
	int m_setPos;                    // the set containing a connected component


	// Methods
public:
	void SetName( string name );
	string GetName();
	void SetLength( double length );
	double GetLength();
	void SetCov( double cov );
	double GetCov();
	void SetID( int id );
	int GetID();
	void AddEdge( PET *p );			// add a pet cluster of this contig
	void SetListPos( list<Contig*>::iterator iter );
	void SetBorderListPos( list<Contig*>::iterator iter );
	bool IsSingleton();			// check if current contig is a singleton
	int GetScaffoldID();
	void SetScaffoldID( int id );
	list<PET*> * GetLeftEdges();
	list<PET*> * GetRightEdges();
	list<Contig*>::iterator GetListPos();
	list<Contig*>::iterator GetBorderListPos();
	bool HasEdge();			// check if current contig has any edge
	bool HasLeftEdge();
	bool HasRightEdge();
	bool IsBorderContig();			// check if current contig is a border contig
	double GetRightDis();
	double GetLeftDis();
	double GetRightDisWithGap();
	double GetLeftDisWithGap();
	void SetLeftDis( double d );
	void SetRightDis( double d );
	void SetLeftDisWithGap( double d );
	void SetRightDisWithGap( double d );
	void SetExtensionOri( int ori );
	int GetExtensionOri();
	void SetInSubgraph();
	void SetNotInSubgraph();
	bool isInSubgraph();
	void SetOri( int ori );
	int GetOri();
	void SetStartPosition( double pos, int gap );
	void SetStartPositionWithGap( double pos );
	double GetStartPosition();
	double GetStartPositionWithGap();
	list< pair<Contig*, int> >::iterator GetUnassignedNodeIter( int type );
	void SetUnassignedNodeIter( int type, list< pair<Contig*, int> >::iterator iter );
	bool IsUnassignedNode( int type );
	void ClearUnassignedNodeIter( int type );			// set iter to end()
	void SetIfInAR( bool ar );		// set if it is in active region
	bool IfInAR();				// check if it is in active region
	bool IfHasHappyDE( int ori );		// check if it has happy dangling edges in right direction
	void SetEndOfList( list< pair<Contig*, int> >::iterator iter );
	// check the number of valid left Edge for certain orientation
	int GetNumOfValidLeftEdges( int ori );
	// generate scaffold String with certain orientation
	string GenScaffoldString( int ori );
	// set scaffold string
	void SetScaffoldString( string sca );
	// get original contig orientation in scaffold
	int GetContigOriInScaffold( string contigName );
	// reverse left and right edges
	void ReverseEdges();
	// get orientation of complex contig
	string GetOriOfFirstContig( int ori );
	string GetNameOfFirstContig();
	bool CheckOldContigOri( Contig *c );
	string GetScaffoldString();

	// gap related
	void SetScaffoldPos( int pos );
	int GetScaffoldPos();
	void SetGap( int gap );
	int GetGap();
	void SetLengthWithGap( double l );
	double GetLengthWithGap();

	// heuristic related
	void SetStep( int s );
	int GetStep();

	// multiple library related
	void AddEdgeMultiLib( PET *p );			// add a pet cluster of this contig, from all libraries
	void RemoveEdgeMultiLib( PET *p, int contigPos );		// remove a pet cluster of this contig, from all libraries
	void Initialize();						// initialize the attributes
	bool HasMultiLibEdge();					// check if this contig has multiple libraries edge
	list<PET*>* GetLeftEdgeMultiLib();		// get the left edge of multiple libraries
	list<PET*>* GetRightEdgeMultiLib();		// get the right edge of multiple libraries


	// dynamically change cluster threshold related functions
	// remove edge p from graph
	void RemoveEdge( PET *p );

	// add 1 to the number of reads mapped to this contig
	void AddOneMappedBases( int bases );
	// calcuate the coverage using number reads mapped on this contig
	void CalculateCov();

	// calculate the distance from the first contig in the scaffold
	void SetDistance( double d, int step );
	double GetDistance();
	void SetOriInExtension( int ori );
	int GetOriInExtension();
	bool IfCalculated();
	void SetIfCalculated( bool c );
	void AddDistance( double dis, int step );

	// separate the scaffolds into independent components
	void SetScafID( int id );
	int GetScafID();

	void InitializeDistance();

	int GetVisitedTime();
	double GetMidPointDistance();        // get the distance of the mid point

	// methods about saving reference genome information
	void SetReferenceName( string name );
	string GetReferenceName();
	void SetReferenceStartPos( int pos );
	int GetReferenceStartPos();
	void SetReferenceEndPos( int pos );
	int GetReferenceEndPos();
	void SetReferenceOri( int ori );
	int GetReferenceOri();
	void SetReferenceIndex( int index );
	int GetReferenceIndex();
	void SetSubgraphReferenceIndex( int index );
	int GetSubgraphReferenceIndex();

	// handle contig type
	void SetContigType( int type );
	int GetContigType();
	void SetIfRepeat( bool type );
	bool IsRepeat();

	// super contigs related methods
	void SplitSuperContig( Contig *previousContig );              // split the super contig and gain the individual information
	bool IsDanglingEdge( Contig *&ARContig, int &disDelta, int &disDeltaWithGap, PET *edge );    // check if this edge with a repeat is still dangling edge or not
	int GetRemovedDis();                 // get the removed distance due to the repeat at the start of supercontig
	list<SUBCONTIG*>* GetSubcontigs();   // get all the sub contigs
	SUBCONTIG* GetLastSubcontigs();
	SUBCONTIG* GetFirstSubContig();
	void SetOriginalRepeat( Contig *repeat );
	Contig* GetOriginalRepeat();

	void AddRepeatEdges( list<PET*> *edges, int type );
	list<PET*>* GetRepeatEdges( int type );

	void SetIfRemove( bool remove );
	bool NeedRemove();

	// get the beginning/end subcontig in this super contig
	SUBCONTIG* GetEndSubContig( int position, bool reverse );

	// delete corresponding repeat edge
	void RemoveOriginalRepetitiveEdge( PET *edge, int type );
	// delete corresponding repeat edge (multiple libraries)
	void RemoveOriginalRepetitiveEdgeMultiLib( PET *edge, int type );

	// traversing graph related methods
	void SetConnectedSetPos( int setPos );
	int GetConnectedSetPos();

	// clear the status of edges
	void ClearStatusOfEdges( list<PET*> *edges );

	void DeleteEdges( list<PET*> *edges );
	
	// delete edges for repetitive contigs (no need to check occurrence of edges, directly delete)
	void DeleteEdgesForRepeat( list<PET*> *edges );
	void DeleteEdgesForRepeat( int ori );

	// create subcontig for repeat
	void CreateSubContigForRepeat();

private:
	string toScaffoldString();
	//inline void DeleteEdges( list<PET*> *edges );
};
