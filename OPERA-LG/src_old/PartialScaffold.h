#pragma once

// the class of partial scaffold, including active region and dangling edges
#include <list>
#include <string>
#include <utility>
#include "Contig.h"
#include "PET.h"

using namespace std;

class PartialScaffold
{
public:
	PartialScaffold(void);
	PartialScaffold( PartialScaffold *p );
	~PartialScaffold(void);

	// Attributes
private:
	int m_numOfContigs;			// number of contigs in active region
	int m_numOfHappyDanglingEdges;			// the number of correct dangling edges
	int m_numOfUnhappyDanglingEdges;		// the number of unhappy dangling edges
	int m_numOfUnhappyEdges;				// the number of total unhappy edges
	bool m_isEmptyScaffold;				// check if it is an natural empty scaffold
	bool m_isMannualEmptyScaffold;		// check if it is a mannually created empty scaffold(break)

	// the changes to its parent partial scaffold
	list<Contig*> *m_removedContigsInAR;		// the removed contigs from active region
	Contig *m_addedContig;				// the added contig to active region
	list<PET*> *m_removedPetsInHappyDE;			// the removed edges from happy dangling edges
	list<PET*> *m_addedPetsInHappyDE;			// the added edges to happy dangling edges
	list<PET*> *m_removedPetsInUnhappyDE;			// the removed edges from unhappy dangling edges
	list<PET*> *m_addedPetsInUnhappyDE;			// the added edges to unhappy dangling edges
	list<PET*> *m_addedUnhappyEdges;			// the added unhappy edges
	list< pair<Contig*, int> > *m_addedUnassignedNodes;		// the added contigs in unassigned nodes
															// int represent the orientation
	list< pair<Contig*, int> > *m_removedUnassignedNodes;	// the removed contigs in unassigned nodes

	string m_ARString;
	string m_DEString;

	PartialScaffold *m_parent;					// the parent partial scaffold of this one

	// check if this scaffold is break or not
	bool m_ifBreak;								// record if this partial scaffold already break

	// record if this scaffold is the first scaffold
	bool m_isFirstScaffold;

	// record the previous unassigned nodes order
	bool m_ifModifyOrderOfUnassignedNodes;
	list< pair<Contig*, int> > *m_unassignedNodes;

	// record the position of the used unassigned node
	int m_posOfUnassignedNode;

	// record if the added contig is repeat
	bool m_addedContigIsRepeat;

	// Methods
public:
	// add a contig to active region
	void AddNode( Contig *c );
	// add happy dangling edge
	void AddAddedHappyDE( PET *p );
	// add unhappy dangling edge
	void AddAddedUnhappyDE( PET *p, Contig *c );
	// set active region string
	void SetARString( string s );
	// set unhappy dangling edges string
	void SetUnhappyDEString( string s );
	// add removed unassigned nodes
	void AddRemovedUnassignedNodes( Contig *c, int ori );
	// add added unassigned nodes
	void AddAddedUnassignedNodes( Contig *c, int ori );
	// set the parent scaffold
	void SetParent( PartialScaffold *p );
	// get the parent scaffold
	PartialScaffold* GetParent();
	// add removed edge from happy dangling edges
	void AddRemovedHappyDE( PET *p );
	// add removed edge from unhappy dangling edges
	void AddRmovedUnhappyDE( PET *p );
	// add unhappy edge number
	void AddUnhappyEdgeNumber( int i );
	// add added unhappy non-dangling edges
	void AddAddedUnhappyEdge( PET *p );
	// add removed contigs in active region
	void AddRemovedContigInAR( Contig *c );
	// check if current scaffold has been breaked
	bool IfBreaked();
	// set current scaffold as being breaked
	void Break();
	// set number of unhappy edges
	void SetNumOfUnhappyEdges( int i );

	// get the attributes
	Contig* GetAddedContigInAR();
	list<Contig*>* GetRemovedContigsInAR();
	list<PET*>* GetAddedHappyDE();
	list<PET*>* GetRemovedHappyDE();
	list<PET*>* GetAddedUnhappyDE();
	list<PET*>* GetRemovedUnhappyDE();
	list<PET*>* GetAddedUnhappyEdges();
	list< pair<Contig*, int> >* GetAddedUnassignedNodes();
	list< pair<Contig*, int> >* GetRemovedUnassignedNodes();
	// get the number of unhappy edges
	int GetNumOfUnhappyEdges();
	int GetNumOfContigsInAR();
	int GetNumOfHappyDE();
	int GetNumOfUnhappyDE();
	string GetScaffoldString();

	// minus 1 from number of active region
	void RemoveOneFromActiveRegion();

	// check and label if it is the mannualy generated empty scaffold
	void SetIfEmptyScaffold( bool empty );
	bool IfEmptyScaffold();
	void SetIfMannualEmptyScaffold( bool empty );
	bool IfMannualEmptyScaffold();

	string* GetARString();
	string* GetDEString();

	// methods related with the first scaffolds
	void SetFirstScaffold();
	bool IsFirstScaffold();

	// about saving the previous unassigned contigs order
	void SavePreviousUnassignedNodesOrder( list< pair<Contig*, int> > *order );
	list< pair<Contig*, int> > *GetPreviousUnassignedNodesOrder();
	bool IfRecalculatdUnassignedNodes();

	void IncreasePosOfUnassignedNode();
	int GetPosOfUnassignedNode();

	void SetTypeOfAddedContig( bool isRepeat );
	bool IfAddedContigIsRepeat();
};

