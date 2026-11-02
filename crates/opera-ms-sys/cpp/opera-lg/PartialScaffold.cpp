#include "PartialScaffold.h"

PartialScaffold::PartialScaffold(void)
{
	m_removedContigsInAR = new list<Contig*>;		
	m_addedContig = NULL;	
	m_removedPetsInHappyDE = new list<PET*>;	
	m_addedPetsInHappyDE = new list<PET*>;	
	m_removedPetsInUnhappyDE = new list<PET*>;	
	m_addedPetsInUnhappyDE = new list<PET*>;
	m_addedUnhappyEdges = new list<PET*>;	
	m_addedUnassignedNodes = new list< pair<Contig*, int> >;
	m_removedUnassignedNodes = new list< pair<Contig*, int> >;
	m_parent = NULL;
	m_numOfUnhappyEdges = 0;
	m_numOfContigs = 0;
	m_numOfHappyDanglingEdges = 0;
	m_numOfUnhappyDanglingEdges = 0;
	m_ifBreak = false;
	m_isEmptyScaffold = false;
	m_isMannualEmptyScaffold = false;

	m_isFirstScaffold = false;
	m_ifModifyOrderOfUnassignedNodes = false;
	m_unassignedNodes = new list< pair<Contig*, int> >;
	
	m_posOfUnassignedNode = 0;

	m_addedContigIsRepeat = false;
}

PartialScaffold::PartialScaffold( PartialScaffold *p ){
	m_removedContigsInAR = new list<Contig*>;		
	m_addedContig = NULL;	
	m_removedPetsInHappyDE = new list<PET*>;	
	m_addedPetsInHappyDE = new list<PET*>;	
	m_removedPetsInUnhappyDE = new list<PET*>;	
	m_addedPetsInUnhappyDE = new list<PET*>;
	m_addedUnhappyEdges = new list<PET*>;	
	m_addedUnassignedNodes = new list< pair<Contig*, int> >;
	m_removedUnassignedNodes = new list< pair<Contig*, int> >;
	m_parent = NULL;
	m_numOfUnhappyEdges = p->GetNumOfUnhappyEdges();
	m_numOfContigs = p->GetNumOfContigsInAR();
	m_numOfHappyDanglingEdges = p->GetNumOfHappyDE();
	m_numOfUnhappyDanglingEdges = p->GetNumOfUnhappyDE();
	m_ifBreak = false;
	m_isEmptyScaffold = false;
	m_isMannualEmptyScaffold = false;
	m_isFirstScaffold = false;
	m_ifModifyOrderOfUnassignedNodes = false;
	m_unassignedNodes = new list< pair<Contig*, int> >;
	m_posOfUnassignedNode = 0;

	m_addedContigIsRepeat = false;
}

PartialScaffold::~PartialScaffold(void)
{
	m_removedContigsInAR->clear();		delete m_removedContigsInAR;		
	m_removedPetsInHappyDE->clear();	delete m_removedPetsInHappyDE;	
	m_addedPetsInHappyDE->clear();		delete m_addedPetsInHappyDE;	
	m_removedPetsInUnhappyDE->clear();  delete m_removedPetsInUnhappyDE;	
	m_addedPetsInUnhappyDE->clear();	delete m_addedPetsInUnhappyDE;
	m_addedUnhappyEdges->clear();		delete m_addedUnhappyEdges;	
	m_addedUnassignedNodes->clear();	delete m_addedUnassignedNodes;
	m_removedUnassignedNodes->clear();	delete m_removedUnassignedNodes;
	m_unassignedNodes->clear();             delete m_unassignedNodes;
}

// add a contig to active region
void PartialScaffold::AddNode( Contig *c ){
	c->SetIfInAR( true );
	m_numOfContigs++;
	m_addedContig = c;
}
	
// add happy dangling edge
void PartialScaffold::AddAddedHappyDE( PET *p ){
	m_addedPetsInHappyDE->push_back( p );
	m_numOfHappyDanglingEdges++;
	p->SetDE( true );
}

// add unhappy dangling edge
// c is the contig which is not used
void PartialScaffold::AddAddedUnhappyDE( PET *p, Contig *c ){
	m_addedPetsInUnhappyDE->push_back( p );
	m_numOfUnhappyDanglingEdges++;
	m_numOfUnhappyEdges++;
	p->SetDE( true );
	p->SetUnhappy( c );
}

// set active region string
void PartialScaffold::SetARString( string s ){
	m_ARString = s;
}

// set unhappy dangling edges string
void PartialScaffold::SetUnhappyDEString( string s ){
	m_DEString = s;
}

// get the number of unhappy edges
int PartialScaffold::GetNumOfUnhappyEdges(){
	return m_numOfUnhappyEdges;
}

// add removed unassigned nodes
void PartialScaffold::AddRemovedUnassignedNodes( Contig *c, int ori ){
	pair<Contig*, int> p( c, ori );
	m_removedUnassignedNodes->push_back( p );
}

// add added unassigned nodes
void PartialScaffold::AddAddedUnassignedNodes( Contig *c, int ori ){
	pair<Contig*, int> p( c, ori );
	m_addedUnassignedNodes->push_front( p );
}

// set the parent scaffold
void PartialScaffold::SetParent( PartialScaffold *p ){
	m_parent = p;
}
	
// get the parent scaffold
PartialScaffold* PartialScaffold::GetParent(){
	return m_parent;
}

// add removed edge from happy dangling edges
void PartialScaffold::AddRemovedHappyDE( PET *p ){

	if( p->GetOriginalStartContig()->GetName() == "178734" && p->GetOriginalEndContig()->GetName() == "178762" ){
		cerr<<"found removed happy de\n";
	}
	m_removedPetsInHappyDE->push_back( p );
	p->SetDE( false );
	if( m_numOfHappyDanglingEdges > 0 )
		m_numOfHappyDanglingEdges--;
}

// add unhappy edge number
void PartialScaffold::AddUnhappyEdgeNumber( int i ){
	m_numOfUnhappyEdges += i;
}

// add removed edge from unhappy dangling edges
void PartialScaffold::AddRmovedUnhappyDE( PET *p ){
	m_removedPetsInUnhappyDE->push_back( p );
	p->SetDE( false );
	if( m_numOfUnhappyDanglingEdges > 0 )
		m_numOfUnhappyDanglingEdges--;
}

// add added unhappy non-dangling edges
void PartialScaffold::AddAddedUnhappyEdge( PET *p ){
	m_addedUnhappyEdges->push_back( p );
	p->SetUnhappy();
}

// add removed contigs in active region
void PartialScaffold::AddRemovedContigInAR( Contig *c ){
	m_removedContigsInAR->push_back( c );
	c->SetIfInAR( false );
	m_numOfContigs--;
	if( m_numOfContigs < 0 )
		m_numOfContigs = 0;
}

Contig* PartialScaffold::GetAddedContigInAR(){
	return m_addedContig;
}

list<Contig*>* PartialScaffold::GetRemovedContigsInAR(){
	return m_removedContigsInAR;
}

list<PET*>* PartialScaffold::GetAddedHappyDE(){
	return m_addedPetsInHappyDE;
}

list<PET*>* PartialScaffold::GetRemovedHappyDE(){
	return m_removedPetsInHappyDE;
}

list<PET*>* PartialScaffold::GetAddedUnhappyDE(){
	return m_addedPetsInUnhappyDE;
}

list<PET*>* PartialScaffold::GetRemovedUnhappyDE(){
	return m_removedPetsInUnhappyDE;
}

list<PET*>* PartialScaffold::GetAddedUnhappyEdges(){
	return m_addedUnhappyEdges;
}

list< pair<Contig*, int> >* PartialScaffold::GetAddedUnassignedNodes(){
	return m_addedUnassignedNodes;
}

list< pair<Contig*, int> >* PartialScaffold::GetRemovedUnassignedNodes(){
	return m_removedUnassignedNodes;
}

// check if current scaffold has been breaked
bool PartialScaffold::IfBreaked(){
	return m_ifBreak;
}

// set current scaffold as being breaked
void PartialScaffold::Break(){
	m_ifBreak = true;
}

// set number of unhappy edges
void PartialScaffold::SetNumOfUnhappyEdges( int i ){
	m_numOfUnhappyEdges = i;
}

int PartialScaffold::GetNumOfContigsInAR(){
	return m_numOfContigs;
}

int PartialScaffold::GetNumOfHappyDE(){
	return m_numOfHappyDanglingEdges;
}

int PartialScaffold::GetNumOfUnhappyDE(){
	return m_numOfUnhappyDanglingEdges;
}

string* PartialScaffold::GetARString(){
	return &m_ARString;
}

string* PartialScaffold::GetDEString(){
	return &m_DEString;
}

string PartialScaffold::GetScaffoldString(){
	return m_ARString + "\t" + m_DEString;
}

// minus 1 from number of active region
void PartialScaffold::RemoveOneFromActiveRegion(){
	this->m_numOfContigs--;
	if( m_numOfContigs < 0 )
		m_numOfContigs = 0;
}

// check and label if it is the mannualy generated empty scaffold
void PartialScaffold::SetIfEmptyScaffold( bool empty ){
	m_isEmptyScaffold = empty;
}

bool PartialScaffold::IfEmptyScaffold(){
	return m_isEmptyScaffold;
}

void PartialScaffold::SetIfMannualEmptyScaffold( bool empty ){
	m_isMannualEmptyScaffold = empty;
}

bool PartialScaffold::IfMannualEmptyScaffold(){
	return m_isMannualEmptyScaffold;
}

// methods related with the first scaffolds
void PartialScaffold::SetFirstScaffold(){
	m_isFirstScaffold = true; 
}


bool PartialScaffold::IsFirstScaffold(){
	return m_isFirstScaffold;
}

// about saving the previous unassigned contigs order
void PartialScaffold::SavePreviousUnassignedNodesOrder( list< pair<Contig*, int> > *order ){
	m_ifModifyOrderOfUnassignedNodes = true;
	m_unassignedNodes->clear();
	for( list< pair<Contig*, int> >::iterator iter = order->begin(); iter != order->end(); iter++ ){
		m_unassignedNodes->push_back( *iter );
	}
}

list< pair<Contig*, int> > * PartialScaffold::GetPreviousUnassignedNodesOrder(){
	return m_unassignedNodes;
}

bool PartialScaffold::IfRecalculatdUnassignedNodes(){
	return m_ifModifyOrderOfUnassignedNodes;
}

void PartialScaffold::IncreasePosOfUnassignedNode(){
	m_posOfUnassignedNode++;
}

int PartialScaffold::GetPosOfUnassignedNode(){
	return m_posOfUnassignedNode;
}

void PartialScaffold::SetTypeOfAddedContig( bool isRepeat ){
	m_addedContigIsRepeat = isRepeat;
}

bool PartialScaffold::IfAddedContigIsRepeat(){
	return m_addedContigIsRepeat = false;
}
