#include "PET.h"

PET::PET(void)
{
	m_originalClusters = new list<PET*>;
	m_repetitiveEdge = NULL;
	m_uniqueEdge = NULL;
	m_originalEdge = NULL;
}

PET::~PET(void)
{
	for( list<PET*>::iterator iter = m_originalClusters->begin(); iter != m_originalClusters->end(); iter++ ){
		delete *iter;
	}
	m_originalClusters->clear();
	delete m_originalClusters;

	if( !m_startContigOriginal->IsRepeat() )
		delete m_startContigOriginal;
	
	if( !m_endContigOriginal->IsRepeat() )
		delete m_endContigOriginal;

	if( m_uniqueEdge != NULL )
		delete m_uniqueEdge;

	if( m_repetitiveEdge != NULL )
		delete m_repetitiveEdge;
}

PET::PET( Contig *c1, int o1, Contig *c2, int o2, int dis, int std, int size ){
	//cerr<<c1->GetName()<<" "<<o1<<" "<<c2->GetName()<<" "<<dis<<" "<<std<<" "<<size<<"\n";
	//cerr<<" "<<o1<<" "<<" "<<dis<<" "<<std<<" "<<size<<"\n";
	m_startContig = c1;
	m_startContigOri = o1;
	m_endContig = c2;
	m_endContigOri = o2;
	m_distance = dis;
	m_distanceWithGap = dis;
	m_std = std;
	m_size = size;
	m_isVisited = false;
	m_isUnhappy = false;		
	m_isInSubgraph = false;
	m_isDE = false;
	m_unusedContig = NULL;

	GenClusterString();
	m_originalClusters = new list<PET*>;

	// set original information
	//cerr<<"set original information\n";
	if( m_startContig->IsRepeat() )
		m_startContigOriginal = m_startContig;
	else
		m_startContigOriginal = new Contig( m_startContig->GetName(), m_startContig->GetLength(), m_startContig->GetCov() );

	if( m_endContig->IsRepeat() )
		m_endContigOriginal = m_endContig;
	else
		m_endContigOriginal = new Contig( m_endContig->GetName(), m_endContig->GetLength(), m_endContig->GetCov() );
	
	//cerr<<"finish setting\n";

	m_startContigOriOriginal = m_startContigOri;
	m_endContigOriOriginal = m_endContigOri;
	m_distanceOriginal = m_distance;

	m_repetitiveEdge = NULL;
	m_uniqueEdge = NULL;
	m_originalEdge = NULL;
}

PET::PET( Contig *c1, int o1, Contig *c2, int o2, int dis, int std, int size, PET *edge ){
	m_startContig = c1;
	m_startContigOri = o1;
	m_endContig = c2;
	m_endContigOri = o2;
	m_distance = dis;
	m_distanceWithGap = dis;
	m_std = std;
	m_size = size;
	m_isVisited = false;
	m_isUnhappy = false;		
	m_isInSubgraph = false;
	m_isDE = false;
	m_unusedContig = NULL;

	GenClusterString();
	m_originalClusters = new list<PET*>;

	// set original information
	if( m_startContig->IsRepeat() )
		m_startContigOriginal = edge->GetOriginalStartContig();
	else
		m_startContigOriginal = new Contig( m_startContig->GetName(), m_startContig->GetLength(), m_startContig->GetCov() );

	if( m_endContig->IsRepeat() )
		m_endContigOriginal = edge->GetOriginalEndContig();
	else
		m_endContigOriginal = new Contig( m_endContig->GetName(), m_endContig->GetLength(), m_endContig->GetCov() );
	
	m_startContigOriOriginal = edge->GetOriginalStartContigOri();
	m_endContigOriOriginal = edge->GetOriginalEndContigOri();
	m_distanceOriginal = edge->GetOriginalDistance();

	m_repetitiveEdge = NULL;
	m_uniqueEdge = NULL;
	m_originalEdge = NULL;
}

// add one cluster to super cluster
void PET::AddOneClusterToSuperCluster( PET *e ){
	m_originalClusters->push_back( e );
}

// set all variables about a cluster
void PET::SetProperties( Contig *c1, int o1, Contig *c2, int o2, int dis, int std, int size ){
	m_startContig = c1;
	m_startContigOri = o1;
	m_endContig = c2;
	m_endContigOri = o2;
	m_distance = dis;
	m_distanceWithGap = dis;
	m_std = std;
	m_size = size;
	m_isVisited = false;
	m_isUnhappy = false;		
	m_isInSubgraph = false;
	m_isDE = false;
	m_unusedContig = NULL;

	GenClusterString();

	// set original information
	if( m_startContig->IsRepeat() )
		m_startContigOriginal = m_startContig;
	else
		m_startContigOriginal = new Contig( m_startContig->GetName(), m_startContig->GetLength(), m_startContig->GetCov() );

	if( m_endContig->IsRepeat() )
		m_endContigOriginal = m_endContig;
	else
		m_endContigOriginal = new Contig( m_endContig->GetName(), m_endContig->GetLength(), m_endContig->GetCov() );
	
	m_startContigOriOriginal = m_startContigOri;
	m_endContigOriOriginal = m_endContigOri;
	m_distanceOriginal = m_distance;
}

int PET::GetSize(){
	return m_size;
}

Contig * PET::GetStartContig(){
	return m_startContig;
}

Contig * PET::GetEndContig(){
	return m_endContig;
}

int PET::GetStartContigOri(){
	return m_startContigOri;
}

int PET::GetEndContigOri(){
	return m_endContigOri;
}

ostream& operator<<( ostream& out, const PET& p ){
	out<<p.m_startContig->GetID()<<"\t"<<p.m_endContigOri<<"\t"
		<<p.m_endContig->GetID()<<"\t"<<p.m_endContigOri<<"\t"
		<<p.m_distance<<"\t"<<p.m_std<<"\t"<<p.m_size;
	return out;
}

// output the string of current pet
string PET::ToString(){
	string result = m_startContig->GetName() + "\t";
	if( m_startContigOri == PLUS ) 
		result += "+\t";
	else
		result += "-\t";
	result += m_endContig->GetName() + "\t";
	if( m_endContigOri == PLUS )
		result += "+\t";
	else
		result += "-\t";
	result += itos( m_distance ) + "\t" + itos( m_std ) + "\t" + itos( m_size );
	return result;
}

// check if current edge is unhappy
bool PET::IfUnhappy(){
	return m_isUnhappy;
}

// check if current edge is in subgraph
bool PET::IfInSubgraph(){
	return m_isInSubgraph;
}

// check if this edge has already been visited in deletion mode of contig
bool PET::IfVisited(){
	return m_isVisited;
}

// visit the edge in deletion mode of contig for the first time
void PET::VisitEdge(){
	m_isVisited = true;
}

// get the positon of contig c: start or end
int PET::GetPositionOfContig( Contig *c ){
	if( m_startContig == c )
		return START;
	else
		return END;
}

// get the orientation of contig c
int PET::GetOrientationOfContig( Contig *c ){
	if( m_startContig == c )
		return m_startContigOri;
	else
		return m_endContigOri;
}

// replace the contig in pos using c
void PET::ReplaceContig( int pos, Contig *c ){
	if( pos == START )
		m_startContig = c;
	else
		m_endContig = c;
}

// set the orientation of contig in certain pos
void PET::SetOri( int pos, int ori ){
	if( pos == START )
		m_startContigOri = ori;
	else
		m_endContigOri = ori;
}

void PET::SetDis( int dis ){
	m_distance = dis; 
}

int PET::GetDis(){
	return m_distance;
}

void PET::SetDisWithGap( int dis ){
	m_distanceWithGap = dis;
}

int PET::GetDisWithGap(){
	return m_distanceWithGap;
}

// label this edge as in subgraph
void PET::SetInSubgraph(){
	m_isInSubgraph = true;
}

// label this edge as not in subgraph
void PET::SetNotInSubgraph(){
	m_isInSubgraph = false;
}

// check if this edge is in subgraph
bool PET::isInSubgraph(){
	return m_isInSubgraph;
}

// get the other contig rather than c
Contig* PET::GetOtherContig( Contig *c ){
	if( m_startContig == c )
		return m_endContig;
	else
		return m_startContig;
}

// label the edge as unhappy
// c is the unused contigs
void PET::SetUnhappy( Contig *c ){
	m_isUnhappy = true;
	m_unusedContig = c;
}

// label the edge as happy
void PET::SetHappy(){
	m_isUnhappy = false;
	m_unusedContig = NULL;
}

int PET::GetID(){
	return m_id;
}

void PET::SetID( int id ){
	m_id = id;
}

// set if it is a dangling edge
void PET::SetDE( bool isDE ){
	m_isDE = isDE;
}

// check if it is a dangling edge
bool PET::IsDE(){
	return m_isDE;
}

// get the standard deviation of this edge
int PET::GetStd() const{
	return m_std;
}

// set the happy dangling edge iterator
void PET::SetHappyDEIter( list<PET*>::iterator iter ){
	m_happyDEIter = iter;
}

// get the happy dangling edge iterator
list<PET*>::iterator PET::GetHappyDEIter(){
	return m_happyDEIter;
}

// set the unhappy dangling edge iterator
void PET::SetUnhappyDEIter( list<PET*>::iterator iter ){
	m_unhappyDEIter = iter;
}

// get the unhappy danglin edge iterator
list<PET*>::iterator PET::GetUnhappyDEIter(){
	return m_unhappyDEIter;
}

// find which contig is in active region
Contig* PET::FindARContig(){
	if( m_startContig->IsRepeat() ){
		if( m_endContig->IfInAR() )
			return m_endContig;
		else
			return m_startContig;
	}

	if( m_endContig->IsRepeat() ){
		if( m_startContig->IfInAR() )
			return m_startContig;
		else
			return m_endContig;
	}

	if( m_startContig->IfInAR() )
		return m_startContig;
	else
		return m_endContig;
}

// set this edge as unhappy
void PET::SetUnhappy(){
	m_isUnhappy = true;
}

// get unused contig
Contig* PET::GetUnusedContig(){
	return m_unusedContig;
}

void PET::SetOriString( string s ){
	m_oriString = s;
}

string PET::GetOriString(){
	return m_oriString;
}

// generate original cluster string
void PET::GenClusterString(){
	m_oriString = m_startContig->GetName() + "\t";
	if( m_startContigOri == PLUS )
		m_oriString += "+\t";
	else
		m_oriString += "-\t";

	m_oriString += m_endContig->GetName() + "\t";
	if( m_endContigOri == PLUS )
		m_oriString += "+\t";
	else
		m_oriString += "-\t";

	m_oriString += itos( m_distance ) + "\t" + itos( m_std ) + "\t" + itos( m_size );
}

// set the edge iterater, edges in all libraries
void PET::SetEdgeIterMultiLib( list<PET*>::iterator iter, int type, int contigPos ){
	m_allLibraryIter[ contigPos ] = iter;
	m_allLibraryEdgeType[ contigPos ] = type;
}

list<PET*>::iterator PET::GetEdgeIterMultiLib( int contigPos ){
	return m_allLibraryIter[ contigPos ];
}

int PET::GetEdgeTypeMultiLib( int contigPos ){
	return m_allLibraryEdgeType[ contigPos ];
}

// add the corresponding library and iterator of edges
void PET::AddIterInLib( PetLibrary *lib, list<PET*>::iterator iter ){
	m_lib = lib;
	m_edgeIterMultiLib = iter;
}

// delete this edge from library list
void PET::DeleteFromLib(){
	m_lib->GetEdges()->erase( m_edgeIterMultiLib );
}

// initialize all properties
void PET::Initialize()
{
	m_isVisited = false;
	m_isUnhappy = false;		
	m_isInSubgraph = false;
	m_isDE = false;
	m_unusedContig = NULL;
}

// calculate anchor region
void PET::AddAnchorString( string s ){
	m_oriString += s;
}

// get all clusters forming this super cluster
list<PET*>* PET::GetClustersOfSuperCluster(){
	return m_originalClusters;
}

void PET::SetRepetitiveEdge( PET *edge ){
	//if( edge == NULL && m_repetitiveEdge != NULL ){
	if( edge != NULL && m_repetitiveEdge != NULL ){
		//cerr<<"Try to delete previous repetitive edge\n";
		delete m_repetitiveEdge;
		//cerr<<"Done\n";
	}

	if( edge == NULL )
		m_repetitiveEdge = NULL;
	else
		m_repetitiveEdge = edge;
}

PET* PET::GetRepetitiveEdge(){
	return m_repetitiveEdge;
}

void PET::SetOriginalEdge( PET *edge ){
	m_originalEdge = edge;
}

PET* PET::GetOriginalEdge(){
	return m_originalEdge;
}

void PET::SetOriginalStartContig( Contig *c ){
	m_startContigOriginal = c;
}

Contig* PET::GetOriginalStartContig(){
	return m_startContigOriginal;
}

void PET::SetOriginalEndContig( Contig *c ){
	m_endContigOriginal = c;
}

Contig* PET::GetOriginalEndContig(){
	return m_endContigOriginal;
}

void PET::SetOriginalStartContigOri( int ori ){
	m_startContigOriOriginal = ori;
}

int PET::GetOriginalStartContigOri(){
	return m_startContigOriOriginal;
}

void PET::SetOriginalEndContigOri( int ori ){
	m_endContigOriOriginal = ori;
}

int PET::GetOriginalEndContigOri(){
	return m_endContigOriOriginal;
}

void PET::SetOriginalDistance( int dis ){
	m_distanceOriginal = dis;
}

int PET::GetOriginalDistance(){
	return m_distanceOriginal;
}


Contig* PET::GetOtherContigOriginal( Contig *c ){
	if( m_startContigOriginal == c )
		return m_endContigOriginal;
	else
		return m_startContigOriginal;
}	

int PET::GetOrientationOfContigOriginal( Contig *c ){
	if( m_startContigOriginal == c )
		return m_startContigOriOriginal;
	else
		return m_endContigOriOriginal;
}

int PET::GetPositionOfContigOriginal( Contig *c ){
	if( m_startContigOriginal == c )
		return START;
	else
		return END;
}

void PET::RemoveFromMultiLib(){
	m_lib->RemoveEdge( m_edgeIterMultiLib );
}

void PET::SetUniqueEdge( PET *edge ){
	//if( edge == NULL && m_uniqueEdge != NULL ){
	if( m_uniqueEdge != NULL ){
		//cerr<<"unique edge:\t"<<m_uniqueEdge->GetStartContig()->GetName()<<"\t"<<m_uniqueEdge->GetEndContig()->GetName()<<endl;
		delete m_uniqueEdge;
	}

	if( edge == NULL )
		m_uniqueEdge = NULL;
	else
		m_uniqueEdge = edge;
}

PET* PET::GetUniqueEdge(){
	return m_uniqueEdge;
}

void PET::SetRepetitiveEdgeOnly( PET *edge ){
	m_repetitiveEdge = edge;
}

void PET::SetUniqueEdgeOnly( PET *edge ){
	m_uniqueEdge = edge;
}
