#include "Contig.h"

Contig::Contig(void)
{
	m_leftEdges = new list<PET*>;
	m_rightEdges = new list<PET*>;
	m_ifInSubgraph = false;
	m_gapSize = 0;
	m_scaffold = toScaffoldString();
	m_isBorderContig = false;
	m_extentionDirection = NEITHER;
	m_emptyList = new list< pair<Contig*, int> >;
	m_isInAR = false;

	m_name = "xxx";
	m_id = -1;

	m_leftEdgesMultiLib = new list<PET*>;
	m_rightEdgesMultiLib = new list<PET*>;
	
	m_numberOfBases = 0;

	m_distance = 0;
	m_ifCalculated = false;
	m_visitedTime = 0;
	m_weightedVisitedTime = 0;

	m_referenceGenomeName = "";            
	m_referenceStartPos = -1;                 
	m_referenceEndPos = -1;                   // the start position on reference genome
	m_referenceOri = -1;                      // the orientation of this contig on reference genome
	m_referenceIndex = -1;
	m_subgraphReferenceIndex = -1;

	m_scaffoldIDSet.clear();
	m_scafIDSet.clear();

	m_removedDis = 0;

	m_isRepeat = false;
	m_needToRemove = false;
}

void Contig::Initialize(){
	m_ifInSubgraph = false;
	m_gapSize = 0;
	m_isBorderContig = false;
	m_extentionDirection = NEITHER;
	m_isInAR = false;
	m_distance = 0;
	m_ifCalculated = false;
	m_visitedTime = 0;
	m_weightedVisitedTime = 0;
	m_subgraphReferenceIndex = -1;

	m_scaffoldIDSet.clear();
	m_scafIDSet.clear();
}

void Contig::InitializeDistance(){
	m_distance = 0;
	m_visitedTime = 0;
	m_ifCalculated = false;
	m_weightedVisitedTime = 0;
}


Contig::Contig( string name, double length, double cov ){
	this->m_name = name;
	this->m_length = length;
	this->m_lengthWithGap = length;
	this->m_cov = cov;
	this->m_ori = PLUS;

	m_leftEdges = new list<PET*>;
	m_rightEdges = new list<PET*>;
	m_ifInSubgraph = false;
	m_gapSize = 0;
	m_scaffold = toScaffoldString();
	m_isBorderContig = false;
	m_extentionDirection = NEITHER;
	m_emptyList = new list< pair<Contig*, int> >;
	m_isInAR = false;

	m_leftEdgesMultiLib = new list<PET*>;
	m_rightEdgesMultiLib = new list<PET*>;

	m_numberOfBases = 0;
	m_distance = 0;
	m_ifCalculated = false;
	m_visitedTime = 0;

	m_referenceGenomeName = "";            
	m_referenceStartPos = -1;                 
	m_referenceEndPos = -1;                   // the start position on reference genome
	m_referenceOri = -1;                      // the orientation of this contig on reference genome
	m_referenceIndex = -1;
	m_subgraphReferenceIndex = -1;

	m_scaffoldIDSet.clear();
	m_scafIDSet.clear();

	m_removedDis = 0;

	m_isRepeat = false;
	m_needToRemove = false;
}

Contig::~Contig(void)
{
	//if( !m_isRepeatCopy )
		DeleteEdges( m_leftEdges );
		//else
		//DeleteEdgesForRepeat( m_leftEdges );
	delete m_leftEdges;

	//if( !m_isRepeatCopy )
		DeleteEdges( m_rightEdges );
		//else
		//DeleteEdgesForRepeat( m_rightEdges );
	delete m_rightEdges;
	delete m_emptyList;

	delete m_leftEdgesMultiLib;
	delete m_rightEdgesMultiLib;

	for( list<SUBCONTIG*>::iterator subIter = m_subcontigs.begin(); subIter != m_subcontigs.end(); subIter++ ){
		delete *subIter;
	}
	m_subcontigs.clear();
}

// add 1 to the number of reads mapped to this contig
void Contig::AddOneMappedBases( int bases ){
	m_numberOfBases += bases;
}

// calcuate the coverage using number reads mapped on this contig
void Contig::CalculateCov(){
	m_cov = m_numberOfBases / m_length;
}

void Contig::DeleteEdges( list<PET*> *edges ){
	list<PET*>::iterator iter = edges->begin();
	while( iter != edges->end() ){
		if( !(*iter)->IfVisited() ){
			(*iter)->VisitEdge();

			if( (*iter)->GetOriginalEdge() != NULL ){
				(*iter)->GetOriginalEdge()->SetRepetitiveEdgeOnly( NULL );
				(*iter)->GetOriginalEdge()->SetUniqueEdgeOnly( NULL );
				delete (*iter);
			}
			else{
				(*iter)->SetRepetitiveEdge( NULL );
				(*iter)->SetUniqueEdge( NULL );
			}

			iter = edges->erase( iter );
		}
		else{
			PET* temp = *iter;
			iter = edges->erase( iter );
			delete temp;
		}
	}
}

// delete edges for repetitive contigs (no need to check occurrence of edges, directly delete)
void Contig::DeleteEdgesForRepeat( list<PET*> *edges ){
	list<PET*>::iterator iter = edges->begin();
	while( iter != edges->end() ){
		(*iter)->GetOriginalEdge()->SetRepetitiveEdge( NULL );
		(*iter)->GetOriginalEdge()->SetUniqueEdge( NULL );
		
		PET *tempEdge = *iter;
		iter = edges->erase( iter );

		delete tempEdge;
	}
}

void Contig::DeleteEdgesForRepeat( int ori ){
	list<PET*> *edges;
	if( ori == LEFT ){
		edges = m_leftEdges;
	}
	else
		edges = m_rightEdges;

	list<PET*>::iterator iter = edges->begin();
	while( iter != edges->end() ){
		//(*iter)->GetOriginalEdge()->SetRepetitiveEdge( NULL );
		//(*iter)->GetOriginalEdge()->SetUniqueEdge( NULL );
		
		PET *tempEdge = *iter;
		iter = edges->erase( iter );

		delete tempEdge;
	}
	
}

// clear the status of edges
void Contig::ClearStatusOfEdges( list<PET*> *edges ){
	list<PET*>::iterator iter = edges->begin();
	while( iter != edges->end() ){
		if( (*iter)->GetOriginalEdge() != NULL ){
			(*iter)->GetOriginalEdge()->SetRepetitiveEdge( NULL );
			(*iter)->GetOriginalEdge()->SetUniqueEdge( NULL );
		}
		(*iter)->SetRepetitiveEdge( NULL );
		(*iter)->SetUniqueEdge( NULL );
		iter++;
	}
}

// set contig name
void Contig::SetName( string name ){
	m_name = name;
}

// get contig name
string Contig::GetName(){
	return m_name;
}

// set contig length
void Contig::SetLength( double length ){
	m_length = length;
}

// get contig length
double Contig::GetLength(){
	return m_length;
}

// set contig coverage
void Contig::SetCov( double cov ){
	m_cov = cov;
}

// get contig coverage
double Contig::GetCov(){
	return m_cov;
}

// set contig id
void Contig::SetID( int id ){
	m_id = id;
}

// get contig ID
int Contig::GetID(){
	return m_id;
}


// add a pet cluster of this contig
void Contig::AddEdge( PET *p ){	
	// check if it is left edge
	if( ( p->GetEndContig() == this && p->GetEndContigOri() == PLUS )
		|| ( p->GetStartContig() == this && p->GetStartContigOri() == MINUS ) ){
			m_leftEdges->push_back( p );
	}
	else{
		m_rightEdges->push_back( p );
	}
}

// add a pet cluster of this contig, from all libraries
void Contig::AddEdgeMultiLib( PET *p ){
	// check if it is left edge
	if( ( p->GetEndContig() == this && p->GetEndContigOri() == PLUS )
		|| ( p->GetStartContig() == this && p->GetStartContigOri() == MINUS ) ){
			m_leftEdgesMultiLib->push_front( p );
			if( p->GetEndContig() == this )
				p->SetEdgeIterMultiLib( m_leftEdgesMultiLib->begin(), LEFT, END );
			else
				p->SetEdgeIterMultiLib( m_leftEdgesMultiLib->begin(), LEFT, START );
	}
	else{
		m_rightEdgesMultiLib->push_front( p );
		if( p->GetEndContig() == this )
			p->SetEdgeIterMultiLib( m_rightEdgesMultiLib->begin(), RIGHT, END );
		else
			p->SetEdgeIterMultiLib( m_rightEdgesMultiLib->begin(), RIGHT, START );
	}
}

// remove a pet cluster of this contig, from all libraries
void Contig::RemoveEdgeMultiLib( PET *p, int contigPos ){
	if( p->GetEdgeTypeMultiLib( contigPos ) == LEFT ){
		m_leftEdgesMultiLib->erase( p->GetEdgeIterMultiLib( contigPos ) );
	}
	else
		m_rightEdgesMultiLib->erase( p->GetEdgeIterMultiLib( contigPos ) );
}

// set the list iterator of current contig
void Contig::SetListPos( list<Contig*>::iterator iter ){
	m_listPos = iter;
}

// set the border list iterator of current contig
void Contig::SetBorderListPos( list<Contig*>::iterator iter ){
	m_borderListPos = iter;
	m_isBorderContig = true;
}

// generate the scaffold string
string Contig::toScaffoldString(){
	char buffer[ 1000 ];
	sprintf( buffer, "%s\tBE\t%.0f\t%d", m_name.c_str(), m_length, m_gapSize );
	string result( buffer );
	return result; 
}

// check if current contig is a singleton
bool Contig::IsSingleton(){
	if( m_leftEdges->empty() && m_rightEdges->empty() )
		return true;
	else
		return false;
}

// get the super contig(scaffold) id for this contig
int Contig::GetScaffoldID(){
	if( !m_isRepeat )
		return m_scaffoldID;
	else
		return *(m_scaffoldIDSet.begin());
}

list<PET*> * Contig::GetLeftEdges(){
	return m_leftEdges;
}

list<PET*> * Contig::GetRightEdges(){
	return m_rightEdges;
}

list<Contig*>::iterator Contig::GetListPos(){
	return m_listPos;
}

list<Contig*>::iterator Contig::GetBorderListPos(){
	return m_borderListPos;
}

// check if current contig has any edge
bool Contig::HasEdge(){
	return (!m_leftEdges->empty() || !m_rightEdges->empty());
}

// check if current contig is a border contig
bool Contig::IsBorderContig(){
	return m_isBorderContig;
}

// get the left distance
double Contig::GetLeftDis(){
	return m_leftDistance;
}

// get the right distance
double Contig::GetRightDis(){
	return m_rightDistance;
}

// set the direction of extension
void Contig::SetExtensionOri( int ori ){
	if( m_extentionDirection == NEITHER || ori == NEITHER )
		m_extentionDirection = ori;
	else if( ori == BOTH )
		m_extentionDirection = BOTH;
	else if( ori != m_extentionDirection )
		m_extentionDirection = BOTH;
}

// get the direction of extension
int Contig::GetExtensionOri(){
	return m_extentionDirection;
}

// label this contig as in subgraph
void Contig::SetInSubgraph(){
	m_ifInSubgraph = true;
}

// check if this contig is in subgraph
bool Contig::isInSubgraph(){
	return m_ifInSubgraph;
}

// check if this contig has left edges
bool Contig::HasLeftEdge(){
	return !m_leftEdges->empty();
}

// check if this contig has right edges
bool Contig::HasRightEdge(){
	return !m_rightEdges->empty();
}

// set orientation in scaffold
void Contig::SetOri( int ori ){
	m_ori = ori;
}
// get orientation in scaffold
int Contig::GetOri(){
	return m_ori;
}

// set start position with and without gap
void Contig::SetStartPosition( double pos, int gap ){
	m_startPosition = pos;
	m_startPositionWithGap = pos + (double)gap;
}

void Contig::SetStartPositionWithGap( double pos ){
	m_startPositionWithGap = pos;
}

void Contig::SetLeftDisWithGap( double d ){
	m_leftDistanceWithGap = d;
}

void Contig::SetRightDisWithGap( double d ){
	m_rightDistanceWithGap = d;
}

double Contig::GetRightDisWithGap(){
	return m_rightDistanceWithGap;
}

double Contig::GetLeftDisWithGap(){
	return m_leftDistanceWithGap;
}


// get start position without gap
double Contig::GetStartPosition(){
	return m_startPosition;
}

// get start position with gap
double Contig::GetStartPositionWithGap(){
	return m_startPositionWithGap;
}

// get the iterator in unassigned node set
list< pair<Contig*, int> >::iterator Contig::GetUnassignedNodeIter( int type ){
	return m_unassignedNodeIter[ type ];
}

// set unassigned node iterator
void Contig::SetUnassignedNodeIter( int type, list< pair<Contig*, int> >::iterator iter ){
	m_unassignedNodeIter[ type ] = iter;
	m_isUnassignedNode[ type ] = true;
}

// check if it is an unassigned node with certain orientation
bool Contig::IsUnassignedNode( int type ){
	return m_isUnassignedNode[ type ];
	/*if( m_unassignedNodeIter[ type ] != m_endOfList )
		return true;
	return false;*/
}

// set iter to end()
void Contig::ClearUnassignedNodeIter( int type ){
	//m_unassignedNodeIter[ type ] = m_endOfList;
	m_isUnassignedNode[ type ] = false;
}

// set if it is in active region
void Contig::SetIfInAR( bool ar ){
	m_isInAR = ar;
}

// check if it is in active region
bool Contig::IfInAR(){
	return m_isInAR;
}		

// check if it has happy dangling edges in right direction
bool Contig::IfHasHappyDE( int ori ){
	list<PET*> *edges;
	if( ori == PLUS )
		edges = m_rightEdges;
	else
		edges = m_leftEdges;

	for( list<PET*>::iterator iter = edges->begin(); iter != edges->end(); iter++){
		/*if( m_name == "178734" ){
			cerr<<"checking: "<<(*iter)->GetOriString()<<endl;
			cerr<<"\tIs de: "<<(*iter)->IsDE()<<"\tIs happy: "<<!(*iter)->IfUnhappy()<<endl;
		}
		*/
		if( (*iter)->IsDE() && !(*iter)->IfUnhappy() && (*iter)->isInSubgraph() )
			return true;
	}

	return false;
}

void Contig::SetEndOfList( list< pair<Contig*, int> >::iterator iter ){
	m_endOfList = iter;
}

// check the number of valid left Edge for certain orientation
// invalid edges include all dangling edges and unhappy edges
int Contig::GetNumOfValidLeftEdges( int ori ){
	int result = 0;
	list<PET*> *edges;
	if( ori == PLUS )
		edges = m_leftEdges;
	else
		edges = m_rightEdges;

	for( list<PET*>::iterator iter = edges->begin(); iter != edges->end(); iter++ ){
		if( (*iter)->isInSubgraph() && !(*iter)->IsDE() && !(*iter)->IfUnhappy() )
			result++;
	}

	return result;
}

// generate scaffold String with certain orientation
string Contig::GenScaffoldString( int ori ){
	string result = "";

	if( ori == PLUS ){
		list<string> *contigString = new list<string>;
		Split( m_scaffold, "\t", contigString );
		result = "";
		list<string>::reverse_iterator iter = contigString->rbegin();
		iter++;
		while( iter != contigString->rend() ){
			result = (*iter) + "\t" + result;
			iter++;
		}
		delete contigString;
		return result;
	}

	// generate the minus scaffold string
	list<string> *lines = new list<string>;
	Split( m_scaffold, "\n", lines );
	bool firstLine = true;
	//cout<<m_scaffold<<endl;
	//cout<<"lines: "<<lines->size()<<endl;
	for( list<string>::reverse_iterator iter = lines->rbegin(); iter != lines->rend(); iter++ ){
		string line = *iter;
		vector<string> *contigString = new vector<string>;
		//cout<<"\tgenerating contig string: "<<line<<endl;
		Split( line, "\t", contigString );
		//cout<<"\tsize of string: "<<contigString->size()<<endl;
		string tempString;
		if( firstLine )
			firstLine = false;
		else{
			result += contigString->at( 3 ) + "\n";
		}

		if( contigString->at( 1 ) == "BE" )
			tempString = contigString->at( 0 ) + "\tEB\t" + contigString->at( 2 ) + "\t";
		else
			tempString = contigString->at( 0 ) + "\tBE\t" + contigString->at( 2 ) + "\t";
		result += tempString;
		delete contigString;
	}

	//result += "0\n";

	delete lines;
	return result;
}

void Contig::SetLeftDis( double d ){
	m_leftDistance = d;
}
	
void Contig::SetRightDis( double d ){
	m_rightDistance = d;
}

void Contig::SetScaffoldID( int id ){
	m_scaffoldID = id;

	if( m_isRepeat ) {
		m_scaffoldIDSet.insert( id );
	}
}

// set scaffold string
void Contig::SetScaffoldString( string sca ){
	m_scaffold = sca;
}

void Contig::SetNotInSubgraph(){
	m_ifInSubgraph = false;
}

// get original contig orientation in scaffold
int Contig::GetContigOriInScaffold( string contigName ){
	list<string> lines;
	Split( m_scaffold, "\n", &lines );
	for( list<string>::iterator iter = lines.begin(); iter != lines.end(); iter++ ){
		vector<string> line;
		Split( *iter, "\t", &line );
		if( line.at( 0 ) == contigName ){
			if( line.at( 1 ) == "BE" )
				return PLUS;
			else
				return MINUS;
		}
	}
	return -1;
}

// reverse left and right edges
void Contig::ReverseEdges(){
	list<PET*> *temp = m_leftEdges;
	m_leftEdges = m_rightEdges;
	m_rightEdges = temp;
}

// get orientation of complex contig
string Contig::GetOriOfFirstContig( int ori ){
	list<string> line;
	Split( m_scaffold, "\t", &line );

	list<string>::iterator iter = line.begin();
	iter++;
	string result = *iter;

	if( ori == PLUS ){
		if( result == "BE" )
			return "BE";
		else
			return "EB";
	}
	else{
		if( result == "BE" )
			return "EB";
		else
			return "BE";
	}
}

string Contig::GetNameOfFirstContig(){
	list<string> line;
	Split( m_scaffold, "\t", &line );

	list<string>::iterator iter = line.begin();
	return *iter;
}

string Contig::GetScaffoldString(){
	return m_scaffold;
}

// check if old contig has the same orientation in new contig
bool Contig::CheckOldContigOri( Contig *c ){
	// get old contig information
	list<string> line;
	Split( c->GetScaffoldString(), "\t", &line );

	list<string>::iterator iter = line.begin();
	string oldContigName = *iter;
	iter++;
	string oldContigOri = *iter;

	// get new contig information
	list<string> lines;
	Split( m_scaffold, "\n", &lines );

	iter = lines.begin();
	while( iter != lines.end() ){
		vector<string> content;
		Split( *iter, "\t", &content );
		if( content.at( 0 ) == oldContigName ){
			if( content.at( 1 ) == oldContigOri )
				return true;
			else
				return false;
		}
		iter++;
	}

	return false;
}

// gap related
// set the id in scaffold
void Contig::SetScaffoldPos( int pos ){
	m_scaffoldPos = pos;
}

// get the id in scaffold
int Contig::GetScaffoldPos(){
	return m_scaffoldPos;
}

// set gap size after this contig
void Contig::SetGap( int gap ){
	m_gapSize = gap;
}

int Contig::GetGap(){
	return m_gapSize;
}

// heuristic related
void Contig::SetStep( int s ){
	m_step = s;
}

int Contig::GetStep(){
	return m_step;
}

// check if this contig has multiple libraries edge
bool Contig::HasMultiLibEdge()
{
	if( m_leftEdgesMultiLib->empty() && m_rightEdgesMultiLib->empty() )
		return false;
	else
		return true;
}

// get the left edge of multiple libraries
list<PET*>* Contig::GetLeftEdgeMultiLib()
{
	return m_leftEdgesMultiLib;
}

// get the right edge of multiple libraries
list<PET*>* Contig::GetRightEdgeMultiLib()
{
	return m_rightEdgesMultiLib;
}

// remove edge p from graph
void Contig::RemoveEdge( PET *p )
{
	m_leftEdges->remove( p );
	m_rightEdges->remove( p );
}

// set the distance of this contig to the first contig in the scaffold
void Contig::SetDistance( double d, int step )
{
	m_distance = d;
	m_ifCalculated = true;
	m_visitedTime = 1;
	m_weightedVisitedTime = 1;
	//m_weightedVisitedTime = 1.0/(double)step;
}

// get the distance of this contig to the first contig in the scaffold
double Contig::GetDistance()
{
	return m_distance;
}

// set the orientation of this contig during extension of finding unassigned nodes
void Contig::SetOriInExtension( int ori )
{
	m_oriInExtension = ori;
}

// get the orientation of this contig during extension of finding unassigned nodes
int Contig::GetOriInExtension()
{
	return m_oriInExtension;
}

// check if the distance of this contig has been calculated yet
bool Contig::IfCalculated(){
	return m_ifCalculated;
}

// set if the contigs have been calculated or not
void Contig::SetIfCalculated( bool c ){
	m_ifCalculated = c;
	if( c == false ){
		m_visitedTime = 0;
	}
}

// check if need to add a new distance to original distance
void Contig::AddDistance( double dis, int step ){
	if( m_type != REPEAT ){
		//if( abs( dis - m_distance ) < 10000 ){
		//m_distance = ( m_distance * m_visitedTime + dis ) / ( m_visitedTime + 1 );
		m_distance = m_distance * m_weightedVisitedTime;
		m_visitedTime++;
		m_weightedVisitedTime += 1.0/(double)step;
		m_distance = ( m_distance + dis * (1.0/(double)step) ) / m_weightedVisitedTime;
		//}
	}
	else{
		// it is a repeat. only remember the smallest position
		if( m_visitedTime == 0 ){
			m_distance = dis;
			m_visitedTime = 1;
			m_weightedVisitedTime = 1;
		}
		else{
			if( dis < m_distance ){
				m_distance = dis;
			}
		}
	}
}

void Contig::SetLengthWithGap( double l ){
	m_lengthWithGap = l;
}

double Contig::GetLengthWithGap(){
	return m_lengthWithGap;
}

void Contig::SetScafID( int id ){
	m_scafID = id;

	if( m_isRepeat ){
		m_scafIDSet.insert( id );
	}
}

int Contig::GetScafID(){
	if( !m_isRepeat )
		return m_scafID;
	else
		return *(m_scafIDSet.begin());
}

int Contig::GetVisitedTime(){
	return m_visitedTime;
}

// get the distance of the mid point
double Contig::GetMidPointDistance(){
	return m_distance + m_length / 2;      // return the position of the middle point
}

// methods about saving reference genome information
void Contig::SetReferenceName( string name ){
	m_referenceGenomeName = name;
}

string Contig::GetReferenceName(){
	return m_referenceGenomeName;
}

void Contig::SetReferenceStartPos( int pos ){
	m_referenceStartPos = pos;
}

int Contig::GetReferenceStartPos(){
	return m_referenceStartPos;
}

void Contig::SetReferenceEndPos( int pos ){
	m_referenceEndPos = pos;
}

int Contig::GetReferenceEndPos(){
	return m_referenceEndPos;
}

void Contig::SetReferenceOri( int ori ){
	m_referenceOri = ori;
}

int Contig::GetReferenceOri(){
	return m_referenceOri;
}
void Contig::SetReferenceIndex( int index ){
	m_referenceIndex = index;
}

int Contig::GetReferenceIndex(){
	return m_referenceIndex;
}

void Contig::SetSubgraphReferenceIndex( int index ){
	m_subgraphReferenceIndex = index;
}

int Contig::GetSubgraphReferenceIndex(){
	return m_subgraphReferenceIndex;
}

void Contig::SetContigType( int type ){
	m_type = type; 
}

int Contig::GetContigType(){
	return m_type;
}


void Contig::SetIfRepeat( bool type ){
	m_isRepeat = type;
}
	
bool Contig::IsRepeat(){
	return m_isRepeat;
}

// split the super contig and gain the individual information
// according to the orientation of current super contig
void Contig::SplitSuperContig( Contig *previousContig ){
	for( list<SUBCONTIG*>::iterator subIter = m_subcontigs.begin(); subIter != m_subcontigs.end(); subIter++ ){
		delete *subIter;
	}
	m_subcontigs.clear();

	// get the last contig of the previous super contig
	string lastContigName;
	int lastContigOri = 0;
	
	if( previousContig != NULL ){
		if( !previousContig->IsRepeat() ){
			SUBCONTIG *lastContig = previousContig->GetLastSubcontigs();
			lastContigName = lastContig->m_contigName;
			lastContigOri = lastContig->m_ori;
		}
		else{
			lastContigName = previousContig->GetName();
			lastContigOri = previousContig->GetOri();
		}
	}

	vector<string> lines;
	Split( m_scaffold, "\n", &lines );
	
	m_removedDis = 0;
	if( m_ori == PLUS ){
		for( int i = 0; i < (int) lines.size(); i++ ){
			vector<string> columns;
			Split( lines.at( i ), "\t", &columns );
			SUBCONTIG *newSubcontig = new SUBCONTIG;
				
			newSubcontig->m_contigName = columns.at( 0 );
			if( columns.at( 1 ) == "BE" )
				newSubcontig->m_ori = PLUS;
			else
				newSubcontig->m_ori = MINUS;
			
			newSubcontig->m_length = atof( columns.at( 2 ).c_str() );
			newSubcontig->m_gapSize = atoi( columns.at( 3 ).c_str() );

			if( previousContig != NULL && i == 0 && lastContigOri == newSubcontig->m_ori && lastContigName == columns.at( 0 ) ){
				// remove the first contig
				m_removedDis = atoi( columns.at( 2 ).c_str() );
				delete newSubcontig;
				continue;
			}
			else			
				m_subcontigs.push_back( newSubcontig );
		}
	}
	else{
		// reverse the scaffolds
		for( int i = lines.size() - 1; i >= 0; i-- ){
			vector<string> columns;
			Split( lines.at( i ), "\t", &columns );
			SUBCONTIG *newSubcontig = new SUBCONTIG;
				
			newSubcontig->m_contigName = columns.at( 0 );
			if( columns.at( 1 ) == "BE" )
				newSubcontig->m_ori = MINUS;
			else
				newSubcontig->m_ori = PLUS;
			
			newSubcontig->m_length = atof( columns.at( 2 ).c_str() );
			newSubcontig->m_gapSize = atoi( columns.at( 3 ).c_str() );

			if( previousContig != NULL 
			    && i == (int) lines.size() - 1 
			    && lastContigOri != newSubcontig->m_ori 
			    && lastContigName == columns.at( 0 ) ){
				// remove the first contig
				m_removedDis = atoi( columns.at( 2 ).c_str() );
				delete newSubcontig;
				continue;
			}
			else
				m_subcontigs.push_back( newSubcontig );
		}
	}
}   

// get the beginning/end subcontig in this super contig
SUBCONTIG* Contig::GetEndSubContig( int position, bool reverse ){
	vector<string> lines;
	Split( m_scaffold, "\n", &lines );
	
	SUBCONTIG *newSubcontig = new SUBCONTIG;
	if( position == START ){
		// get the first contig
		vector<string> columns;
		Split( lines.at( 0 ), "\t", &columns );
				
		newSubcontig->m_contigName = columns.at( 0 );
		if( !reverse ){
			if( columns.at( 1 ) == "BE" )
				newSubcontig->m_ori = PLUS;
			else
				newSubcontig->m_ori = MINUS;
		}
		else{
			if( columns.at( 1 ) == "BE" )
				newSubcontig->m_ori = MINUS;
			else
				newSubcontig->m_ori = PLUS;
		}
			
		newSubcontig->m_length = atof( columns.at( 2 ).c_str() );
	}
	else{
		// get the last subcontig
		vector<string> columns;
		Split( lines.at( lines.size() - 1 ), "\t", &columns );
				
		newSubcontig->m_contigName = columns.at( 0 );
		if( !reverse ){
			if( columns.at( 1 ) == "BE" )
				newSubcontig->m_ori = PLUS;
			else
				newSubcontig->m_ori = MINUS;
		}
		else{
			if( columns.at( 1 ) == "BE" )
				newSubcontig->m_ori = MINUS;
			else
				newSubcontig->m_ori = PLUS;
		}
		
		newSubcontig->m_length = atof( columns.at( 2 ).c_str() );
	}

	return newSubcontig;
}

SUBCONTIG* Contig::GetLastSubcontigs(){
	return *m_subcontigs.rbegin();
}

SUBCONTIG* Contig::GetFirstSubContig(){
	return *m_subcontigs.begin();
}

// check if this edge with a repeat can be satisfied or not
bool Contig::IsDanglingEdge( Contig *&ARContig, int &disDelta, int &disDeltaWithGap, PET *edge ){
	// get the repetitive contig and corresponding orientation
	Contig *cInAR = edge->FindARContig();
	int cInAROri = cInAR->GetOri();
	int cInAROriInPet = edge->GetOrientationOfContig( cInAR );

	// c is a repeat
	Contig *c = edge->GetOtherContig( cInAR );
	//cout<<c->GetName()<<endl;    
	int cOri;
	if( cInAROri == cInAROriInPet )
		cOri = edge->GetOrientationOfContig( c );
	else
		cOri = GetOppositeOri( edge->GetOrientationOfContig( c ) );

	
	// get the upperbound of the edge
	int dis = 0;
	int disWithGap = 0;
	int upperbound = edge->GetDis() + Configure::STD_TIMES * edge->GetStd();

	// check if the edge can be satisfied or not
	for( list<SUBCONTIG*>::iterator subIter = m_subcontigs.begin(); subIter != m_subcontigs.end(); subIter++ ){
		if( m_removedDis > 0 && subIter == m_subcontigs.begin() ){
			// don't consider the first contig if it is the same as the tail of the previous contig in active region
			continue;
		}
		
		if( dis > upperbound ){
			// do not find, this edge is still happy
			return true;
		}
		
		if( c->GetName() == (*subIter)->m_contigName && cOri == (*subIter)->m_ori ){
			// find the contig, create a unique edge
			//CreateUniqueEdge( cInAR, cInAROri, this, this->m_ori, dis, edge );
			ARContig = cInAR;
			disDelta = dis;
			disDeltaWithGap = disWithGap;
			return false;
		}
		
		dis += (*subIter)->m_length;
		disWithGap += (*subIter)->m_length;
		disWithGap += (*subIter)->m_gapSize;
	}

	return true;
}    

void Contig::SetOriginalRepeat( Contig *repeat ){
	m_originalRepeat = repeat;
}

Contig* Contig::GetOriginalRepeat(){
	return m_originalRepeat;
}

// get the removed distance due to the repeat at the start of supercontig
int Contig::GetRemovedDis(){
	return m_removedDis;
}            

// get all the sub contigs
list<SUBCONTIG*>* Contig::GetSubcontigs(){
	return &m_subcontigs;
}

void Contig::AddRepeatEdges( list<PET*> *edges, int type ){
	if( type == LEFT ){
		m_repeatLeftEdges->clear();
		m_repeatLeftEdges->insert(m_repeatLeftEdges->end(), edges->begin(), edges->end() );
	}
	else{
		m_repeatRightEdges->clear();
		m_repeatRightEdges->insert(m_repeatRightEdges->end(), edges->begin(), edges->end() );
	}
}

list<PET*>* Contig::GetRepeatEdges( int type ){
	if( type == LEFT ){
		return m_repeatLeftEdges;
	}
	else{
		return m_repeatRightEdges;
	}
}

void Contig::SetIfRemove( bool remove ){
	m_needToRemove = remove;
}

bool Contig::NeedRemove(){
	return m_needToRemove;
}

// delete corresponding repeat edge
void Contig::RemoveOriginalRepetitiveEdge( PET *edge, int type ){
	if( type == LEFT ){
		// check left edges
		for( list<PET*>::iterator edgeIter = m_leftEdges->begin(); edgeIter != m_leftEdges->end(); edgeIter++ ){
			if( edge == *edgeIter ){
				//PET *temp = *edgeIter;
				edgeIter = m_leftEdges->erase( edgeIter );
				//delete temp;
				break;
			}
		}
	}
	else{
		// check right edges
		for( list<PET*>::iterator edgeIter = m_rightEdges->begin(); edgeIter != m_rightEdges->end(); edgeIter++ ){
			if( edge == *edgeIter ){
				//PET *temp = *edgeIter;
				edgeIter = m_rightEdges->erase( edgeIter );
				//delete temp;
				break;
			}
		}
	}

	
}

// delete corresponding repeat edge
void Contig::RemoveOriginalRepetitiveEdgeMultiLib( PET *edge, int type ){
	if( type == LEFT ){
		// check left edges
		for( list<PET*>::iterator edgeIter = m_leftEdgesMultiLib->begin(); edgeIter != m_leftEdgesMultiLib->end(); edgeIter++ ){
			if( edge == *edgeIter ){
				//PET *temp = *edgeIter;
				edgeIter = m_leftEdgesMultiLib->erase( edgeIter );
				//delete temp;
				break;
			}
		}
	}
	else{
		// check right edges
		for( list<PET*>::iterator edgeIter = m_rightEdgesMultiLib->begin(); edgeIter != m_rightEdgesMultiLib->end(); edgeIter++ ){
			if( edge == *edgeIter ){
				//PET *temp = *edgeIter;
				edgeIter = m_rightEdgesMultiLib->erase( edgeIter );
				//delete temp;
				break;
			}
		}
	}

	
}

// traversing graph related methods
void Contig::SetConnectedSetPos( int setPos ){
	m_setPos = setPos;
}

int Contig::GetConnectedSetPos(){
	return m_setPos;
}


// create subcontig for repeat
void Contig::CreateSubContigForRepeat(){
	SUBCONTIG *newSubcontig = new SUBCONTIG;
				
	newSubcontig->m_contigName = m_name;
	newSubcontig->m_ori = m_ori;
	newSubcontig->m_length = m_length;
	newSubcontig->m_gapSize = m_gapSize;
	m_subcontigs.push_back( newSubcontig );
}

