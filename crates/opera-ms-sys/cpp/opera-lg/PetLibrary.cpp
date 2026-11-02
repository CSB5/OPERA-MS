#include "PetLibrary.h"


PetLibrary::PetLibrary(string name)
{
	m_fileName = name;
	m_singlePetsMap = new map<pair<int, int>, multiset<SinglePet*, lessDistance>*>();
	m_edges = new list<PET*>();
	//m_threshold = -1;
	m_possibilityOfDistance = new vector<double>();
	m_mapInfo = new list<mappingInfo*>();
	m_mean = 0;
	m_std = 0;
	m_ori = 0;
	m_fileNameWithoutPath = "";
	m_minDistance = 0;
	m_maxDistance = 0;
	m_readLength = 0;

}


PetLibrary::~PetLibrary(void)
{
	delete m_singlePetsMap;
	m_edges->clear();
	delete m_edges;

	m_possibilityOfDistance->clear();
	delete m_possibilityOfDistance;

	for( list<mappingInfo*>::iterator iter = m_mapInfo->begin(); iter != m_mapInfo->end(); iter++ ){
		delete *iter;
	}
	m_mapInfo->clear();
	delete m_mapInfo;
}

void PetLibrary::SetMean( int mean ){
	m_mean = mean;
}

/*void PetLibrary::SetThreshold( int thres )
{
	m_threshold = thres;
}

int PetLibrary::GetThreshold( ){
	return m_threshold;
}
*/
	
void PetLibrary::SetStd( int std ){
	m_std = std;
}

void PetLibrary::SetOri( int ori ){
	m_ori = ori;
}

int PetLibrary::GetOri(){
	return m_ori;
}

int PetLibrary::GetMean(){
	return m_mean;
}

int PetLibrary::GetStd(){
	return m_std;
}

// check if a contig pair exist in map
void PetLibrary::InsertPair( int c1, int c2, SinglePet *pet ){
	pair<int, int> newPair( c1, c2 );
	map<pair<int, int>, multiset<SinglePet*, lessDistance>*>::iterator mapIter = m_singlePetsMap->find( newPair );
	if( mapIter == m_singlePetsMap->end() ){
		// new pair, create a new set
		multiset<SinglePet*, lessDistance> *newSet = new multiset<SinglePet*, lessDistance>;
		newSet->insert( pet );
		m_singlePetsMap->insert( pair<pair<int, int>, multiset<SinglePet*, lessDistance>*>( newPair, newSet ) );
	}
	else{
		// insert the single pet to previous set
		mapIter->second->insert( pet );
	}
}

void PetLibrary::SetNameWithoutPath( string name ){
	m_fileNameWithoutPath = name;
}

string PetLibrary::GetNameWithoutPath(){
	return m_fileNameWithoutPath;
}

// get the pair map
map<pair<int, int>, multiset<SinglePet*, lessDistance>*>* PetLibrary::GetSinglePetsMap(){
	return m_singlePetsMap;
}

// add a cluster
void PetLibrary::AddCluster( PET *pet ){
	m_edges->push_front( pet );
	pet->AddIterInLib( this, m_edges->begin() );
}

// get the pet clusers
list<PET*>* PetLibrary::GetEdges(){
	return m_edges;
}

// get file name
string PetLibrary::GetFileName()
{
	return this->m_fileName;
}

// initialize the possibility array of all possible distance
void PetLibrary::InitPossibility( int min, int max )
{
	m_minDistance = min;
	m_maxDistance = max;

	for( int i = 0; i <= max; i++ ){
		m_possibilityOfDistance->push_back( 0 );
	}
}

// add one occurrance of a certain distance
void PetLibrary::AddOccuOfDistance( int dis )
{
	(*m_possibilityOfDistance)[ dis ]++;
}

// Calculate possibility
void PetLibrary::CalculatePossibility()
{
	double sum = 0;
	for( int i = 0; i <= m_maxDistance; i++ ){
		sum += m_possibilityOfDistance->at( i );
	}

	for( int i = 0; i <= m_maxDistance; i++ ){
		(*m_possibilityOfDistance)[ i ] = (*m_possibilityOfDistance)[ i ] / sum;
		//cout<<i<<"\t"<<(*m_possibilityOfDistance)[ i ]<<endl;
	}
}

// get the vector of possibility
vector<double>* PetLibrary::GetPossibilities()
{
	return m_possibilityOfDistance;
}

int PetLibrary::GetMinDis()
{
	return m_minDistance;
}


int PetLibrary::GetMaxDis()
{
	return m_maxDistance;
}

void PetLibrary::SetReadLength( int l )
{
	m_readLength = l;
}

int PetLibrary::GetReadLength()
{
	return m_readLength;
}

// add one mapping information 
void PetLibrary::AddMapInfo( Contig *c, int pos, string ori, int readLength )
{
	mappingInfo *newInfo = new mappingInfo;
	newInfo->contig = c;
	newInfo->position = pos;
	newInfo->orientation = ori;
	newInfo->readLength = readLength;

	m_mapInfo->push_back( newInfo );
}  

// get the mapping information
list<mappingInfo*>* PetLibrary::GetMappingInfo()
{
	return m_mapInfo;
}

void PetLibrary::RemoveEdge( list<PET*>::iterator iter ){
	m_edges->erase( iter );
}
