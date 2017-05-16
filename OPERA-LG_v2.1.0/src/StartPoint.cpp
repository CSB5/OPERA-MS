#include "StartPoint.h"

StartPoint::StartPoint(void)
{
}

StartPoint::StartPoint( Contig *c, int ori ){
	m_contig = c;
	m_ori = ori;

	list<PET*> *edges;

	if( m_ori == PLUS )
		edges = c->GetLeftEdges();
	else
		edges = c->GetRightEdges();

	m_numOfUnhappyEdges = 0;

	for( list<PET*>::iterator iter = edges->begin(); iter != edges->end(); iter++ ){
		if( (*iter)->isInSubgraph() )
			m_numOfUnhappyEdges++;
	}
}

StartPoint::~StartPoint(void)
{
}

int StartPoint::GetNumOfUnhappyEdges() const{
	return m_numOfUnhappyEdges;
}

Contig * StartPoint::GetContig() const{
	return m_contig;
}

int StartPoint::GetOri() const{
	return m_ori;
}

void StartPoint::SetNumOfUnhappyEdges( int n ){
	m_numOfUnhappyEdges = n;
}
