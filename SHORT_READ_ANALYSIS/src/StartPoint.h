#pragma once

#include "Contig.h"
#include "Configure.h"

class StartPoint
{
public:
	StartPoint(void);
	StartPoint( Contig *c, int ori );
	~StartPoint(void);

	// Attributes
private:
	Contig *m_contig;
	int m_numOfUnhappyEdges;
	int m_ori;

	// Methods
public:
	int GetNumOfUnhappyEdges() const;
	Contig * GetContig() const;
	int GetOri() const;
	void SetNumOfUnhappyEdges( int n );
};
