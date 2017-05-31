#include "LibInfo.h"


LibInfo::LibInfo(void)
{
	m_fileName = "";
	m_mean = -1;
	m_std = -1;
	m_threshold = 5;
	m_ori = IN;
	m_mapType = BOWTIE;
}


LibInfo::~LibInfo(void)
{
}

string LibInfo::GetFileName()
{
	return m_fileName;
}

int LibInfo::GetMean()
{ 
	return m_mean;
}
	
int LibInfo::GetStd()
{
	return m_std;
}
	
int LibInfo::GetThreshold()
{
	return m_threshold;
}

int LibInfo::GetOri()
{
	return m_ori;
}

int LibInfo::GetMapType()
{
	return m_mapType;
}

void LibInfo::SetFileName( string fileName )
{
	m_fileName = fileName;
}
	
void LibInfo::SetMean( int mean )
{
	m_mean = mean;
}

void LibInfo::SetStd( int std )
{
	m_std = std;
}
	
void LibInfo::SetThreshold( int threshold )
{
	m_threshold = threshold;
}

void LibInfo::SetOri( int ori )
{
	m_ori = ori;
}

void LibInfo::SetMapType( int mapType )
{
	m_mapType = mapType;
}
