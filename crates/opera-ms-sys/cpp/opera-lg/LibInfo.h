#pragma once

#include <string>
class Configure;
#include "Configure.h"
using namespace std;

class LibInfo
{
public:
	LibInfo(void);
	virtual ~LibInfo(void);

private:
	string m_fileName;
	int m_mean;
	int m_std;
	int m_threshold;
	int m_ori;
	int m_mapType;

public:
	string GetFileName();
	int GetMean();
	int GetStd();
	int GetThreshold();
	int GetOri();
	int GetMapType();
	void SetFileName( string fileName );
	void SetMean( int mean );
	void SetStd( int std );
	void SetThreshold( int threshold );
	void SetOri( int ori );
	void SetMapType( int mapType );
};

