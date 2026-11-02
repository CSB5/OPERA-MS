#pragma once

// record the scaffold results after dealing with a subgraph

#include <string>

using namespace std;

class ScaffoldResult
{
public:
	ScaffoldResult(void);
	~ScaffoldResult(void);

	// Methods
public:
	int GetID();
	void SetID( int id );
	double GetLength();
	string GetScaffoldString();
	void SetScaffoldString( string s );
	void SetLength( double l );
	void SetCov( double cov );
	double GetCov();
	void SetLengthWithGap( double l );
	double GetLengthWithGap();

	// Attributes
private:
	string m_scaffold;
	int m_id;			// the id in result array
	double m_length;		// the length of the new scaffold
	double m_cov;			// the average coverage of new scaffold
	double m_lengthWithGap;         // the length of the new scaffold with gap
};
