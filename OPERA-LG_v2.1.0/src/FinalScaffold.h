#pragma once

#include <string>

using namespace std;

class FinalScaffold
{
public:
	FinalScaffold(void);
	FinalScaffold( string name, double length, string scaffold );
	FinalScaffold( string name, double length, string scaffold, double cov );
	virtual ~FinalScaffold(void);

	// Attributes
private:
	string m_name;					// name of this scaffold
	double m_length;					// legnth of this scaffold, without gap sizes
	string m_scaffold;				// the scaffold contents if it is a scaffold
	double m_cov;					// the coverage of final scaffolds

	// Methods
public:
	double GetLength() const;
	string GetName() const;
	bool IsScaffold(); 
	string GetScaffold() const;
	double GetCov() const;
};


