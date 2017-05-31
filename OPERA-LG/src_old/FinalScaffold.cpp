#include "FinalScaffold.h"


FinalScaffold::FinalScaffold(void)
{
}

FinalScaffold::FinalScaffold( string name, double length, string scaffold )
{
	m_name = name;
	m_length = length;
	m_scaffold = scaffold;
}

FinalScaffold::FinalScaffold( string name, double length, string scaffold, double cov )
{
	m_name = name;
	m_length = length;
	m_scaffold = scaffold;
	m_cov = cov;
}

FinalScaffold::~FinalScaffold(void)
{
}

double FinalScaffold::GetLength() const
{
	return m_length;
}

string FinalScaffold::GetName() const
{
	return m_name;
}

bool FinalScaffold::IsScaffold()
{
	if( m_scaffold == "" )
		return false;
	else
		return true;
}

string FinalScaffold::GetScaffold() const
{
	return m_scaffold;
}

double FinalScaffold::GetCov() const
{
	return m_cov;
}


