#include "ScaffoldResult.h"

ScaffoldResult::ScaffoldResult(void)
{
	m_id = -1;
}

ScaffoldResult::~ScaffoldResult(void)
{
}

int ScaffoldResult::GetID(){
	return m_id;
}

void ScaffoldResult::SetID( int id ){
	m_id = id;
}


double ScaffoldResult::GetLength(){
	return m_length;
}

string ScaffoldResult::GetScaffoldString(){
	return m_scaffold;
}

void ScaffoldResult::SetScaffoldString( string s ){
	m_scaffold = s;
}

void ScaffoldResult::SetLength( double l ){
	m_length = l;
}

void ScaffoldResult::SetCov( double cov ){
	m_cov = cov;
}
	
double ScaffoldResult::GetCov(){
	return m_cov;
}

void ScaffoldResult::SetLengthWithGap( double l ){
	m_lengthWithGap = l;
}

double ScaffoldResult::GetLengthWithGap(){
	return m_lengthWithGap;
}
