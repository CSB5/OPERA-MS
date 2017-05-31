#pragma once

#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "CommonFunction.h"

using namespace std;

struct Gap{
	int m_realGapSize;
	int m_estimatedGapSize;
	double m_newMu;
};

class GapCorrecter
{
public:
	GapCorrecter(void);
	GapCorrecter( int max, int min, int contigLength, int mean, int std, int readLength, vector<double> *pos );
	~GapCorrecter(void);

	// Attributes
private:
	vector<Gap*> *m_gapList;
	int m_contigLength;                   // the summation of the length of two contigs
	int m_maxGapSize;                     // the maximum gap size
	int m_mean;                           // the mean of the library
	int m_std;                            // the standard deviation of the library
	vector<double> *m_possibility;        // the possibility of each distance in the library
	int m_minGapSize;                     // the minimun gap size
	double m_minEstimatedGapSize;         // the minimum estimated gap size
	double m_maxEstimatedGapSize;         // the maximum estimated gap size
	int m_readLength;                     // the read length of this library

	bool m_useNormalDistribution;          // label if use normal distribution instead of empirical distribution


	// Methods
private:
	void Init();
	double CalculatePossibility( int value );   // calculate the possibility of a certain value

public:
	void CalculateTable();                // construct the relationship between real gap size and estimated gap size
	int PrintTable( string fileName );                     // print the constructed table
	int GetRealGapSize( double estimatedGapSize );         // get the corresponding real gap size of a given estimated gap size
};
