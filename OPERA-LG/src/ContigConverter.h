#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <list>
#include "Configure.h"
#include "Contig.h"
#include "Graph.h"
#include "CommonFunction.h"
#include <ext/hash_map>
#include <stdexcept> 


using namespace std;
using namespace __gnu_cxx; 


class ContigConverter
{
public:
	ContigConverter(void);
	~ContigConverter(void);

	// Methods
public:
	// Convert contig file
	int ConvertContigFile( string fileName, Graph *graph, list<PetLibrary*> *libs );

private:
	// analyze velvet file
	int AnalyzeVelvet( ifstream *contigReader, vector<Contig*> *contigs );
	// analyze soapdenovo file
	int AnalyzeSoapDenovo( ifstream *contigReader, vector<Contig*> *contigs );
	// analyze opera scaffold file
	int AnalyzeOperaScaf( ifstream *contigReader, vector<Contig*> *contigs );
	// analyze soapdenovo file
	int AnalyzeSoapDenovoContig( ifstream *contigReader, vector<Contig*> *contigs );
	// analyze normal fasta file
	int AnalyzeNormalFasta( ifstream *contigReader, vector<Contig*> *contigs );
	// filter repeat
	void FilterRepeat( vector<Contig*> *contigs, Graph *graph );
	// analyze opera's contig file
	int AnalyzeStatistic( ifstream *contigReader, vector<Contig*> *contigs );
	// remove small contigs only
	void RemoveSmallContigs( vector<Contig*> *contigs, Graph *graph );
	// print list of contigs
	int PrintContigs( list<Contig*> *contigs, string fileName );
	inline void DeleteContigs( list<Contig*> *contigs );

	// calculate coverage using mapping information
	int CalculateCovUsingMapping( list<PetLibrary*> *libs );

	// check if two lines represent a pair of reads
	bool IsPair( string firstRead, string secondRead );
	// calculate the distance of two reads in the same contig
	int CalculateReadDisOnSameContig( vector<string> *firstColumn, vector<string> *secondColumn );
	// calculate mean and std using Interquartile Range
	int IQR( vector<double> *dis, double &mean, double &std, string libName, double &percentageOfOutliers, double &numberOfUsedPairs, 
		 PetLibrary *currentLibrary );

	// Attributes
private:
	// read contig file into a vector
	vector<Contig*> *myContigs;
	list<Contig*> *m_repeatContigs;
	list<Contig*> *m_smallContigs;
	hash_map<const char*, int, hash<const char*>, eqName> *m_contigNameHashMap;

	string m_libString;
};
