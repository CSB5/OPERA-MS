#include "AbstractDist.h"

AbstractDist::AbstractDist(int input_kmerlen)
{
	kmerlen = input_kmerlen;
	init();
}

AbstractDist::~AbstractDist()
{
	delete(kmap);
	if (rank1 != NULL)
	{
		free(rank1);
	}
	if (rank2 != NULL)
	{
		free(rank2);
	}
	if (tmprank != NULL)
	{
		free(tmprank);
	}
}

void AbstractDist::init()
{
	kmap = new kmerMap(kmerlen, true);
	entrynum = kmap->getEntryNum();
	rank1 = (int*)malloc(sizeof(int) * entrynum);
	rank2 = (int*)malloc(sizeof(int) * entrynum);
	tmprank = (int*)malloc(sizeof(int) * entrynum);
	normalize = false;
}

void AbstractDist::setNormalization(bool input)
{
	normalize = input;
}

void AbstractDist::resetDist()
{
	if (rank1 != NULL)
	{
		memset(rank1, '\0', sizeof(int) * entrynum);
	}
	if (rank2 != NULL)
	{
		memset(rank2, '\0', sizeof(int) * entrynum);
	}
	percent_N = 0;
}

void AbstractDist::computeRank(double *profile, int *rank)
{
	int i;
	// Initialize rank
	for (i = 0; i < entrynum; i++)
	{
		tmprank[i] = i;
	}
	qs.input(profile, tmprank, entrynum);
	qs.sort_all();
	for (i = 0; i < entrynum; i++)
	{
		rank[(int)(tmprank[i])] = i;
	}
}

double AbstractDist::getNpercentage()
{
	return percent_N;
}
