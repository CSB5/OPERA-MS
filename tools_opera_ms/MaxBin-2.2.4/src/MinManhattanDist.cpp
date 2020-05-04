#include "MinManhattanDist.h"
#include <math.h>

void MinManhattanDist::init()
{
	delete(kmap);
	kmap = new kmerMap(kmerlen, false);
	entrynum = kmap->getEntryNum();
	free(rank1);
	free(rank2);
	free(tmprank);
	rank1 = NULL;
	rank2 = NULL;
	tmprank = NULL;
}

double MinManhattanDist::getDist(const char *inputseq1, const char *inputseq2)
{
	int i, j;
	double f;

	resetDist();

	// Seq1
	pro1 = new Profiler(kmerlen, inputseq1, kmap);
	profile1 = pro1->getProfile();
	percent_N = pro1->getPercentN();

	// Seq2
	pro2 = new Profiler(kmerlen, inputseq2, kmap);
	profile2 = pro2->getProfile();
	if (percent_N < pro2->getPercentN())
	{
		percent_N = pro2->getPercentN();
	}

	// Calculate Ratio/Manhattan Distance
	f = 0;
	for (i = 0; i < entrynum; i++)
	{
		j = kmap->getReverseMapping(i);
		if (i > j)
		{
			if (profile1[i] <= profile1[j])
			{
				profile1[j] = 0;
			}
			else
			{
				profile1[i] = profile1[j];
				profile1[j] = 0;
			}
			if (profile2[i] <= profile2[j])
			{
				profile2[j] = 0;
			}
			else
			{
				profile2[i] = profile2[j];
				profile2[j] = 0;
			}
		}
		else if (i == j)
		{
			profile1[i] = 0;
			profile1[j] = 0;
			profile2[i] = 0;
			profile2[j] = 0;
		}
	}
	f = 0;
	for (i = 0; i < entrynum; i++)
	{
		if (i > j)
		{
			if (profile1[i] > profile2[i])
			{
				f = f + (double)(profile1[i] - profile2[i]);
			}
			else
			{
				f = f + (double)(profile2[i] - profile1[i]);
			}
		}
	}
	if (f >= 1)
	{
		printf("%s\n%s\n%f\n", inputseq1, inputseq2, f);
		exit(-1);
	}

	delete(pro1);
	delete(pro2);
	return(f);
}

double MinManhattanDist::getDist(double *input_pro1, double *input_pro2)
{
	int i, j;
	double f, *p1, *p2;

	p1 = (double*)malloc(sizeof(double) * entrynum);
	memcpy(p1, input_pro1, sizeof(double) * entrynum);
	p2 = (double*)malloc(sizeof(double) * entrynum);
	memcpy(p2, input_pro2, sizeof(double) * entrynum);

	resetDist();

	// Calculate Ratio/Manhattan Distance
	f = 0;
	for (i = 0; i < entrynum; i++)
	{
		j = kmap->getReverseMapping(i);
		if (i > j)
		{
			if (p1[i] <= p1[j])
			{
				p1[j] = 0;
			}
			else
			{
				p1[i] = p1[j];
				p1[j] = 0;
			}
			if (p2[i] <= p2[j])
			{
				p2[j] = 0;
			}
			else
			{
				p2[i] = p2[j];
				p2[j] = 0;
			}
		}
		else if (i == j)
		{
			p1[i] = 0;
			p1[j] = 0;
			p2[i] = 0;
			p2[j] = 0;
		}
	}
	f = 0;
	for (i = 0; i < entrynum; i++)
	{
		if (i > j)
		{
			if (p1[i] > p2[i])
			{
				f = f + (double)(p1[i] - p2[i]);
			}
			else
			{
				f = f + (double)(p2[i] - p1[i]);
			}
		}
	}

	free(p1);
	free(p2);
	return(f);
}

