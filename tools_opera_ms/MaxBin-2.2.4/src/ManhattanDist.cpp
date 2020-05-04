#include "ManhattanDist.h"

double ManhattanDist::getDist(const char *inputseq1, const char *inputseq2)
{
	int i;
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

	// Calculate Manhattan Distance
	f = 0;
	for (i = 0; i < entrynum; i++)
	{
		if (profile1[i] >= profile2[i])
		{
			f = f + (double)(profile1[i] - profile2[i]);
		}
		else
		{
			f = f + (double)(profile2[i] - profile1[i]);
		}
	}

	if (normalize == true)
	{
		f = f / 2;
	}

	delete(pro1);
	delete(pro2);
	return(f);
}

double ManhattanDist::getDist(double *input_pro1, double *input_pro2)
{
	int i;
	double f, *p1, *p2;

	p1 = (double*)malloc(sizeof(double) * entrynum);
	memcpy(p1, input_pro1, sizeof(double) * entrynum);
	p2 = (double*)malloc(sizeof(double) * entrynum);
	memcpy(p2, input_pro2, sizeof(double) * entrynum);

	resetDist();

	// Calculate Manhattan Distance
	f = 0;
	for (i = 0; i < entrynum; i++)
	{
		if (p1[i] >= p2[i])
		{
			f = f + (double)(p1[i] - p2[i]);
		}
		else
		{
			f = f + (double)(p2[i] - p1[i]);
		}
	}

	if (normalize == true)
	{
		f = f / 2;
	}

	free(p1);
	free(p2);
	return(f);
}
