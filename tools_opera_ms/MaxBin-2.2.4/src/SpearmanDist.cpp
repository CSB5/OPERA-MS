#include "SpearmanDist.h"

double SpearmanDist::getDist(const char *inputseq1, const char *inputseq2)
{
	int i;
	double j;

	resetDist();

	// Seq1
	pro1 = new Profiler(kmerlen, inputseq1, kmap);
	profile1 = pro1->getProfile();
	percent_N = pro1->getPercentN();
	computeRank(profile1, rank1);

	// Seq2
	pro2 = new Profiler(kmerlen, inputseq2, kmap);
	profile2 = pro2->getProfile();
	if (percent_N < pro2->getPercentN())
	{
		percent_N = pro2->getPercentN();
	}
	//printf("%s\n%s\nN\%\: %f\n", inputseq1, inputseq2, percent_N);
	computeRank(profile2, rank2);

	// The percentage of N is not taken into consideration so far.
	// Could be considered in the future...

	// Calculate Spearman Footrule Distance
	j = 0;
	for (i = 0; i < entrynum; i++)
	{
		if (rank1[i] >= rank2[i])
		{
			j = j + (rank1[i] - rank2[i]);
		}
		else
		{
			j = j + (rank2[i] - rank1[i]);
		}
	}

	if (normalize == true)
	{
		j = j / (entrynum * (entrynum + 1));
	}

	delete(pro1);
	delete(pro2);
	return(j);
}

double SpearmanDist::getDist(double *input_pro1, double *input_pro2)
{
	int i;
	double j, *p1, *p2;

	p1 = (double*)malloc(sizeof(double) * entrynum);
	memcpy(p1, input_pro1, sizeof(double) * entrynum);
	p2 = (double*)malloc(sizeof(double) * entrynum);
	memcpy(p2, input_pro2, sizeof(double) * entrynum);

	resetDist();

	// Seq1
	computeRank(p1, rank1);

	// Seq2
	computeRank(p2, rank2);

	// Calculate Spearman Footrule Distance
	j = 0;
	for (i = 0; i < entrynum; i++)
	{
		if (rank1[i] >= rank2[i])
		{
			j = j + (rank1[i] - rank2[i]);
		}
		else
		{
			j = j + (rank2[i] - rank1[i]);
		}
	}

	if (normalize == true)
	{
		j = j / (entrynum * (entrynum + 1));
	}

	free(p1);
	free(p2);
	return(j);
}

