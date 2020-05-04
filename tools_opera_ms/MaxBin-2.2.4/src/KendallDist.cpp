#include "KendallDist.h"

double KendallDist::getDist(const char *inputseq1, const char *inputseq2)
{
	int i, j, k, tot;
	double f;

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
	computeRank(profile2, rank2);

	// Calculate Kendall's Distance
	k = 0;
	tot = (entrynum * (entrynum - 1)) / 2;
	for (i = 0; i < entrynum; i++)
	{
		for (j = i + 1; j < entrynum; j++)
		{
			if (((rank1[i] > rank1[j]) && (rank2[i] > rank2[j])) || ((rank1[i] < rank1[j]) && (rank2[i] < rank2[j])))
			{
				k++;
			}
		}
	}
	f = (double)(tot - k) / (double)tot;

	delete(pro1);
	delete(pro2);
	return(f);
}
