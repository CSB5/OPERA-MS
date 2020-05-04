#include "EManager.h"
#include <map>

#define MAX_ABUND_FILE 1024 // This definition also appears in EManager.cpp

int main(int argc, char *argv[])
{
	int i;
	char *inputfasta = NULL, *seed_file = NULL, *out = NULL;
	char **abund_file = NULL;
	int threadnum = 1;
	int abund_count;
	int maxem = -1;
	bool free_out = false;
	bool isverbose = false;
	unsigned int min_length = 0;
	double prob_threshold = 0;

	abund_file = (char**)malloc(sizeof(char*) * MAX_ABUND_FILE);
	memset(abund_file, '\0', sizeof(char*) * MAX_ABUND_FILE);
	abund_count = 0;

	// Handle parameters
	for (i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-fasta") == 0)
		{
			i++;
			inputfasta = argv[i];
		}
		else if (strcmp(argv[i], "-out") == 0)
		{
			i++;
			out = argv[i];
		}
		else if (strcmp(argv[i], "-seed") == 0)
		{
			i++;
			seed_file = argv[i];
		}
		else if (strncmp(argv[i], "-abund", 6) == 0)
		{
			i++;
			abund_file[abund_count] = argv[i];
			abund_count++;
		}
		else if (strcmp(argv[i], "-max_run") == 0)
		{
			i++;
			maxem = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-min_contig_length") == 0)
		{
			i++;
			min_length = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-prob_threshold") == 0)
		{
			i++;
			prob_threshold = (double)atof(argv[i]);
		}
		else if (strcmp(argv[i], "-verbose") == 0)
		{
			isverbose = true;
		}
		else if (strcmp(argv[i], "-thread") == 0)
		{
			i++;
			threadnum = atoi(argv[i]);
		}
		else
		{
			printf("Unrecognized token: %s\n", argv[i]);
			printf("Usage: MaxBin -fasta (input fasta) -seed (seed file) -abund (abundance file) [-abund2 (abundfile) -abund3 (abundfile) -abund4 ... ] [-out (output file)] [-max_run (max EM iteration; default 50)] [-thread (thread num; default 1)] [-prob_threshold (probability threshold)]\n");
			exit(-1);
		}
	}

	if (inputfasta == NULL || seed_file == NULL || abund_file == NULL)
	{
		printf("Usage: MaxBin -fasta (input fasta) -seed (seed file) -abund (abundance file) [-abund2 (abundfile) -abund3 (abundfile) -abund4 ... ] [-out (output file)] [-max_run (max EM iteration; default 50)] [-thread (thread num; default 1)] [-prob_threshold (probability threshold)]\n");
		exit(-1);
	}
	if (out == NULL)
	{
		out = (char*)malloc(sizeof(char) * (strlen(inputfasta) + 5));
		memset(out, '\0', sizeof(char) * (strlen(inputfasta) + 5));
		sprintf(out, "%s.out", inputfasta);
		free_out = true;
	}
	if (abund_count == 0)
	{
		printf("Please input at least one abund file.\n");
		printf("Usage: MaxBin -fasta (input fasta) -seed (seed file) -abund (abundance file) [-abund2 (abundfile) -abund3 (abundfile) -abund4 ... ] [-out (output file)] [-max_run (max EM iteration; default 50)] [-thread (thread num; default 1)] [-prob_threshold (probability threshold)]\n");
		exit(-1);
	}

	EManager em(inputfasta, NULL, out);
	em.setAbundNum(abund_count);
	for (i = 0; i < abund_count; i++)
	{
		em.addAbund(abund_file[i]);
	}
	em.setVerbose(isverbose);
	if (threadnum > abund_count)
	{
		threadnum = abund_count;
	}
	em.setThreadNum(threadnum);
	if (maxem != -1)
	{
		em.setMaxEM(maxem);
	}
	if (min_length > 0)
	{
		em.setMinLength(min_length);
	}
	if (prob_threshold > 0)
	{
		em.setProbThreshold(prob_threshold);
	}
	i = em.run(seed_file);
	if (i == -1)
	{
		printf("Please input at least one abundance files. Program exit.\n");
	}

	free(abund_file);
	if (free_out == true)
	{
		free(out);
	}
}

