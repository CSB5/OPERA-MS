#include "EManager.h"
#include "SpearmanDist.h"
#include "ManhattanDist.h"
#include "EucDist.h"
#include "AbundanceLoader.h"
#include <fstream>
/*
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
*/
#include <math.h>

//#define USE_TWO_DIST
#define MAX_ABUND_FILE 1024 // This definition also appears in main.cpp
#define pi 3.141592653589793

EManager::EManager(char *input_fs, char *input_abundance, char *input_output)
{
	logtime_start();
	fastafile = input_fs;
	outputfile = input_output;
	init();
	seq = new fastaReader(input_fs);
	seqnum = seq->getNum();
	if (input_abundance != NULL)
	{
		setAbundNum(1);
		addAbund(input_abundance);
	}
}

EManager::~EManager()
{
	int i, j;
	for (i = 0; i < ab_num; i++)
	{
		//delete(ab_loader[i]);
		free(seq_abundance[i]);
		free(ab_name[i]);
		
	}
	//free(ab_loader);
	free(seq_abundance);
	if (seed_profile != NULL)
	{
		for (i = 0; i < seed_num; i++)
		{
			delete(seed_profile[i]);
#ifdef USE_TWO_DIST
			delete(seed_profile2[i]);
#endif
			//free(seed_header[i]);
			if (seed_count[i] > 0)
			{
				free(bin_name[i]);
			}
		}
		free(seed_profile);
#ifdef USE_TWO_DIST
		free(seed_profile2);
#endif
		free(seed_header);
		free(seed_count);
		free(bin_name);
		j = seq->getNum();
		for (i = 0; i < j; i++)
		{
			delete(seq_profile[i]);
#ifdef USE_TWO_DIST
			delete(seq_profile2[i]);
#endif
			free(seq_prob[i]);
		}
		free(seq_profile);
		free(is_profile_N);
#ifdef USE_TWO_DIST
		free(seq_profile2);
		free(is_profile2_N);
#endif
		free(seq_prob);
		for (i = 0; i < seed_num; i++)
		{
			free(seed_abundance[i]);
		}
		free(seed_abundance);
		free(seq_bin);
		free(is_estimated);
	}

	delete(kmap);
#ifdef USE_TWO_DIST
	delete(kmap2);
#endif
	if (intranormaldistr != NULL)
	{
		delete(intranormaldistr);
	}
	if (internormaldistr != NULL)
	{
		delete(internormaldistr);
	}
#ifdef USE_TWO_DIST
	if (intranormaldistr2 != NULL)
	{
		delete(intranormaldistr2);
	}
	if (internormaldistr2 != NULL)
	{
		delete(internormaldistr2);
	}
#endif
	delete(mylog);
	delete(seq);
}

void EManager::init()
{
	ab_num = 0;
	ab_curr = 0;
	//ab_loader = (AbundanceLoader**)malloc(sizeof(AbundanceLoader*) * MAX_ABUND_FILE);
	//memset(ab_loader, '\0', sizeof(AbundanceLoader*) * MAX_ABUND_FILE);
	ab_name = (char**)malloc(sizeof(char*) * MAX_ABUND_FILE);
	memset(ab_name, '\0', sizeof(char*) * MAX_ABUND_FILE);
	seed_profile = NULL;
	seq_profile = NULL;
	is_profile_N = NULL;
#ifdef USE_TWO_DIST
	seed_profile2 = NULL;
	seq_profile2 = NULL;
	is_profile2_N = NULL;
#endif
	seed_abundance = NULL;
	seq_abundance = NULL;
	seq_prob = NULL;
	seq_bin = NULL;
	seed_count = NULL;
	bin_name = NULL;
	is_estimated = NULL;
	intranormaldistr = NULL;
	internormaldistr = NULL;
	intranormaldistr2 = NULL;
	internormaldistr2 = NULL;
	dist_prob = NULL;
	abund_prob = NULL;
#ifdef USE_TWO_DIST
	dist_prob2 = NULL;
#endif
	//poissondistr = NULL;
	kmer_len = 4;
#ifdef USE_TWO_DIST
	kmer_len2 = 4;
#endif
	kmap = new kmerMap(kmer_len, true);
#ifdef USE_TWO_DIST
	kmap2 = new kmerMap(kmer_len2, true);
#endif
	min_seq_length = 1000;
	max_EM = 50;
	VERY_SMALL_DOUBLE = 1e-20;
	MIN_PROB_THRESHOLD = 0.5;
	MAX_DIST_ESTIMATE = 100000;
	FASTA_LINE = 70;
	STABLE_BIN_COUNT = 5;
	threadnum = 1;

	// Setup log file
	sprintf(str, "%s.log", outputfile);
	mylog = new logger(str);
	mylog->setVerbose(false);
}

void EManager::setVerbose(bool is_verbose)
{
	mylog->setVerbose(is_verbose);
}

void EManager::setThreadNum(int num)
{
	threadnum = num;
}

void EManager::setMaxEM(int maxem)
{
	max_EM = maxem;
}

void EManager::setMinLength(unsigned int min_length)
{
	min_seq_length = min_length;
	sprintf(str, "Minimum contig length set to %d.\n", min_seq_length);
	mylog->writelog(str, true);
}

void EManager::setProbThreshold(double prob_threshold)
{
	if (prob_threshold >= 0)
	{
		MIN_PROB_THRESHOLD = prob_threshold;
	}
	else
	{
		sprintf(str, "Cannot set probability threshold to %f. Stick to probability threshold %Lf.\n", prob_threshold, MIN_PROB_THRESHOLD);
		mylog->writelog(str, true);
	}
}


void EManager::setAbundNum(int num)
{
	if (num >= MAX_ABUND_FILE)
	{
		sprintf(str,"Abundance profile number exceeds maximum allowed number. Program stop.\n");
		mylog->writelog(str, true);
		exit(-1);
	}
	ab_num = num;
}

void EManager::addAbund(char *input_abundfile)
{
	AbundanceLoader temp_abund(input_abundfile);
	int i;
	/*
	if (ab_num >= MAX_ABUND_FILE)
	{
		sprintf(str, "Abundance profile number exceeds maximum allowed number. Program stop.\n");
		mylog->writelog(str, true);
		exit(-1);
	}
	*/
	if (temp_abund.is_parse_success() == false)
	{
		sprintf(str, "Failed to get abudance information from input file [%s].\nThe error lies in the following line:\n===\n", input_abundfile);
		mylog->writelog(str, true);
		mylog->writelog(temp_abund.get_parse_error(), true);
		sprintf(str, "\n===\nPlease check your input abundance files.\n");
		mylog->writelog(str, true);
		exit(-1);
	}

	if (seq_abundance == NULL)
	{
		seq_abundance = (long double**)malloc(sizeof(long double*) * ab_num);
		memset(seq_abundance, '\0', sizeof(long double*) * ab_num);
	}
	seq_abundance[ab_curr] = (long double*)malloc(sizeof(long double) * seqnum);
	memset(seq_abundance[ab_curr], '\0', sizeof(long double) * seqnum);
	for (i = 0; i < seqnum; i++)
	{
		seq_abundance[ab_curr][i] = temp_abund.getAbundance(seq->getHeaderByNum(i));
		if (seq_abundance[ab_curr][i] == -1)
		{
			sprintf(str, "Failed to get Abundance information for contig [%s] in file [%s]\n", seq->getHeaderByNum(i), input_abundfile);
			mylog->writelog(str, true);
			exit(-1);
		}
	}

	ab_name[ab_curr] = (char*)malloc(strlen(input_abundfile) + 10);
	memset(ab_name[ab_curr], '\0', strlen(input_abundfile) + 10);
	memcpy(ab_name[ab_curr], input_abundfile, strlen(input_abundfile));

	ab_curr++;

	/*
	ab_loader[ab_num] = new AbundanceLoader(input_abundfile);
	ab_name[ab_num] = (char*)malloc(strlen(input_abundfile) + 10);
	memset(ab_name[ab_num], '\0', strlen(input_abundfile) + 10);
	memcpy(ab_name[ab_num], input_abundfile, strlen(input_abundfile));
	ab_num++;
	*/
}

int EManager::run(char *seedfile)
{
	if (ab_num == 0)
	{
		return(-1);
	}

	// Open seed file and look for the sequences used as seeds. One sequence per line
	fstream *fs = new fstream(seedfile, ios::in);

	sprintf(str, "Reading seed list...\n");
	mylog->writelog(str, true);

	seq->resetMark();
	seed_num = 0;
	while (!fs->eof())
	{
		memset(str, '\0', 1024);
		fs->getline(str, 1023);
		while (str[strlen(str) - 1] == '\n' || str[strlen(str) - 1] == '\r')
		{
			str[strlen(str) - 1] = '\0';
		}
		if (str[0] == '\0')
		{
			continue;
		}
		seq->setMark(str);
		seed_num++;

		mylog->writelog("\t", false);
		mylog->writelog(str, false);
		mylog->writelog("\n", false);

		if (fs->eof())
		{
			break;
		}
		else if (fs->fail())
		{
			fs->clear();
		}
	}
	fs->close();
	delete(fs);

	if (seed_num > 1)
	{
		init_EM();
		//sprintf(str, "Test run...\n");
		//mylog->writelog(str, true);
		//run_EM(10);
		// Only classify sequences >= threshold. For eliminating bins without sequences
		//classify(MIN_PROB_THRESHOLD, min_seq_length);
		//filter_seed();
		run_EM(max_EM);
		classify(MIN_PROB_THRESHOLD, min_seq_length);
		write_result();
		logtime_end();
	}

	return(0);
}

void EManager::init_EM()
{
	int i, j, k;
	//long double d;
	// Initialize EM data structures

	sprintf(str, "Looking for seeds in sequences.\n");
	mylog->writelog(str, true);

	j = 0;
	for (i = 0; i < seqnum; i++)
	{
		if (seq->isMark(i) == true)
		{
			j++;
		}
	}
	if (j != seed_num)
	{
		sprintf(str, "Cannot find all seeds in contigs. Found %i seeds instead of %i in seed file.\nPlease check if the seed files.\n", j, seed_num);
		mylog->writelog(str, true);
		exit(-1);
	}
	seed_num = j;

	seq_profile = (Profiler**)malloc(sizeof(Profiler*) * seqnum);
	memset(seq_profile, '\0', sizeof(Profiler*) * seqnum);
	seed_profile = (Profiler**)malloc(sizeof(Profiler*) * seed_num);
	memset(seed_profile, '\0', sizeof(Profiler*) * seed_num);
	is_profile_N = (bool*)malloc(sizeof(bool) * seqnum);
	memset(is_profile_N, '\0', sizeof(bool) * seqnum);
#ifdef USE_TWO_DIST
	seq_profile2 = (Profiler**)malloc(sizeof(Profiler*) * seqnum);
	memset(seq_profile2, '\0', sizeof(Profiler*) * seqnum);
	seed_profile2 = (Profiler**)malloc(sizeof(Profiler*) * seed_num);
	memset(seed_profile2, '\0', sizeof(Profiler*) * seed_num);
	is_profile2_N = (bool*)malloc(sizeof(bool) * seqnum);
	memset(is_profile2_N, '\0', sizeof(bool) * seqnum);
#endif
	seed_abundance = (long double**)malloc(sizeof(long double*) * seed_num);
	memset(seed_abundance, '\0', sizeof(long double*) * seed_num);
	for (i = 0; i < seed_num; i++)
	{
		seed_abundance[i] = (long double*)malloc(sizeof(long double) * ab_num);
		memset(seed_abundance[i], '\0', sizeof(long double) * ab_num);
	}
	seed_header = (char**)malloc(sizeof(char*) * seed_num);
	memset(seed_header, '\0', sizeof(char*) * seed_num);
	seed_count = (int*)malloc(sizeof(int) * seed_num);
	memset(seed_count, '\0', sizeof(int) * seed_num);
	seq_bin = (int*)malloc(sizeof(int) * seqnum);
	memset(seq_bin, '\0', sizeof(int) * seqnum);
	is_estimated = (bool*)malloc(sizeof(bool) * seqnum);
	memset(is_estimated, '\0', sizeof(bool) * seqnum);

	j = 0;
	for (i = 0; i < seqnum; i++)
	{
		seq_profile[i] = new Profiler(kmer_len, seq->getSeqByNum(i), kmap);
		if (seq_profile[i]->getPercentN() == 1)
		{
			is_profile_N[i] = true;
		}
#ifdef USE_TWO_DIST
		seq_profile2[i] = new Profiler(kmer_len2, seq->getSeqByNum(i), kmap2);
		if (seq_profile2[i]->getPercentN() == 1)
		{
			is_profile2_N[i] = true;
		}
#endif
		//d = ab_loader->getAbundance(seq->getHeaderByNum(i));
		//if (seq->isMark(i) == true && d > 0)
		if (seq->isMark(i) == true)
		{
			seed_profile[j] = new Profiler(kmer_len, seq->getSeqByNum(i), kmap);
#ifdef USE_TWO_DIST
			seed_profile2[j] = new Profiler(kmer_len2, seq->getSeqByNum(i), kmap2);
#endif
			for (k = 0; k < ab_num; k++)
			{
				//seed_abundance[j][k] = ab_loader[k]->getAbundance(seq->getHeaderByNum(i));
				seed_abundance[j][k] = seq_abundance[k][i];
			}
			seed_header[j] = seq->getHeaderByNum(i);
			sprintf(str, "\t%s", seed_header[j]);
			mylog->writelog(str, true);
			for (k = 0; k < ab_num; k++)
			{
				sprintf(str, " [%Lf]", seed_abundance[j][k]);
				mylog->writelog(str, true);
			}
			sprintf(str, "\n");
			mylog->writelog(str, true);
			j++;
		}
		else
		{
			seq->isMark(i) == false;
		}
	}
	estimate_normaldistr();

	if (j != seed_num)
	{
		seed_num = j;
	}
	sprintf(str, "Get %d seeds.\n", seed_num);
	mylog->writelog(str, true);

	seq_prob = (long double**)malloc(sizeof(long double*) * seqnum);
	for (i = 0; i < seqnum; i++)
	{
		seq_prob[i] = (long double*)malloc(sizeof(long double) * seed_num);
		memset(seq_prob[i], '\0', sizeof(long double) * seed_num);
	}
}

void EManager::estimate_normaldistr()
{
	long double mean, std, mean2, std2;
#ifdef USE_TWO_DIST
	long double mean_2, std_2, mean2_2, std2_2;
#endif

	// Pre-defined mean and standard deviation found from simulation test on 3181 IMG bacterial and archaeal genomes

	// Euclidean distance
	// 4-mer
	mean = 0;
	std = 0.01037897 / 2;
	//mean = 0.01505219 / 2;
	//std = 0.01037897 / 2;
	//mean = 0;
	//std = 0.01037897;
	mean2 = 0.0676654;
	std2 = 0.03419337;
#ifdef USE_TWO_DIST
	// 5-mer
	mean_2 = 0.01228021 / 2;
	std_2 = 0.008654644 / 2;
	//mean_2 = 0;
	//std_2 = 0.008654644;
	mean2_2 = 0.04346744;
	std2_2 = 0.02088791;
	// 6-mer
	//mean = 0.01032181;
	//std = 0.007806874;
	//mean = 0.02672312;
	//std = 0.01257017;

	// Spearman footrule distance
	// 4-mer
	//mean = 0.05849619;
	//std = 0.03473218;
	//mean2 = 0.2523695;
	//std2 = 0.113352;
	mean_2 = 0.05849619 / 2;
	std_2 = 0.03473218 / 2;
	mean2_2 = 0.2523695;
	std2_2 = 0.113352;

	// Manhattan distance
	// 4-mer
	//mean = 0.06337266;
	//std = 0.03575501;
	//mean2 = 0.2886038;
	//std2 = 0.1432431;
#endif

	
	// intranormaldistr = new normal(mean, std);
	// internormaldistr = new normal(mean2, std2);
	intranormaldistr = new NormalDistribution(mean, std);
	internormaldistr = new NormalDistribution(mean2, std2);
#ifdef USE_TWO_DIST
	// intranormaldistr2 = new normal(mean_2, std_2);
	// internormaldistr2 = new normal(mean2_2, std2_2);
	intranormaldistr2 = new NormalDistribution(mean_2, std_2);
	internormaldistr2 = new NormalDistribution(mean2_2, std2_2);
#endif

	// Estimate mean and variance for the distances in the input dataset.
	/*
	int i, n1, n2;
	long double *dist;
	long double d, mean, std;
	EucDist edist(kmer_len);
	//SpearmanDist edist(kmer_len);
	edist.setNormalization(true);

	// Boost random number generator
	boost::mt19937 rng;
	boost::uniform_int<> drawer(0, seq->getNum() - 1);
	boost::variate_generator<boost::mt19937, boost::uniform_int<> > dice(rng, drawer);
	sprintf(str, "Estimating the parameters for distance probability function.\nSampling %d pairs of sequences...\n", MAX_DIST_ESTIMATE);
	mylog->writelog(str, true);

	dist = (long double*)malloc(sizeof(long double) * MAX_DIST_ESTIMATE);
	memset(dist, '\0', sizeof(long double) * MAX_DIST_ESTIMATE);

	mean = 0;
	for (i = 0; i < MAX_DIST_ESTIMATE; i++)
	{
		n1 = 0;
		n2 = 0;
		while (n1 == n2)
		{
			n1 = dice();
			n2 = dice();
		}
		d = edist.getDist(seq_profile[n1]->getProfile(), seq_profile[n2]->getProfile());
		dist[i] = d;
		mean = mean + d;
	}

	mean = mean / MAX_DIST_ESTIMATE; // Get mean value
	std = 0;
	for (i = 0; i < MAX_DIST_ESTIMATE; i++)
	{
		std = std + pow(mean - dist[i], 2);
	}
	std = sqrt(std / (MAX_DIST_ESTIMATE - 1));

	sprintf(str, "Mean = %Lf, Standard deviation = %Lf\n", mean, std);
	mylog->writelog(str, false);

	free(dist);

	// Create normal distribution
	intranormaldistr = new normal(mean, std);
	*/
}

void EManager::run_EM(int run_time)
{
	int run, i, j, k, p, diff_count, tempbin, stable_count;
	long double sum, d, max;
	int maxcount;
#ifdef USE_TWO_DIST
	long double d2;
#endif
	EucDist edist(kmer_len);
	edist.setNormalization(true);
#ifdef USE_TWO_DIST
	//EucDist edist2(kmer_len2);
	//SpearmanDist edist(kmer_len);
	SpearmanDist edist2(kmer_len2);
	//ManhattanDist edist(kmer_len);
	edist2.setNormalization(true);
#endif

	dist_prob = (long double*)malloc(sizeof(long double) * seed_num);
	memset(dist_prob, '\0', sizeof(long double) * seed_num);
#ifdef USE_TWO_DIST
	dist_prob2 = (long double*)malloc(sizeof(long double) * seed_num);
	memset(dist_prob2, '\0', sizeof(long double) * seed_num);
#endif
	abund_prob = (long double**)malloc(sizeof(long double*) * seed_num);
	memset(abund_prob, '\0', sizeof(long double*) * seed_num);
	for (i = 0; i < seed_num; i++)
	{
		abund_prob[i] = (long double*)malloc(sizeof(long double) * ab_num);
		memset(abund_prob[i], '\0', sizeof(long double) * ab_num);
	}

	sprintf(str, "\nStart EM process.\n");
	mylog->writelog(str, true);

	// Allocate memory for threads
	threadarr = (std::thread**)malloc(sizeof(thread*) * threadnum);
	memset(threadarr, '\0', sizeof(thread*) * threadnum);

	// Execute EM algorithm
	stable_count = 0;
	/*
	poissondistr = (poisson****)malloc(sizeof(poisson***) * seed_num);
	memset(poissondistr, '\0', sizeof(poisson***) * seed_num);
	for (i = 0; i < seed_num; i++)
	{
		poissondistr[i] = (poisson***)malloc(sizeof(poisson**) * ab_num);
		memset(poissondistr[i], '\0', sizeof(poisson**) * ab_num);
		for (j = 0; j < ab_num; j++)
		{
			poissondistr[i][j] = (poisson**)malloc(sizeof(poisson*) * threadnum);
			memset(poissondistr[i][j], '\0', sizeof(poisson*) * threadnum);
		}
	}
	*/
	for (run = 0; run < run_time; run++)
	{
		sprintf(str, "Iteration %d\n", run + 1);
		mylog->writelog(str, true);

		// Expectation
		// For each run, calculate the distance between the sequences and seeds and calculate probability
		for (i = 0; i < seed_num; i++)
		{
			//sprintf(str, "Abundance info [%s] [%d]: %Lf\n", seed_header[i], i + 1, seed_abundance[i]);
			sprintf(str, "Abundance info [%s] [%d]:", seed_header[i], i + 1);
			mylog->writelog(str, false);
			for (j = 0; j < ab_num; j++)
			{
				sprintf(str, " %Lf", seed_abundance[i][j]);
				mylog->writelog(str, false);
			}
			sprintf(str, "\n");
			mylog->writelog(str, false);
			/*
			for (j = 0; j < ab_num; j++)
			{
				for (k = 0; k < threadnum; k++)
				{
					poissondistr[i][j][k] = new poisson(seed_abundance[i][j]);
				}
			}
			*/
		}
		diff_count = 0;
		for (i = 0; i < seqnum; i++)
		{
			if (seq->getSeqLenByNum(i) >= min_seq_length && is_profile_N[i] == false)
			{
				// Calculate tetramer probability
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
					d = edist.getDist((double*)seq_profile[i]->getProfile(), (double*)seed_profile[j]->getProfile());
					dist_prob[j] = get_prob_dist(d);
					if (dist_prob[j] < VERY_SMALL_DOUBLE)
					{
						dist_prob[j] = VERY_SMALL_DOUBLE;
					}
					sum = sum + dist_prob[j];
				}
				for (j = 0; j < seed_num; j++)
				{
					dist_prob[j] = dist_prob[j] / sum;
				}
#ifdef USE_TWO_DIST
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
					d = edist2.getDist((double*)seq_profile2[i]->getProfile(), (double*)seed_profile2[j]->getProfile());
					dist_prob2[j] = get_prob_dist2(d);
					if (dist_prob2[j] < VERY_SMALL_DOUBLE)
					{
						dist_prob2[j] = VERY_SMALL_DOUBLE;
					}
					sum = sum + dist_prob2[j];
				}
				for (j = 0; j < seed_num; j++)
				{
					dist_prob2[j] = dist_prob2[j] / sum;
				}
#endif

				// Insert thread function here
				k = 0;
				while (k < ab_num)
				{
					maxcount = k + threadnum;
					if (maxcount > ab_num)
					{
						maxcount = ab_num;
					}
					for (p = k; p < maxcount; p++)
					{
						//threadarr[p - k] = new std::thread(&EManager::threadfunc_E, this, ab_loader[p]->getAbundance(seq->getHeaderByNum(i)), p, p - k);
						threadarr[p - k] = new std::thread(&EManager::threadfunc_E, this, seq_abundance[p][i], p, p - k);
					}
					for (p = k; p < maxcount; p++)
					{
						threadarr[p - k]->join();
						delete(threadarr[p - k]);
					}
					k = maxcount;
				}

/* Put into thread function threadfunc_E
				for (k = 0; k < ab_num; k++)
				{
					sum = 0;
					for (j = 0; j < seed_num; j++)
					{
						//abund_prob[j][k] = get_prob_abund(ab_loader[k]->getAbundance(seq->getHeaderByNum(i)), poissondistr[j][k][0]);
						abund_prob[j][k] = get_prob_abund(seq_abundance[k][i], poissondistr[j][k][0]);
						if (abund_prob[j][k] < VERY_SMALL_DOUBLE)
						{
							abund_prob[j][k] = VERY_SMALL_DOUBLE;
						}
						sum += abund_prob[j][k];
					}
					for (j = 0; j < seed_num; j++)
					{
						abund_prob[j][k] = abund_prob[j][k] / sum;
					}
				}
*/
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
#ifdef USE_TWO_DIST
					seq_prob[i][j] = dist_prob[j] * dist_prob2[j];
#else
					seq_prob[i][j] = dist_prob[j];
#endif
/*
if (seq_prob[i][j] != seq_prob[i][j])
{
	sprintf(str, "dist_prob[%d][%d] is nan.\n", i, j);
	mylog->writelog(str, false);
	exit(-1);
}
*/
					for (k = 0; k < ab_num; k++)
					{
						seq_prob[i][j] = seq_prob[i][j] * abund_prob[j][k];
/*
if (abund_prob[j][k] != abund_prob[j][k])
{
	sprintf(str, "%d: abund_prob[%d][%d] is nan.\n", i, j, k);
	mylog->writelog(str, false);
	exit(-1);
}
*/
					}
					sum = sum + seq_prob[i][j];
				}
				/*
				if (sum < VERY_SMALL_DOUBLE)
				{
					sum = 1;
				}
				*/
				max = 0;
				tempbin = -1;
				sprintf(str, "[%s]", seq->getHeaderByNum(i));
				mylog->writelog(str, false);
				for (j = 0; j < seed_num; j++)
				{
					seq_prob[i][j] = seq_prob[i][j] / sum;
					if (max < seq_prob[i][j])
					{
						max = seq_prob[i][j];
						tempbin = j;
					}
					sprintf(str, "\t(%d)%Lf", j + 1, seq_prob[i][j]);
					mylog->writelog(str, false);
				}
				if (seq_bin[i] != tempbin)
				{
					diff_count++;
					seq_bin[i] = tempbin;
				}
				mylog->writelog("\n", false);
				is_estimated[i] = true;
			}
		}

		if (diff_count == 0)
		{
			stable_count++;
		}
		else
		{
			stable_count = 0;
		}

		// Maximization
		i = 0;
		while (i < seed_num)
		{
			maxcount = i + threadnum;
			if (maxcount > seed_num)
			{
				maxcount = seed_num;
			}
			for (j = i; j < maxcount; j++)
			{
				threadarr[j - i] = new std::thread(&EManager::threadfunc_M, this, j);
			}
			for (j = i; j < maxcount; j++)
			{
				threadarr[j - i]->join();
				delete(threadarr[j - i]);
			}
			i = maxcount;
		}

/* Put into the thread function threadfunc_M
		for (i = 0; i < seed_num; i++)
		{
			// Thread function

			for (j = 0; j < ab_num; j++)
			{
				seed_abundance[i][j] = 0;
			}
			seed_profile[i]->reset();
#ifdef USE_TWO_DIST
			seed_profile2[i]->reset();
#endif
			d = 0;
			for (j = 0; j < seqnum; j++)
			{
				if (seq->getSeqLenByNum(j) >= min_seq_length && seq_prob[j][i] == seq_prob[j][i] && is_profile_N[j] == false)
				{
					// Update seed abundances
					for (k = 0; k < ab_num; k++)
					{
						//seed_abundance[i][k] = seed_abundance[i][k] + (ab_loader[k]->getAbundance(seq->getHeaderByNum(j)) * seq->getSeqLenByNum(j) * seq_prob[j][i]);
						seed_abundance[i][k] = seed_abundance[i][k] + (seq_abundance[k][j] * seq->getSeqLenByNum(j) * seq_prob[j][i]);
					}
					d = d + seq->getSeqLenByNum(j) * seq_prob[j][i];
					// Update seed profiles
					seed_profile[i]->addProfile(seq_profile[j], (long double)seq->getSeqLenByNum(j) * seq_prob[j][i]);
#ifdef USE_TWO_DIST
					seed_profile2[i]->addProfile(seq_profile2[j], (long double)seq->getSeqLenByNum(j) * seq_prob[j][i]);
#endif
				}
			}
			seed_profile[i]->calcProfile();
#ifdef USE_TWO_DIST
			seed_profile2[i]->calcProfile();
#endif
			for (k = 0 ; k < ab_num; k++)
			{
				seed_abundance[i][k] = seed_abundance[i][k] / d;
			}
		}
*/
		/*
		for (i = 0; i < seed_num; i++)
		{
			for (j = 0; j < ab_num; j++)
			{
				for (k = 0; k < threadnum; k++)
				{
					delete(poissondistr[i][j][k]);
				}
			}
		}
		*/

		if (stable_count >= STABLE_BIN_COUNT)
		{
			break;
		}
	}
	/*
	for (i = 0; i < seed_num; i++)
	{
		for (j = 0; j < ab_num; j++)
		{
			free(poissondistr[i][j]);
		}
		free(poissondistr[i]);
	}
	free(poissondistr);
	*/
	sprintf(str, "\nEM finishes successfully.\n");
	mylog->writelog(str, true);

	/* Write out dist information -- not used anymore
	   Also the get_probability function is outdated and needs to be re-written if this function is to revive.

	long double **seed_prob;
	sprintf(str, "%s.dist", outputfile);
	fstream *dist = new fstream(str, ios::out);
	poissondistr = (poisson***)malloc(sizeof(poisson**) * seed_num);
	memset(poissondistr, '\0', sizeof(poisson**) * seed_num);
	for (i = 0; i < seed_num; i++)
	{
		poissondistr[i] = (poisson**)malloc(sizeof(poisson*) * ab_num);
		memset(poissondistr[i], '\0', sizeof(poisson*) * ab_num);
	}
	seed_prob = (long double**)malloc(sizeof(long double*) * seed_num);
	memset(seed_prob, '\0', sizeof(long double*) * seed_num);
	for (i = 0; i < seed_num; i++)
	{
		for (j = 0; j < ab_num; j++)
		{
			poissondistr[i][j] = new poisson(seed_abundance[i][j]);
		}
		seed_prob[i] = (long double*)malloc(sizeof(long double) * seed_num);
		memset(seed_prob[i], '\0', sizeof(long double) * seed_num);
	}
	for (i = 0; i < seed_num; i++)
	{
		sprintf(str, "\tBin%3d", i + 1);
		dist->write(str, strlen(str));
	}
	dist->write("\n", 1);
	for (i = 0; i < seed_num; i++)
	{
		sum = 0;
		for (j = 0; j < seed_num; j++)
		{
			d = edist.getDist(seed_profile[i]->getProfile(), seed_profile[j]->getProfile());
#ifdef USE_TWO_DIST
			d2 = edist2.getDist(seed_profile2[i]->getProfile(), seed_profile2[j]->getProfile());
			seed_prob[i][j] = get_probability(d, d2, seed_abundance[j], poissondistr[i]);
#else
			seed_prob[i][j] = get_probability(d, seed_abundance[j], poissondistr[i]);
#endif
			sum = sum + seed_prob[i][j];
		}
		sprintf(str, "Bin%3d", i + 1);
		dist->write(str, strlen(str));
		for (j = 0; j < seed_num; j++)
		{
			seed_prob[i][j] = seed_prob[i][j] / sum;
			// Special handling of NAN (caused when sum is too small...)
			if (seed_prob[i][j] != seed_prob[i][j])
			{
				seed_prob[i][j] = 0;
			}
			sprintf(str, "\t%Lf", seed_prob[i][j]);
			dist->write(str, strlen(str));
		}
		dist->write("\n", 1);
	}
	dist->close();
	delete(dist);

	for (i = 0; i < seed_num; i++)
	{
		for (j = 0; j < ab_num; j++)
		{
			if (poissondistr[i][j] != NULL)
			{
				delete(poissondistr[i][j]);
			}
		}
		free(poissondistr[i]);
		free(seed_prob[i]);
	}
	free(seed_prob);
	free(poissondistr);
	*/

	free(dist_prob);
#ifdef USE_TWO_DIST
	free(dist_prob2);
#endif
	for (i = 0; i < seed_num; i++)
	{
		free(abund_prob[i]);
	}
	free(abund_prob);

	free(threadarr);
}

bool EManager::classify(long double min_prob, unsigned int min_seqlen)
{
	int i, j, k, p;
	long double max, sum, d;
	bool ret;

	EucDist edist(kmer_len);
	edist.setNormalization(true);
#ifdef USE_TWO_DIST
	//EucDist edist2(kmer_len2);
	//SpearmanDist edist(kmer_len);
	SpearmanDist edist2(kmer_len2);
	//ManhattanDist edist(kmer_len);
	edist2.setNormalization(true);
#endif

	// Allocate memory for threads
	threadarr = (std::thread**)malloc(sizeof(thread*) * threadnum);
	memset(threadarr, '\0', sizeof(thread*) * threadnum);

	dist_prob = (long double*)malloc(sizeof(long double) * seed_num);
	memset(dist_prob, '\0', sizeof(long double) * seed_num);
#ifdef USE_TWO_DIST
	dist_prob2 = (long double*)malloc(sizeof(long double) * seed_num);
	memset(dist_prob2, '\0', sizeof(long double) * seed_num);
#endif
	abund_prob = (long double**)malloc(sizeof(long double*) * seed_num);
	memset(abund_prob, '\0', sizeof(long double*) * seed_num);
	for (i = 0; i < seed_num; i++)
	{
		abund_prob[i] = (long double*)malloc(sizeof(long double) * ab_num);
		memset(abund_prob[i], '\0', sizeof(long double) * ab_num);
	}

	sprintf(str, "\nClassifying sequences based on the EM result.\nMinimum probability for binning: %0.2Lf\n", MIN_PROB_THRESHOLD);
	mylog->writelog(str, true);

	/*
	poissondistr = (poisson****)malloc(sizeof(poisson***) * seed_num);
	memset(poissondistr, '\0', sizeof(poisson***) * seed_num);
	for (i = 0; i < seed_num; i++)
	{
		poissondistr[i] = (poisson***)malloc(sizeof(poisson**) * ab_num);
		memset(poissondistr[i], '\0', sizeof(poisson**) * ab_num);
		for (j = 0; j < ab_num; j++)
		{
			poissondistr[i][j] = (poisson**)malloc(sizeof(poisson*) * threadnum);
			memset(poissondistr[i][j], '\0', sizeof(poisson*) * threadnum);
			for (k = 0; k < threadnum; k++)
			{
				poissondistr[i][j][k] = new poisson(seed_abundance[i][j]);
			}
		}
	}
	*/

	for (i = 0; i < seqnum; i++)
	{
		if (seq->getSeqLenByNum(i) >= min_seqlen && is_profile_N[i] == false)
		{
			if (is_estimated[i] == false)
			{
				// Test of separating abundance probability and tetramer probability
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
					d = edist.getDist((double*)seq_profile[i]->getProfile(), (double*)seed_profile[j]->getProfile());
					dist_prob[j] = get_prob_dist(d);
					sum = sum + dist_prob[j];
				}
				for (j = 0; j < seed_num; j++)
				{
					dist_prob[j] = dist_prob[j] / sum;
				}
#ifdef USE_TWO_DIST
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
					d = edist2.getDist((double*)seq_profile2[i]->getProfile(), (double*)seed_profile2[j]->getProfile());
					dist_prob2[j] = get_prob_dist2(d);
					sum = sum + dist_prob2[j];
				}
				for (j = 0; j < seed_num; j++)
				{
					dist_prob2[j] = dist_prob2[j] / sum;
				}
#endif

				// Insert thread function here
				// threadfunc(abund, k, threadcount);
				int maxcount;
				k = 0;
				while (k < ab_num)
				{
					maxcount = k + threadnum;
					if (maxcount > ab_num)
					{
						maxcount = ab_num;
					}
					for (p = k; p < maxcount; p++)
					{
						//threadarr[p - k] = new std::thread(&EManager::threadfunc_E, this, ab_loader[p]->getAbundance(seq->getHeaderByNum(i)), p, p - k);
						threadarr[p - k] = new std::thread(&EManager::threadfunc_E, this, seq_abundance[p][i], p, p - k);
					}
					for (p = k; p < maxcount; p++)
					{
						threadarr[p - k]->join();
						delete(threadarr[p - k]);
					}
					k = maxcount;
				}
/*
				for (k = 0; k < ab_num; k++)
				{
					sum = 0;
					for (j = 0; j < seed_num; j++)
					{
						abund_prob[j][k] = get_prob_abund(seq_abundance[k][i], poissondistr[j][k][0]);
						sum += abund_prob[j][k];
					}
					for (j = 0; j < seed_num; j++)
					{
						abund_prob[j][k] = abund_prob[j][k] / sum;
					}
				}
*/
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
#ifdef USE_TWO_DIST
					seq_prob[i][j] = dist_prob[j] * dist_prob2[j];
#else
					seq_prob[i][j] = dist_prob[j];
#endif
					for (k = 0; k < ab_num; k++)
					{
						seq_prob[i][j] = seq_prob[i][j] * abund_prob[j][k];
					}
					sum = sum + seq_prob[i][j];
				}
				sprintf(str, "[%s]", seq->getHeaderByNum(i));
				mylog->writelog(str, false);
				for (j = 0; j < seed_num; j++)
				{
					seq_prob[i][j] = seq_prob[i][j] / sum;
					// Special handling of NAN (caused when sum is too small...)
					sprintf(str, "\t(%d)%Lf", j + 1, seq_prob[i][j]);
					mylog->writelog(str, false);
					if (seq_prob[i][j] != seq_prob[i][j])
					{
						seq_prob[i][j] = 0;
					}
				}
				mylog->writelog("\n", false);
			}

			max = 0;
			for (j = 0; j < seed_num; j++)
			{
				if (seq_prob[i][j] > max)
				{
					max = seq_prob[i][j];
					seq_bin[i] = j;
				}
			}
			if (max <= min_prob)
			{
				seq_bin[i] = -1;
			}
			else
			{
				seed_count[seq_bin[i]]++;
			}
			if (seq_bin[i] != -1)
			{
				sprintf(str, "seq [%s] assigned to bin [%d]\n", seq->getHeaderByNum(i), seq_bin[i]);
				mylog->writelog(str, false);
			}

			/* The writing of prob_dist file is disabled. MaxBin will no longer produce this file.
			   If the writing of prob_dist file is to be resumed, get_prob_abund, poissondistr, and ab_loader need to be re-designed.

				sprintf(str, "%s.prob_dist", outputfile);
				fstream *distf = new fstream(str, ios::out);
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
					d = edist.getDist(seq_profile[i]->getProfile(), seed_profile[j]->getProfile());
					seq_prob[i][j] = get_prob_dist(d);
					sum = sum + seq_prob[i][j];
				}
				sprintf(str, "[%s]", seq->getHeaderByNum(i));
				distf->write(str, strlen(str));
				for (j = 0; j < seed_num; j++)
				{
					seq_prob[i][j] = seq_prob[i][j] / sum;
					// Special handling of NAN (caused when sum is too small...)
					sprintf(str, "\t(%d)%Lf", j + 1, seq_prob[i][j]);
					distf->write(str, strlen(str));
					if (seq_prob[i][j] != seq_prob[i][j])
					{
						seq_prob[i][j] = 0;
					}
				}
				distf->write("\n", 1);
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
					seq_prob[i][j] = get_prob_abund(ab_loader->getAbundance(seq->getHeaderByNum(i)), poissondistr[j]);
					sum = sum + seq_prob[i][j];
				}
				for (j = 0; j < seed_num; j++)
				{
					seq_prob[i][j] = seq_prob[i][j] / sum;
					// Special handling of NAN (caused when sum is too small...)
					sprintf(str, "\t(%d)%Lf", j + 1, seq_prob[i][j]);
					distf->write(str, strlen(str));
					if (seq_prob[i][j] != seq_prob[i][j])
					{
						seq_prob[i][j] = 0;
					}
				}
				distf->write("\n", 1);

				distf->close();
				delete(distf);
			*/
		}
		else
		{
			seq_bin[i] = -1;
		}
	}
	/*
	for (i = 0; i < seed_num; i++)
	{
		for (j = 0; j < ab_num; j++)
		{
			for (k = 0; k < threadnum; k++)
			{
				if (poissondistr[i][j][k] != NULL)
				{
					delete(poissondistr[i][j][k]);
				}
			}
			free(poissondistr[i][j]);
		}
		free(poissondistr[i]);
	}
	free(poissondistr);
	*/

	free(dist_prob);
#ifdef USE_TWO_DIST
	free(dist_prob2);
#endif
	for (i = 0; i < seed_num; i++)
	{
		free(abund_prob[i]);
	}
	free(abund_prob);

	free(threadarr);

	ret = false;
	for (i = 0; i < seed_num; i++)
	{
		if (seed_count[i] == 0)
		{
			ret = true;
			break;
		}
	}
	return(ret);
}

void EManager::filter_seed()
{
	// seed_profile
	// seed_abundance
	// seed_header
	// seed_num
	int i, j, k;
	int temp_num;
	Profiler **temp_profile;
	long double **temp_abundance;
	char **temp_header;
	int *temp_count;

	// Get the number of bins that are not empty
	temp_num = 0;
	for (i = 0; i < seed_num; i++)
	{
		if (seed_count[i] > 0)
		{
			temp_num++;
		}
	}
	sprintf(str, "Filter out %d bins with no results during test run.\n", seed_num - temp_num);
	mylog->writelog(str, true);

	temp_profile = (Profiler**)malloc(sizeof(Profiler*) * temp_num);
	memset(temp_profile, '\0', sizeof(Profiler*) * temp_num);
	temp_abundance = (long double**)malloc(sizeof(long double*) * temp_num);
	memset(temp_abundance, '\0', sizeof(long double*) * temp_num);
	for (i = 0; i < seed_num; i++)
	{
		temp_abundance[i] = (long double*)malloc(sizeof(long double) * ab_num);
		memset(temp_abundance[i], '\0', sizeof(long double) * ab_num);
	}
	temp_header = (char**)malloc(sizeof(char*) * temp_num);
	memset(temp_header, '\0', sizeof(char*) * temp_num);
	temp_count = (int*)malloc(sizeof(int) * temp_num);
	memset(temp_count, '\0', sizeof(int) * temp_num);

	j = 0;
	for (i = 0; i < seed_num; i++)
	{
		if (seed_count[i] > 0)
		{
			temp_profile[j] = seed_profile[i];
			for (k = 0; k < ab_num; k++)
			{
				temp_abundance[j][k] = seed_abundance[i][k];
			}
			temp_header[j] = seed_header[i];
			temp_count[j] = seed_count[i];
			j++;
		}
		else
		{
			delete(seed_profile[i]);
		}
	}
	free(seed_profile);
	for (i = 0; i < seed_num; i++)
	{
		free(seed_abundance[i]);
	}
	free(seed_abundance);
	free(seed_header);
	free(seed_count);
	seed_profile = temp_profile;
	seed_abundance = temp_abundance;
	seed_header = temp_header;
	seed_num = temp_num;
	seed_count = temp_count;

	for (i = 0; i < seqnum; i++)
	{
		free(seq_prob[i]);
	}
	free(seq_prob);

	seq_prob = (long double**)malloc(sizeof(long double*) * seqnum);
	for (i = 0; i < seqnum; i++)
	{
		seq_prob[i] = (long double*)malloc(sizeof(long double) * seed_num);
		memset(seq_prob[i], '\0', sizeof(long double) * seed_num);
	}
}

void EManager::write_result()
{
	int i, j, k;
	fstream **fs, *un;

	// Name the bins -- ignore bins with no sequences.
	bin_name = (char**)malloc(sizeof(char*) * seed_num);
	memset(bin_name, '\0', sizeof(char*) * seed_num);
	j = 1;
	k = 0;
	for (i = 0; i < seed_num; i++)
	{
		if (seed_count[i] > 0)
		{
			bin_name[i] = (char*)malloc(sizeof(char) * (strlen(outputfile) + 20));
			memset(bin_name[i], '\0', sizeof(char) * (strlen(outputfile) + 20));
			//if (seed_num < 1000)
			//{
			//	sprintf(bin_name[i], "%s.%03d.fasta", outputfile, j);
			//}
			//else
			//{
				sprintf(bin_name[i], "%s.%04d.fasta", outputfile, j);
			//}
			j++;
		}
		else
		{
			k++;
		}
	}
	sprintf(str, "Ignoring %d bins without any sequences.\n", k);
	mylog->writelog(str, true);

	// Write summary report into report page
	sprintf(str, "%s.summary", outputfile);
	un = new fstream(str, ios::out);
	for(i = 0; i < seed_num; i++)
	{
		if (seed_count[i] > 0)
		{
			sprintf(str, "Bin [%s]", bin_name[i]);
			un->write(str, strlen(str));
			for (j = 0; j < ab_num; j++)
			{
				sprintf(str, "\t%0.4Lf", seed_abundance[i][j]);
				un->write(str, strlen(str));
			}
			un->write("\n", 1);
			for (j = 0; j < seqnum; j++)
			{
				if (seq_bin[j] == i)
				{
					sprintf(str, "\t%s", seq->getHeaderByNum(j));
					un->write(str, strlen(str));
					for (k = 0; k < ab_num; k++)
					{
						//sprintf(str, "\t%0.4f", ab_loader[k]->getAbundance(seq->getHeaderByNum(j)));
						sprintf(str, "\t%0.4Lf", seq_abundance[k][j]);
						un->write(str, strlen(str));
					}
					un->write("\n", 1);
				}
			}
			un->write("\n", 1);
		}
	}
	sprintf(str, "\nBins without any sequences:\n");
	un->write(str, strlen(str));
	j = 1;
	for(i = 0; i < seed_num; i++)
	{
		if (seed_count[i] == 0)
		{
			sprintf(str, "%d:", j);
			un->write(str, strlen(str));
			for (j = 0; j < ab_num; j++)
			{
				sprintf(str, " (%0.4Lf)", seed_abundance[i][j]);
				un->write(str, strlen(str));
			}
			un->write("\n", 1);
			j++;
		}
	}
	un->close();
	delete(un);

	fs = (fstream**)malloc(sizeof(fstream*) * seed_num);
	memset(fs, '\0', sizeof(fstream*) * seed_num);
	for (i = 0; i < seed_num; i++)
	{
		if (seed_count[i] > 0)
		{
			fs[i] = new fstream(bin_name[i], ios::out);
		}
	}
	sprintf(str, "%s.noclass", outputfile);
	un = new fstream(str, ios::out);

	j = 0;
	for (i = 0; i < seqnum; i++)
	{
		if (seq_bin[i] != -1)
		{
			sprintf(str, ">%s\n", seq->getHeaderByNum(i));
			fs[seq_bin[i]]->write(str, strlen(str));
			write_fasta(seq->getSeqByNum(i), fs[seq_bin[i]]);
		}
		else
		{
			j++;
			sprintf(str, ">%s\n", seq->getHeaderByNum(i));
			un->write(str, strlen(str));
			write_fasta(seq->getSeqByNum(i), un);
		}
	}

	for (i = 0; i < seed_num; i++)
	{
		if (seed_count[i] > 0)
		{
			fs[i]->close();
			delete(fs[i]);
		}
	}
	free(fs);
	un->close();
	delete(un);

	sprintf(str, "Number of unclassified sequences: %d (%2.2f%%)\n", j, ((float)j / (float)(seq->getNum())) * 100);
	mylog->writelog(str, true);
}

void EManager::write_fasta(char *seq, fstream *fs)
{
	int len, i;
	char *c = seq;
	len = strlen(seq);
	i = 0;
	while(i < len)
	{
		if (i + FASTA_LINE >= len)
		{
			i = len;
			fs->write(c, strlen(c));
			fs->write("\n", 1);
		}
		else
		{
			i = i + FASTA_LINE;
			fs->write(c, FASTA_LINE);
			fs->write("\n", 1);
			c = c + FASTA_LINE;
		}
	}
}

/*
long double EManager::get_probability(long double distance, long double curr_abund, poisson *poisson_distr)
{
	long double ret = get_prob_dist(distance) * get_prob_abund(curr_abund, poisson_distr);
	return ret;
}

long double EManager::get_probability(long double distance1, long double distance2, long double curr_abund, poisson *poisson_distr)
{
	long double ret = get_prob_dist(distance1) * get_prob_dist(distance2) * get_prob_abund(curr_abund, poisson_distr);
	return ret;
}
*/

long double EManager::get_prob_dist(long double distance)
{
	//long double d_intra = pdf(*intranormaldistr, distance);
	//long double d_inter = pdf(*internormaldistr, distance);
	long double d_intra = intranormaldistr->prob(distance);
	long double d_inter = internormaldistr->prob(distance);
	long double ret = d_intra / (d_inter + d_intra);
	//long double ret = pow(d_intra / (d_inter + d_intra), 2);
	return ret;
}

/* Using hardwired distance distribution -- not performing well
long double EManager::get_prob_dist(long double distance)
{
	// Initialize hardwired distribution
	long double dist_prob[650] = {0.0457,0.0381,0.0422,0.0402,0.0438,0.0391,0.0411,0.0372,0.0372,0.0331,0.0339,0.0334,0.0298,0.0290,0.0265,0.0262,0.0249,0.0218,0.0215,0.0207,0.0189,0.0171,0.0186,0.0149,0.0160,0.0142,0.0123,0.0113,0.0106,0.0106,0.0106,0.0092,0.0078,0.0092,0.0083,0.0065,0.0067,0.0064,0.0065,0.0062,0.0058,0.0051,0.0050,0.0043,0.0045,0.0042,0.0035,0.0033,0.0036,0.0029,0.0029,0.0029,0.0028,0.0024,0.0026,0.0024,0.0023,0.0023,0.0017,0.0020,0.0018,0.0016,0.0018,0.0016,0.0012,0.0013,0.0015,0.0015,0.0011,0.0009,0.0011,0.0008,0.0010,0.0009,0.0010,0.0010,0.0007,0.0007,0.0008,0.0008,0.0006,0.0007,0.0007,0.0006,0.0006,0.0005,0.0006,0.0005,0.0005,0.0006,0.0004,0.0004,0.0005,0.0003,0.0005,0.0004,0.0004,0.0002,0.0003,0.0002,0.0002,0.0003,0.0002,0.0001,0.0003,0.0002,0.0002,0.0002,0.0001,0.0002,0.0002,0.0002,0.0002,0.0002,0.0001,0.0002,0.0002,0.0002,0.0001,0.0002,0.0002,0.0002,0.0001,0,0.0001,0.0002,0.0002,0.0002,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0,0,0.0001,0.0001,0.0001,0.0001,0,0.0001,0,0,0,0,0,0,0.0001,0.0001,0.0001,0,0,0,0,0,0.0001,0,0,0,0,0,0.0001,0,0,0,0,0.0001,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0001,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0001,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0001,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0001,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	long double interval = 0.001;

	int i = distance / interval;
	if (i >= 650)
	{
		return 0;
	}
	else
	{
		return dist_prob[i], 2;
	}
}
*/

long double EManager::get_prob_dist2(long double distance)
{
	// long double d_intra = pdf(*intranormaldistr2, distance);
	// long double d_inter = pdf(*internormaldistr2, distance);
	long double d_intra = intranormaldistr2->prob(distance);
	long double d_inter = internormaldistr2->prob(distance);
	long double ret = d_intra / (d_inter + d_intra);
	return ret;
}

long double EManager::get_prob_abund(long double curr_abund, long double lambda)
{
	long double ret;
	long double l;
	if (lambda == 0 || curr_abund == 0)
	{
		ret = 0;
	}
	else
	{
		// This probability function is adapted from http://www.masaers.com/2013/10/08/Implementing-Poisson-pmf.html
		l = log(lambda);
		ret = exp((curr_abund * l) - lgamma(curr_abund + 1.0) - lambda);
	}
	return ret;
}

/*
long double EManager::get_prob_abund(long double curr_abund, poisson *poisson_distr)
{
	long double ret;
	if (poisson_distr != NULL && curr_abund > 0)
	{
		//ret = sqrt(pdf(*poisson_distr, curr_abund));
		ret = pdf(*poisson_distr, curr_abund);
	}
	else
	{
		ret = 0;
	}
	return ret;
}
*/

void EManager::logtime_start()
{
	start_t = time(NULL);
}

void EManager::logtime_end()
{
	end_t = time(NULL);
	mylog->writetime(end_t - start_t);
}

void EManager::threadfunc_E(double abund, int k, int t_count)
{
	long double sum;
	int j;
	sum = 0;
	for (j = 0; j < seed_num; j++)
	{
		//abund_prob[j][k] = get_prob_abund(abund, poissondistr[j][k][t_count]);
		abund_prob[j][k] = get_prob_abund(abund, seed_abundance[j][k]);
		if (abund_prob[j][k] < VERY_SMALL_DOUBLE)
		{
			abund_prob[j][k] = VERY_SMALL_DOUBLE;
		}
		sum += abund_prob[j][k];
	}
	for (j = 0; j < seed_num; j++)
	{
		abund_prob[j][k] = abund_prob[j][k] / sum;
	}
}


void EManager::threadfunc_M(int i)
{
	long double d;
	int j, k;
	unsigned int len;

	for (j = 0; j < ab_num; j++)
	{
		seed_abundance[i][j] = 0;
	}
	seed_profile[i]->reset();
#ifdef USE_TWO_DIST
	seed_profile2[i]->reset();
#endif
	d = 0;
	for (j = 0; j < seqnum; j++)
	{
		thread_mutex.lock();
		len = seq->getSeqLenByNum(j);
		thread_mutex.unlock();

		if (len >= min_seq_length && seq_prob[j][i] == seq_prob[j][i] && is_profile_N[j] == false)
		{
			// Update seed abundances
			for (k = 0; k < ab_num; k++)
			{
				seed_abundance[i][k] = seed_abundance[i][k] + (seq_abundance[k][j] * len * seq_prob[j][i]);
			}
			d = d + len * seq_prob[j][i];
			// Update seed profiles
			seed_profile[i]->addProfile(seq_profile[j], len * seq_prob[j][i]);
#ifdef USE_TWO_DIST
			seed_profile2[i]->addProfile(seq_profile2[j], len * seq_prob[j][i]);
#endif
		}
	}
	seed_profile[i]->calcProfile();
#ifdef USE_TWO_DIST
	seed_profile2[i]->calcProfile();
#endif
	for (k = 0 ; k < ab_num; k++)
	{
		seed_abundance[i][k] = seed_abundance[i][k] / d;
	}
}

