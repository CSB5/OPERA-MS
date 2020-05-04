#ifndef __EMANAGER_H__
#define __EMANAGER_H__

#include "global_inc.h"
#include "fastaReader.h"
#include "EucDist.h"
#include "kmerMap.h"
#include "Profiler.h"
#include "logger.h"
#include "NormalDistribution.h"
#include <time.h>
#include <thread>
#include <mutex>
/*
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/poisson.hpp>
using boost::math::normal;
using boost::math::poisson;
*/

class EManager
{
	public:
		EManager(char *input_fs, char *input_abundance, char *input_output);
		~EManager();
		void setAbundNum(int num);
		void addAbund(char *input_abundfile);
		int run(char *seedfile);
		void setVerbose(bool is_verbose);
		void setMaxEM(int maxem);
		void setThreadNum(int num);
		void setMinLength(unsigned int min_length);
		void setProbThreshold(double prob_threshold);

	private:
		// Variables
		char *fastafile;
		char *outputfile;
		int kmer_len, kmer_len2;
		kmerMap *kmap, *kmap2;
		fastaReader *seq;
		int ab_num;
		int ab_curr;
		//AbundanceLoader **ab_loader;
		char **ab_name;
		int max_EM;
		/*
		normal *intranormaldistr, *internormaldistr;
		normal *intranormaldistr2, *internormaldistr2;
		poisson ****poissondistr;
		*/
		NormalDistribution *intranormaldistr, *internormaldistr;
		NormalDistribution *intranormaldistr2, *internormaldistr2;
		long double *dist_prob, **abund_prob;
#ifdef USE_TWO_DIST
		long double *dist_prob2;
#endif
		logger *mylog;
		time_t start_t, end_t;

		unsigned int min_seq_length;
		int seqnum;
		int seed_num;
		Profiler **seed_profile;
		Profiler **seq_profile;
		Profiler **seed_profile2;
		Profiler **seq_profile2;
		bool *is_profile_N;
		bool *is_profile2_N;
		long double **seq_abundance;
		long double **seed_abundance;
		char **seed_header;
		long double **seq_prob;
		int *seq_bin;
		int *seed_count;
		char **bin_name;
		bool *is_estimated;

		int threadnum;
		std::thread **threadarr;
		std::mutex thread_mutex;

		char str[1024];
		long double MIN_PROB_THRESHOLD;
		int MAX_DIST_ESTIMATE;
		int FASTA_LINE;
		long double VERY_SMALL_DOUBLE;
		int STABLE_BIN_COUNT;

		// Functions
		void init();
		void estimate_normaldistr();
		void init_EM();
		void run_EM(int run_time);
		bool classify(long double min_prob, unsigned int min_seqlen);
		void filter_seed();
		void write_result();
		void write_fasta(char *seq, fstream *fs);
		//long double get_probability(long double distance, long double curr_abund, poisson *poisson_distr);
		//long double get_probability(long double distance1, long double distance2, long double curr_abund, poisson *poisson_distr);
		long double get_prob_dist(long double distance);
		long double get_prob_dist2(long double distance);
		//long double get_prob_abund(long double curr_abund, poisson *poisson_distr);
		long double get_prob_abund(long double curr_abund, long double lambda);
		void logtime_start();
		void logtime_end();
		void threadfunc_E(double abund, int k, int t_count);
		void threadfunc_M(int i);
};


#endif
