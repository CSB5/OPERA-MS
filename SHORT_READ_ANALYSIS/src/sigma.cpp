#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>


#include "sigma.h"

#include "contig_reader.h"
#include "mapping_reader.h"
#include "edge_reader.h"
#include "cluster.h"
#include "cluster_graph.h"
#include "probability_distribution.h"

//std::string Sigma::contigs_file_type;
int Sigma::AR_TYPE;
int Sigma::R_ESTIMATION_TYPE;
int Sigma::COMPUTE_SCORE;
int Sigma::SPLITTING_PENALTY;
int Sigma::USE_WINDOW;

std::string Sigma::contigs_file;
std::vector<std::string> Sigma::mapping_files;
std::vector<std::string> Sigma::edges_files;

std::string Sigma::sigma_contigs_file;
std::string Sigma::SAMDIR;

std::string Sigma::output_dir;
std::vector<std::string> Sigma::skipped_filter_edges_files;
std::vector<std::string> Sigma::skipped_cluster_edges_files;
std::vector<std::string> Sigma::filtered_edges_files;
std::string Sigma::clusters_file;


int Sigma::num_samples;

int Sigma::contig_len_thr;
int Sigma::contig_edge_len;
int Sigma::contig_window_len;
int Sigma::kmer_size;
long int Sigma::total_assembly_size;
long int Sigma::total_assembly_nb_contig;

std::string Sigma::pdist_type;

double Sigma::R_VALUE;

void Sigma::readConfigFile(char* config_file) {
	ParamsMap params;

	std::ifstream config_fp(config_file, std::ifstream::in);

	if (config_fp.is_open()) {
		std::string line;

		while (!config_fp.eof()) {
			std::getline(config_fp, line);

			if (line.empty()) continue;

			std::size_t comment_pos = line.find_first_of("#");

			if (comment_pos != std::string::npos) {
				line.erase(comment_pos, line.size() - comment_pos);
			}

			std::size_t eq_pos = line.find_first_of("=");

			if (eq_pos == std::string::npos) continue;

			std::size_t key_s_pos = line.find_first_not_of(" \t\n");
			std::size_t key_e_pos = line.find_last_not_of(" \t\n", eq_pos - 1);

			std::size_t value_s_pos = line.find_first_not_of(" \t\n", eq_pos + 1);
			std::size_t value_e_pos = line.find_last_not_of(" \t\n");

			std::string key = line.substr(key_s_pos, key_e_pos - key_s_pos + 1);
			std::string value = line.substr(value_s_pos, value_e_pos - value_s_pos + 1);

			params.insert(std::make_pair(key, value));
		}

		config_fp.close();
	} else {
		if(strcmp(config_file, "help") == 0){
			fprintf(stderr,"Usage:\n  bin/sigma <config-file>\n");
			exit(0);
		}
		else{
			fprintf(stderr, "Error opening file: %s\n", config_file);
			exit(EXIT_FAILURE);
		}
	}

	configure(&params);
}

void Sigma::configure(ParamsMap* params) {
	//contigs_file_type = getStringValue(params, std::string("contigs_file_type"));

	contigs_file = getStringValue(params, std::string("contigs_file"));
	mapping_files = getVectorValue(params, std::string("mapping_files"));
	edges_files = getVectorValue(params, std::string("edges_files"));
	kmer_size = getIntValue(params, std::string("kmer_size"));

	
	
	sigma_contigs_file = getStringValue(params, std::string("sigma_contigs_file"));

	output_dir = getStringValue(params, std::string("output_dir"));


	SAMDIR = getStringValue(params, std::string("samtools_dir"));
	
	for (auto it = edges_files.begin(); it != edges_files.end(); ++it) {
		std::string file_path = *it;

		std::size_t slash_pos = file_path.find_last_of('/');

		if (slash_pos == std::string::npos) {
			slash_pos = 0;
		} else {
			slash_pos++;
		}

		std::string file_name = file_path.substr(slash_pos, file_path.size() - slash_pos);

		skipped_filter_edges_files.push_back(output_dir + "/skipped_preprocess_" + file_name);
		skipped_cluster_edges_files.push_back(output_dir + "/skipped_hierarchical_" + file_name);
		filtered_edges_files.push_back(output_dir + "/filtered_" + file_name);
	}

	clusters_file = output_dir + "/clusters";

	num_samples = (int) mapping_files.size();

	contig_len_thr = getIntValue(params, std::string("contig_len_thr"));
	contig_edge_len = getIntValue(params, std::string("contig_edge_len"));
	contig_window_len = getIntValue(params, std::string("contig_window_len"));

	
	
	if (contig_len_thr == -1) contig_len_thr = 500;
	if (contig_edge_len == -1) contig_edge_len = 100;
	if (contig_window_len == -1) contig_window_len = 300;

	pdist_type = getStringValue(params, std::string("pdist_type"));

	if (pdist_type == "-") pdist_type = std::string("NegativeBinomial");

	
	
	AR_TYPE = getIntValue(params, std::string("AR_TYPE"));
	R_VALUE = getDoubleValue(params, std::string("R_VALUE"));
	COMPUTE_SCORE = getIntValue(params, std::string("COMPUTE_SCORE"));
	SPLITTING_PENALTY = getIntValue(params, std::string("SPLITTING_PENALTY"));
	USE_WINDOW = getIntValue(params, std::string("USE_WINDOW"));
	if(AR_TYPE == -1) AR_TYPE = 1;
	if(COMPUTE_SCORE == -1) COMPUTE_SCORE = 1;
	if(SPLITTING_PENALTY == -1) SPLITTING_PENALTY = 0;
	if(USE_WINDOW == -1) USE_WINDOW = 1;
	//R = getDoubleValue(params, std::string("R"));
}

int Sigma::getIntValue(ParamsMap* params, std::string key) {
	auto it = params->find(key);

	if (it != params->end())
		return atoi((*it).second.c_str());

	return -1;
}

double Sigma::getDoubleValue(ParamsMap* params, std::string key) {
	auto it = params->find(key);

	if (it != params->end())
		return atof((*it).second.c_str());

	return -1.0;
}

std::string Sigma::getStringValue(ParamsMap* params, std::string key) {
	auto it = params->find(key);

	if (it != params->end())
		return (*it).second;

	return std::string("-");
}

std::vector<std::string> Sigma::getVectorValue(ParamsMap* params, std::string key) {
	std::vector<std::string> vector;

	auto it = params->find(key);

	if (it != params->end()) {
		std::stringstream ss((*it).second);
		std::string element;

		while (std::getline(ss, element, ',')) {
			vector.push_back(element);
		}
	}

	return vector;
}



/* sigma constructor!!! */
Sigma::Sigma(char* config_file) {
	
	readConfigFile(config_file);
	init();

}

Sigma::~Sigma() {
	
	delete contigs;

}


void Sigma::init() {
	
	contigs = new ContigMap();

	ContigReader* contig_reader;
	contig_reader = new AllReader();

	fprintf(stderr, "Loading contigs from %s...\n", Sigma::contigs_file.c_str());
	Sigma::total_assembly_size = contig_reader->read(Sigma::contigs_file.c_str(), contigs);
	
	delete contig_reader;
}

//Only work for single sample
void Sigma::set_contig_coverage(const char* qname, int flag, const char* contig_id, int read_pos, int read_length){

	if(qname[strlen(qname) -1] != '2' ){
			
		auto it = contigs->find(contig_id);
		//Unmapped read or contig is too small
		if (it != contigs->end()){

			//Read in reverse complement
			if(flag == 16) {
				read_pos = read_pos + read_length - 1;
			}

			sContig* contig = (*it).second;

			--read_pos; // POS is 1-based

			if (read_pos >= contig->left_edge() && read_pos <= contig->right_edge()) {
				if (Sigma::contig_window_len > 0) {
					const int window_index = (read_pos - contig->left_edge()) / Sigma::contig_window_len;
					contig->read_counts()[0][window_index]++;
				} else {
					contig->read_counts()[0][0]++;
					
				}
				contig->sum_read_counts()[0]++;
			}
		}
	}
}


void Sigma::save_contigs(){
	ContigIO::save_contigs(contigs, Sigma::sigma_contigs_file.c_str());
}
