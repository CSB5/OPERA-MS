#include <cstdlib>
#include <cstdio>
#include <ctime>

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>


#include "sigma.h"

#include "contig_reader.h"
#include "mapping_reader.h"
#include "edge_reader.h"
#include "contig.h"
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
		fprintf(stderr, "Error opening file: %s\n", config_file);
		exit(EXIT_FAILURE);
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

int main(int argc, char** argv) {
	time_t start, finish;

	if (argc > 4) {
		fprintf(stderr, "Usage: ./sigma config_file [opt tree file]\n");
		exit(EXIT_FAILURE);
	}

	Sigma::readConfigFile(argv[1]);

	ContigMap contigs;

	if (Sigma::num_samples == 0) {
		fprintf(stderr, "Loading contig information from %s...\n", Sigma::sigma_contigs_file.c_str());
		time(&start);
		ContigIO::load_contigs(Sigma::sigma_contigs_file.c_str(), &contigs);
		time(&finish);
		fprintf(stderr, "DONE! %.4f sec\n\n", difftime(finish, start));


		/*Get total assembly size*/

		/* Check if it is already computed */
 	 	struct stat buffer;
 	 	std::string assembly_size_name = Sigma::output_dir + "/assembly_size.dat";   
  		if (stat (assembly_size_name.c_str(), &buffer) == 0){   //file exists
  			FILE* assembly_size_file = fopen(assembly_size_name.c_str(), "r");
			if (assembly_size_file != NULL) {
				if (fscanf(assembly_size_file, "%ld\n", &Sigma::total_assembly_size) != 1){
					fprintf(stderr, "Error with reading assembly size from file: %s\n", assembly_size_name.c_str());
					exit(EXIT_FAILURE);
				}
				if (fscanf(assembly_size_file, "%ld\n", &Sigma::total_assembly_nb_contig) != 1){
					fprintf(stderr, "Error total contig number with reading assembly size from file: %s\n", assembly_size_name.c_str());
					exit(EXIT_FAILURE);
				}
				

			}
			fclose(assembly_size_file);

  		} else {	// file does not exist, calculate total assembly size and write it into file
			ContigReader* contig_reader;

			//if (Sigma::contigs_file_type == "SOAPdenovo") {
			//	contig_reader = new SOAPdenovoReader();
			//} else if (Sigma::contigs_file_type == "Velvet") {
			//	contig_reader = new VelvetReader();
			//} else if(Sigma::contigs_file_type == "RAY") {
			//	contig_reader = new RAYReader();
			//}else {
			//	fprintf(stderr, "Unknown contigs_file_type: %s\n", Sigma::contigs_file_type.c_str());
			//	exit(EXIT_FAILURE);
			//}

            contig_reader = new AllReader();
			//Sigma::total_assembly_size = 
			//The static varailbes Sigma::total_assembly_size and Sigma::total_assembly_nb_contig are initialize in this fuction

			contig_reader->get_assembly_size(Sigma::contigs_file.c_str());
			delete contig_reader;
		}

	} else {
		ContigReader* contig_reader;
        /** We used to care about what format the contigs used to be in.**/

		//if (Sigma::contigs_file_type == "SOAPdenovo") {
		//	contig_reader = new SOAPdenovoReader();
		//} else if (Sigma::contigs_file_type == "Velvet") {
		//	contig_reader = new VelvetReader();
		//} else if(Sigma::contigs_file_type == "RAY") {
		//	contig_reader = new RAYReader();
		//}else {
		//	fprintf(stderr, "Unknown contigs_file_type: %s\n", Sigma::contigs_file_type.c_str());
		//	exit(EXIT_FAILURE);
		//}

        contig_reader = new AllReader();

		fprintf(stderr, "Loading contigs from %s...\n", Sigma::contigs_file.c_str());
		time(&start);
		Sigma::total_assembly_size = contig_reader->read(Sigma::contigs_file.c_str(), &contigs);
		time(&finish);
		fprintf(stderr, "DONE! %.4f sec\n\n", difftime(finish, start));

		delete contig_reader;

		MappingReader* mapping_reader = new SAMReader();
		
		for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
			fprintf(stderr, "Loading mapping from %s...\n", Sigma::mapping_files[sample_index].c_str());
			time(&start);
			
			//need to read the mapping files multiple times no save from now
			std::string fake_sam = Sigma::mapping_files[sample_index] + ".sam";
			std::string cmd = "mkfifo " + fake_sam;
			fprintf(stderr, " *** EXEC %s\n", cmd.c_str());
			if( system (cmd.c_str()) );
			
			//run samtools to get the fake sam
			//cmd = "samtools view " + Sigma::mapping_files[sample_index] + " | head -n 100000 " + " > " + fake_sam + "&";
			cmd = "samtools view " + Sigma::mapping_files[sample_index] + " > " + fake_sam + "&";
			fprintf(stderr, " *** EXEC %s\n", cmd.c_str());
			if( system (cmd.c_str()) );

			//Read the fake sam file
			fprintf(stderr, "Loading mapping from %s...\n", fake_sam.c_str());

			mapping_reader->read(fake_sam.c_str(), sample_index, &contigs);
			//mapping_reader->read(Sigma::mapping_files[sample_index].c_str(), sample_index, &contigs);
			time(&finish);
			
			cmd = "rm " + fake_sam;
			fprintf(stderr, " *** EXEC %s\n", cmd.c_str());
			if( system (cmd.c_str()) );

			fprintf(stderr, "DONE! %.4f sec\n\n", difftime(finish, start));//std::cin.get();
		}

		delete mapping_reader;

		if (Sigma::sigma_contigs_file != "-") {
			fprintf(stderr, "Saving contig information to %s...\n", Sigma::sigma_contigs_file.c_str());
			time(&start);
			ContigIO::save_contigs(&contigs, Sigma::sigma_contigs_file.c_str());
			time(&finish);
			fprintf(stderr, "DONE! %.4f sec\n\n", difftime(finish, start));
		}
	}

	fprintf(stderr, "Number of contigs: %ld\n\n", contigs.size());
	fprintf(stderr, "Total number of contigs: %ld\n\n", Sigma::total_assembly_nb_contig);
	
	EdgeReader* edge_reader = new OperaBundleReader();

	EdgeSet edges_set;

	for (int bundle_index = 0; bundle_index < (int) Sigma::edges_files.size(); ++bundle_index) {
		fprintf(stderr, "Loading edges from %s...\n", Sigma::edges_files[bundle_index].c_str());
		fprintf(stderr, "Saving edges skipped after filtering to %s...\n", Sigma::skipped_filter_edges_files[bundle_index].c_str());
		time(&start);
		edge_reader->read(Sigma::edges_files[bundle_index].c_str(), &contigs, &edges_set, Sigma::skipped_filter_edges_files[bundle_index].c_str());
		time(&finish);
		fprintf(stderr, "DONE! %.4f sec\n\n", difftime(finish, start));
	}

	EdgeQueue edges;

	for (auto it = edges_set.begin(); it != edges_set.end(); ++it) {
		Edge edge = *it;

		edge.computeDistance();

		edges.push(edge);
	}

	edges_set.clear();

	fprintf(stderr, "Number of edges: %ld\n\n", edges.size());

	
	ProbabilityDistribution* prob_dist;

	if (Sigma::pdist_type == "Poisson") {
		prob_dist = new PoissonDistribution();
	} else if (Sigma::pdist_type == "NegativeBinomial") {
		if (Sigma::R_VALUE <= 0.0) {
			Sigma::R_VALUE = compute_R(&contigs);
			//prob_dist = new NegativeBinomialDistribution(compute_R(&contigs));
		} 
		prob_dist = new NegativeBinomialDistribution(Sigma::R_VALUE);
		
	} else {
		fprintf(stderr, "Unknown pdist_type: %s\n", Sigma::pdist_type.c_str());
		exit(EXIT_FAILURE);
	}

	fprintf(stderr, "Generating cluster graph...\n");
	time(&start);
	
	ClusterGraph* graph;
	if(argc == 2){
		graph = new ClusterGraph(&contigs, &edges);
	}
	else{
		fprintf(stderr, " *** Use guide tree\n");
		graph = new ClusterGraph(&contigs, argv[2]);
	}
	time(&finish);
	fprintf(stderr, "DONE! %.4f sec\n\n", difftime(finish, start));
	fprintf(stderr, "Number of trees: %ld\n\n", graph->roots()->size());
	

	if(Sigma::COMPUTE_SCORE == 1){
		fprintf(stderr, "Computing scores...\n");
		time(&start);
		graph->computeScores(prob_dist);
		time(&finish);
		fprintf(stderr, "DONE! %.4f sec\n\n", difftime(finish, start));
	}
	
	delete prob_dist;
	
	if(Sigma::COMPUTE_SCORE == 1){
		fprintf(stderr, "Computing models...\n");
		time(&start);
		graph->computeModels();
		time(&finish);
		fprintf(stderr, "DONE! %.4f sec\n\n", difftime(finish, start));
	}
	
	for (int bundle_index = 0; bundle_index < (int) Sigma::edges_files.size(); ++bundle_index) {
		fprintf(stderr, "Saving filtered edges to %s...\n", Sigma::filtered_edges_files[bundle_index].c_str());
		fprintf(stderr, "Saving edges removed in clustering to %s...\n", Sigma::skipped_cluster_edges_files[bundle_index].c_str());
		time(&start);
		edge_reader->filter(Sigma::edges_files[bundle_index].c_str(), &contigs, Sigma::filtered_edges_files[bundle_index].c_str(), Sigma::skipped_cluster_edges_files[bundle_index].c_str());
		time(&finish);
		fprintf(stderr, "DONE! %.4f sec\n\n", difftime(finish, start));
	}

	fprintf(stderr, "Saving clusters to %s...\n", Sigma::clusters_file.c_str());
	time(&start);
	graph->saveClusters(Sigma::clusters_file.c_str());
	time(&finish);
	fprintf(stderr, "DONE! %.4f sec\n\n", difftime(finish, start));

	delete edge_reader;
	delete graph;

	return 0;
}


///mnt/pnsg10_projects/bertrandd/sigma/script/run_mummer_large_ref.pl /mnt/pnsg10_projects/bertrandd/sigma/HMP/DATA/HMREFG/all_seqs_bact.fa /mnt/pnsg10_projects/bertrandd/sigma/HMP/DATA/HMREFG/MUMMER_SPLIT/ OPERA_MS/output_dir/debug_output/backup/scaffoldSeq_500k.fasta OPERA_MS/output_dir/debug_output/backup/MMUMMER_SPLIT_500000 OPERA_MS/output_dir/debug_output/backup/MAPPING OPERA_MS/output_dir/debug_output/backup/MAPPING/HMP.coords 15

