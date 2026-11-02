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
#include "scontig.h"
#include "cluster.h"
#include "cluster_graph.h"
#include "probability_distribution.h"

/* sigma main program (the class itself is in sigma.cpp) */

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
			
			cmd = Sigma::SAMDIR + "/samtools view " + Sigma::mapping_files[sample_index] + " > " + fake_sam + "&";
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


