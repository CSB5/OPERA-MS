#include <cstdlib>
#include <cstdio>
#include <iostream>
#include "edge_reader.h"

#include "sigma.h"

EdgeReader::~EdgeReader() {}


OperaBundleReader::OperaBundleReader() {}

void OperaBundleReader::read(const char* edges_file, const ContigMap* contigs, EdgeSet* edges, const char* skipped_edges_file) {
	char id1[256], id2[256], line[1024], line_lib[255];
	int size, stdev, distance, edge_bundle_size_thr;
	double lib_mean = 3000, lib_stdev = 0;

	FILE* edges_fp = fopen(edges_file, "r");
	int arrival_rate_ratio_th = 0;
	
	if (edges_fp != NULL) {
		FILE* skipped_edges_fp = fopen(skipped_edges_file, "w");

		if (skipped_edges_fp == NULL) {
			fprintf(stderr, "Error opening file: %s\n", skipped_edges_file);
			exit(EXIT_FAILURE);
		}


		/* Obtain library cut off:
		*  if librabry mean + std <= 3k
		*	 edge_cut_off = 2
		*  else edge_cut_off = 5
		*/
		std::string edges_file_s(edges_file);
    	std::string edge_log(edges_file_s.substr(0,edges_file_s.find_last_of("/")+1) + "lib.txt");

    	FILE* edge_log_fp = fopen(edge_log.c_str(), "r");
    	if (edge_log_fp == NULL) {
    		fprintf(stderr, "Error opening file: %s\n", edge_log.c_str());
    		exit(EXIT_FAILURE);
    	}

    	while (fgets(line_lib, 255, edge_log_fp))
    	{
        	if (std::string(line_lib).find("Mean length") != std::string::npos) {
            	lib_mean = ::atof(std::string(line_lib).substr(31).c_str());
        	} else if (std::string(line_lib).find("Standard deviation") != std::string::npos) {
            	lib_stdev = ::atof(std::string(line_lib).substr(38).c_str());
        	}
    	}	
    	fclose(edge_log_fp);

    	if ((lib_mean + lib_stdev) <= 3000) {
    		edge_bundle_size_thr = 2;
    	} else {
    		//edge_bundle_size_thr = 5;
		edge_bundle_size_thr = 2;
    	}
	//fprintf(stderr, "\n\nLibrary %s has cutoff %d\n\n", edge_log.c_str(), edge_bundle_size_thr);

		while (!feof(edges_fp)) {
			if( fscanf(edges_fp, "%[^\n]\n", line) );

			// [ID1]\t[ORIENTATION1]\t[ID2]\t[ORIENTATION2]\t[DISTANCE]\t[STDEV]\t[SIZE]\n
			if (sscanf(line, "%s\t%*c\t%s\t%*c\t%d\t%d\t%d%*[^\n]", id1, id2, &distance, &stdev, &size) == 5) {
				auto it1 = contigs->find(id1);
				auto it2 = contigs->find(id2);

				/* Edge size has to be larger than edge_bundle_size_thr and the following needs to be fulfilled
				   if distance < 0 : 6*stdev + kmer_size + distance > 0 (opera condition) because
				   overlap has to be smaller than kmer size  */
				if (it1 != contigs->end() && it2 != contigs->end() && (size >= edge_bundle_size_thr) 
				       &&	((distance >= 0) || ((distance < 0) && (6 * stdev + Sigma::kmer_size + distance > 0)))) {
					Contig* contig1 = (*it1).second;
					Contig* contig2 = (*it2).second;
					//The cutoff
					//Quick ack need to be updated
					double ar1 = contig1->sum_read_counts()[0] / (double)contig1->length();
					double ar2 = contig2->sum_read_counts()[0] / (double)contig2->length();
					//
					if (contig1 != contig2 &&
					    (arrival_rate_ratio_th == 0 || !(ar1 > ar2 * arrival_rate_ratio_th || ar2 > ar1 * arrival_rate_ratio_th))
					    ) {
//                        fprintf(stderr, "EDGE READ SUCCESS : %s\n", line);
						edges->insert(Edge(contig1, contig2));
					}
					else{
						fprintf(skipped_edges_fp, "%s\t%d\n", line, 7);	
					}
				} else {
					int skip_flag = 0;
					if (!(it1 != contigs->end() && it2 != contigs->end())) { skip_flag += 1; }
					if (!(size >= edge_bundle_size_thr)) { skip_flag += 3; }
					if (!((distance >= 0) || ((distance < 0) && (6 * stdev + Sigma::kmer_size + distance > 0)))) { skip_flag += 5; }
					fprintf(skipped_edges_fp, "%s\t%d\n", line, skip_flag);
				}
			}

            else{
//                fprintf(stderr, "EDGE NOT READ SUCCESSFULLY: %s\n", line);
            }
		}

		fclose(skipped_edges_fp);

		fclose(edges_fp);
	} else {
		fprintf(stderr, "Error opening file: %s\n", edges_file);
		exit(EXIT_FAILURE);
	}
}


void OperaBundleReader::filter(const char* edges_file, const ContigMap* contigs, const char* filtered_edges_file, const char* skipped_edges_file) {
	char id1[256], id2[256], line[1024], line_lib[255];
	int size, stdev, distance, edge_bundle_size_thr;
	double lib_mean = 3000, lib_stdev = 0;

	FILE* edges_fp = fopen(edges_file, "r");

	if (edges_fp != NULL) {
		FILE* filtered_edges_fp = fopen(filtered_edges_file, "w");

		if (filtered_edges_fp == NULL) {
			fprintf(stderr, "Error opening file: %s\n", filtered_edges_file);
			exit(EXIT_FAILURE);
		}

		FILE* skipped_edges_fp = fopen(skipped_edges_file, "w");
		if (skipped_edges_fp == NULL) {
			fprintf(stderr, "Error opening file: %s\n", skipped_edges_file);
			exit(EXIT_FAILURE);
		}


		/* Obtain library cut off:
		*  if librabry mean + std <= 3k
		*	 edge_cut_off = 2
		*  else edge_cut_off = 5
		*/
		std::string edges_file_s(edges_file);
    	std::string edge_log(edges_file_s.substr(0,edges_file_s.find_last_of("/")+1) + "lib.txt");

    	FILE* edge_log_fp = fopen(edge_log.c_str(), "r");
    	if (edge_log_fp == NULL) {
    		fprintf(stderr, "Error opening file: %s\n", edge_log.c_str());
    		exit(EXIT_FAILURE);
    	}

    	while (fgets(line_lib, 255, edge_log_fp))
    	{
        	if (std::string(line_lib).find("Mean length") != std::string::npos) {
            	lib_mean = ::atof(std::string(line_lib).substr(31).c_str());
        	} else if (std::string(line_lib).find("Standard deviation") != std::string::npos) {
            	lib_stdev = ::atof(std::string(line_lib).substr(38).c_str());
        	}
    	}	
    	fclose(edge_log_fp);
    	
    	if ((lib_mean + lib_stdev) <= 3000) {
    		edge_bundle_size_thr = 2;
    	} else {
    		edge_bundle_size_thr = 2;
    	}



		while (!feof(edges_fp)) {
			if( fscanf(edges_fp, "%[^\n]\n", line) );

			// [ID1]\t[ORIENTATION1]\t[ID2]\t[ORIENTATION2]\t[DISTANCE]\t[STDEV]\t[SIZE]\n

			if (sscanf(line, "%s\t%*c\t%s\t%*c\t%d\t%d\t%d%*[^\n]", id1, id2, &distance, &stdev, &size) == 5) {
				auto it1 = contigs->find(id1);
				auto it2 = contigs->find(id2);

				if (it1 != contigs->end() && it2 != contigs->end() && (size >= edge_bundle_size_thr) &&
					((distance >= 0) || ((distance < 0) && (6 * stdev + Sigma::kmer_size + distance > 0)))) {
					if (((*it1).second->cluster() == (*it2).second->cluster())) {
						fprintf(filtered_edges_fp, "%s\n", line);
					} else {
						fprintf(skipped_edges_fp, "%s\n", line);
					}
					
				}
			}
		}

		fclose(filtered_edges_fp);
		fclose(skipped_edges_fp);
		fclose(edges_fp);
	} else {
		fprintf(stderr, "Error opening file: %s\n", edges_file);
		exit(EXIT_FAILURE);
	}
}
