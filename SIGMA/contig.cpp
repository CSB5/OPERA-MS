#include <cstdlib>

#include <algorithm>
#include <vector>
#include <map>

#include <iostream>
#include <fstream>

#include <math.h>

#include <random>

#include "contig.h"

#include "sigma.h"



Contig::Contig(std::string id, int length) : id_(id), length_(length) {
	if (Sigma::contig_window_len > 0) {
		num_windows_ = (length_ - 2 * Sigma::contig_edge_len) / Sigma::contig_window_len;

		const int remainder = length_ - num_windows_ * Sigma::contig_window_len;
		left_edge_ = remainder / 2;
		right_edge_ = length_ - 1 - (remainder - left_edge_);
	} else {
		num_windows_ = 1;

		left_edge_ = Sigma::contig_edge_len;
		right_edge_ = length_ - 1 - Sigma::contig_edge_len;
	}

	modified_length_ = right_edge_ - left_edge_ + 1;

	initReadCounts();

	cluster_ = NULL;
}

Contig::Contig(std::string id, int length, int left_edge, int right_edge, int num_windows) :
	id_(id), length_(length),
	left_edge_(left_edge), right_edge_(right_edge), num_windows_(num_windows) {
	modified_length_ = right_edge_ - left_edge_ + 1;

	initReadCounts();

	cluster_ = NULL;
}

void Contig::initReadCounts() {
	sum_read_counts_ = new int[Sigma::num_samples];

	for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
		sum_read_counts_[sample_index] = 0;
	}

	read_counts_ = new int*[Sigma::num_samples];

	for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
		read_counts_[sample_index] = new int[num_windows_];

		for (int window_index = 0; window_index < num_windows_; ++window_index) {
			read_counts_[sample_index][window_index] = 0;
		}
	}
}

Contig::~Contig() {
	delete[] sum_read_counts_;

	for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
		delete[] read_counts_[sample_index];
	}

	delete[] read_counts_;
}

std::string Contig::id() const { return id_; }
int Contig::length() const { return length_; }

int Contig::modified_length() const { return modified_length_; }
int Contig::left_edge() const { return left_edge_; }
int Contig::right_edge() const { return right_edge_; }
int Contig::num_windows() const { return num_windows_; }

int* Contig::sum_read_counts() const { return sum_read_counts_; }
int** Contig::read_counts() const { return read_counts_; }

Cluster* Contig::cluster() const { return cluster_; }

void Contig::set_cluster(Cluster* cluster) { cluster_ = cluster; }


void ContigIO::save_contigs(const ContigMap* contigs, const char* sigma_contigs_file_path) {
	FILE* sigma_contigs_fp = fopen(sigma_contigs_file_path, "w");

	if (sigma_contigs_fp != NULL) {
		fprintf(sigma_contigs_fp, "%d %d %d %d\n", Sigma::num_samples, Sigma::contig_len_thr, Sigma::contig_edge_len, Sigma::contig_window_len);

		for (auto it = contigs->begin(); it != contigs->end(); ++it) {
			Contig* contig = (*it).second;

			fprintf(sigma_contigs_fp, "%s\t%d\t%d\t%d\t%d\n",
					contig->id().c_str(), contig->length(), contig->left_edge(), contig->right_edge(), contig->num_windows());

			for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
				fprintf(sigma_contigs_fp, "%d\n", contig->sum_read_counts()[sample_index]);

				for (int window_index = 0; window_index < contig->num_windows(); ++window_index) {
					fprintf(sigma_contigs_fp, "%d ", contig->read_counts()[sample_index][window_index]);
				}

				fprintf(sigma_contigs_fp, "\n");
			}
		}

		fclose(sigma_contigs_fp);
	} else {
		fprintf(stderr, "Error opening file: %s\n", sigma_contigs_file_path);
		exit(EXIT_FAILURE);
	}
}

void ContigIO::load_contigs(const char* sigma_contigs_file_path, ContigMap* contigs) {
	FILE* sigma_contigs_fp = fopen(sigma_contigs_file_path, "r");

	if (sigma_contigs_fp != NULL) {
		if( fscanf(sigma_contigs_fp, "%d %d %d %d\n", &Sigma::num_samples, &Sigma::contig_len_thr, &Sigma::contig_edge_len, &Sigma::contig_window_len) );

		while (!feof(sigma_contigs_fp)) {
			char id[256];
			int length, left_edge, right_edge, num_windows;

			if( fscanf(sigma_contigs_fp, "%s\t%d\t%d\t%d\t%d\n", id, &length, &left_edge, &right_edge, &num_windows) );

			Contig* contig = new Contig(id, length, left_edge, right_edge, num_windows);

			for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
				if( fscanf(sigma_contigs_fp, "%d\n", &contig->sum_read_counts()[sample_index]) );

				for (int window_index = 0; window_index < contig->num_windows(); ++window_index) {
					if( fscanf(sigma_contigs_fp, "%d ", &contig->read_counts()[sample_index][window_index]) );
				}

				if( fscanf(sigma_contigs_fp, "\n") );
			}

			contigs->insert(std::make_pair(id, contig));
		}

		fclose(sigma_contigs_fp);
	} else {
		fprintf(stderr, "Error opening file: %s\n", sigma_contigs_file_path);
		exit(EXIT_FAILURE);
	}
}





//Conservative R computation that remove 20% of the outlier in the assembly
//This method only works for single sample assembly
double compute_R(ContigMap* contigs) {
	
	double numerator = 0;
	double denominator = 0;
	
	
	double current_R = 0;
	double max_R = 0.001;
	//min_R should be positif
	double min_R = 1;
	
	std::vector<double*>* contig_info = new std::vector<double*>();
	int sample_index;

	std::string r_estimation_dist_file = Sigma::output_dir + "/r_estimate_distribution.dat";   
	FILE* r_estimate = fopen(r_estimation_dist_file.c_str(), "w");
	
	for (auto it = contigs->begin(); it != contigs->end(); ++it) {
		Contig* contig = (*it).second;

		if (contig->length() < 10000) continue;

		//fprintf(stderr, " *** contig length %d",contig->length());
		sample_index = 0;
		//for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
		int* read_counts = contig->read_counts()[sample_index];

		double* c_info = new double[3];
		
		//Mean
		double mean = 0.0;
		for (int window_index = 0; window_index < contig->num_windows(); ++window_index) {
			mean += read_counts[window_index];
		}
		mean /= contig->num_windows();
		
		
		//Variance
		double variance = 0.0;
		for (int window_index = 0; window_index < contig->num_windows(); ++window_index) {
			variance += (read_counts[window_index] - mean) * (read_counts[window_index] - mean);
		}
		variance /= contig->num_windows();

		if(variance - mean > 0){//Check for the division by zero
			c_info[0] = mean;
			c_info[1] = variance;
			//fprintf(stderr, "%f\t%f\n", mean, variance);
			//Need to check for the division by zero
			current_R = pow(mean, 2) / ( variance - mean);
			c_info[2] = current_R;
			contig_info->push_back(c_info);
			fprintf(r_estimate, "%s\t%f\t%f\t%f\n", contig->id().c_str(), mean, variance, current_R);
		}
	       
	}
	fclose(r_estimate);
	//fprintf(stderr, "*** Size:\n");
	//fprintf(stderr, "*** Size: %i\n", (int)contig_info->size());
	double least_square_R;
	if(contig_info->size() != 0){
		
		//Sort the contig info vector
		std::sort (contig_info->begin(), contig_info->end(), sort_contig_R);
	
	
		//Exclude the top/last 10% of the contig R values
		int excluded_contig = (int) (contig_info->size() / 10);
		std::vector<double*>::iterator it_start = contig_info->begin() + excluded_contig;
		std::vector<double*>::iterator it_last = contig_info->end() - excluded_contig;
		min_R = (*it_start)[2];
	
		double mean,variance;
		for (std::vector<double*>::iterator it = it_start; it != it_last; ++it){
			mean = (*it)[0];
			variance = (*it)[1];
			numerator += pow(mean, 4);
			denominator += pow(mean, 2) * variance - pow(mean, 3);
		}

	
		least_square_R = numerator / denominator;
	
	}
	
	else{
		fprintf(stderr, " *** WARNING no contigs larger than 10kb for r estimation");
		min_R = 20;
		least_square_R = 50;
		max_R = 150;
		fprintf(stderr, " *** Use predefined R values: min_R %f least square %f  max_R %f\n", min_R, least_square_R, max_R);
	}
	double res = -1;
	res = min_R;
	
	std::string r_estimation_file = Sigma::output_dir + "/r_estimate_value.dat";   
	r_estimate = fopen(r_estimation_file.c_str(), "w");
	fprintf(r_estimate, "%f\t%f\n", min_R, least_square_R);
	fclose(r_estimate);
	
	fprintf(stderr, " *** Global R: min_R %f least square %f  max_R %f\n *** RES: %f\n", min_R, least_square_R, max_R, res);
	
	//delete the vector
	for (std::vector<double*>::iterator it = contig_info->begin(); it != contig_info->end(); ++it){
		delete[] (*it);
	}
	
	delete contig_info;
		
	return res;
}

bool sort_contig_R (double* t, double* g) {
	return (t[2] < g[2]);
}

//This function work for multi smple but may be very sensitive to contigs outlying R value  
double compute_R_multi(ContigMap* contigs) {
	
	double numerator = 0;
	double denominator = 0;

	
	double current_R = 0;
	double max_R = 0.001;
	//min_R should be positif
	bool init_bin = false;
	double min_R = 1;
	
	for (auto it = contigs->begin(); it != contigs->end(); ++it) {
		Contig* contig = (*it).second;

		if (contig->length() < 10000) continue;

		//fprintf(stderr, " *** contig length %d",contig->length());
		
		for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
			int* read_counts = contig->read_counts()[sample_index];
			
			double mean = 0.0;

			for (int window_index = 0; window_index < contig->num_windows(); ++window_index) {
				mean += read_counts[window_index];
			}
			
			mean /= contig->num_windows();

			double variance = 0.0;
			
			for (int window_index = 0; window_index < contig->num_windows(); ++window_index) {
				variance += (read_counts[window_index] - mean) * (read_counts[window_index] - mean);
			}

			variance /= contig->num_windows();
            
			//fprintf(stderr, "%f\t%f\n", mean, variance);
			
			if(denominator == 0){
				min_R = current_R;
			}
			
			if(variance - mean > 0){
				numerator += pow(mean, 4);
				denominator += pow(mean, 2) * variance - pow(mean, 3);
			}
			
			//Need to check for the division by zero
			current_R = pow(mean, 2) / ( variance - mean);
			
			if(current_R > max_R){
				max_R = current_R;
			}
			if(current_R > 0 && 
			   (!init_bin || (init_bin && current_R  < min_R))){
				min_R = current_R;
				init_bin = true;
			}
		}
	}
	

	double least_square_R = numerator / denominator;
	

	double res = -1;
	/*if(Sigma::R_ESTIMATION_TYPE == 1){
		res = least_square_R;
	}
	else if(Sigma::R_ESTIMATION_TYPE == 2){
		res = max_R;
	}
	else if(Sigma::R_ESTIMATION_TYPE == 3){
		res = min_R;
	}

	if(Sigma::USE_WINDOW == 0){
		res = res / Sigma::contig_window_len;//need to be fixed
		}*/
	
	res = min_R;
	
	std::string r_estimation_file = Sigma::output_dir + "/r_estimate_value.dat";   
	FILE* r_estimate = fopen(r_estimation_file.c_str(), "w");
	fprintf(r_estimate, "%f\t%f\n", min_R, least_square_R);
	fclose(r_estimate);
	
	fprintf(stderr, " *** Global R: min_R %f least square %f  max_R %f\n *** RES: %f\n", min_R, least_square_R, max_R, res);
	
	
	return res;
}

//?// New additions from here.

int sample_dist( std::map< int, int > dist, int x_val ){

    int sample_val;

    if( dist.find( x_val ) != dist.end() ){
        sample_val = dist[ x_val ];
    }
    else {
        int l = x_val-1;
        int r = x_val+1;
        while( l != -1 ){
            if( dist.find(l) != dist.end() ) {
                break;
            }
            l--;
        }
        while( r < dist.rbegin()->first + 1 ){
            if( dist.find(r) != dist.end() ) {
                break;
            }
            r++;
        }

        if( l != -1 && r < dist.rbegin()->first + 1){
            sample_val =  (dist[l] + dist[r])/2;
        }
        else if( l != -1) {
            sample_val =  dist[l];
        }
        else if( r < dist.rbegin()->first + 1 ){
            sample_val =  dist[r];
        }
        else{
            sample_val = -1;
        }
    }
    return sample_val;
}
