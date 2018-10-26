#include <cstdlib>

#include "cluster.h"

#include "sigma.h"

Cluster::Cluster() {
	num_contigs_ = -1;
	//contigs_ = new Contig*[num_contigs_];
	//contigs_[0] = contig;

	length_ = 0;

	sum_read_counts_ = new int[Sigma::num_samples];
	mean_arrival_rates_ = new double[Sigma::num_samples];
	median_arrival_rates_ = new double[Sigma::num_samples];

	
	child1_ = NULL;
	child2_ = NULL;

	score_ = 0.0;
	model_score_ = 0.0;
	modeled_ = false;
	connected_ = false;

}

Cluster::Cluster(Contig* contig) {
	num_contigs_ = 1;
	contigs_ = new Contig*[num_contigs_];

	contigs_[0] = contig;

	length_ = contig->modified_length();

	sum_read_counts_ = new int[Sigma::num_samples];
	mean_arrival_rates_ = new double[Sigma::num_samples];
	median_arrival_rates_ = new double[Sigma::num_samples];

	for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
		sum_read_counts_[sample_index] = contig->sum_read_counts()[sample_index];
		mean_arrival_rates_[sample_index] = sum_read_counts_[sample_index] / (double) length_;
	}
	compute_median_arrival_rate();
	
	
	child1_ = NULL;
	child2_ = NULL;

	score_ = 0.0;
	model_score_ = 0.0;
	modeled_ = false;
	connected_ = false;

	contig->set_cluster(this);
}

Cluster::Cluster(Cluster* child1, Cluster* child2) {
	num_contigs_ = child1->num_contigs_ + child2->num_contigs_;
	contigs_ = new Contig*[num_contigs_];
	

	for (int contig_index = 0; contig_index < child1->num_contigs(); ++contig_index) {
		contigs_[contig_index] = child1->contigs_[contig_index];
	}

	int contig_index_offset = child1->num_contigs();

	for (int contig_index = 0; contig_index < child2->num_contigs(); ++contig_index) {
		contigs_[contig_index_offset + contig_index] = child2->contigs_[contig_index];
	}
	
	sum_read_counts_ = new int[Sigma::num_samples];
	mean_arrival_rates_ = new double[Sigma::num_samples];
	median_arrival_rates_ = new double[Sigma::num_samples];
	

	update_cluster_length_read_count(child1, child2);

	child1_ = child1;
	child2_ = child2;

	score_ = 0.0;
	model_score_ = 0.0;
	modeled_ = false;
	connected_ = false;

	for (int contig_index = 0; contig_index < num_contigs_; ++contig_index) {
		contigs_[contig_index]->set_cluster(this);
	}

	delete[] child1->contigs_;
	delete[] child2->contigs_;
}

void Cluster::update_cluster_length_read_count(Cluster* child1, Cluster* child2){
	length_ = child1->length_ + child2->length_;
	
	if(num_contigs_ == -1){
		num_contigs_ = child1->num_contigs_ + child2->num_contigs_;
	}

	for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
		sum_read_counts_[sample_index] = child1->sum_read_counts_[sample_index] + child2->sum_read_counts_[sample_index];
		mean_arrival_rates_[sample_index] = sum_read_counts_[sample_index] / (double) length_;
	}
	
	//fprintf(stdout, " *** update_cluster_length_read_count: length %d\tsum_read_counts_ %d\tmean_arrival_rates_  %f\n", length_, sum_read_counts_[0], mean_arrival_rates_[0]);

	//compute_median_arrival_rate();
}


Cluster::~Cluster() {
	delete[] sum_read_counts_;
	delete[] mean_arrival_rates_;
	delete[] median_arrival_rates_;
}

Contig** Cluster::contigs() const { return contigs_; }
void Cluster::set_contigs(Contig** contigs) { contigs_ = contigs; }
void Cluster::set_contigs(std::vector<Contig*> *v_contig){
	contigs_ = new Contig*[v_contig->size()];
	int cmp = 0;
	for (auto it = v_contig->begin(); it != v_contig->end(); ++it) {
		contigs_[cmp] = *it; 
		cmp++;
	}
}

int Cluster::num_contigs() const { return num_contigs_; }
int Cluster::length() const { return length_; }
int* Cluster::sum_read_counts() const { return sum_read_counts_; }
//
double* Cluster::arrival_rates() const { 
	double* res = NULL;
	if(Sigma::AR_TYPE == 1){
		res = mean_arrival_rates_; 
	}
	else if(Sigma::AR_TYPE == 2){
		res = median_arrival_rates_; 
	}
	return res;
}

void Cluster::compute_median_arrival_rate(){
	//To compute the median read count work only if use 1 sample
	std::vector<int>* contig_bin_read_count = new std::vector<int>();
	double median;
	Contig* contig;
	
	for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
		for (int contig_index = 0; contig_index  < num_contigs(); ++contig_index) {
			contig = contigs_[contig_index];
			for (int window_index = 0; window_index < contig->num_windows(); ++window_index) {
				contig_bin_read_count->push_back(contig->read_counts()[sample_index][window_index]);
			}
		}
		
		size_t size = contig_bin_read_count->size();
		
		std::sort(contig_bin_read_count->begin(), contig_bin_read_count->end());
		
		if (size  % 2 == 0){
			median = ((*contig_bin_read_count)[size / 2 - 1] + (*contig_bin_read_count)[size / 2]) / 2;
		}
		else{
			median = (*contig_bin_read_count)[size / 2];
		}
		median_arrival_rates_[sample_index] = median/Sigma::contig_window_len;
		contig_bin_read_count->clear();
	}
	delete contig_bin_read_count;
}

double* Cluster::mean_arrival_rates()const { return mean_arrival_rates_;}
double* Cluster::median_arrival_rates()const { return median_arrival_rates_;}

Cluster* Cluster::child1() const { return child1_; }
Cluster* Cluster::child2() const { return child2_; }

double Cluster::score() const { return score_; }
double Cluster::model_score() const { return model_score_; }
bool Cluster::modeled() const { return modeled_; }
bool Cluster::connected() const { return connected_; }

void Cluster::set_score(double score) { score_ = score; }
void Cluster::set_model_score(double model_score) { model_score_ = model_score; }
void Cluster::set_modeled(bool modeled) { modeled_ = modeled; }
void Cluster::set_connected(bool connected) { connected_ = connected; }

void Cluster::set_cluster(Cluster* c1, Cluster* c2) {
	child1_ = c1;
	child2_ = c2;
}

void Cluster::set_ID(std::string* ID){cluster_ID = ID;}
std::string* Cluster::get_ID()const{return cluster_ID;}


//void Cluster::set_iterator_un_set(std::unordered_set<Cluster*>::iterator it){ it_un_set_ = it; }
//std::unordered_set<Cluster*>::iterator Cluster::get_iterator_un_set()const { return it_un_set_; }

