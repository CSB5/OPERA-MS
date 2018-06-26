#include <cstdlib>

#include "edge.h"

#include "sigma.h"
#include "cluster.h"

Edge::Edge(Contig* contig1, Contig* contig2) :
		contig1_(contig1), contig2_(contig2) {}

void Edge::computeDistance() {
	//const double arr_rate_zero_thr = 1e-6;
	
	distance_ = 0.0;
	double curr_distance = 0;
	
	//int num_non_zero_samples = Sigma::num_samples;
	for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
		double arr_rate1 = contig1_->sum_read_counts()[sample_index] / (double) (contig1_->modified_length());
		double arr_rate2 = contig2_->sum_read_counts()[sample_index] / (double) (contig2_->modified_length());
		
		/*if (arr_rate1 < arr_rate_zero_thr && arr_rate2 < arr_rate_zero_thr) {
			num_non_zero_samples--;
			} else*/
		if (arr_rate1 < arr_rate2) {
			curr_distance = (arr_rate2 - arr_rate1) / arr_rate2;
		} else {
			curr_distance = (arr_rate1 - arr_rate2) / arr_rate1;
		}

		if(curr_distance > distance_){
			distance_ = curr_distance;
		}
	}

	//distance_ /= num_non_zero_samples;
}

void Edge::computeDistanceCluster() {
	const double arr_rate_zero_thr = 1e-6;

	distance_ = 0.0;

	int num_non_zero_samples = Sigma::num_samples;
	for (int sample_index = 0; sample_index < Sigma::num_samples; ++sample_index) {
		double arr_rate1 = contig1_->cluster()->mean_arrival_rates()[sample_index];
		double arr_rate2 = contig2_->cluster()->mean_arrival_rates()[sample_index];
		
		if (arr_rate1 < arr_rate_zero_thr && arr_rate2 < arr_rate_zero_thr) {
			num_non_zero_samples--;
		} else if (arr_rate1 < arr_rate2) {
			distance_ += (arr_rate2 - arr_rate1) / arr_rate2;
		} else {
			distance_ += (arr_rate1 - arr_rate2) / arr_rate1;
		}
	}

	distance_ /= num_non_zero_samples;
}

Contig* Edge::contig1() const { return contig1_; }
Contig* Edge::contig2() const { return contig2_; }
double Edge::distance() const { return distance_; }


bool operator==(const Edge& edge1, const Edge& edge2) {
	return ((edge1.contig1_->id() == edge2.contig1_->id() && edge1.contig2_->id() == edge2.contig2_->id()) ||
	(edge1.contig1_->id() == edge2.contig2_->id() && edge1.contig2_->id() == edge2.contig1_->id()));
}


bool EdgeComparator::operator()(const Edge& edge1, const Edge& edge2) const {
	return edge1.distance() > edge2.distance();
}


size_t EdgeHash::operator()(const Edge& edge) const {
	if (edge.contig1()->id() < edge.contig2()->id()) {
		return string_hash_(edge.contig1()->id() + " " + edge.contig2()->id());
	} else {
		return string_hash_(edge.contig2()->id() + " " + edge.contig1()->id());
	}
}
