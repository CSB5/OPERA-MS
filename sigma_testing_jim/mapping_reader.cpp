#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <math.h>
#include <iomanip>
#include <iostream>
using namespace std;

#include "mapping_reader.h"

#include "sigma.h"
#include "string.h"

MappingReader::~MappingReader() {}


SAMReader::SAMReader() {}

void SAMReader::read(const char* mapping_file, int sample_index, ContigMap* contigs) {
	if (mapping_file[0] == '-') {
		readInputStream(stdin, sample_index, contigs);
	} else {
		FILE* mapping_fp = fopen(mapping_file, "r");

		if (mapping_fp != NULL) {
			readInputStream(mapping_fp, sample_index, contigs);
		} else {
			fprintf(stderr, "Error opening file: %s\n", mapping_file);
			exit(EXIT_FAILURE);
		}

		fclose(mapping_fp);
	}

	for (auto it = contigs->begin(); it != contigs->end(); ++it) {
		Contig* contig = (*it).second;

		for (int window_index = 0; window_index < contig->num_windows(); ++window_index) {
			contig->sum_read_counts()[sample_index] += contig->read_counts()[sample_index][window_index];
		}
	}
}

void SAMReader::readInputStream(FILE* mapping_fp, int sample_index, ContigMap* contigs) {
	char contig_id[256];
	int read_pos;
	int flag;
	char qname[256];
	char rdsequence[256];

	while (!feof(mapping_fp)) {
		// QNAME\tFLAG\tRNAME\tPOS\tMAPQ\tCIGAR\tRNEXT\tPNEXT\tTLEN\tSEQ\tQUAL\n
		
		// if (fscanf(mapping_fp, "%*[^\t]\t%*[^\t]\t%[^\t]\t%d\t%*[^\n]\n", contig_id, &read_pos) == 2) {
		if (fscanf(mapping_fp, "%[^\t]\t%d\t%[^\t]\t%d\t%*[^\t]\t%*[^\t]\t%*[^\t]\t%*[^\t]\t%*[^\t]\t%[^\t]\t%*[^\n]\n", qname, &flag, contig_id, &read_pos, rdsequence) == 5) {
			auto it = contigs->find(contig_id);
			
			//Unmapped read or contig is too small
			if (it == contigs->end()) continue;

			//CHECK only 1st read to get fragment count
			////float nearest = roundf(qname * 10); 	// / 10;
			//std::ostringstream ss;
			//ss << std::setprecision( 1 ) << std::fixed << nearest;
			//std::string rd = ss.str();
			//std::string rd2 = ".2";
			//fprintf(stderr, "qname %s\n", qname);
			if(qname[strlen(qname) -1] == '2' ) continue;
			//if(strstr(rd.c_str(),rd2.c_str())) continue;
			//if(qnamenum == 2) continue;
			//fprintf(stderr, "passed qname %s\n", qname);cin.get();
			
			//if FLAG=16 (reverse), add rdsequence to read_pos to obtain right end of read 1
			if(flag == 16) {
				read_pos = read_pos + static_cast<int>(strlen(rdsequence)) - 1;
				//read_pos = read_pos - 1;
			}

			Contig* contig = (*it).second;

			--read_pos; // POS is 1-based

			if (read_pos >= contig->left_edge() && read_pos <= contig->right_edge()) {
				if (Sigma::contig_window_len > 0) {
					const int window_index = (read_pos - contig->left_edge()) / Sigma::contig_window_len;

					contig->read_counts()[sample_index][window_index]++;
				} else {
					contig->read_counts()[sample_index][0]++;
				}
			}
		} else {
			fprintf(stderr, "Invalid SAM file\n");
			exit(EXIT_FAILURE);
		}
	}
}
