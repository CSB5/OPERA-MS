#! /usr/bin/python

# Script for remapping the contigs to the assembly
# Need to supply paf mapping file along with the assembly and the contig file


import os
import sys
from collections import defaultdict

# Function to return the complementary base
def complement(base) :
	if base == "A" :
		return "T"
	elif base == "T" :
		return "A"
	elif base == "G" :
		return "C"
	elif base == "C" :
		return "G"

# Function to provide the reverse complement of the input sequence
def reverse_complement(seq) :
	# Reverse slice inverts the sequence
	rev = seq[::-1]
	rev_comp = ""
	for i in range(0,len(rev)) :
		rev_comp += complement(rev[i])

	return rev_comp


# Function to identify the contigs in the mapping and duplicates name in case of reats contigs
def check_valid_map(line, valid_contig_map, invalid_contig_set, contig_mapping_ID, contig_cmp_id) :
	
	contig_name = line.split()[0]
	contig_len = int(line.split()[1])
	contig_start = int(line.split()[2])
	contig_end = int(line.split()[3])


	map_len = abs(contig_start - contig_end)
	map_fraction = float(map_len) / float(contig_len)
	
	#if map_fraction > 0.99 :
        if map_fraction > 0 :
		if contig_name in valid_contig_map :
			del valid_contig_map[contig_name]
			invalid_contig_set.add(contig_name)
			print >> sys.stderr, "Invalid Repeat contig :	" + contig_name + "_" + str(contig_len) + "_" + str(map_len)
			return False
		elif contig_name in invalid_contig_set :
			return False
		else :
			contig_orientation = line.split()[4]
			scaffold_start = int(line.split()[7]) #- contig_start
			scaffold_end = int(line.split()[8]) #+ (contig_len - contig_end)

                        #To get multiple mapping of the same contig
                        if contig_name not in contig_mapping_ID:
                                contig_mapping_ID[contig_name] = []

                        contig_ID = contig_name + "_XXX_" + str(contig_cmp_id)
                        contig_mapping_ID[contig_name].append(contig_ID)
                        #
                        #
                        valid_contig_map[contig_ID] = {"contig_start":contig_start,"contig_end":contig_end,"contig_len":contig_len,"contig_orientation":contig_orientation,"scaffold_start":scaffold_start,"scaffold_end":scaffold_end}
			return True

	else :
		print >> sys.stderr, "Invalid Short map :	" + contig_name + "_" + str(contig_len) + "_" + str(map_len)
		return False

# Function returns dictionary with all the contig sequences in the + orientation
def get_sequences(contig_path, valid_contig_map, contig_mapping_ID) :
        #print >> sys.stderr, " ************** Get sequence "
	contig_seq = {}
	extract_flag = 0
	contig_name = ""

	with open(contig_path,'r') as contig_file :

		for line in contig_file :

			if extract_flag == 1 :
				#
                                seq = line.rstrip()
                                #print >> sys.stderr, " Contig name "  + contig_name
                                #if valid_contig_map[contig_name]["contig_orientation"] == "-" :
                                for contig_m_id in contig_mapping_ID[contig_name]:
                                        #
                                        c_seq = seq
				        if valid_contig_map[contig_m_id]["contig_orientation"] == "-" :
					        c_seq = reverse_complement(seq)
				                #contig_seq[contig_name] = seq
                                        #print >> sys.stderr, " Contig ID for seq "  + contig_m_id
                                        contig_seq[contig_m_id] = c_seq
				contig_name = ""
				extract_flag = 0

			else :
				if line[1:].rstrip() in contig_mapping_ID : #valid_contig_map :
					contig_name = line[1:].rstrip()
                                        
					extract_flag = 1
        #print >> sys.stderr, " ************** End gett sequence "                
	return contig_seq

# Function to populate and return a sorted list of tuples(scaffold_start,contig_name) of the start sites of the valid contigs in the scaffold
def get_start_sites(valid_contig_map, flag_extend_contig) :
	start_sites_map = {}
        contig_start_pos = 0
	for contig_name in valid_contig_map :
		# Adjust the start site for mismatch at the beginning of contig
                if(flag_extend_contig) :
		        contig_start_pos = valid_contig_map[contig_name]["scaffold_start"] - valid_contig_map[contig_name]["contig_start"]
                else :
                        contig_start_pos = valid_contig_map[contig_name]["scaffold_start"]
                start_sites_map[contig_start_pos] = contig_name

	return sorted(start_sites_map.items())


# Function to replace the contigs in the sequence and print it to file
def get_final_assembly(initial_assembly_seq, valid_contig_map, contig_seq,contig_start_sites, flag_extend_contig) :
	pos = 0
	initial_assembly = list(initial_assembly_seq)
	#print >> sys.stderr, initial_assembly
	#print len(initial_assembly)
	#print "HERE"
	final_assembly = ""
	# Contigs are sorted in order of starting point
        next_contig_end = 0
        contig_cut_end = 0
	for next_contig_start,next_contig_name in contig_start_sites :

                if(flag_extend_contig) :
		        next_contig_end = valid_contig_map[next_contig_name]["scaffold_end"] + valid_contig_map[next_contig_name]["contig_len"] - valid_contig_map[next_contig_name]["contig_end"]
                else :
                        next_contig_end = valid_contig_map[next_contig_name]["scaffold_end"]
                		
		# Debug
		#print    next_contig_name + "\t" + str(next_contig_start) + "\t" + str(next_contig_end) + "\t" + str(valid_contig_map[next_contig_name]["contig_len"])
		
		# Print all the gap contigs
		for i in range(pos, next_contig_start) :
			final_assembly += initial_assembly[i]
                        
		if(pos < next_contig_start) :
			pos = next_contig_start

		# Cut the current contig in case of overlap with previous contig.
                #0 in case of non overlap
                if(flag_extend_contig) :
                        contig_cut = pos - int(next_contig_start)
                else :
                        contig_cut = pos - int(next_contig_start) + valid_contig_map[next_contig_name]["contig_start"]
		# Adjust coords for unmapped end segments
		
		# Get the sequence of the contig which is to be inserted
		insert_seq = contig_seq[next_contig_name]

                if(flag_extend_contig) :
                        contig_cut_end = valid_contig_map[next_contig_name]["contig_len"]
                else :
                        contig_cut_end = valid_contig_map[next_contig_name]["contig_end"]
                                                                
		for i in range(contig_cut, contig_cut_end) :
			final_assembly += insert_seq[i]

                #????????
		if(pos < next_contig_end) :
			pos = next_contig_end# + 1

	for i in range(pos, len(initial_assembly)) :
		final_assembly += initial_assembly[i]

	return  final_assembly



def main(args) :

	paf_path = args[0]
	initial_assembly_path = args[1]
	contig_path = args[2]

        flag_extend_contig = 0
        
	valid_contig_map = {}
        contig_mapping_ID = {}
	invalid_contig_set = set([])
	contig_seq = {}
	initial_assembly = ""

	# File handler for writing the final output file
	final_assembly_file = open(".".join(initial_assembly_path.split(".")[:-1])+"_remapped.fasta",'w')
	
        ID_cmp = 1
	with open(paf_path,'r') as  paf_file :
		for line in paf_file :
			temp_discard = check_valid_map(line, valid_contig_map, invalid_contig_set, contig_mapping_ID, ID_cmp)
                        ID_cmp += 1 

	contig_seq = get_sequences(contig_path, valid_contig_map, contig_mapping_ID)

	contig_start_sites = get_start_sites(valid_contig_map, flag_extend_contig)

	with open(initial_assembly_path,'r') as initial_assembly_file :
		for line in initial_assembly_file :
			
			if line.startswith(">") :
				
				if initial_assembly != "" :
					final_assembly = get_final_assembly(initial_assembly, valid_contig_map, contig_seq,contig_start_sites, flag_extend_contig)
					print >> final_assembly_file, final_assembly
				
				print >> final_assembly_file, line.rstrip()
				initial_assembly = ""
			
			else :
				initial_assembly += line.rstrip()		
				##?## Need to update for more than one gap filled scaffold

		if initial_assembly != "" :
			final_assembly = get_final_assembly(initial_assembly, valid_contig_map, contig_seq,contig_start_sites, flag_extend_contig)
			print >> final_assembly_file, final_assembly



if __name__ == "__main__" :
	main(sys.argv[1:])
	exit(0)




