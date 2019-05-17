#!/usr/bin/python

import os
import sys
from collections import defaultdict
from sets import Set

def main(args) :

        #Indenty the contigs fromn the scaffold
	contig_set = Set([])
	read_set = Set([])
        c1 = ""
        c2 = ""
	with open(args[2], 'r') as scaf_file :
            c1 = ""
            c2 = ""
	    for line in scaf_file :
                c2 = line.split()[0]
		contig_set.add(c1 + ":" + c2)
                c1 = c2

	print >> sys.stderr, "Finished reading names to retain"
	print >> sys.stderr, "Size of name list is " + str(len(contig_set))

        #Indenty read in a scaffold edge
	with open(args[0],'r') as edge_file :
		for line in edge_file :
                    c1 = line.split()[0]
                    c2 = line.split()[2]
		    if (c1 + ":" + c2) in contig_set or (c2 + ":" + c1) in contig_set:
			read_set.update(Set(["@" + name for name in line.split()[4].split(";")]))


        #Read the single read file
        with open(args[1],'r') as single_file :
                for line in single_file :
                        #print >> sys.stderr, line.split()[0]
                        read_set.update(Set(["@" + line.split()[0]]))

        #Identify the read sequence in the long read file
	line_count = 0
	with open(args[3],'r') as src_file :
		for line in src_file :

			if line_count == 0:
				if line.split()[0] in read_set :
					line_count = 1
                                        #print >> sys.stderr, line.split()[0]
					read_set.remove(line.split()[0])

			if line_count != 0 :
				print line.rstrip()
				line_count += 1
				line_count %= 5

	if len(read_set) != 0 :
		print >> sys.stderr, "Some issue"
		print >> sys.stderr, "Size of name list is " + str(len(read_set))


	##?##
	#Need to combine this with Filter_By_Name.py

if __name__ == '__main__' :                                                                                                                                                                                 
        main(sys.argv[1:])                                                                                                                                                                                  
        exit(0)    



