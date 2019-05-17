#! /mnt/software/unstowable/anaconda/bin/python
import sys
import os
import argparse
import subprocess

def generating_mapping(kraken_report):
    kraken_dict = {}
    with open (kraken_report, "r") as fp:
        for line in fp:
            line = line.split("\t")
            if line[3] != 'S' and line[3] != '-':
                continue
            
            species = line[-1].split()
            try:
                species = species[0] + "_" + species[1]
            except:
                print("not a species?? : {}".format(species))
                continue               
            kraken_dict[line[4]] = [species]
    return kraken_dict
   

def kraken_binning(kraken_dict, kraken_out, output_file, abund_dict):

    cmd = "cat {} | cut -f1-4".format(kraken_out)
    kraken_result = subprocess.check_output(cmd, shell=True)

    #print(kraken_result)
    
    kraken_result = kraken_result.strip()
    
    kraken_result = kraken_result.split("\n")
    contigs_bin_dict = {}

    with open (output_file, "w") as fp:
        for item in kraken_result:
            tax_id = item.split("\t")[2]
            contigs = item.split("\t")[1]
            
            if tax_id in kraken_dict:            
                species = kraken_dict[tax_id][0]
                fp.write(contigs + "\t" + species + "\n")
                #output_file = ["{}/BINNING/{}.fasta".format(output_folder, species), "{}/ABUND/sample_{}_bin_{}.abund1".format(output_folder, sample, species)]
                #contigs_bin_dict[contigs] = output_file            
        
    return contigs_bin_dict
    
    
def main(args):
    kraken_dict = generating_mapping(args.kraken_report)

    abund_dict = {}
    #with open(args.abund, "r") as fp:
    #    for line in fp:
    #        abund_dict[line.split("\t")[0]] = line
                     
    contigs_bin_dict = kraken_binning(kraken_dict, args.kraken_out, args.output, abund_dict)

    

if __name__ == "__main__":   
    parser = argparse.ArgumentParser("kraken binning")

    mandatory = parser.add_argument_group("mandatory arguments")
    
    mandatory.add_argument("-k", "--kraken_out",
                           required=True,
                           help="kraken out file")
    mandatory.add_argument("-r", "--kraken_report",
                           required=True,
                           help="kraken out report file")
    mandatory.add_argument("-c", "--contig",
                           #required=True,
                           help="assembled contigs fasta file")
    mandatory.add_argument("-o", "--output",
                           required=True,
                           help="output file")
    mandatory.add_argument("-a", "--abund",
                           #required=True,
                           help="abundance file")
    args=parser.parse_args()
    
    main(args)




