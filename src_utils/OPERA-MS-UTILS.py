import gzip
import sys
import os
import argparse
import subprocess
import config_parser
import re
from binner_analysis import *
from genfunc import *
from checkm_analysis import *
from novel_analysis import *
from kraken_analysis import *
   
# default to look for polish contig, if polished contig not exist, use unpolish.
# if user want to always use unpolish contig, specifiy flag --unploshed SHOULD WE ADD THAT ?

#For problem due to gz file that are uncompress before running opera-ms
def get_long_read_file(read_file):
    res_file = read_file
    if not os.path.isfile(res_file):
        res_file = res_file + ".gz"
        if not os.path.isfile(res_file):
            exit("Long read file not found : " + read_file + " or " + res_file)
    return res_file

#def run_hybrid_binning(contig_file, short_read1, short_read2, assembly_dir, sample_name, nb_thread):
            
def download_utils_db(db_type):
    try:
        cmd = "mkdir {}/../utils_db".format(util_dir)
        run_exe(cmd, True)
    except:
        pass
    
    #mash_db = "/home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/OPERA-MS-DEV/OPERA-MS/genomeDB_Sketch.msh";
    if db_type == "read-concordance":
        if  not os.path.exists(("{}/../utils_db/kraken_db/hash.k2d".format(util_dir))):
        
            cmd = "wget -nc -P {}/../utils_db/ ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz".format(util_dir)
            run_exe(cmd, True)
            cmd = "tar -xvzf {}/../utils_db/minikraken_8GB_202003.tgz -C {}/../utils_db".format(util_dir, util_dir)
            run_exe(cmd, True)
            cmd = cmd = "mv {}/../utils_db/minikraken_8GB_20200312 {}/../utils_db/kraken_db".format(util_dir, util_dir)
            run_exe(cmd, True)
            cmd = "rm {}/../utils_db/minikraken_8GB_202003.tgz".format(util_dir)
            run_exe(cmd, True)
    if db_type == "novel-species":
        cmd = "wget --no-check-certificate -nc -O {}/../utils_db/GTDB.msh https://ndownloader.figshare.com/files/22471499".format(util_dir)
        run_exe(cmd, True)
        cmd = "wget --no-check-certificate -nc -O {}/../utils_db/MAG.msh https://ndownloader.figshare.com/files/22471505".format(util_dir)
        run_exe(cmd, True)
    
def read_taxonomy_file(genomes_dir, taxonomy, db_genome_dir, genome_list, genome_length):
    OUT_LIST = open(genome_list, "w")
    OUT_LENGTH = open(genome_length, "w")
    FILE = open(taxonomy, "r")
    line_list = {}
    tax_info = {}
    genome = ""
    species_name = ""

    DIRECTORY_MAX_NB_GENOME = 4000
    dir_id = 1
    current_db_genomes_dir = db_genome_dir + "_" + str(dir_id)
    create_dir(current_db_genomes_dir)
    nb_genomes_in_dir = 0
    for line in FILE:
        #print line
        line_list = (line.rstrip('\n')).split("\t")
        genome = line_list[0]
        tax_info = line_list[1].split(";")
        species_name = ""
        for t in tax_info:
            current_tax = t.split("__")
            if current_tax[0] == "s":
                species_name = current_tax[1].replace(" ", "_")
                
        if species_name == "":
            exit("Malformed taxonomy file: " + taxonomy + "\n" + line)
        else:
            #
            if nb_genomes_in_dir == DIRECTORY_MAX_NB_GENOME:
                dir_id += 1
                current_db_genomes_dir = db_genome_dir + "_" + str(dir_id)
                nb_genomes_in_dir = 0
                create_dir(current_db_genomes_dir)
            #copy the file in the opera-ms-db directory
            novel_genome_name = current_db_genomes_dir + "/" + species_name + "__" + genome  #Need to fix this !!!
            #Check for gzip file
            run_exe("cp {}/{} {}".format(genomes_dir, genome, novel_genome_name), True)
            OUT_LIST.write(novel_genome_name + "\n")
            #Write the length
            OUT_LENGTH.write("{}\t{}\n".format(novel_genome_name, "\t".join([str(x) for x in compute_genome_length(novel_genome_name)])))
            nb_genomes_in_dir += 1
    OUT_LENGTH.close()
    OUT_LIST.close()
    FILE.close()

def compute_genome_length(genome):
    res = [0, 0]
    with gzip.open(genome, "r") as FILE:
        for line in FILE:
            if not (line[0] == ">"):
                res[1] += (len(line)-1)
            else:
                res[0] += 1
    return res
        
def create_mash_sketch(genome_list, out_file, nb_thread):
    run_exe("{}/mash sketch -o {} -p {} -l {}".format(util_dir, out_file, nb_thread, genome_list), True) #other potential parameters -k -s
    
def opera_ms_db(genomes_dir, taxonomy, db_name, nb_thread):
    #Create the output database directory
    create_dir(db_name)
    genome_db = db_name + "/genomes"
        
    #Read the taxonomy file
    #The genome name will be renamed according to the taxonomy file: SPECIES_NAME__GENOME_NAME
    genome_list = db_name + "/genomes_list.txt"
    genome_size = db_name + "/genomes_length.txt"
    read_taxonomy_file(genomes_dir, taxonomy, genome_db, genome_list, genome_size)
            
    #Create mash sketch
    create_mash_sketch(genome_list, db_name+"/genomes.msh", nb_thread)
            
    #Create 8Gb kraken db    


def run_circular_sequence_identification(assembly_dir):
    scaffold_file = assembly_dir + "/intermediate_files/opera_long_read/scaffolds.scaf"
    edge_files_dir = assembly_dir + "/intermediate_files/read_mapping/"
    ana_dir = assembly_dir + "/circular_sequence"
    contig_info_file = assembly_dir + "/contig_info.txt"
    contig_file = get_contig_file(assembly_dir)
    create_dir(ana_dir)
    run_exe(util_dir + "/../bin/detect_circular_scaffold.pl " + contig_file + " " + scaffold_file + " " + edge_files_dir + " " + contig_info_file + " " + ana_dir, True)

def check_software(cmd, tool):
    try:
        run_exe(cmd, False)
        print("{}    functioning".format(tool))
    except Exception as e:
        print("{}    *** NOT FUNCTIONING ***".format(tool))
        
def check_installation():

    cmd = "{}/metabat -h 2> /dev/null".format(util_dir)    
    check_software(cmd, "MetaBAT2")

    cmd = "{}/kraken2/kraken2 -v > /dev/null ".format(util_dir)
    check_software(cmd, "Kraken2")
        
    cmd = "{}/checkm -h > /dev/null".format(util_dir)    
    check_software(cmd, "CheckM")
    
    cmd = "perl {}/maxbin2/run_MaxBin.pl -h > /dev/null".format(util_dir)
    check_software(cmd, "MaxBin2")
    
    
def main(args):
    
    command = args.command
    nb_thread = 0
    if command == "read-concordance" or command == "circular-sequence" or command == "binning" or command == "bin-evaluation" or command == "mash" :
        
        #Parse the config file
        config_dict = read_opera_ms_config_file(args.config)
                        
        #Set the number of thread        
        nb_thread = args.thread

        if command == "binning":
            bin_method = args.binner
            
            if args.sample_name == None:                
                sample_name = os.path.basename(os.path.normpath(config_dict["OUTPUT_DIR"]))                
            else:
                sample_name = args.sample_name
            #sample_name = config_dict["OUTPUT_DIR"].split("/")[-1]
            run_binner(bin_method, sample_name, config_dict["OUTPUT_DIR"], config_dict["ILLUMINA_READ_1"], config_dict["ILLUMINA_READ_2"], nb_thread)
            
        elif command == "bin-evaluation":
            checkm_analysis(config_dict["OUTPUT_DIR"], args.binner, nb_thread, args.high_qual_mags, args.medium_qual_mags)

        elif command == "read-concordance":
       
            abundance_threshold = args.abundance_threshold
            
            run_kraken2(config_dict["OUTPUT_DIR"], config_dict["ILLUMINA_READ_1"], config_dict["ILLUMINA_READ_2"], get_long_read_file(config_dict["LONG_READ"]), nb_thread, float(abundance_threshold))
                
        elif command == "circular-sequence":
            run_circular_sequence_identification(config_dict["OUTPUT_DIR"])
            
    #Command without config file
    else:
        if command == "novel-species":
            run_novel_species_analysis(args.known_species, args.taxonomy_database, args.thread, args.configs, args.binner, args.mags_qual, args.out)
            
        if command == "opera-ms-db":
            opera_ms_db(args.genomes_dir, args.taxonomy, args.db_name, args.thread)
            
        elif command == "utils-db":            
            download_utils_db(args.dbtype)            

        elif command == "check-dependency":
            check_installation()

def run_binner(binner, sample_name, assembly_dir, short_read1, short_read2, nb_thread):
    
    contig_file = get_contig_file(assembly_dir)
    
    if binner == "maxbin2":
        run_maxbin2(contig_file, short_read1, assembly_dir, sample_name, nb_thread)
        
    if binner == "metabat2":
        run_metabat2(contig_file, short_read1, short_read2, assembly_dir, sample_name, nb_thread)
        
    if binner == "hybrid":
        run_hybrid_binning(contig_file, short_read1, short_read2, assembly_dir, sample_name, nb_thread)
        

if __name__ == "__main__":   

    parser = argparse.ArgumentParser()
    
    #group = parser.add_mutually_exclusive_group()
    #The type of software
    subparsers = parser.add_subparsers(help='commands', dest='command')

    #this
    mandatory = parser.add_argument_group("mandatory arguments")

    #opera-db
    opera_db_parser = subparsers.add_parser('opera-ms-db', help='Generate a custom OPERA-MS genome databasee')
    mandatory = opera_db_parser.add_argument_group("mandatory arguments")
    opera_db_parser._action_groups[-1].add_argument("-g", "--genomes-dir",  required=True, help='Directory that contains genome files')
    opera_db_parser._action_groups[-1].add_argument("-x", "--taxonomy",  required=True, help='Species name of each genomes')
    opera_db_parser._action_groups[-1].add_argument("-d", "--db-name",  required=True, help='Database name')
    opera_db_parser.add_argument("-t", "--thread", help='Number of threads [default: 2]',  default=2, type = int)

    #kraken
    kraken_parser = subparsers.add_parser('read-concordance', parents=[config_parser.parser], help='Compute the abundance profile correlation between long and short-reads')
    kraken_parser.add_argument("-a", "--abundance-threshold", default=0.1, help="Lower percentage abundance threshold [default: 0.1]", type=float)
    kraken_parser.add_argument("-t", "--thread", help='Number of threads [default: 2]', default=2, type = int)
    
    #binner
    binner_parser = subparsers.add_parser('binning', parents=[config_parser.parser], help='Streamline binning of OPERA-MS assembled contigs using MetaBAT2 or MaxBin2')
    binner_parser.add_argument("-b", "--binner",  required=False, default = "metabat2", choices=["maxbin2", "metabat2"], help='binning method [default: metabat2]' )
    binner_parser.add_argument("-s", "--sample-name",  required=False, help='Sample name [default: OPERA-MS output folder]')
    binner_parser.add_argument("-t", "--thread", help='Number of threads [default: 2]', default=2, type = int)
       
    #checkm
    checkm_parser = subparsers.add_parser('bin-evaluation', parents=[config_parser.parser], help='Streamline bin evaluation using CheckM')
    checkm_parser.add_argument("-b", "--binner",  required=False, default = "metabat2", choices=["maxbin2", "metabat2", "opera_ms_clusters"], help = "Bins for evaluation [default: MetaBat2]")
    checkm_parser.add_argument("-H", "--high-qual-mags",  default="90,5", help = 'High quality bins completeness and contamination thresholds [default: 90,5]', type=str)
    checkm_parser.add_argument("-M", "--medium-qual-mags",  default="50,10", help = 'Medium quality bins completeness and contamination thresholds [default: 50,10]', type=str)
    checkm_parser.add_argument("-t", "--thread", help='Number of threads [default: 2]', default=2, type = int)
    #circular identification
    #circular_sequence_parser = subparsers.add_parser('circular-sequence', parents=[config_parser.parser], help='Identify circular sequences')
    
    #novel species
    novel_species_parser = subparsers.add_parser('novel-species', help='Identification of OPERA-MS MAGs most closely related species and identification of novel species')
    mandatory = novel_species_parser.add_argument_group("mandatory arguments")
    novel_species_parser._action_groups[-1].add_argument("-o", "--out",  required=True, help='Output directory')    
    #
    novel_species_parser.add_argument("-k", "--known-species",  required=False, default = "{}/../utils_db/MAG.msh".format(util_dir), help=argparse.SUPPRESS)#'Mash sketch of known species reference genomes (default utils_db/small_newgut_segata.msh)')
    novel_species_parser.add_argument("-x", "--taxonomy-database",  required=False, default = "{}/../utils_db/GTDB.msh".format(util_dir), help=argparse.SUPPRESS) #'Mash sketch of reference genomes with taxonomy info (default utils_db/genomes)')    
    novel_species_parser.add_argument("-b", "--binner",  required=False, default = "metabat2",  choices=["maxbin2", "metabat2", "opera_ms_clusters"], help='Bins used for the analysis [default: MetaBat2]')
    novel_species_parser.add_argument('configs', metavar='C', nargs='+', help='Path to OPERA-MS configuration file(s)')
    #
    novel_species_parser.add_argument("-q", "--mags-qual", help='Quality of the MAGS used [default: high]', choices=["high", "medium"], default="high")
    novel_species_parser.add_argument("-c", "--cluster-threshold", help='Maximum distance at which 2 genomes are considered of the same species [default: 0.05]', default=0.05, type = float)
    #
    novel_species_parser.add_argument("-t", "--thread", help='Number of threads [default: 2]', default=2, type = int)
       
    #utils-db
    utils_db_parser = subparsers.add_parser('utils-db', help=' Download the database required by the utils command')
    mandatory = utils_db_parser.add_argument_group("mandatory arguments")
    utils_db_parser._action_groups[-1].add_argument("-db", "--dbtype", choices = ["read-concordance", "novel-species"], required=True, help='read-concordance or novel species analysis database')
    
    #check if the utils sofware are functional in the current system
    check_install_parser = subparsers.add_parser('check-dependency', help='Check which OPERA-MS-UTILS software are functional in the current system')   
    
    args=parser.parse_args()    
    
    #print(args.checkm)#print(args.metabat2)
    main(args)
