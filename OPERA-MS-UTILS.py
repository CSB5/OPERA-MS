import sys
import os
import argparse
import subprocess
import config_parser

util_dir = "/home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/OPERA-MS-DEV/OPERA-MS/utils"
"""
qsub -terse -m a -M \$USER_PRINCIPAL_NAME -cwd -V -l mem_free=10G,h_rt=12:0:0 -pe OpenMP 16 -N CFSAN28_MAXBIN -b y perl /mnt/projects/bertrandd/opera_lg/SOFTWARE/MAXBIN/MaxBin-2.2.4/run_MaxBin.pl -min_contig_length 500 -contig ANALYSIS/ASSEMBLY/28h/contig.fasta -reads Illumina_reads/R1_28h_mergedrep.fastq.gz -out ANALYSIS/ASSEMBLY/24h/MAXBIN2/28h_bin -thread 16

/mnt/software/unstowable/miniconda3-4.6.14/envs/metabat2-v2.12.1/bin/metabat -i ANALYSIS/ASSEMBLY/40h/contig.fasta -o test_meta/40h/metabin
"""

# default to look for polish contig, if polished contig not exist, use unpolish.
# if user want to always use unpolish contig, specifiy flag --unploshed

def create_dir(new_dir):
    try:
        os.mkdir(new_dir)
    except:
        pass

    
def run_exe(cmd):
    run = 1
    print(cmd)
    if run:
        return subprocess.check_output(cmd, shell=True)

    
def run_binner(binner, sample_name, assembly_dir, short_read1, short_read2, nb_thread):
    
    contig_file = "{}/contigs.fasta".format(assembly_dir)
    
    if os.path.exists('{}/contigs.polished.fasta'.format(assembly_dir)):
        contig_file = "{}/contigs.polished.fasta".format(assembly_dir)

    if binner == "maxbin2":
        create_dir("assembly_dir/maxbin2")
        run_exe("perl {}/run_MaxBin.pl -min_contig_length 500 -contig {} -reads {} -out {}/maxbin2/{}_bin -thread {}".format(util_dir, contig_file, short_read1, assembly_dir, sample_name, nb_thread))
        
    if binner == "metabat2":
        run_metabat2(contig_file, short_read1, short_read2, assembly_dir, sample_name, nb_thread)
        

def run_metabat2(contig_file, short_read1, short_read2, assembly_dir, sample_name, nb_thread):

    out_dir = assembly_dir + "/metabat2"
    
    #Check if a mapping exists to the unpolished data
    bam_file = assembly_dir + "/intermediate_files/polished_assembly/contigs.fa.bam"

    #Run the mapping
    if not os.path.isfile(bam_file):
        # bwa ref index
        run_exe("{}/bwa index {} > {}/bwa_index.out 2> {}/bwa_index.err".format(util_dir, contig_file, assembly_dir, assembly_dir))

        # bwa mapping
        run_exe("{}/bwa mem -t {} {} {} {} | {}/samtools view -Sub - | {}/samtools sort - {} > {}/bwa.out 2>> {}/bwa.err".format(util_dir, nb_thread, contig_file, short_read1, short_read2, util_dir, util_dir, contig_file, assembly_dir, assembly_dir))

        # bwa bam index
        run_exe("{}/samtools index >> {}/bwa.out 2>> {}/bwa.err".format(util_dir, bam_file, assembly_dir, assembly_dir))


    #Create the output directory
    create_dir(out_dir)
    
    # metabat depth file
    run_exe("{}/jgi_summarize_bam_contig_depths --outputDepth {}/output_depth {}".format(util_dir, out_dir, bam_file))
    
    # run metabat
    run_exe("{}/metabat --unbinned -i {} -a {}/output_depth -o {}/{}_bin -t {}".format(util_dir, contig_file, out_dir, out_dir, sample_name, nb_thread))


def run_checkm(assembly_dir, binner, nb_thread):
    
    suffix = "fa"
    if binner == "maxbin2":
        suffix="fasta"
        
    bin_dir = assembly_dir + "/" + binner
    checkm_dir = bin_dir + "/checkm"
    create_dir(checkm_dir)

    eval_file = bin_dir + "/eval.dat"
    if not os.path.exists(eval_file):
        run_exe("{}/checkm lineage_wf -t {}  -x {} {}/ {}".format(util_dir, nb_thread, suffix, bin_dir, checkm_dir))
        run_exe("{}/checkm qa -o 2 {}/lineage.ms {} > {}".format(util_dir, checkm_dir, checkm_dir, eval_file))


def run_kraken2(assembly_dir, read1, read2, long_read, nb_thread, abundance_threshold):
    kraken_db = util_dir + "/../utils_db/minikraken2_v2_8GB"
    
    out_dir = assembly_dir + "/kraken2"
    create_dir(out_dir)
    #Short read
    out_file = out_dir + "/short_read.out"
    out_file_report = out_file + ".report"
    compress_format = "--gzip-compressed"
    if not os.path.exists(out_file_report):
        run_exe(util_dir + "/kraken2" + " --db " +  kraken_db + " --threads " + nb_thread + " --paired " + compress_format + " --output " + out_file + " --report " + out_file_report + " " + read1 + " " + read2)
    #
    #long read
    out_file = out_dir + "/long_read.out"
    out_file_report = out_file + ".report"
    compress_format = ""
    if not os.path.exists(out_file_report):
        run_exe(util_dir + "/kraken2" + " --db " +  kraken_db + " --threads " + nb_thread + " " + compress_format + " --output " + out_file + " --report " + out_file_report + " " + long_read)
    
    #Comapare the abundance profile
    compare_abundance_profile(out_dir, out_dir + "/short_read.out.report", out_dir + "/long_read.out.report", "S", abundance_threshold)
    compare_abundance_profile(out_dir, out_dir + "/short_read.out.report", out_dir + "/long_read.out.report", "G", abundance_threshold)
                
def compare_abundance_profile(out_dir, short_read_profile, long_read_profile, tax_level, abundance_threshold):
    tax_abundance_comparison = {}
    read_profile(short_read_profile, tax_level, abundance_threshold, tax_abundance_comparison, 0)
    read_profile(long_read_profile, tax_level, abundance_threshold, tax_abundance_comparison, 1)
    #
    OUT = open(out_dir + "/" + tax_level + "_abundance_comparison.txt", "w")
    OUT.write("Tax_name\tShort_read_abundance\tLong_read_abundance\n")
    for tax in tax_abundance_comparison:
        OUT.write("{}\t{}\t{}\n".format(tax, tax_abundance_comparison[tax][0], tax_abundance_comparison[tax][1]))

def read_profile(profile, tax_level, abundance_threshold, abundance_comparison, col):
    FILE = open(profile, "r")
    for line in FILE:
        #print line
        line_list = line.split("\t")
        tax = line_list[3]
        abundance = float(line_list[0].strip())
        tax_name = ((line_list[5].lstrip()).rstrip()).replace(" ", "_")
        if tax == tax_level and abundance > abundance_threshold:
            #print tax + " |" + abundance + "| |" + tax_name + "|\n"
            if tax_name not in abundance_comparison:
                abundance_comparison[tax_name] = {0:0, 1:0}
            abundance_comparison[tax_name][col] = abundance
    
def run_mash(all_binning_method, assembly_dir):
    mash_db = "/home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/OPERA-MS-DEV/OPERA-MS/genomeDB_Sketch.msh";
    mash_exe = "/home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/OPERA-MS-DEV/OPERA-MS/utils/mash"
    all_binning_method_array = all_binning_method.split(",")
    for binning_method in all_binning_method_array:
        binning_dir = assembly_dir + "/" + binning_method
        cmd = mash_exe +  " dist -p 1 -d 0.2 " + mash_db + " " + binning_dir + "/*fa* > " + binning_dir + "/mash.dist"
        run_exe(cmd)
        
def download_utils_db():
    mash_db = "/home/bertrandd/PROJECT_LINK/OPERA_LG/META_GENOMIC_HYBRID_ASSEMBLY/OPERA-MS-DEV/OPERA-MS/genomeDB_Sketch.msh";
    kraken_db = "/mnt/genomeDB/misc/softwareDB/kraken2/standard-20190108";


    
def main(args):

    config_file = args.config
    command = args.command

    #Parse the config file
    config_dict = {}
    with open(config_file, "r") as fp:
        for line in fp:
            line = line.split("#")[0]
            line = line.split()
            # if empty line
            if not line:
                continue
            config_dict[line[0]] = " ".join(line[1:])

    #Set the number of thread
    nb_thread = config_dict["NUM_PROCESSOR"]
    if args.thread != None:
        nb_thread = args.thread
            
    if command == "maxbin2" or command == "metabat2":
        sample_name = args.sample_name
        if sample_name == None:
            sample_name = config_dict["OUTPUT_DIR"].split("/")[-1]
        run_binner(command, sample_name, config_dict["OUTPUT_DIR"], config_dict["ILLUMINA_READ_1"], config_dict["ILLUMINA_READ_2"], nb_thread)

    elif command == "checkm":
        run_checkm(config_dict["OUTPUT_DIR"], args.binner, nb_thread)

    elif command == "kraken2":
        abundance_threshold = 0.1
        if args.abundance_threshold != None:
            abundance_threshold = args.abundance_threshold
        run_kraken2(config_dict["OUTPUT_DIR"], config_dict["ILLUMINA_READ_1"], config_dict["ILLUMINA_READ_2"], config_dict["LONG_READ"], nb_thread, float(abundance_threshold))

    elif command == "mash":
        run_mash(mash, config_dict["OUTPUT_DIR"])

    elif command == "novel-species":
        print("TO DO")

    elif command == "opera-ms-db":
        print("TO DO")
        
    elif command == "utils-db":
        print("TO DO")

    elif command == "check-install":
        print("TO DO")
        
if __name__ == "__main__":   

    parser = argparse.ArgumentParser()
    
    #group = parser.add_mutually_exclusive_group()
    #The type of software
    subparsers = parser.add_subparsers(help='commands', dest='command')
    #maxbin2
    maxbin_parser = subparsers.add_parser('maxbin2', parents=[config_parser.parser], help='Run MaxBin2')
    maxbin_parser.add_argument("-s", "--sample-name", help="Sample name [default output directory]")
    
    #metabat2
    metabat_parser = subparsers.add_parser('metabat2', parents=[config_parser.parser], help='Run MetaBAT2')
    metabat_parser.add_argument("-s", "--sample-name", help="Sample name [default output directory]")
    
    #kraken
    kraken_parser = subparsers.add_parser('kraken2', parents=[config_parser.parser], help='Run Kraken2 on the illumina reads and nanopore reads and compare the abundance profiles')
    kraken_parser.add_argument("-a", "--abundance-threshold", help="Lower percentage abundance threshold [default 0.1]")
    
    #chechm
    checkm_parser = subparsers.add_parser('checkm', parents=[config_parser.parser], help='Run CheckM on a set of bins')
    checkm_parser._action_groups[-1].add_argument("-b", "--binner",  required=True, choices=["maxbin2", "metabat2", "opera_ms_cluster"])
        
    #mash
    mash_parser = subparsers.add_parser('mash', parents=[config_parser.parser], help='Run Mash')
    
    #novel species
    novel_species_parser = subparsers.add_parser('novel-species', parents=[config_parser.parser], help='Run novel species identification')
    
    #opera-db
    opera_db_parser = subparsers.add_parser('opera-ms-db', help='Create a OPERA-MS genome database')
    
    #utils-db
    utils_db_parser = subparsers.add_parser('utils_db', help='Download all the data base required by the utils software')

    check_install_parser = subparsers.add_parser('check_install', help='Check which OPERA-MS-UTILS software are functional in the current system')
    
    args=parser.parse_args()
    print(args)
    
    #print(args.checkm)#print(args.metabat2)
    main(args)
