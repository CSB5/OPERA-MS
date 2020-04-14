import os
import subprocess

util_dir = os.path.dirname(os.path.realpath(__file__)) + "/../tools_opera_ms_utils/"

def create_dir(new_dir):
    try:
        os.mkdir(new_dir)
    except:
        pass
    

def get_binner_file_type(binner):
    file_type = "fasta"
    if binner == "metabat2":
        file_type ="fa"
    return file_type


def run_exe(cmd, print_out):
    run = 1
    if print_out:
        print(cmd)
        
    if run:
        return subprocess.check_output(cmd, shell=True)
    

def read_opera_ms_config_file(config_file):
    config_dict = {}
    with open(config_file, "r") as fp:
        for line in fp:
            line = line.split("#")[0]
            line = line.split()
            # if empty line
            if not line:
                continue
            config_dict[line[0]] = " ".join(line[1:])
    return config_dict


def get_contig_file(assembly_dir):
    contig_file = "{}/contigs.fasta".format(assembly_dir)
    if os.path.exists('{}/contigs.polished.fasta'.format(assembly_dir)):
        contig_file = "{}/contigs.polished.fasta".format(assembly_dir)
    return contig_file
