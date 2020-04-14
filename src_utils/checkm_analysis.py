from genfunc import *
import re


def checkm_analysis(assembly_dir, binner, nb_thread, high_qual_mags, medium_qual_mags):

    #Identify the correct suffix for check analysis
    suffix = get_binner_file_type(binner)
   
    #Identify the bins directory
    binner_dir = assembly_dir + "/" + binner
    bin_dir = binner_dir + "/all"

    if os.path.exists(bin_dir):
        #Create the required directories
        checkm_dir = binner_dir + "/checkm"
        create_dir(checkm_dir)
        
        #Run the checkm analysis
        eval_file = binner_dir + "/eval.dat"
        run_checkm(bin_dir, checkm_dir, suffix, nb_thread, eval_file)
    
        #Get the high and medium quality mags directories and file description
        identify_high_medium_bin(binner_dir, eval_file, high_qual_mags, medium_qual_mags, suffix)
    else:
        print("Checkm analysis can not performed. Please run OPERA-MS-UTILS.py {} first.\n".format(binner))

        
def identify_high_medium_bin(binner_dir, eval_file, high_qual_mags, medium_qual_mags, suffix):
    high_qual_threshold = [float(item) for item in high_qual_mags.split(',')]
    medium_qual_threshold = [float(item) for item in medium_qual_mags.split(',')]
    #
    high_dir = binner_dir + "/high_quality"
    medium_dir = binner_dir + "/medium_quality"

    try:
        run_exe("rm -r " + high_dir + " " + medium_dir, True)
    except:
        pass
    
    create_dir(high_dir)
    create_dir(medium_dir)
    #
    print(high_qual_threshold)
    FILE = open(eval_file, "r")
    OUT = open(binner_dir + "/bin_info.txt", "w")
    bin_id = bin_status = bin_out_dir = ""
    bin_size, bin_nb_contig, bin_longest_contig, bin_n50, bin_contamination, bin_completness = 0, 0, 0, 0, 0, 0
    OUT.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("ID", "Status", "Completness", "Contamination",  "Size", "#_contigs", "Contig_N50", "Longest_contig"))
    for line in FILE:
        #print line
        line_list = re.split(" \s+", line)
        if len(line_list) != 1 and line_list[1] != "Bin Id":
            bin_id = line_list[1]
            bin_status = "LOW"
            bin_completness = float(line_list[6])
            bin_contamination = float(line_list[7])

            if bin_completness >= high_qual_threshold[0] and bin_contamination <= high_qual_threshold[1]:
                bin_status = "HIGH"
                bin_out_dir = high_dir
            elif bin_completness >= medium_qual_threshold[0] and bin_contamination <= medium_qual_threshold[1]:
                bin_status = "MEDIUM"
                bin_out_dir = medium_dir
            #
            bin_size = line_list[9]
            bin_nb_contig = line_list[11]
            bin_n50 = line_list[14]
            bin_longest_contig = line_list[18]
            #Get the stats
            OUT.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(bin_id, bin_status, bin_completness, bin_contamination,  bin_size, bin_nb_contig, bin_n50, bin_longest_contig))

            if bin_status == "HIGH" or bin_status == "MEDIUM":
                run_exe("ln -s ../all/{}.{} {}/{}.{}".format(bin_id, suffix, bin_out_dir, bin_id, suffix), True)
            
    FILE.close()
    OUT.close()
  
        
def run_checkm(bin_dir, checkm_dir, suffix, nb_thread, eval_file):
    #Run the checkm analysis
    if not os.path.exists(eval_file):
        run_exe("{}/checkm lineage_wf -t {}  -x {} {}/ {}".format(util_dir, nb_thread, suffix, bin_dir, checkm_dir), True)
        run_exe("{}/checkm qa -o 2 {}/lineage.ms {} > {}".format(util_dir, checkm_dir, checkm_dir, eval_file), True)
    else:
        print("CheckM result detected skip analysis")
