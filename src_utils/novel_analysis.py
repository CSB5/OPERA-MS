from genfunc import *
import os


def is_good_quality(bin_status, mag_selected_qual):
    res = False 
    if bin_status.lower() == "high" or bin_status.lower() == mag_selected_qual:#If high quality or medium and mag_qual is medium the bin is used
        res = True
    return res


def get_sample_mags(conf, binner, mag_selected_qual, mag_out, mag_info):
    config_dict = read_opera_ms_config_file(conf)
    binner_dir = config_dict["OUTPUT_DIR"] + "/" + binner

    try :
        print(binner_dir + "/bin_info.txt")
        FILE = open(binner_dir + "/bin_info.txt", "r")
       
        for line in FILE:
            info = (line.rstrip('\n')).split("\t")
            bin_id = info[0]
            bin_status = info[1]
            if is_good_quality(bin_status, mag_selected_qual):
                mag_info[bin_id] = {"INFO":list(info),
                                    "BEST_HIT_SPECIES_TAX":"NA",
                                    "BEST_HIT_GENOME_TAX":"NA",
                                    "BEST_HIT_DIST_TAX":1,
                                    #
                                    "BEST_HIT_SPECIES_NOVEL":"NA",
                                    "BEST_HIT_GENOME_NOVEL":"NA",
                                    "BEST_HIT_DIST_NOVEL":1,
                                    #
                                    "STATUS":"KNOWN",
                                    "NOVEL_SPECIES_ID":"NA"}
                mag_link = mag_out + "/" + bin_id + ".fasta"
                if not os.path.exists(mag_link):
                    
                    run_exe("ln -s " +
                            binner_dir + "/" + bin_status.lower() + "_quality/" + bin_id + "." + get_binner_file_type(binner) +
                            " " +
                            mag_link, True)
                    
    except Exception as e:
        print(e)
        print("WARNING binning not completed in {}. Sample excluded from analysis\n".format(binner_dir))
        pass

    
def get_mag_id(mag_path):
    return ((mag_path.split("/"))[-1]).replace(".fasta", "")


def get_ref_genome_name(ref_genome_path):
    return ((ref_genome_path.split("/"))[-1]).split("__")


def update_best_hit_status(mag_info, file_info, ana_type):
    #print(file_info + "\t" + ana_type)
    FILE = open(file_info, "r")
    for line in FILE:
        (ref_genome_path, mag_path, distance, pvalue, kmer_info) = (line.rstrip('\n')).split('\t')
        mag_id = get_mag_id(mag_path)
        #
        if ana_type == "TAX":
            ref_genome_name = get_ref_genome_name(ref_genome_path)
            mag_info[mag_id]["BEST_HIT_SPECIES" + "_" + ana_type] = ref_genome_name[0]
            mag_info[mag_id]["BEST_HIT_GENOME" + "_" + ana_type] = ref_genome_name[1]
        else:
            mag_info[mag_id]["BEST_HIT_GENOME" + "_" + ana_type] = (ref_genome_path.split("/"))[-1]
        mag_info[mag_id]["BEST_HIT_DIST" + "_" + ana_type] = float(distance)
        

def update_novel_status(mag_info, file_info):
    FILE = open(file_info, "r")
    FILE.readline() #skip the header
    for line in FILE:
        (mag_id, cluster_id) = (line.rstrip('\n')).split('\t')
        mag_id = mag_id
        if mag_info[mag_id]["BEST_HIT_DIST_NOVEL"] > 0.05 and mag_info[mag_id]["BEST_HIT_DIST_TAX"] > 0.05:#Only bin for which both the best hit in the taxonomy file and in the genome file are novel are considered novel
            mag_info[mag_id]["STATUS"] = "NOVEL"
            mag_info[mag_id]["NOVEL_SPECIES_ID"] = "NS_" + cluster_id
            
        else:
            # if either analysis is not present in the db, give NO_REF as its reference
            if mag_info[mag_id]["BEST_HIT_DIST_NOVEL"] == 1.0:
                mag_info[mag_id]["BEST_HIT_GENOME_NOVEL"] = "NO_REF"
            elif mag_info[mag_id]["BEST_HIT_DIST_TAX"] == 1.0:
                mag_info[mag_id]["BEST_HIT_GENOME_TAX"] = "NO_REF"


def run_novel_species_analysis(ref_known_species, ref_taxonomy_species, nb_thread, configs, binner, mag_qual, out):
    create_dir(out)
    mag_out=out + "/MAGs"
    #Get the required MAGS for each samples creat a soft link in out_mags and collect the MAGs info
    create_dir(mag_out)
    mag_info = {}
    print("***Parsing all config files\n")
    for conf in configs:
        get_sample_mags(conf, binner, mag_qual, mag_out, mag_info)
    
    #Run the novel species analysis
    out_analysis_novel = out+"/novel_species_analysis"

    print("\n***running novel species analysis...\n")
    
    if not os.path.exists(out_analysis_novel+"/novel/cluster.tsv"):
        os.system("python {}/../src_utils/novel_species_identification.py -t {} -s 10000 --ref {} --query {} --file_type {} --outdir {}".format(util_dir, nb_thread, ref_known_species, mag_out, "fasta", out_analysis_novel))

    out_analysis_taxonomy = out+"/species_taxonomy_analysis"
    print("\n***running taxonomy analysis...\n")
    
    if not os.path.exists(out_analysis_taxonomy+"/novel/cluster.tsv"):
        os.system("python {}/../src_utils/novel_species_identification.py -t {} -s 1000 --ref {} --query {} --file_type {} --outdir {}".format(util_dir, nb_thread, ref_taxonomy_species, mag_out, "fasta", out_analysis_taxonomy))
        
    #Update the statistics about the analysis
    print("\n***update statistics for taxonomy analysis...\n")
    update_best_hit_status(mag_info, out_analysis_taxonomy+"/query_ref/min_hit.dat", "TAX")
    #
    print("\n***update statistics for novel analysis...\n")
    update_best_hit_status(mag_info, out_analysis_novel+"/query_ref/min_hit.dat", "NOVEL")
    print("\n***update novel status...\n")
    update_novel_status(mag_info, out_analysis_novel+"/novel/cluster.tsv")

    #Print the final info file
    OUT = open(out+"/MAGs_info.txt", "w")
    print("\n***writing output file...\n")
    OUT.write("MAG_ID\tQuality\tCompleteness\tContamination\tSize\t#_contigs\tContig_N50\tLongest_contig\tClosest_GTDB_species\tClosest_GTDB_genome\tDistance_to_closest_GTDB_genome\tClosest_MAG\tDistance_to_closest_MAG\tNovel_status\tNovel_species_ID\n")
    for mag in sorted(mag_info, key=lambda x: float(mag_info[x]["INFO"][2]), reverse=True):
        OUT.write(
            "\t".join(mag_info[mag]["INFO"]) + "\t" +
            #
            mag_info[mag]["BEST_HIT_SPECIES_TAX"] + "\t" +
            mag_info[mag]["BEST_HIT_GENOME_TAX"] + "\t" +
            str(mag_info[mag]["BEST_HIT_DIST_TAX"]) + "\t" +
            #
            mag_info[mag]["BEST_HIT_GENOME_NOVEL"] + "\t" +
            str(mag_info[mag]["BEST_HIT_DIST_NOVEL"]) + "\t" +
            #
            mag_info[mag]["STATUS"] + "\t" +
            mag_info[mag]["NOVEL_SPECIES_ID"] + "\n")
                       
