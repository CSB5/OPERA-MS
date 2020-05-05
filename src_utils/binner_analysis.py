import os
from genfunc import *

#util_dir = os.path.dirname(os.path.realpath(__file__)) + "/../tools_opera_ms/"

def run_maxbin2(contig_file, short_read1, assembly_dir, sample_name, nb_thread):
    out_dir = assembly_dir + "/maxbin2"
    create_dir(out_dir)
    out_dir = assembly_dir + "/maxbin2/all"
    create_dir(out_dir)
    run_exe("perl {}/maxbin/run_MaxBin.pl -min_contig_length 500 -contig {} -reads {} -out {}/{}_bin -thread {}".format(util_dir, contig_file, short_read1, out_dir, sample_name, nb_thread), True)


def run_metabat2(contig_file, short_read1, short_read2, assembly_dir, sample_name, nb_thread):

    out_dir = assembly_dir + "/metabat2"
    
    #Check if a mapping exists to the unpolished data
    bam_file = assembly_dir + "/intermediate_files/polished_assembly/contigs.fa.bam"
    
    #Run the mapping
    if not os.path.isfile(bam_file):
        # bwa ref index
        run_exe("{}/bwa index {} > {}/bwa_index.out 2> {}/bwa_index.err".format(util_dir, contig_file, assembly_dir, assembly_dir), True)

        # bwa mapping
        run_exe("{}/bwa mem -t {} {} {} {} | {}/samtools view -Sub - | {}/samtools sort - {} > {}/bwa.out 2>> {}/bwa.err".format(util_dir, nb_thread, contig_file, short_read1, short_read2, util_dir, util_dir, contig_file, assembly_dir, assembly_dir), True)

        # bwa bam index
        bam_file = contig_file + ".bam"
        run_exe("{}/samtools index {} >> {}/bwa.out 2>> {}/bwa.err".format(util_dir, bam_file, assembly_dir, assembly_dir), True)

    #Create the output directory
    create_dir(out_dir)
    out_dir = assembly_dir + "/metabat2/all"
    create_dir(out_dir)
    
    # metabat depth file
    run_exe("{}/jgi_summarize_bam_contig_depths --outputDepth {}/output_depth {}".format(util_dir, out_dir, bam_file), True)    
    # run metabat
    run_exe("{}/metabat --unbinned -i {} -a {}/output_depth -o {}/{}_bin -t {}".format(util_dir, contig_file, out_dir, out_dir, sample_name, nb_thread), True)



