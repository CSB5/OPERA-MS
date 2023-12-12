# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-12-12 10:15:51
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-12-12 16:29:35

import concurrent.futures
import argparse
import gzip
import os, shutil

try:
    import tqdm
    USED_TQDM = True
except ImportError:
    USED_TQDM = False

def process(args):
    print ('Read taxonomy file')
    with gzip.open(args.taxonomy, 'rt') as f:
        taxinfo = dict(
            (line.strip().split('\t')) 
            for line in f
        )

    # Removal of the GB and RS left part
    taxinfo = {seqID[3:]: seqtax for seqID, seqtax in taxinfo.items()}

    print ('Check taxonomic information')
    ignored = set()
    for seqID, seqtax in taxinfo.items():
        taxons = seqtax.split(';')
        species = [taxon[3:] for taxon in taxons if taxon[0] == 's']

        if not species:
            print (f'Unable to find a species name for seqID {seqID}, sequence will be ignore')
            print (f'Taxon info for this seqID: {seqtax}')
            ignored.add(seqID)

        if len(species) > 1:
            raise Exception('Found more than one species for seqID {seqID}.\
                Taxon info for this seqID: {seqtax}')

        taxinfo[seqID] = species[0].replace(' ', '_')

    print ('Read genome file')
    with gzip.open(args.genomes_reps, 'rt') as f:
        fnames = sorted(
            line.strip() for line in f
            )

    print ('Check concordance')
    for fname in fnames:
        # removal of `_genomic.fna.gz`
        bname = id_from_path(fname)
        if bname not in taxinfo:
            raise Exception(f'Unable to find seqID {bname} in taxonomic file')

    print (f'Define genome size and seq numbers. Number of threads: {args.threads}')
    seqinfos = multi_threads_seqinfos(fnames) if args.threads > 1 else single_thread_seqinfos(fnames)

    print ('Create the output files')
    os.makedirs(args.outdir, exist_ok=True)

    local_paths = {}
    fun = shutil.move if args.move else os.symlink

    for idx, fname in tqdm.tqdm(enumerate(fnames)):
        # I keep the exact same formating than what has been done previously
        # I don't want to break any weird basename parsing in OPERA-MS
        gdir = 'genomes_' + str((idx // 4000) + 1)
        mgdir = os.path.join(args.outdir, gdir)
        if not os.path.isdir(mgdir): os.makedirs(mgdir)

        gid = id_from_path(fname)
        basename = taxinfo[gid] + '__' + gid + '_genomic.fna.gz'
        gpath = os.path.join(mgdir, basename)
        
        local_paths[fname] = os.path.join(gdir, basename)
        if os.path.isfile(gpath): os.remove(gpath)
        fun(os.path.realpath(fname), gpath)

    # The genome length is also used by OPERA-MS
    # And need a path relative to the opera-ms executable (I assume)
    outfile = os.path.join(args.outdir, 'genomes_length.txt')
    with open(outfile, 'w') as f:
        for fname, nfname in local_paths.items():
            gid = id_from_path(fname)
            count, size = seqinfos[gid]
            nfname = os.path.join('OPERA-MS-DB', nfname)
            f.write(f'{nfname}\t{count}\t{size}\n')

def single_thread_seqinfos(fnames):
    seqinfos = {}
    
    if USED_TQDM:
        iterator = tqdm.tqdm(fnames)
    else:
        iterator = fnames

    for fname in iterator:
        seqres = fasta_info(fname)
        seqinfos.update(seqres)

    return seqinfos

def multi_threads_seqinfos(fnames):
    seqinfos = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        if USED_TQDM:
            iterator = tqdm.tqdm(executor.map(fasta_info, fnames), total=len(fnames))
        else:
            iterator = executor.map(fasta_info, fnames)

        for seqres in iterator:
            seqinfos.update(seqres)

    return seqinfos

def id_from_path(fname):
    return os.path.basename(fname)[:-15]

def fasta_info(fname):
    count = size = 0
    with gzip.open(fname, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
                continue
            else:
                size += len(line.strip())   

    gid = id_from_path(fname)
    return {gid: (count, size)}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('genomes_reps', help='list of genomes', type=str)
    parser.add_argument('taxonomy', help='taxonomy file (must be the concatenation of arc and bac)', type=str)
    parser.add_argument('--outdir', help='output directory', type=str, default='operams_db')
    parser.add_argument('--threads', help='number of threads', type=int, default=1)
    parser.add_argument('--move', help='move instead of making symlink', action="store_true")

    args = parser.parse_args()
    process(args)

if __name__ == '__main__':
    main()


    