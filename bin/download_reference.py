#!/usr/bin/env python
'''
=============================================================
HEADER
=============================================================
INSTITUTION: BU-ISCIII
AUTHOR: Guillermo J. Gorines Cordero
MAIL: guillermo.gorines@urjc.es
VERSION: 0.1
CREATED: Early 2022
REVISED: 18-2-2022
DESCRIPTION: 

INPUT:
    -FILE: file containing the ranking of references from kmerfinder.
           created by the script find_common_references
    -REFERENCE: file with the NCBI reference list
    -OUTDIR: name of the output dir

OUTPUT:
OPTIONS:
USAGE:
    python download_reference.py -file [FILE] -reference [REFERENCE] -out_dir [OUTDIR]
REQUIREMENTS:
    -Python >= 3.6
DISCLAIMER:
TO DO: 
================================================================
END_OF_HEADER
================================================================
'''

import sys
import argparse
import os

import wget

def parse_args(args=None):
    Description = 'download the reference files (fna, faa, gff) from the reference NCBI file.'
    Epilog = """Usage example: python download_reference.py -file <file with the references created by find_common_reference> -reference <file from the NCBI with all bacterial references> -out_dir <output directory>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('-file', help="File containing the ranking of references from kmerfinder.")
    parser.add_argument('-reference', help="File containing the paths to bacterial references.")
    parser.add_argument('-out_dir', help="Output directory.")
    
    return parser.parse_args(args)

def download_references (file, reference, out_dir):
    '''
    Downloads the top reference from the NCBI database 
    '''

    reference_ends = ['_genomic.fna.gz','_protein.faa.gz', '_genomic.gff.gz']
    
    # extract the most common reference from file
    with open(file) as infile:
        infile = infile.readlines()
        infile = [item.replace("\n","").split("\t") for item in infile]
        top_reference = infile[0][0]

    # create the outdir (do nothing if already there)
    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass

    # open the reference and find the reference
    with open(reference) as infile:
        infile = infile.readlines()
        infile = [item.replace("\n","").split("\t") for item in infile]
        infile = [row for row in infile if row[0] in top_reference]
    

    # get url and reference file
    url = infile[0][19]  

    for r_end in reference_ends:
        
        out_file = out_dir + "/" + top_reference + r_end
        file_url = url + '/' + top_reference + r_end
        
        print(out_file)
        print(file_url)

        wget.download(file_url, out_file)
    
    return

def main(args=None):
    args = parse_args(args)
    download_references(args.file, args.reference, args.out_dir)

if __name__ == '__main__':
    sys.exit(main())
