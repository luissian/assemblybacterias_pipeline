#!/usr/bin/env python
'''
if [ $# -eq 2 ]
    then
        assembly_summary_file='/processing_Data/bioinformatics/references/bacteria/latest_db/assembly_summary_bacteria.txt'
        bacterial_id=$1
        output_dir=$2

        nucleotide_end='_genomic.fna.gz'
        protein_end='_protein.faa.gz'
        gff_end='_genomic.gff.gz'

        ftp_path=$(cat $assembly_summary_file | grep "$bacterial_id" | cut -f20)
        ftp_path=$ftp_path/$bacterial_id

        cd $output_dir

        wget $ftp_path$nucleotide_end
        wget $ftp_path$protein_end
        wget $ftp_path$gff_end
        gunzip *.gz
    else
        echo "You must supply:"
        echo "As first argumen: Reference Genome RefSeq ID (ex.: GCF_000191485.1_ASM19148v1)"
        echo "As second argumen: Absolute path to copy reference (ex.: /processing_Data/bioinformatics/services_and_colaborations/CNM/bacteriologia/{service_name}/REFERENCES/)"
fi
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
        top_reference = infile[0]
    
    # create the outdir
    os.mkdir(out_dir)

    # open the reference and find the reference
    with open(reference) as infile:
        infile = infile.readlines()
        infile = infile.replace("\n","")split("\t")
        infile = [row for row in infile if row[] == top_reference]
    
    # get url
        url = infile[]  
    

    for r_end in reference_ends:
        
        out_file = os.file.join(out_dir,reference_file + r_end)
        file_url = url + '/' + reference_file + r_end
        
        wget.download(file_url, out_file)
    
    return

def main(args=None):
    args = parse_args(args)
    download_references(args.file, args.reference, args.out_dir)

if __name__ == '__main__':
    sys.exit(main())
