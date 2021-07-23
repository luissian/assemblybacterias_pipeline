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
import urllib.request

def parse_args(args=None):
    Description = 'download the reference files (fna, faa, gff) from the reference ncbi file.'
    Epilog = """Example usage: python download_reference.py -url <ftp path> -out_dir <output directory>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('-url', help="url path to download the query reference")
    parser.add_argument('-out_dir', help="Output directory.")
    return parser.parse_args(args)

def download_references (url, out_dir):
    reference_ends = ['_genomic.fna.gz','_protein.faa.gz', '_genomic.gff.gz']
    os.mkdir(out_dir)
    reference_file = url.split('/')[-1]
    for r_end in reference_ends:
        out_file = os.path.join(out_dir,reference_file + r_end)
        urllib.request.urlretrieve(url + '/' + reference_file + r_end , out_file)
        # urllib.request.urlretrieve("http://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/236/605/GCA_003236605.1_ASM323660v1/GCA_003236605.1_ASM323660v1_genomic.fna.gz", "GCA_003236605.1_genomic.fna.gz")
    return

def main(args=None):
    args = parse_args(args)
    download_references(args.url, args.out_dir)



if __name__ == '__main__':
    sys.exit(main())
