#!/usr/bin/env python

import os
import sys
import errno
import argparse

def parse_args(args=None):
    Description = 'Fetch kmerfinder result files and get the most used reference.'
    Epilog = """Example usage: python find_common_reference.py -d <input directory> <output file>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('-d', help="Input directory.")
    parser.add_argument('-o', help="Output file.")
    return parser.parse_args(args)


def group_references(kmer_result_dir, out_file):
    '''
    '''
    reference_assembly = {}
    for k_file in os.listdir(kmer_result_dir):
        with open (os.path.join(kmer_result_dir, k_file), 'r') as fh:
            file_lines = fh.readlines()
        heading = file_lines[0].split('\t')
        first_line = file_lines[1].split('\t')
        index_assembly  = heading.index('# Assembly')
        reference = first_line[index_assembly]
        if reference not in reference_assembly:
            index_description = heading.index('Description')
            reference_assembly[reference] = [0, first_line[index_description]]
        reference_assembly[reference][0] += 1
    bigger_value = 0
    order_reference = dict(sorted(reference_assembly.items(), key = lambda x: x[1][0], reverse = True))
    with open (out_file, 'w') as f_out:
        for key, value in order_reference.items():
            f_out.write(key + '\t' + str(value[0]) + '\t' + value[1] + '\n')
    return

def main(args=None):
    args = parse_args(args)
    group_references(args.d,args.o)



if __name__ == '__main__':
    sys.exit(main())
