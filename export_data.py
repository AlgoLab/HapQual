#!/usr/bin/env python3
# coding: utf-8
"""
Copyright (C) 2018 Simone Ciccolella

Distributed under the MIT license.

You should have received a copy of the MIT license along with this program.
"""

import argparse, os, sys, errno
import re, pickle
import logging

def generate_argparser():
    parser = argparse.ArgumentParser(description='hapqual v2 - propotype', add_help=True)
    parser.add_argument('--mat', action='store', type=str, required=True,
                        help='path of the mat file.')
    parser.add_argument('--phasedvcf', action='store', type=str, required=True,
                        help='path of whatshap modified output')
    parser.add_argument('-o', '--outdir', action='store', type=str, required=True,
                        help='output directory.')
    return parser

'''
The function export 4 dict from the input file, structured as follows:
['position_genotype', 'pblocks', 'pos_read', 'pos_read_nucl']

position_genotype := {
    POSITION: {
        GENOTYPE: 'N1/N2'
        READS: [
            [NUCL, GEN],
            ...
        ]
    },
    ....
}

pblocks := {
    PHASING_BLOCK : [
        [POSITION, PHASING],
        ....
    ],
    .....
}

pos_read := {
        POSITION : set(READ1, READ2, ...),
        .....
}

pos_read_nucl := {
    [POSITION, READ] : [NUCL, GEN],
    ....
}
'''

def export_phased_only(phased_vcf, mat_file):
    snp_positions = dict()
    from collections import defaultdict

    phasing_blocks = defaultdict(list)

    logging.info('Reading VCF file...')
    # TODO: this should be modified to be more general
    pattern_vcf = re.compile(r'^(?P<chrom>.+?)\s+(?P<pos>\d+)\s+(?P<id>.+?)\s+'
                        r'(?P<ref>\w+)\s+(?P<alt>\w+)\s+(?P<qual>.+?)\s+'
                        r'(?P<filter>.+?)\s+(?P<info>.+?)\s+(?P<format>.+?)\s+'
                        # genotype for phased only
                        r'(?P<phasing>\d\|\d):.+:(?P<pblock>.+?)$'
                        )
    with open(phased_vcf, 'r') as fin:
        for line in fin:
            # print(line.split()[-1])
            match = pattern_vcf.match(line)
            if match:
                snp_positions[int(match.group('pos'))] = '%s/%s' % (match.group('ref'), match.group('alt'))
                phasing_blocks[match.group('pblock')].append([int(match.group('pos')), match.group('phasing')])
                # if match.group('pblock') != 'PATMAT':
                #     print(match.group('pblock'))
                #     print(line.split()[-1])
                #     print(line)


    logging.info('Parsed %d SNP positions and found %d phased blocks' % (len(snp_positions), len(phasing_blocks)))

    pattern_mat = re.compile(r'^(?P<read>.+?)\s+(?P<pos>\d+)\s+'
                            r'(?P<allele>\w+?)\s+(?P<genotype>\d)\s+'
                            r'(?P<phred>\d+).*$')
    
    reads = defaultdict(list)
    positions_to_reads = defaultdict(set)
    repos_to_nucl = dict()
    logging.info('Reading MAT file...')
    read_count = 0
    with open(mat_file, 'r') as fin:
        for line in fin:
            match = pattern_mat.search(line)
            if match:
                # TODO: change this when mat creator will be 1-based
                pos = int(match.group('pos')) + 1
                
                if pos in snp_positions:
                    reads[pos].append([match.group('allele'), match.group('genotype')])
                    read_count += 1

                    positions_to_reads[pos].add(match.group('read'))
                    repos_to_nucl[pos, match.group('read')] = [match.group('allele'), match.group('genotype')]

    logging.info('Parsed %d reads at the phased loci.' % read_count)

    # print(positions_to_reads)

    json_export = dict()

    logging.info('Creating the genotype export...')
    json_genotype = dict()
    for pos in snp_positions:
        json_genotype[pos] = {
            'Genotype': snp_positions[pos],
            'Reads': reads[pos]
        }

    logging.info('Sorting blocks...')
    for block in phasing_blocks:
        phasing_blocks[block] = sorted(phasing_blocks[block], key=lambda x: x[0])

    # from collections import namedtuple
    # ExportData = namedtuple('Data', ['position_genotype', 'pblocks', 'pos_read', 'pos_read_nucl'])
    edata = [json_genotype, phasing_blocks, positions_to_reads, repos_to_nucl]

    return edata


def main(args):
    export = export_phased_only(args.phasedvcf, args.mat)

    out_file = os.path.splitext(os.path.basename(args.phasedvcf))[0]
    pickle.dump(export, open('%s/%s.p' % (args.outdir, out_file), 'wb+'))
    
    logging.info('Done.')


def setup(args):
    try:
        os.makedirs(args.outdir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(args.outdir):
            pass
        else:
            raise

    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(message)s',
                        datefmt='%d-%b-%y %H:%M:%S')
    
    main(args)

if __name__ == '__main__':
    parser = generate_argparser()
    args = parser.parse_args()
    setup(args)
