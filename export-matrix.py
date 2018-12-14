#!/usr/bin/env python3

description = '''

   export a SNP/fragment matrix from a BAM/VCF pair, to be used
   downstream by export_data.py --- requires setup.py

'''

import sys
import os, subprocess
import argparse

#
# preamble
#----------------------------------------------------------------------

dir = os.path.dirname(os.path.realpath(__file__))
whatshap = '{}/whatshap/venv/bin/whatshap'.format(dir)

#
# functions
#----------------------------------------------------------------------

# auxiliary function to add arguments to the argument parser
def add_arguments(parser) :

    arg = parser.add_argument

    # positional arguments
    arg('vcf', metavar = 'VCF',
        help = 'VCF file with the phased variants')
    arg('bam', metavar = 'BAM',
        help = 'BAM file of sequencing reads from which the phasing is based')

    # output directory
    arg('--outdir', '-o', metavar = 'OUTPUT', action = 'store',
	type = str, required = True,
        help = 'output VCF file with phase information')

    # realignment
    arg('--reference', '-r', metavar = 'FASTA',
        help = 'reference file, for detecting alleles in realignment mode')


# use whatshap to read in snp/fragment matrix from a bam file
def read_bam(vcf, bam, reference, out_dir) :

    # setup realignment mode
    realignment = ''
    if reference :
        realignment = '--reference {}'.format(reference)

    # run whatshap
    flags = '--ignore-read-groups' # should always have in this case
    rawreal = 'realigned' if realignment else 'raw'
    reads = '{}/wh.{}.reads'.format(out_dir, rawreal)
    matrix = '{}/wh.{}.mat'.format(out_dir, rawreal)
    log = '{}.log'.format(matrix)
    subprocess.run('''

  {} phase {} {} --output-read-list {} {} {}

    '''.format(whatshap, flags, realignment, reads, vcf, bam).split(),
                   stdout = open(matrix,'w'),
                   stderr = open(log,'w'))

    # return matrix file (location)
    return matrix


#
# main
#----------------------------------------------------------------------
def main(argv = sys.argv[1:]) :

    # parse arguments
    parser = argparse.ArgumentParser(description = 'export-matrix', add_help = True)
    add_arguments(parser)
    args = parser.parse_args(argv)

    # assert the existence of whatshap and hapchat core phasing algorithm
    assert os.path.exists(whatshap), 'WhatsHap not found, please run setup.sh'

    # set up dir for intermediate output
    out_dir = args.outdir
    subprocess.run('mkdir -p {}'.format(out_dir).split())

    # read vcf/bam to obtain the matrix
    matrix = read_bam(args.vcf, args.bam, args.reference, out_dir)
    print('matrix file can be found in:', matrix, file = sys.stderr)


if __name__ == '__main__' :
    main()
