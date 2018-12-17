'''

Perform the steps of the example specified in README.md.  Re-run this
workflow simply with:

   snakemake -p

from the shell

note: requires basic unix tools such as: sed, awk, zcat, wget, git,
python tools such as: python3, virtualenv,
and: samtools

'''

import pysam

#
# master rule
#----------------------------------------------------------------------
rule master :
	input :
		'measures/chr1__genotype_plot.pdf'


#
# obtain the VCF file
#----------------------------------------------------------------------
_ftp_ = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/' # parameters for VCF (customizable)
_vcf_ = _ftp_ + 'release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz'

# rename the chromosome field so that it maches that of the BAM file
rule rename_chr :
        output : 'chr1_.vcf'
	input : 'chr1.vcf'

	message : 'renaming chromosome field of {input} to obtain {output}'

	shell : '''

   sed '/^#/ s/##contig=<ID=/##contig=<ID=chr/' {input} | \
      sed '/^[^#]/ s/^/chr/' > {output} '''

# extract chromosome (1)
rule extract_chr :
        output : 'chr{chromosome,[0-9]}.vcf'
	input : 'chrs.vcf'

	message : 'extracting chromosome {wildcards.chromosome} from {input} to obtain {output}'

	shell : '''

   awk '/^#/ || ($1 == "{wildcards.chromosome}")' {input} > {output} '''

# unzip the gzipped VCF
rule unzip :
        output : 'chrs.vcf'
	input : 'chrs.vcf.gz'

	message : 'unzip {input} to obtain {output}'

	shell : 'zcat {input} > {output}'

# download the (gzipped) VCF
rule download_vcf :
	output : 'chrs.vcf.gz'

	message : 'download VCF: {_vcf_} into {output}'

	shell : 'wget -O {output} {_vcf_}'


#
# obtain BAM file
#----------------------------------------------------------------------
_bam_ = _ftp_ + 'data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam'

# index BAM file with samtools
rule index_bam :
        output : '{filename}.bam.bai'
	input : '{filename}.bam'

	message : 'index {input} with samtools to obtain {output}'

	shell : '''

   samtools index {input}
   touch {output} '''

# download the BAM file
rule download_bam :
        output : 'reads.bam'

	message : 'download BAM: {_bam_} into {output}'

	shell : 'wget -0 {output} {_bam_}'


#
# obtain the matrix file, using WhatsHap
#----------------------------------------------------------------------
_branch_ = 'matrix' # parameters for whatshap (customizable)
_virtualenv_ = 'venv'
_whatshap_ = 'whatshap/' + _virtualenv_ + '/bin/whatshap'
_flags_ = '--ignore-read-groups'
_output_ = 'output'

# run whatshap on the pair we have obtained above
#
# note: that this performs the same task as export-matrix.py, in that
# this is the same as calling:
#
#    python3 export-matrix.py --outdir output chr1_.vcf reads.bam
#
# from the shell
rule run_whatshap :
	output : _output_ + '/wh.raw.mat'
	input :
		prgm = _whatshap_,
		vcf = 'chr1_.vcf',
		bam = 'reads.bam',
		bai = 'reads.bam.bai',

	params :
		flags = _flags_,

	log :
		reads = _output_ + '/wh.raw.reads',
		log = _output_ + '/wh.raw.mat.log'

	message : 'run {input.prgm} on {input.vcf}, {input.bam} to obtain {output}'

	shell : '''

   {input.prgm} phase \
      {params.flags} --output-read-list {log.reads} \
         {input.vcf} {input.bam} \
            > {output} 2> {log.log} '''

# obtain and install whatshap
#
# note: that this performs the same task as setup.sh, in that this is
# the same as calling:
#
#    bash setup.sh
#
# from the shell

# install whatshap (into a virtual environment)
rule install_whatshap :
	output : 'whatshap/{venv}/bin/whatshap'
        input : 'whatshap/{venv}/bin/pytest'

	message : 'installing whatshap into {output} ...'

	shell : '''

   cd whatshap
   {wildcards.venv}/bin/pip3 install -e .[dev] '''

# create a virtual environment
rule create_virtualenv :
        output : 'whatshap/{venv}/bin/pytest'
	input : 'whatshap/README.md'
	message : 'creating virtual environment whatshap/{wildcards.venv}'

	shell : '''

   cd whatshap
   virtualenv -p python3 {wildcards.venv}
   {wildcards.venv}/bin/pip3 install networkx Cython nose pytest tox sphinx '''

# clone branch of tweaked version of WhatsHap
rule download_whatshap :
	output : 'whatshap/README.md'

	params :
		branch = _branch_

	message : 'downloading whatshap ...'

	shell : '''

   git clone \
      --branch {params.branch} \
      https://bitbucket.org/whatshap/whatshap '''


#
# compile the matrix data for quality measures computation
#----------------------------------------------------------------------
_compiler_ = 'export_data.py' # parameters for export_data.py
_compiled_ = 'compiled'

# compile matrix data with export_data.py
rule compile_matrix :
        output : _compiled_ + '/{vcf}.p'
	input :
        	prgm = _compiler_,
		matrix = 'output/wh.raw.mat',
		vcf = '{vcf}.vcf'

	message : 'compile {input.matrix}, {input.vcf} with {input.prgm} to obtain {output}'

	shell : '''

   python3 {input.prgm} --mat {input.matrix} --phasedvcf {input.vcf} -o {_compiled_} '''


#
# obtain quality measures from this compiled data
#----------------------------------------------------------------------
_measurer_ = 'quality.py' # parameters for quality.py
_measures_ = 'measures'
_errorrate_ = 0.15
_techcov_ = 40

# compute measures with quality.py
rule measure_quality :
	output : _measures_ + '/{vcf}_genotype_plot.pdf'
	input : 
		prgm = _measurer_,
		data = _compiled_ + '/{vcf}.p'

	params :
		epsilon = _errorrate_,
		techcov = _techcov_

	message : '''

   measure quality of {input.data} with {input.prgm} to produce {output},
   according to error rate {params.epsilon} and coverage {params.techcov} '''

	shell : '''

   python3 {input.prgm} --epsilon {params.epsilon} --techcov {params.techcov} \
      --data {input.data} --outdir {_measures_}
   touch {output} '''
