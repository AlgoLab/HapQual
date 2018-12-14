# HapQual
Computing measures of quality of a haplotype phasing with respect to sequencing reads

Here, we propose a way to asses the quality of the phasing of a set of
_single nucleotide variants_ (SNVs), as phased using a reads-based
phaser (or any method) with respect to a set of sequencing reads,
_e.g._, those reads which were used to obtain the phasing.  This is
described in a document in `documents` directory of this repository,
which we will hereby refer to as the _document_.  The codes in this
directory compute these measures according to this _document_.  The
best way to illustrate these codes is with the following full example


# Usage --- via an example

suppose we wish to study PacBio chromosome 1 of the individual
NA12878.


## obtain VCF file

download the (gzipped) VCF :

    wget -O chrs.vcf.gz ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.vcf.gz

unzip it :

    zcat chrs.vcf.gz > chrs.vcf

extract chromosome 1 :

    awk '/^#/ || ($1 == "1")' chrs.vcf > chr1.vcf

rename the chromosome field so that it matches that of the BAM file
that we will download in the next step :

    sed '/^#/ s/##contig=<ID=/##contig=<ID=chr/' chr1.vcf | sed '/^[^#]/ s/^/chr/' > chr1_.vcf


## obtain BAM file

download the BAM file :

    wget -O reads.bam ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NA12878_PacBio_MtSinai/sorted_final_merged.bam

index it with samtools :

    samtools index reads.bam

to create `reads.bam.bai`


## obtain matrix file, using WhatsHap

we tweaked WhatsHap to output a (type of) SNV/fragment matrix from
VCF/BAM pair.  First we need to set this up, download and install the
branch of WhatsHap

    bash setup.sh

now we run this on the pair we have obtained above

    python3 export-matrix.py --outdir output chr1_.vcf reads.bam

the output will be in `output/wh.raw.mat`


## compile the matrix data for quality measures computation

here we compile the matrix data, given a VCF file, so that we can
measure the quality according to different prior hypotheses (which can
be decided later).  Here we compile :

    python3 export_data.py --mat output/wh.raw.mat --phasedvcf chr1_.vcf -o compiled

the output (a binary file) will appear in `compiled/chr1_.p`.  From
this file we will obtain the quality measures of the phasing in the
next step


## obtain quality measures from this compiled data

here we compute quality measures (genotype, phasing, etc.) from this
compiled data according to various priors which may be changed.  This
is why we have a compiled file, so that we may quickly re-compute
these for (a variety of) other priors.  We compute the quality with

    python3 quality.py --data compiled/chr1_.p --epsilon 0.15 --outdir measures --techcov 40

from our compiled data in `compiled/chr1_.p`, which is PacBio data, so
we assume a prior `--epsilon` (error rate) of `0.15` (15% error rate).
The average coverage of this dataset is ~40, hence the `--techcov 40`.
All other parameters (used to compute the various priors, see the
_document_) have the default values.  The results will appear (as
plots) in the `measures` directory:

    chr1__genotype_plot.pdf

corresponds to the genotype quality measures with

* the distribution of the coverage over the SNV sites
* the _maximum a posteriori_ (MAP) distribution of theta
* the distribution of the p-value the MAP distribution (over a null
  model)

as detailed in section 1 of the _document_

    chr1__simulations_plot.pdf

corresponds to several simulations run to predict priors as detailed
in section 1 of the _document_

    chr1__phasing_plot.pdf

contains

* the distribution of reads which support or oppose the given phasing
  (or neither, which is coverage - supporting - opposing)
* and the log of the phasing quality

as detailed in section 2 of the _document_

    chr1__haplotype_plot.pdf

contains the distribution of the haplotype quality as detailed in
section 3 of the _document_
