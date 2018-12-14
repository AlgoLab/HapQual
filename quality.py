#!/usr/bin/env python3
# coding: utf-8

# Ignore seaborn future warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import argparse
import numpy as np
import json
import scipy.stats as stats
import pickle
import os, sys, errno
import logging
import itertools

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return zip(a, b)

def generate_argparser():
    parser = argparse.ArgumentParser(description='hapqual v2 - propotype', add_help=True)
    parser.add_argument('-d', '--data', action='store', type=str, required=True,
                        help='path of the exported data.')
    parser.add_argument('-e', '--epsilon', action='store', type=float, required=True,
                        help='technology error')
    parser.add_argument('-o', '--outdir', action='store', type=str, required=True,
                        help='output directory.')
    parser.add_argument('--techcov', action='store', type=int, required=True,
                        help='technology average coverage')

    # parser.add_argument('--techcovstd', action='store', type=int, default=2,
    #                     help='technology coverage standard deviation [Default=2]')
    parser.add_argument('--priorcov', action='store', type=int, default=10,
                        help='coverage to use to obtain prior [Default=10]')
    parser.add_argument('--nullsim', action='store', type=int, default=1000,
                        help='number of simulations to calculate H0 std-dev [Default=1000]')
    # parser.add_argument('--disable-dumps', action='store_true', required=False,
    #                     help='disable the use of dumps')
    parser.add_argument('--disable-plots', action='store_true', required=False,
                        help='disable the drawing of plots')
    parser.add_argument('--force-0cov', action='store_true', required=False,
                        help='force the calculation on 0-coverage SNPs')
    return parser

def import_data(path):
    data = pickle.load(open(path, 'rb'))
    return data

def simulate_h0_stdv(error, tot_simulations, tech_cov_mean, tech_cov_stdv, basename):
    stdvs = list()
    coverages = list()
    for i in range(tot_simulations):
        sample_coverage = np.random.normal(tech_cov_mean, tech_cov_stdv)
        stdvs.append(np.sqrt(sample_coverage*(1-error)*error))
        coverages.append(sample_coverage)

    if PLOTS:
        f, ax = plt.subplots(2, 1, figsize=(4,4))
        sns.distplot(coverages, ax=ax[0])
        ax[0].set(title="Coverage simulations")
        
        sns.distplot(stdvs, ax=ax[1], color='orange')
        ax[1].set(title="std-dev simulations")
        plt.tight_layout()

        f.savefig('%s_simulations_plot.pdf' % basename)
    
    simulated_stdv = stats.mode(stdvs)[0][0]
    logging.info('Estimation of H0 std-dev for %d simulations: %.5f' %
            (tot_simulations, simulated_stdv))
    return simulated_stdv

def simulate_priors(sim_cov, error):
    beta_correct = int(sim_cov * error)
    beta_error = sim_cov - beta_correct

    logging.info('Prior estimation for coverage %d: [Correct: %d, Error: %d]' %
            (sim_cov, beta_correct, beta_error))
    return (beta_correct, beta_error)

def calc_genotype_ztest(null_mean, null_stdv, sample_mean, sample_size):
    z_score = float(null_mean - sample_mean) / \
                (null_stdv / np.sqrt(sample_size))
    
    p_value = stats.norm.cdf(z_score, null_mean, null_stdv)
    # p_value = stats.norm.cdf(z_score)
    
    return p_value


def calc_genotype_quality(data_reads, priors, null_mean, null_stdv, basename):
    genotype_MAPs = dict()
    genotype_quality = dict()
    beta_correct, beta_error = priors
    MAPs = list()
    p_values = list()
    coverages = list()

    ZERO_COV = 0

    for position in data_reads:
        genotype = data_reads[position]['Genotype']
        reads = data_reads[position]['Reads']
        coverage = len(reads)
        
        alpha_correct = 0
        alpha_error = 0
        
        if coverage == 0:
            ZERO_COV += 1
        if (not FORCE_ZERO_COV) and coverage == 0:    
            #TODO: change these
            genotype_MAPs[position] = [1, 0, 0]
            genotype_quality[position] = 1
        else:
            coverages.append(coverage)

            for read in reads:
                if read[0] in genotype:
                    alpha_correct += 1
                else:
                    alpha_error += 1
            theta_map = float(beta_correct + alpha_correct - 1) / \
                        (beta_correct + alpha_correct + beta_error + alpha_error)

            p_value = calc_genotype_ztest(null_mean, null_stdv, theta_map, coverage)
        
            genotype_MAPs[position] = [theta_map, alpha_correct, alpha_error]
            genotype_quality[position] = p_value

            if PLOTS:
                MAPs.append(theta_map)
                p_values.append(p_value)
    
    if PLOTS:
        f, ax = plt.subplots(3, 1, figsize=(6, 6))
        sns.distplot(coverages, color='purple', ax=ax[0])
        ax[0].set(title="SNP coverage distribution")

        sns.distplot(MAPs, color='green', ax=ax[1])
        ax[1].set(title="theta_MAP distribution")

        sns.distplot(p_values, color='blue', ax=ax[2])
        ax[2].set(title="p-value distribution")
        ax[2].axvline(x=0.05, color='red', linestyle='--')
        plt.tight_layout()

        f.savefig('%s_genotype_plot.pdf' % basename)
    

    if ZERO_COV > 0:
        logging.warning('WARNING: A total of %d SNPs have 0 coverage' % ZERO_COV)

    return (genotype_MAPs, genotype_quality)

def is_phasing_supported(phasing1, phasing2, gen1, gen2):
    if gen1 == phasing1[0] and gen2 == phasing2[0]:
        return True
    elif gen1 == phasing1[1] and gen2 == phasing2[1]:
        return True
    else:
        return False

def calc_phasing_prob(n, k, error, MAPs=None):
    if MAPs is None:
        prob = k * np.log(1 - error*(2 - error)) + \
                (n-k) * np.log(error*(2-error))
        return prob
    else:
        logging.critical('Not implemented yet.')
        return 1

def calc_haplotype_quality(pblocks, pos_read, pos_read_nucl, error, genotype_quality, basename, MAPs=None):
    haplotypes_qualities = dict()
    
    coverages = list()
    supporting = list()
    opposing = list()
    phasing_qualities = list()

    ZERO_COV = 0
    
    bcount = 0
    for block in pblocks:
        bcount += 1
        logging.debug('Computing block %d of %d' % (bcount, len(pblocks)))

        for p1, p2 in pairwise(pblocks[block]):
            pos1, phas1 = p1
            pos2, phas2 = p2

            phas1 = phas1.split('|')
            phas2 = phas2.split('|')

            bridging_reads = pos_read[pos1] & pos_read[pos2]
            n = len(bridging_reads)

            if n == 0:
                ZERO_COV += 1
            if (not FORCE_ZERO_COV) and n == 0:    
                #TODO: change these
                log_quality = 0
                haplotypes_qualities[pos1, pos2] = 0
            else:
                coverages.append(n)
                
                k = 0
                for read in bridging_reads:
                    _, gen_p1 = pos_read_nucl[pos1, read]
                    _, gen_p2 = pos_read_nucl[pos2, read]

                    if is_phasing_supported(phas1, phas2, gen_p1, gen_p2):
                        k += 1

                log_quality = calc_phasing_prob(n, k, error)
                
                haplotypes_qualities[pos1, pos2] = np.log(genotype_quality[pos1]) + \
                                                    np.log(genotype_quality[pos1]) + \
                                                    log_quality

                if PLOTS:
                    supporting.append(k)
                    opposing.append(n-k)
                    phasing_qualities.append(log_quality)
    
    if PLOTS:
        f, ax = plt.subplots(2, 1, figsize=(6, 5))
        sns.distplot(coverages, ax=ax[0], hist=False,
                kde_kws={"color": "blue", "label": "Coverage"}
        )
        sns.distplot(supporting, ax=ax[0], hist=False,
                kde_kws={"color": "green", "label": "Supporting"}
        )
        sns.distplot(opposing, ax=ax[0], hist=False,
                kde_kws={"color": "red", "label": "Opposing"}
        )
        ax[0].set(title="phasing distributions")


        sns.distplot(phasing_qualities, ax=ax[1], color='orange')
        ax[1].set(title="phasing log-quality")
        plt.tight_layout()
        f.savefig('%s_phasing_plot.pdf' % basename)
    
    if ZERO_COV > 0:
        logging.warning('WARNING: A total of %d pairs of SNPs have 0 bridging reads' % ZERO_COV)

    return haplotypes_qualities

def get_final_quality(hap_quality, basename):
    logging.info('Haplotype quality (log) = %f over a total of %d phased loci.' % \
                (sum(hap_quality.values()), len(hap_quality)))

    if PLOTS:
        f, ax = plt.subplots(1, 1, figsize=(6, 3))
        
        if FORCE_ZERO_COV:
            hap_qual = list(hap_quality.values())
        else:
            hap_qual = [x for x in list(hap_quality.values()) if x != 0]
        sns.distplot(hap_qual, ax=ax, color='green')
        ax.set(title="haplotype log-quality")
        
        plt.tight_layout()
        f.savefig('%s_haplotype_plot.pdf' % basename)

def main(args):
    basename = '%s/%s' % (args.outdir,
                os.path.splitext(os.path.basename(args.data))[0])
    logging.info('Loading exported file...')
    position_genotype, pblocks, pos_read, pos_read_nucl = import_data(args.data)
    logging.info('Imported %d loci.' % len(position_genotype))

    # TODO: add it as input parameters
    tech_stdv = 2

    logging.info('Simulating priors...')
    h0_stdv = simulate_h0_stdv(args.epsilon, args.nullsim, args.techcov, 
                        tech_stdv, basename)
    priors = simulate_priors(args.priorcov, 1 - args.epsilon)
    
    logging.info('Computing genotype quality...')
    genotype_MAPs, genotype_quality = calc_genotype_quality(position_genotype, 
                                priors, args.epsilon, h0_stdv, basename)


    logging.info('Computing phasing quality...')

    hap_quality = calc_haplotype_quality(pblocks, pos_read, pos_read_nucl, 
                args.epsilon, genotype_quality, basename)

    get_final_quality(hap_quality, basename)
    
    
def setup(args):
    global PLOTS
    PLOTS = not args.disable_plots

    global FORCE_ZERO_COV
    FORCE_ZERO_COV = args.force_0cov

    try:
        os.makedirs(args.outdir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(args.outdir):
            pass
        else:
            raise

    verbose = False
    if verbose:
        logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s - %(message)s',
                        datefmt='%H:%M:%S')
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(message)s',
                            datefmt='%H:%M:%S')
        
    if PLOTS:
        global sns
        global plt
        global matplotlib

        import matplotlib
        matplotlib.use('Agg')
        import seaborn as sns
        import matplotlib.pyplot as plt

        sns.set_context("paper")
        sns.set_style("ticks")
        
        main(args)
    else:
        main(args)

if __name__ == '__main__':
    parser = generate_argparser()
    args = parser.parse_args()
    setup(args)
