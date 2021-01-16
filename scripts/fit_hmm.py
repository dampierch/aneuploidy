#!/usr/bin/env python3

## README
## aim: model allele fractions at heterozygous sites in normal and tumor samples
## method: use hidden markov model as implemented in pomegranate
## concept:
##   -- given a series of allele fraction states (one series per chrom),
##   assign probabilities for diploid or haploid states
##   -- regions not in diploid or haploid states could be called aneuploid
## implementation:
##   -- read all chromosomal data before doing analysis
##   -- censor bad sites and estimate LOH mean
##   -- for each chrom, extract and transform allele fractions, store in array
##   -- for each array, fit an HMM and set probabilities
##   -- store state probabilities for all sites in tsv file
## usage example:
##   fit_alleles3ac.py TCGA-*_cnts2R.tsv.gz > TCGA-*_hmmfit.tsv
## input: TCGA-*_cnts2R.tsv
## subject_id    chrom    pos    t_type    maj    maj_cnt    min    min_cnt    maj_fract
## TCGA-3L-AA1B  chr1     69511  norm      A      358        G      312        0.53433
## TCGA-3L-AA1B  chr1     69511  tumor     A      346        G      304        0.53231
## TCGA-3L-AA1B  chr1     826965 norm      A      81         G      144        0.36000
## TCGA-3L-AA1B  chr1     826965 tumor     A      111        G      131        0.45868


import os
import sys
import gzip
import argparse
import numpy as np  ## must use pip3 install --user numpy==1.16.1
from pomegranate import State, NormalDistribution, HiddenMarkovModel


mod_name = 'fit_hmm'


def read_counts(target):
    '''
    -- reads in heterozygous site count data and filters sites to exclude those
    that do not exist in both normal and tumor samples
    -- also extracts major allele fraction of normal sample after site filter to
    use for cleaning data
    '''
    data_miss = []
    data_norm = []
    data_tumor = []
    fields = False
    with gzip.open(target, 'rt') as infile:
        for inline in infile:
            ## capture header
            if not fields:
                infields = inline.strip('\n').split('\t')
                fields = True
                continue
            ## collect data
            data = dict(zip(infields, inline.strip('\n').split('\t')))
            if data['maj_fract'] == 'NA':
                data_miss.append(data)
                continue
            if data['t_type'] == 'norm':
                data_norm.append(data)
            else:
                data_tumor.append(data)
    d_norm = {'::'.join([i['chrom'], i['pos']]): i for i in data_norm}
    d_tumor = {'::'.join([i['chrom'], i['pos']]): i for i in data_tumor}
    l_norm = [d_norm[k] for k in d_norm if k in d_tumor]
    l_tumor = [d_tumor[k] for k in d_tumor if k in d_norm]
    l_norm_miss = [d_norm[k] for k in d_norm if k not in d_tumor]
    l_norm_tum = [d_tumor[k] for k in d_tumor if k not in d_norm]
    data_miss = data_miss + l_norm_miss + l_norm_tum
    ## for censoring
    fract_norm = [float(i['maj_fract']) for i in l_norm]
    return data_miss, l_norm, l_tumor, fract_norm


def set_range(data, range=1.6, median=0.5):
    '''
    parameters:
    -- range: inspired by R boxplot(), range parameter determines how far out
    from observed interquartile range values will be accepted
    -- median: expected median of data
    '''
    min, q1, med, q3, max = fivenum(data)
    r = abs(q3 - q1) * range
    r_min = median - r
    r_max = median + r
    return r_min, r_max


def init_model(start_dip, stay_state, mean_eu, sd_eu, mean_loh):

    ## define distributions
    d_eu = NormalDistribution(mean_eu, sd_eu)  ## euploid enriched at 0
    d_loh = NormalDistribution(mean_loh, sd_eu)  ## loss of heterozygosity enriched at 1
    d_aneu = NormalDistribution(mean_loh / 2.0, sd_eu * 1.4)  ## aneuploid enriched at 1

    ## define states
    s_eu = State(d_eu, name='EU')  ## enriched at 0
    s_loh = State(d_loh, name='LOH')  ## enriched at 1
    s_aneu = State(d_aneu, name='ANEU')  ## enriched at 1

    ## define model and pass in states
    model = HiddenMarkovModel()
    model.add_states(s_eu, s_loh, s_aneu)

    ## define transition matrix (state a, state b, probability)
    model.add_transition(model.start, s_eu, start_dip)
    model.add_transition(model.start, s_loh, 1.0 - start_dip - 0.1)
    model.add_transition(model.start, s_aneu, 0.1)

    model.add_transition(s_eu, s_eu, stay_state)
    model.add_transition(s_eu, s_loh, 1.0 - 4 * stay_state / 5 - 0.001)
    model.add_transition(s_eu, s_aneu, 1.0 - stay_state / 5 - 0.001)
    model.add_transition(s_eu, model.end, 0.002)

    model.add_transition(s_loh, s_loh, stay_state)
    model.add_transition(s_loh, s_eu, 1.0 - 4 * stay_state / 5 - 0.001)
    model.add_transition(s_loh, s_aneu, 1.0 - stay_state / 5 - 0.001)
    model.add_transition(s_loh, model.end, 0.002)

    model.add_transition(s_aneu, s_aneu, stay_state)
    model.add_transition(s_aneu, s_eu, 1.0 - stay_state / 2 - 0.001)
    model.add_transition(s_aneu, s_loh, 1.0 - stay_state / 2 - 0.001)
    model.add_transition(s_aneu, model.end, 0.002)

    ## finalize internal structure
    model.bake()
    ## only train transitions, not emissions
    model.freeze_distributions()

    return model


def fit_hmm(allele_freqs, start_dip, stay_state, mean_eu, sd_eu, mean_loh):
    model = init_model(start_dip, stay_state, mean_eu, sd_eu, mean_loh)
    model.fit([allele_freqs], algorithm='baum-welch')

    ## check
    [norm_stateProbs, norm_transFreqs] = model.forward_backward(allele_freqs)
#   print "norm_stateProbs:\n",norm_stateProbs
#   print "norm_transFreqs:\n",norm_transFreqs[0:100]

    post_states = model.predict(allele_freqs, algorithm='map')
    post_probs = model.predict_proba(allele_freqs)
    return post_states, post_probs, model


def map_states(model):
    post_states = []
    state_pos = {}
    for ix, state in enumerate(model.states):
        post_states.append(state.name)
        state_pos[state.name] = ix
    return state_pos, post_states


def fivenum(data):
    ## five number summary
    return np.percentile(data, [0, 25, 50, 75, 100], interpolation='midpoint')


def main():

    ## provide command line
    print("# " + ' '.join(sys.argv))

    sqrt_sqrt2 = np.sqrt(np.sqrt(2))
    curr_chrom = ''

    target = ''.join([args.count_path, args.subject_id, '_cnts2R.tsv.gz'])
    data_miss, data_norm, data_tumor, fract_norm = read_counts(target)
    r_min, r_max = set_range(fract_norm)

## double check that chrom_data['norm'] and chrom_data['tumor'] are in phase

## do statistics on chrom_data['norm'] to identify outliers, make a new
## chrom_data_norm, chrom_data_tumor, allele_freqs_norm, allele_freqs_tumor
## for good sites


def clean_data(data_norm, data_tumor, r_min, r_max):
    '''
    -- make sure all site pairs are in same order in normal and tumor
    -- remote outlier sites based on range criteria
    '''
    ## all the row data
    chrom_data_norm = []
    chrom_data_tumor = []

    ## censored unscaled fractions
    cen_allele_fracts_norm = []
    cen_allele_fracts_tumor = []

    ## scaled censored fractions
    sc_allele_fracts_norm = []
    sc_allele_fracts_tumor = []

    normal_cnts = 0
    tumor_cnts = 0

    for idx, data_n in enumerate(data_norm):
        data_t = data_tumor[idx]
        if (int(n_data['pos']) != int(t_data['pos'])):
            sys.stderr.write("*** Phase error -- norm: %s %s :: tumor %s %s\n"%(n_data['chrom'],n_data['pos'],t_data['chrom'],t_data['pos']))
            exit(1)

        n_fract = float(n_data['fract'])

        if (iqr_min < n_fract < iqr_max):
            chrom_data_norm.append(n_data)
            chrom_data_tumor.append(t_data)

            normal_cnts += int(n_data['ref_cnt'])+int(n_data['alt_cnt'])
            tumor_cnts += int(t_data['ref_cnt'])+int(t_data['alt_cnt'])

            t_fract = float(t_data['fract'])

            cen_allele_fracts_norm.append(n_fract)
            cen_allele_fracts_tumor.append(t_fract)

            scaled_n_fract = abs(n_fract - 0.5) * 2.0
            scaled_t_fract = abs(t_fract - 0.5) * 2.0

            sc_allele_fracts_norm.append(scaled_n_fract)
            sc_allele_fracts_tumor.append(scaled_t_fract)

    ## calculate std.dev of cen_allele_fracts
    tumor_fa_score = np.std(cen_allele_fracts_tumor)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fit 3 state HMM classifier')
    parser.add_argument(
        '--subject_id', help='subject_id to be analyzed',
        action='store', dest='subject_id',
        default=False
    )
    parser.add_argument(
        '--count_path', help='path to intercalated counts',
        action='store', dest='count_path',
        default='/scratch/chd5n/aneuploidy/hetsites-data/'
    )
    parser.add_argument(
        '--out_path', help='path to output directory',
        action='store', dest='out_path',
        default='/scratch/chd5n/aneuploidy/hmm-fit/'
    )
    args = parser.parse_args()
    main(args)
else:
    print('functions loaded for', mod_name)
