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
##   -- store state probabilities for all sites in R data object

## input: TCGA-*_cnts2R.tsv
## subject_id    chrom    pos    t_type    maj    maj_cnt    min    min_cnt    maj_fract
## TCGA-3L-AA1B  chr1     69511  norm      A      358        G      312        0.53433
## TCGA-3L-AA1B  chr1     69511  tumor     A      346        G      304        0.53231
## TCGA-3L-AA1B  chr1     826965 norm      A      81         G      144        0.36000
## TCGA-3L-AA1B  chr1     826965 tumor     A      111        G      131        0.45868


import fileinput
import sys
import numpy as np  ## must use pip3 install --user numpy==1.16.1
from pomegranate import State, NormalDistribution, HiddenMarkovModel


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
    '''
    five number summary
    '''
    return np.percentile(data, [0, 25, 50, 75, 100], interpolation='midpoint')


def main():

    ## provide command line
    print("# " + ' '.join(sys.argv))

    sqrt_sqrt2 = np.sqrt(np.sqrt(2))

    curr_chrom = ''


THIS SHOULD BE A SEPARATE FUNCTION TO LOAD DATA
    norm_data_raw = []
    norm_freq_raw = []
    tumor_data_raw = []

    have_fields = False
    for line in fileinput.input():

        ## copy previous comments
        if line[0] == '#':
            print(line, end='')

        if not have_fields:
            field_names = line.strip('\n').split('\t')
            have_fields = True

            continue

        r_data = dict(zip(field_names, line.strip('\n').split('\t')))

        if r_data['fract'] == 'NA':
            sys.stderr.write("*** %s missing data: %s" % (fileinput.filename(),line))
            continue

        if r_data['t_type'] == 'norm':
            norm_data_raw.append(r_data)
            ## used for censoring
            norm_freq_raw.append(float(r_data['fract']))
        else:
            tumor_data_raw.append(r_data)


## done reading data

## double check that chrom_data['norm'] and chrom_data['tumor'] are in phase

## do statistics on chrom_data['norm'] to identify outliers, make a new
## chrom_data_norm, chrom_data_tumor, allele_freqs_norm, allele_freqs_tumor
## for good sites

# need python implementation of this
# fv5 <- fivenum(wide_thet_data$fract_norm)
# iqr_x = abs(fv5[4]-fv5[2])*1.6
# c_iqr_min = fv5[3] - iqr_x
# c_iqr_max = fv5[3] + iqr_x
# wide_het_data$is.good = wide_het_data$fract_norm >= c_iqr_min & wide_het_data$fract_norm<= c_iqr_max

THIS SHOULD BE A SEPARATE FUNCTION TO CENSOR DATA
    n_min, n_q1, n_med, n_q3, n_max = fivenum(norm_freq_raw)

    iqr_x = abs(n_q3 - n_q1) * 1.6
    iqr_min = 0.5 - iqr_x
    iqr_max = 0.5 + iqr_x

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

    for ix, n_data in enumerate(norm_data_raw):
        t_data = tumor_data_raw[ix]
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
