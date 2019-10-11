#!/usr/bin/env python3

## hetcnts_2R.py :: input subject_[normal/tumor]_hetcnts.bed :: output TSV

    # from within het_counter.sh
    # usage: hetcnts_2R.py --sample_only ${subject_id}

    # this script takes bed files with simple allele counts at heterozygous
    # sites from normal and tumor bed files, then writes allele counts to file
    # for R to process

    # fields: chr pos t_status ref ref_cnt alt alt_cnt fract

    # this script ensures no loss of phasing (same location in both tumor/normal)


import sys
import argparse


data_home = '/scratch/chd5n/aneuploidy/'
anno_home = data_home + 'raw-data/annotations/'
seq_home = data_home + 'raw-data/sequencing/'
crunch_path = seq_home + 'crunch/'
in_fields = ['chrom', 'pos0', 'pos', 'an_type', 'score', 'strand', 'info_str']
out_fields = ['subject_id', 'chrom', 'pos', 't_type', 'maj', 'maj_cnt', 'min', 'min_cnt', 'fract']


def parse_info(info_str):
    '''
    parses the info string to return major and minor allele (simple) depth at
    given position
    '''
    (af_str, ad_str, sd_str) = info_str.split(';')
    (maj_str, min_str) = sd_str.split('::')[1].split(',')
    return (maj_str.split('=') + min_str.split('='))


def same_pos(a_data, b_data):
    '''
    determines whether locations are truly the same in each tumor-normal pair
    of input files
    '''
    if a_data['chrom'] == b_data['chrom'] and int(a_data['pos']) == int(b_data['pos']):
        return 0
    elif a_data['chrom'] == b_data['chrom']:
        if int(a_data['pos']) < int(b_data['pos']):
            return -1
        else:
            return 1
    else:
        if (a_data['chrom'] < b_data['chrom']):
            return -1
        else:
            return 1


def print_het_cnts(out_f, data, out_fields):
    '''
    assigns observed counts to their respective allele, calculates ratio, then
    writes information to file
    '''
    allele_data = parse_info(data['info_str'])
    data['maj'] = allele_data[0]
    data['maj_cnt'] = allele_data[1]
    data['min'] = allele_data[2]
    data['min_cnt'] = allele_data[3]
    maj_flt = float(data['maj_cnt'])
    min_flt = float(data['min_cnt'])
    if (maj_flt < 1 and min_flt < 1):
        sys.stderr.write('%s - %s - division by zero: %s:%s\n' % (data['subject_id'], data['t_type'], data['chrom'], data['pos']))
        data['fract'] = 'NA'
    else:
        data['fract'] = '%.5f' % (maj_flt/(maj_flt + min_flt))
    out_f.write('\t'.join([data[x] for x in out_fields])+'\n')


def main(args):
    ## prepare normal data as list of dicts
    norm_data = []
    with open(args.normal_hetcnts_bed, 'r') as in_f:
        for in_line in in_f:
            data = dict(list(zip(in_fields,in_line.strip('\n').split('\t'))))
            data['t_type'] = 'norm'
            data['subject_id'] = args.subject_id
            norm_data.append(data)
    ## prepare tumor data as list of dicts
    tumor_data = []
    with open(args.tumor_hetcnts_bed, 'r') as in_f:
        for in_line in in_f:
            data = dict(list(zip(in_fields,in_line.strip('\n').split('\t'))))
            data['t_type'] = 'tumor'
            data['subject_id'] = args.subject_id
            tumor_data.append(data)
    ## setup output
    with open(crunch_path + args.subject_id + '_cnts2R.tsv', 'w') as out_f:
        out_f.write('\t'.join(out_fields)+'\n')
        ## confirm that both sets of data have same coordinates
        tx = 0
        for norm_d in norm_data:
            tumor_d = tumor_data[tx]
            rel_pos = same_pos(norm_d, tumor_d)
            if rel_pos == 0:
                print_het_cnts(out_f, norm_d, out_fields)
                print_het_cnts(out_f, tumor_d, out_fields)
                tx = tx + 1
            elif rel_pos < 0:  # norm_d is in front of tumor_d
                sys.stderr.write('*** phase error: tumor missing data\n')
                print_het_cnts(sys.stderr, norm_d, out_fields)
                print_het_cnts(sys.stderr, tumor_d, out_fields)
                print_het_cnts(out_f, norm_d, out_fields)
                continue
            else:
                sys.stderr.write('*** phase error: normal missing data\n')
                print_het_cnts(sys.stderr, norm_d, out_fields)
                print_het_cnts(sys.stderr, tumor_d, out_fields)
                print_het_cnts(out_f, tumor_d, out_fields)
                print_het_cnts(out_f, norm_d, out_fields)
                print_het_cnts(out_f, tumor_data[tx+1], out_fields)
                tx = tx + 2
                continue


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='transform .bed file of allele frequencies in TSV for R')
    parser.add_argument('--subject_id',help='subject_id to be analyzed',action='store',dest='subject_id',default=False)
    parser.add_argument('--normal_hetcnts_bed',help='path to normal hetcnts bed file',action='store',dest='normal_hetcnts_bed',default=False)
    parser.add_argument('--tumor_hetcnts_bed',help='path to tumor hetcnts bed file',action='store',dest='tumor_hetcnts_bed',default=False)
    args = parser.parse_args()

    main(args)






######### TESTING

# parser=argparse.ArgumentParser()
# parser.add_argument('--sample_only',help='infer sample_normal, sample_tumor names',action='store_true',default=False)
# parser.add_argument('in_files', metavar='FILE', nargs='*')
# args=parser.parse_args('--sample_only TCGA-T9-A92H'.split())
#
# args = parser.parse_args('--subject TCGA-T9-A92H --normal_hetcnts_bed /scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/TCGA-T9-A92H_normal_hetcnts.bed --tumor_hetcnts_bed /scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/TCGA-T9-A92H_tumor_hetcnts.bed'.split())
#
# norm_data = []
# with open(args.normal_hetcnts_bed, 'r') as in_f:
#     for in_line in in_f:
#         data = dict(list(zip(in_fields,in_line.strip('\n').split('\t'))))
#         data['t_type'] = 'norm'
#         data['subject_id'] = args.subject_id
#         norm_data.append(data)
#
# tumor_data = []
# with open(args.tumor_hetcnts_bed, 'r') as in_f:
#     for in_line in in_f:
#         data = dict(list(zip(in_fields,in_line.strip('\n').split('\t'))))
#         data['t_type'] = 'tumor'
#         data['subject_id'] = args.subject_id
#         tumor_data.append(data)
