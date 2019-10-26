#!/usr/bin/env python3

## hetcnts_2R.py :: input subject_[normal/tumor]_hetcnts.bed :: output TSV

    # from within het_counter.sh
    # usage: hetcnts_2R.py --sample_only ${subject_id} --normal_hetcnts_bed ${crunch_dir}${subject_id}_normal_hetcnts.bed --tumor_hetcnts_bed ${crunch_dir}${subject_id}_tumor_hetcnts.bed 2> subject_R_missing.err

    # this script takes bed files with simple allele counts at heterozygous
    # sites from normal and tumor bed files, then writes allele counts to file
    # for R to process; ensures no loss of hetsites even if some are not in both
    # samples


import sys
import argparse


data_home = '/scratch/chd5n/aneuploidy/'
anno_home = data_home + 'raw-data/annotations/'
seq_home = data_home + 'raw-data/sequencing/'
crunch_path = seq_home + 'crunch/'
in_fields = ['chrom', 'pos0', 'pos', 'an_type', 'score', 'strand', 'info_str']
out_fields = ['subject_id', 'chrom', 'pos', 't_type', 'maj', 'maj_cnt', 'min', 'min_cnt', 'maj_fract']


def get_hetcnts_list(file_name, tissue_type):
    data_list = []
    with open(file_name,'r') as in_f:
        for in_line in in_f:
            data = dict(list(zip(in_fields,in_line.strip('\n').split('\t'))))
            data['t_type'] = tissue_type
            data['subject_id'] = args.subject_id
            data_list.append(data)
    return data_list


def parse_info(info_str):
    '''
    parses the info string to return major and minor allele (simple) depth at
    given position
    '''
    af_str, ad_str, sd_str = info_str.split(';')
    maj_str, min_str = sd_str.split('::')[1].split(',')
    return maj_str.split('=') + min_str.split('=')


def get_rel_pos(norm_dict, tumor_dict):
    '''
    determines relative position of currently loaded hetsites for each
    tumor-normal pair of inputs
    '''
    if norm_dict['chrom'] == tumor_dict['chrom'] and int(norm_dict['pos']) == int(tumor_dict['pos']):
        return 0
    elif norm_dict['chrom'] == tumor_dict['chrom']:
        if int(norm_dict['pos']) < int(tumor_dict['pos']):
            return -1
        else:
            return 1
    else:
        if norm_dict['chrom'] < tumor_dict['chrom']:
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
    if maj_flt < 1 and min_flt < 1:
        sys.stderr.write('%s - %s - division by zero: %s:%s\n' % (data['subject_id'], data['t_type'], data['chrom'], data['pos']))
        data['maj_fract'] = 'NA'
    else:
        data['maj_fract'] = '%.5f' % (maj_flt/(maj_flt + min_flt))
    out_f.write('\t'.join([data[x] for x in out_fields])+'\n')


def main(args):
    ## prepare normal data as list of dicts
    norm_data = get_hetcnts_list(args.normal_hetcnts_bed, 'norm')
    ## prepare tumor data as list of dicts
    tumor_data = get_hetcnts_list(args.tumor_hetcnts_bed, 'tumor')
    ## setup output
    with open(crunch_path + args.subject_id + '_cnts2R.tsv', 'w') as out_f:
        out_f.write('\t'.join(out_fields)+'\n')
        ## confirm that both sets of data have same coordinates
        tumor_counter = 0
        for norm_dict in norm_data:
            tumor_dict = tumor_data[tumor_counter]
            rel_pos = get_rel_pos(norm_dict, tumor_dict)
            if rel_pos == 0:
                ## when hetsite is same in normal and tumor,
                ## print allele fractions for both
                print_het_cnts(out_f, norm_dict, out_fields)
                print_het_cnts(out_f, tumor_dict, out_fields)
                tumor_counter = tumor_counter + 1
            elif rel_pos < 0:
                ## when hetsite in normal is to left of tumor,
                ## print err, then print allel fraction for normal
                sys.stderr.write('*** phase error: tumor missing data\n')
                print_het_cnts(sys.stderr, norm_dict, out_fields)
                print_het_cnts(sys.stderr, tumor_dict, out_fields)
                print_het_cnts(out_f, norm_dict, out_fields)
                ## keep iterating over norm_data until norm pos
                ## catches up with tumor
                continue
            else:
                ## when hetsite in normal is to right of tumor,
                ## print err, then print allele fraction for lagging
                ## tumor, current normal, and next tumor, then adv
                ## tumor counter; maybe tumor will catch up
                sys.stderr.write('*** phase error: normal missing data\n')
                print_het_cnts(sys.stderr, norm_dict, out_fields)
                print_het_cnts(sys.stderr, tumor_dict, out_fields)
                print_het_cnts(out_f, tumor_dict, out_fields)
                print_het_cnts(out_f, norm_dict, out_fields)
                tumor_counter = tumor_counter + 1
                tumor_dict = tumor_data[tumor_counter]
                print_het_cnts(out_f, tumor_dict, out_fields)
                tumor_counter = tumor_counter + 1


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='transform .bed file of allele frequencies in TSV for R')
    parser.add_argument('--subject_id',help='subject_id to be analyzed',action='store',dest='subject_id',default=False)
    parser.add_argument('--normal_hetcnts_bed',help='path to normal hetcnts bed file',action='store',dest='normal_hetcnts_bed',default=False)
    parser.add_argument('--tumor_hetcnts_bed',help='path to tumor hetcnts bed file',action='store',dest='tumor_hetcnts_bed',default=False)
    args = parser.parse_args()

    main(args)






######### TESTING

# parser = argparse.ArgumentParser()
# parser.add_argument('--subject_id',help='subject_id to be analyzed',action='store',dest='subject_id',default=False)
# parser.add_argument('--normal_hetcnts_bed',help='path to normal hetcnts bed file',action='store',dest='normal_hetcnts_bed',default=False)
# parser.add_argument('--tumor_hetcnts_bed',help='path to tumor hetcnts bed file',action='store',dest='tumor_hetcnts_bed',default=False)
# args = parser.parse_args('--subject TCGA-T9-A92H --normal_hetcnts_bed /scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/TCGA-T9-A92H_normal_hetcnts.bed --tumor_hetcnts_bed /scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/TCGA-T9-A92H_tumor_hetcnts.bed'.split())
#
# norm_data = get_hetcnts_list(args.normal_hetcnts_bed, 'norm')
# tumor_data = get_hetcnts_list(args.tumor_hetcnts_bed, 'tumor')

## notes
    # in theory, there should never be tumor data for which normal is missing,
    # since hetsites are chosen based on hetsites in normal; if that hetsite was
    # in normal, it should be in normal again
    # HOWEVER
    # in practice, depth and quality thresholds may cause the loss of sites in
    # the normal samples
