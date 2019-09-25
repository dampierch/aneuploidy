#!/usr/bin/env python3

## parse_gdc_info.py

    # this script takes tab-delimited output from gdc_write_anno.py (pheno.tsv)
    # and extracts file_id, file_name, subject_id, tissue_type and reformats
    # to match WP *.file_info file

    # *.file_info file fields: subject_id, tissue_type, file_name, file_id
    # importantly, *.file_info is sorted by subject_id and tissue_type

import fileinput
import re
import pandas as pd

data_home = '/scratch/chd5n/aneuploidy/'
anno_path = data_home + 'raw-data/annotations/'
input_file = anno_path + 'pheno.tsv'
output_file = anno_path + 'coad-read.file_info'
error_file = anno_path + 'coad-read.sub_errors'
normal = ['Blood Derived Normal', 'Solid Tissue Normal']
tumor = 'Primary Tumor'
print_fields = ['subject_id', 'tissue_type', 'file_name', 'file_id']

sub_set = {}
sub_list = []

with fileinput.input(files=(input_file)) as in_f:
    for d_line in in_f:
        if in_f.isfirstline():
            header_fields = d_line.strip().split('\t')
            continue
        data = dict(list(zip(header_fields, d_line.strip().split('\t'))))
        if data['subject_id'] in sub_set:
            if data['tissue_type'] in normal:
                sub_set[data['subject_id']]['normal'].append(data)
            elif data['tissue_type'] == tumor:
                sub_set[data['subject_id']]['tumor'].append(data)
            else:
                pass
        else:
            sub_list.append(data['subject_id'])
            sub_set[data['subject_id']] = { 'id': data['subject_id'], 'normal': [], 'tumor':[] }
            if data['tissue_type'] in normal:
                sub_set[data['subject_id']]['normal'].append(data)
            elif data['tissue_type'] == tumor:
                sub_set[data['subject_id']]['tumor'].append(data)
            else:
                pass

# we want only subjects with at least one normal and one tumor sample
# we are indifferent between blood derived normal and solid tissue normal
# we want highest total sequences

problems = []
with open(output_file, 'w+') as out_f:
    for sub in sub_list:
        if len(sub_set[sub]['normal']) > 0 and len(sub_set[sub]['tumor']) > 0:
            df_norm = pd.DataFrame(sub_set[sub]['normal'])
            df_tum = pd.DataFrame(sub_set[sub]['tumor'])
            norm_select = df_norm['total_seq']==df_norm['total_seq'].max()
            tum_select = df_tum['total_seq']==df_tum['total_seq'].max()
            norm = df_norm.loc[norm_select, print_fields].values.tolist()[0]
            norm[1] = 'normal'
            tum = df_tum.loc[tum_select, print_fields].values.tolist()[0]
            tum[1] = 'tumor'
            out_f.write('\t'.join(norm) + '\n')
            out_f.write('\t'.join(tum) + '\n')
        else:
            problems.append(sub)

with open(error_file, 'w+') as out_f:
    for sub in problems:
        out_f.write(sub + '\n')
