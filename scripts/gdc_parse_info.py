#!/usr/bin/env python3


## gdc_parse_info.py

    ## this script takes tab-delimited output from gdc_write_anno.py (pheno.tsv)
    ## and extracts file_id, file_name, subject_id, tissue_type and reformats
    ## to match WP *.file_info file

    ## *.file_info file fields: subject_id, tissue_type, file_name, file_id
    ## importantly, *.file_info is sorted by subject_id and tissue_type


import pandas as pd


def get_subject_info(input_file):
    '''
    read pheno file info into dict of dicts with keys set to subject_ids and
    values set to dicts of info for files with separate list of dicts for
    normal and tumor samples for each subject
    '''
    normal_labs = ['Blood Derived Normal', 'Solid Tissue Normal']
    tumor_lab = 'Primary Tumor'
    sub_dict = {}
    with open(input_file,'r') as in_f:
        line_num = 0
        for d_line in in_f:
            if line_num == 0:
                header_fields = d_line.strip().split('\t')
                line_num = line_num + 1
                continue
            data = dict(list(zip(header_fields,d_line.strip().split('\t'))))
            if data['subject_id'] in sub_dict:
                if data['tissue_type'] in normal_labs:
                    sub_dict[data['subject_id']]['normal'].append(data)
                elif data['tissue_type'] == tumor_lab:
                    sub_dict[data['subject_id']]['tumor'].append(data)
                else:
                    pass
            else:
                sub_dict[data['subject_id']] = { 'id': data['subject_id'], 'normal': [], 'tumor': [] }
                if data['tissue_type'] in normal_labs:
                    sub_dict[data['subject_id']]['normal'].append(data)
                elif data['tissue_type'] == tumor_lab:
                    sub_dict[data['subject_id']]['tumor'].append(data)
                else:
                    pass
            line_num = line_num + 1
    return sub_dict


def parse_info(output_file,sub_dict):
    '''
    keep only subjects with at least one normal and one tumor sample;
    indifferent between blood derived normal and solid tissue normal;
    want samples with highest total sequences for subjects with multiple samples
    '''
    problems = []
    print_fields = ['subject_id', 'tissue_type', 'file_name', 'file_id']
    with open(output_file, 'w+') as out_f:
        for sub in sub_dict.keys():
            if len(sub_dict[sub]['normal']) > 0 and len(sub_dict[sub]['tumor']) > 0:
                df_norm = pd.DataFrame(sub_dict[sub]['normal'])
                df_tum = pd.DataFrame(sub_dict[sub]['tumor'])
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
    return problems


def write_subject_errors(error_file,problems):
    with open(error_file, 'w+') as out_f:
        for sub in problems:
            out_f.write(sub + '\n')


def main():
    anno_path = '/scratch/chd5n/aneuploidy/raw-data/annotations/'
    input_file = anno_path + 'pheno.tsv'
    output_file = anno_path + 'coad-read.file_info'
    error_file = anno_path + 'coad-read.sub_errors'
    sub_dict = get_subject_info(input_file)
    problems = parse_info(output_file,sub_dict)
    write_subject_errors(error_file,problems)


main()




# ## for testing
# file_info = pd.read_csv('/scratch/chd5n/aneuploidy/raw-data/annotations/coad-read.file_info',sep='\t',names=['subject_id','tissue_type','file_name','file_id'])
# for sub in sub_dict.keys():
#     if len(sub_dict[sub]['normal']) > 0 and len(sub_dict[sub]['tumor']) > 0:
#         df_norm = pd.DataFrame(sub_dict[sub]['normal'])
#         df_tum = pd.DataFrame(sub_dict[sub]['tumor'])
#         norm_select = df_norm['total_seq']==df_norm['total_seq'].max()
#         tum_select = df_tum['total_seq']==df_tum['total_seq'].max()
#         norm = df_norm.loc[norm_select, print_fields].values.tolist()[0]
#         tum = df_tum.loc[tum_select, print_fields].values.tolist()[0]
#         fi_ss = file_info.loc[file_info['subject_id']==sub,:]
#         if list(norm[-1]==fi_ss.loc[fi_ss['tissue_type']=='normal','file_id'])[0] and list(tum[-1]==fi_ss.loc[fi_ss['tissue_type']=='tumor','file_id'])[0]:
#             pass
#         else:
#             print('problem with',sub)
