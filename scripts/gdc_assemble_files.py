#!/usr/bin/env python3

## gdc_assemble_files.py

    ## this script takes tab-delim output from gdc_parse_info.py (*.file_info)
    ## and extracts subject_id, tissue_type, file_name, file_id and reformats
    ## into *.file_set; also moves files to data crunch directory

    ## *.file_set has three fields: subject_id, file_name_normal, file_name_tumor
    ## file='/home/chd5n/projects/aneuploidy/scripts/gdc_assemble_files.py'
    ## exec(open(file).read())


import os
import subprocess
import argparse


def get_file_info(input_file):
    '''
    read file info into dict of dicts with keys set to subject_ids and
    values set to dicts of info for files with separate dicts for
    normal and tumor samples for each subject

    i.e. for each subject_id as key, create values of normal and tumor; then
    for each value, create dict of file name and file uuid; note that dict
    is modified indirectly with modifications to values
    '''
    subject_dict = {}
    with open(input_file) as in_f:
        for line in in_f:
            if line[0] == '#':
                continue
            subject_id, tissue_type, file_name, file_id = line.strip('\n').split('\t')
            if subject_id in subject_dict.keys():
                sub = subject_dict[subject_id]
            else:
                subject_dict[subject_id] = {}
                sub = subject_dict[subject_id]
            if tissue_type in sub:
                sub[tissue_type].append({'name': file_name, 'uuid': file_id}) ## this code modifies subject_dict
            else:
                sub[tissue_type] = [{'name': file_name, 'uuid': file_id}] ## this code modifies subject_dict
    return subject_dict


def check_download(sub,download_path):
    '''
    check to see if files exist
    '''
    for tissue_type in ('normal','tumor'):
        if tissue_type in sub:
            flag = 'have_' + tissue_type
            for file in sub[tissue_type]:
                if os.path.exists(download_path + file['uuid']):
                    if os.path.exists('%s/%s' % (download_path + file['uuid'], file['name'])):
                        if flag not in sub:
                            sub[flag] = file
    return sub


def assemble_files(download_path,sub,crunch_path):
    '''
    move newly downloaded files to crunch directory for analysis;
    must include bai for pysam pileup
    '''
    cmd_n_bam = 'mv %s/%s %s' % (download_path + sub['have_normal']['uuid'], sub['have_normal']['name'], crunch_path)
    cmd_n_bai = 'mv %s/%s %s' % (download_path + sub['have_normal']['uuid'], sub['have_normal']['name'].split('.')[0] + '.bai', crunch_path)
    cmd_t_bam = 'mv %s/%s %s' % (download_path + sub['have_tumor']['uuid'], sub['have_tumor']['name'], crunch_path)
    cmd_t_bai = 'mv %s/%s %s' % (download_path + sub['have_tumor']['uuid'], sub['have_tumor']['name'].split('.')[0] + '.bai', crunch_path)
    subprocess.call(cmd_n_bam, shell=True)
    subprocess.call(cmd_n_bai, shell=True)
    subprocess.call(cmd_t_bam, shell=True)
    subprocess.call(cmd_t_bai, shell=True)


def process_file_info(output_file,subject_dict,download_path,crunch_path):
    '''
    check for complete pairs and existence of downloaded files; pairs checked
    in gdc_parse_info.py too; if pairs complete and files exist, write output
    as *.file_set and assemble files for analysis
    '''
    problems = []
    with open(output_file, 'w+') as out_f:
        for subject_id in subject_dict.keys():
            sub = subject_dict[subject_id]
            sub = check_download(sub,download_path)
            if 'have_normal' in sub and 'have_tumor' in sub:
                out_f.write('\t'.join([subject_id, sub['have_normal']['name'], sub['have_tumor']['name']]) + '\n')
                assemble_files(download_path,sub,crunch_path)
            else:
                if 'have_normal' not in sub:
                    problems.append('%s missing normal\n' % (subject_id))
                if 'have_tumor' not in sub:
                    problems.append('%s missing tumor\n' % (subject_id))
    return problems


def record_errors(error_file,problems):
    '''
    list subjects that fail pair test or download confirmation test
    '''
    with open(error_file, 'w+') as out_f:
        for sub in problems:
            out_f.write(sub)


def set_current(original,revised):
    cmd = ' '.join(['cp',original,revised])
    subprocess.call(cmd, shell=True)


def main(args):
    dt = args.dt
    set_num = args.set_num
    data_home = '/scratch/chd5n/aneuploidy/'
    anno_home = data_home + 'raw-data/annotations/'
    seq_home = data_home + 'raw-data/sequencing/'
    download_path = seq_home + 'current_set/'
    crunch_path = seq_home + 'crunch/'
    input_file = anno_home + 'coad-read.file_info'
    output_file = anno_home + '_'.join(['coad-read_set', set_num, dt + '.file_set'])
    error_file = anno_home + '_'.join(['coad-read_set', set_num, dt + '.file_err'])
    subject_dict = get_file_info(input_file)
    problems = process_file_info(output_file,subject_dict,download_path,crunch_path)
    record_errors(error_file,problems)
    output_current = anno_home + 'coad-read_current.file_set'
    error_current = anno_home + 'coad-read_current.file_err'
    set_current(output_file,output_current)
    set_current(error_file,error_current)


parser = argparse.ArgumentParser(description='Assemble downloaded files')
parser.add_argument('--datetime',help='date and time of execution',action='store',dest='dt')
parser.add_argument('--set_num',help='sample set',action='store',dest='set_num')
args = parser.parse_args()
main(args)




## to move files back

# for subject_id in subject_list:
#     sub = subject_set[subject_id]
#     for tissue_type in ('normal','tumor'):
#         if tissue_type in sub:
#             flag = 'have_' + tissue_type
#             for file in sub[tissue_type]:
#                 if os.path.exists(crunch_path + file['name']):
#                     if flag not in sub:
#                         sub[flag] = file
#     if 'have_normal' in sub and 'have_tumor' in sub:
#         cmd_n = 'mv %s %s/' % (crunch_path + sub['have_normal']['name'], download_path + sub['have_normal']['uuid'])
#         cmd_t = 'mv %s %s/' % (crunch_path + sub['have_tumor']['name'], download_path + sub['have_tumor']['uuid'])
#         subprocess.call(cmd_n, shell=True)
#         subprocess.call(cmd_t, shell=True)


## to move bai files in case miss that part

# for subject_id in subject_list:
#     sub = subject_set[subject_id]
#     for tissue_type in ('normal','tumor'):
#         if tissue_type in sub:
#             flag = 'have_' + tissue_type
#             for file in sub[tissue_type]:
#                 if os.path.exists(crunch_path + file['name']):
#                     if flag not in sub:
#                         sub[flag] = file
#     if 'have_normal' in sub and 'have_tumor' in sub:
#         cmd_n_bai = 'mv %s/%s %s' % (download_path + sub['have_normal']['uuid'], sub['have_normal']['name'].split('.')[0] + '.bai', crunch_path)
#         cmd_t_bai = 'mv %s/%s %s' % (download_path + sub['have_tumor']['uuid'], sub['have_tumor']['name'].split('.')[0] + '.bai', crunch_path)
#         subprocess.call(cmd_n_bai, shell=True)
#         subprocess.call(cmd_t_bai, shell=True)
