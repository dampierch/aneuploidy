#!/usr/bin/env python3

## assemble_gdc_files.py

    # this script takes tab-delim output from parse_gdc_info.py (*.file_info)
    # and extracts subject_id, tissue_type, file_name, file_id and reformats
    # into *.file_set and also moves files to data crunch directory

    # *.file_set has three fields: subject_id, file_name_normal, file_name_tumor
    # file='/home/chd5n/projects/aneuploidy/scripts/assemble_gdc_files.py'
    # exec(open(file).read())

import fileinput
import os
import subprocess
import datetime

dt = str(datetime.date.today())
data_home = '/scratch/chd5n/aneuploidy/'
anno_home = data_home + 'raw-data/annotations/'
seq_home = data_home + 'raw-data/sequencing/'

input_file = anno_home + 'coad-read.file_info'
download_path = seq_home + 'first_10_pairs/'
crunch_path = seq_home + 'crunch/'

output_file = anno_home + 'coad-read_' + dt + '.file_set'
error_file = anno_home + 'coad-read_' + dt + '.file_errors'

subject_list = []
subject_set = {}

## subject_set
  ## for each subject_id as key, create values of normal and tumor
  ## for each value, create dict of file name and file uuid
  ## note that dict is modified indirectly with modifications to values

with fileinput.input(files=(input_file)) as in_f:
    for line in in_f:
        if (line[0] == '#'):
            continue
        (subject_id, tissue_type, file_name, file_id) = line.strip('\n').split('\t')
        if subject_id in subject_set:
            sub = subject_set[subject_id]
        else:
            subject_set[subject_id] = {}
            sub = subject_set[subject_id]
            subject_list.append(subject_id)
        if tissue_type in sub:
            sub[tissue_type].append({'name': file_name, 'uuid': file_id}) ## this code modifies subject_set
        else:
            sub[tissue_type] = [{'name': file_name, 'uuid': file_id}] ## this code modifies subject_set

## now check for pairs and existence of files
  ## the pairs were checked by the modified parse_gdc_info script but no harm in
  ## checking again
  ## if pairs and files ok, then write .file_set

problems = []
with open(output_file, 'w+') as out_f:
    for subject_id in subject_list:
        sub = subject_set[subject_id]
        ## check to see if files exist
        for tissue_type in ('normal','tumor'):
            if tissue_type in sub:
                flag = 'have_' + tissue_type
                for file in sub[tissue_type]:
                    if os.path.exists(download_path + file['uuid']):
                        if os.path.exists('%s/%s' % (download_path + file['uuid'], file['name'])):
                            if flag not in sub:
                                sub[flag] = file
        ## make sure both halves of pair are available
        ## if so, move them to data crunch path
        if 'have_normal' in sub and 'have_tumor' in sub:
            out_f.write('\t'.join([subject_id, sub['have_normal']['name'], sub['have_tumor']['name']]) + '\n')
            cmd_n = 'mv %s/%s %s' % (download_path + sub['have_normal']['uuid'], sub['have_normal']['name'], crunch_path)
            cmd_t = 'mv %s/%s %s' % (download_path + sub['have_tumor']['uuid'], sub['have_tumor']['name'], crunch_path)
            subprocess.call(cmd_n, shell=True)
            subprocess.call(cmd_t, shell=True)
        else:
            if 'have_normal' not in sub:
                problems.append('%s missing normal\n' % (subject_id))
            if 'have_tumor' not in sub:
                problems.append('%s missing tumor\n' % (subject_id))

with open(error_file, 'w+') as out_f:
    for sub in problems:
        out_f.write(sub)



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
