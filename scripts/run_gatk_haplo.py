#!/usr/bin/env python3

## run_gatk_haplo.py

    # this script takes tab-delim output from gdc_assemble_files.py (*.file_set)
    # and extracts normal file_name, then submits each normal raw bam to
    # gatk_haplo.bash via sbatch subprocess

    # *.file_set fields: subject_id, file_name normal, file_name tumor


import os
import subprocess
import argparse


def sbatch_command(script_home,input_file,crunch_path):
    os.chdir(script_home)
    sub_list = []
    with open(input_file,'r') as in_f:
        for in_line in in_f:
            subject_id, file_name_norm, file_name_tum = in_line.strip('\n').split('\t')
            output_setting = '--output=gatk_haplo_' + subject_id + '.out'
            cmd = ' '.join(['sbatch', output_setting, 'gatk_haplo.sh', crunch_path + file_name_norm])
            subprocess.call(cmd, shell=True)
            sub_list.append(subject_id)
    return sub_list


def write_status(output_file,sub_list):
    with open(output_file,'w+') as out_f:
        out_f.write('\n'.join(sub_list))


def main(input_file):
    home = os.environ['HOME']
    script_home = home + '/projects/aneuploidy/scripts/'
    data_home = '/scratch/chd5n/aneuploidy/'
    seq_home = data_home + 'raw-data/sequencing/'
    crunch_path = seq_home + 'crunch/'
    sub_list = sbatch_command(script_home,input_file,crunch_path)
    output_file = seq_home + 'latest_gatk_haplo.txt'
    write_status(output_file,sub_list)


parser = argparse.ArgumentParser(description='SLURM calls for GATK HaplotypeCaller')
parser.add_argument('--in_file',help='file_set for variant calling',action='store',dest='input_file')
args = parser.parse_args()
main(args.input_file)
