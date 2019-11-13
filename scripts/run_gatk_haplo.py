#!/usr/bin/env python3

## run_gatk_haplo.py

    # this script takes tab-delim output from assemble_gdc_files.py (*.file_set)
    # and extracts normal file_name, then submits each normal raw bam to
    # gatk_haplo.bash via sbatch subprocess

    # *.file_set fields: subject_id, file_name normal, file_name tumor

import os
import fileinput
import subprocess

home = os.environ['HOME']
script_home = home + '/projects/aneuploidy/scripts'
data_home = '/scratch/chd5n/aneuploidy/'
anno_home = data_home + 'raw-data/annotations/'
seq_home = data_home + 'raw-data/sequencing/'
crunch_path = seq_home + 'crunch/'

input_file = anno_home + 'coad-read_2019-09-26.file_set'

os.chdir(script_home)
with open(input_file,'r') as in_f:
    for in_line in in_f:
        (subject_id, file_name_norm, file_name_tum) = in_line.strip('\n').split('\t')
        output_setting = '--output=gatk_haplo_' + subject_id + '.out'
        cmd = ' '.join(['sbatch', output_setting, 'gatk_haplo.sh', crunch_path + file_name_norm])
        subprocess.call(cmd, shell=True)
