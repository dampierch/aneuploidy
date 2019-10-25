#!/usr/bin/env python3

## run_plotter.py

    # this script sends either name of tab-delim output file from
    # assemble_gdc_files.py (*.file_set) or single subject ID to plotting
    # scripts via Rscript subprocess

    # *.file_set fields: subject_id, file_name normal, file_name tumor

import os
import subprocess

home = os.environ['HOME']
script_home = home + '/projects/aneuploidy/scripts'
data_home = '/scratch/chd5n/aneuploidy/'
anno_home = data_home + 'raw-data/annotations/'
seq_home = data_home + 'raw-data/sequencing/'
crunch_path = seq_home + 'crunch/'

file_set = 'coad-read_2019-09-26.file_set'
single_subject = False
subject_id = 'TCGA-AF-3400'

os.chdir(script_home)
if single_subject:
    cmd = ' '.join(['module load gcc/7.1.0 openmpi/3.1.4 R/3.6.0; Rscript', 'plot_haplo_chrom.R', '--args', subject_id])
else:
    cmd = ' '.join(['module load gcc/7.1.0 openmpi/3.1.4 R/3.6.0; Rscript', 'plot_haplo_chrom.R', '--args', file_set])
subprocess.call(cmd, shell=True)
