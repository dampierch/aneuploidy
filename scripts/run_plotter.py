#!/usr/bin/env python3

## run_plotter.py

    # this script sends either name of tab-delim output file from
    # assemble_gdc_files.py (*.file_set) or single subject ID to plotting
    # scripts via Rscript subprocess

    # *.file_set fields: subject_id, file_name normal, file_name tumor


import os
import subprocess
import argparse


def shell_cmd(script_home,single_subject,subject_id,file_set):
    os.chdir(script_home)
    if single_subject:
        cmd = ' '.join(['module load gcc/7.1.0 openmpi/3.1.4 R/3.6.0; Rscript', 'plot_haplo_chrom.R', '--args', subject_id])
    else:
        cmd = ' '.join(['module load gcc/7.1.0 openmpi/3.1.4 R/3.6.0; Rscript', 'plot_haplo_chrom.R', '--args', file_set])
    subprocess.call(cmd, shell=True)


def copy_current(original,revised):
    cmd = ' '.join(['cp',original,revised])
    subprocess.call(cmd, shell=True)


def main(args):
    file_set = args.input_file
    single_sub = args.single_sub
    subject_id = args.subject_id
    dt = args.dt
    set_num = args.set_num
    home = os.environ['HOME']
    script_home = home + '/projects/aneuploidy/scripts/'
    shell_cmd(script_home,single_sub,subject_id,file_set)
    plot_dir = '/scratch/chd5n/aneuploidy/results/plots/'
    fn_before = plot_dir + 'coad-read_current_hetcnts_plot.pdf'
    fn_after = plot_dir + '_'.join(['coad-read_set', set_num, dt, 'hetcnts_plot.pdf'])
    copy_current(fn_before,fn_after)


parser = argparse.ArgumentParser(description='shell calls for heterozygous site plots in R')
parser.add_argument('--in_file',help='file_set for hetsite ploting',action='store',dest='input_file')
parser.add_argument('--single_sub',help='boolean for whether to make plot for single subject',default=False,action='store',dest='single_sub')
parser.add_argument('--subject_id',help='if single_sub True, use this subject_id to make single subject plot',action='store',dest='subject_id')
parser.add_argument('--datetime',help='date and time of execution',action='store',dest='dt')
parser.add_argument('--set_num',help='sample set',action='store',dest='set_num')
args = parser.parse_args()
main(args)
