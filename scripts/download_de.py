#!/usr/bin/env python3

## download_de.py :: input samples for rna-seq analysis :: output count matrices
    ## usage: nohup python download_de.py --set_name {set_name} > ../logs/download_de_set_{set_name}_{dt}.out 2> ../logs/download_de_set_{set_name}_{dt}.err
    ## this script downloads files from gdc with gdc-client using file uuid;
    ## executes download on front end


import subprocess
from datetime import datetime
import os
import argparse
import glob


def get_uuids(in_file):
    '''
    read all uuids into list
    '''
    with open(in_file,'r') as in_f:
        uuids = []
        header = True
        for in_line in in_f:
            if header:
                header = False
                continue
            uuids.append(in_line.strip('\n').split('\t')[-3])
    return uuids


def download_files(uuids,gdc_home,dest):
    '''
    download files specified in uuids on frontend
    '''
    print(' '.join(['download attempt start',str(datetime.now())]))
    new_file_count = 0
    old_file_count = 0
    for uuid in uuids:
        if len(glob.glob(dest+uuid+'/*.htseq.counts')) == 1:  ## if already have downloaded and decompressed, skip
            old_file_count = old_file_count + 1
            continue
        cmd = ' '.join([gdc_home + 'gdc-client download',uuid,'-d',dest])
        subprocess.call(cmd, shell=True)
        # print(cmd)  ## for testing
        new_file_count = new_file_count + 1
        print(uuid + ' attempted')
    print(str(old_file_count) + ' files already downloaded and decompressed')
    print(str(new_file_count) + ' files attempted')
    print(' '.join(['download attempt end',str(datetime.now())]))


def write_downloaded(downloaded_file,uuids,set_name):
    '''
    append to record the uuid set downloaded
    '''
    with open(downloaded_file,'w+') as out_file:
        line1 = ' '.join(['most recent set downloaded:',set_name]) + '\n'
        line2 = ' '.join(uuids) + '\n'
        out_file.write(line1+line2)


def main(set_name):
    '''
    parse uuids and download files
    '''
    ## set variables
    home = os.environ['HOME']
    gdc_home = home + '/Apps/gdc-client/'
    aneuploidy_home = '/scratch/chd5n/aneuploidy/'
    anno_home = aneuploidy_home + 'raw-data/annotations/'
    count_home = aneuploidy_home + 'raw-data/counts/'
    in_file = anno_home + 'rna_set_' + set_name + '.tsv'
    dest = count_home
    downloaded_file = count_home + 'latest_download.txt'
    ## take actions
    uuids = get_uuids(in_file)
    download_files(uuids,gdc_home,dest)
    write_downloaded(downloaded_file,uuids,set_name)


parser = argparse.ArgumentParser(description='Download set of UUIDs')
parser.add_argument('--set_name',help='uuid set name',action='store',dest='set_name')
args = parser.parse_args()
main(args.set_name)
