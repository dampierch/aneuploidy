#!/usr/bin/env python3


  ## this script downloads files from gdc with gdc-client using file uuid;
  ## executes download on front end


import subprocess
from datetime import datetime
import os
import argparse


def get_uuids(in_file):
    '''
    read all uuids into list
    '''
    with open(in_file,'r') as in_f:
        uuids = []
        for in_line in in_f:
            uuids.append(in_line.strip('\n').split('\t')[-1])
    return uuids


def parse_uuids(uuids,set_size):
    '''
    break uuids into manageable sets for download
    '''
    uuid_sets = []
    uuid_set = []
    for i in range(0, len(uuids)):
        if i == 0:
            uuid_set.append(uuids[i])
        elif i % set_size == 0:
            uuid_sets.append(uuid_set)
            uuid_set = [uuids[i]]
        elif i == (len(uuids)-1):
            uuid_set.append(uuids[i])
            uuid_sets.append(uuid_set)
        else:
            uuid_set.append(uuids[i])
    return uuid_sets


def download_files(uuid_sets,set_num,gdc_home,dest,token):
    '''
    download files specified in uuid_set on frontend
    '''
    print(' '.join(['download attempt start',str(datetime.now())]))
    file_count = 0
    for uuid in uuid_sets[set_num]:
        cmd = ' '.join([gdc_home + 'gdc-client download',uuid,'-d',dest,'-t',token])
        subprocess.call(cmd, shell=True)
        file_count = file_count + 1
        print(uuid + ' attempted')
    print(str(file_count) + ' files attempted')
    print(' '.join(['download attempt end',str(datetime.now())]))


def write_downloaded(downloaded_file,uuid_sets,set_num):
    '''
    append to record the uuid set downloaded
    '''
    with open(downloaded_file,'w+') as out_file:
        line1 = ' '.join(['most recent set number downloaded:',str(set_num)]) + '\n'
        line2 = ' '.join(uuid_sets[set_num]) + '\n'
        out_file.write(line1+line2)


def main(set_num):
    '''
    parse uuids and download selected set
    '''
    ## set variables
    home = os.environ['HOME']
    gdc_home = home + '/Apps/gdc-client/'
    aneuploidy_home = '/scratch/chd5n/aneuploidy/'
    anno_home = aneuploidy_home + 'raw-data/annotations/'
    seq_home = aneuploidy_home + 'raw-data/sequencing/'
    in_file = anno_home + 'coad-read.file_info'
    dest = seq_home + 'current_set/'
    token = seq_home + 'gdc-user-token.2019-11-20T23_21_41.960Z.txt'
    set_size = 100
    set_num = int(set_num)
    downloaded_file = seq_home + 'latest_download.txt'
    ## take actions
    uuids = get_uuids(in_file)
    uuid_sets = parse_uuids(uuids,set_size)
    download_files(uuid_sets,set_num,gdc_home,dest,token)
    write_downloaded(downloaded_file,uuid_sets,set_num)


parser = argparse.ArgumentParser(description='Download set of UUIDs')
parser.add_argument('--set_num',help='uuid set number',action='store',dest='set_num')
args = parser.parse_args()
main(args.set_num)
