#!/usr/bin/env python3

## store_hetsite_data.py

    ## this script takes selected output files from hetsite counter pipeline
    ## and stores them in hetsites-data directory for future use


import subprocess
import argparse
import glob


def get_filenames(input_file,crunch_path):
    '''
    read subject_ids into list and get selected filenames with glob; for
    hetsites.bed and vcf files, need to use _gdc_realn prefix
    '''
    subject_ids = []
    prefixes = []
    with open(input_file) as in_f:
        for line in in_f:
            if line[0] == '#':
                continue
            subject_id, file_name_norm, file_name_tum = line.strip('\n').split('\t')
            subject_ids.append(subject_id)
            prefixes.append(file_name_norm.split('.')[0])
            prefixes.append(file_name_tum.split('.')[0])
    hts_beds = []
    vcfs = []
    trash = []
    cts_beds = []
    rcnts = []
    for pre in prefixes:
        hts_fn = glob.glob(crunch_path + pre + '*hetsites.bed')
        vcf_fn = glob.glob(crunch_path + pre + '*snp.indel.vcf')
        hts_beds.append(hts_fn)
        vcfs.append(vcf_fn)
        trs_fn_list = []
        for suf in ['bai','bam','idx']:
            trs_fn_list = trs_fn_list + glob.glob(crunch_path + pre + '*' + suf)
        trash.append(trs_fn_list)
    for id in subject_ids:
        bed_fn_list = []
        for suf in ['errcnts.bed','hetcnts.bed']:
            bed_fn = glob.glob(crunch_path + id + '*' + suf)
            bed_fn_list = bed_fn_list + bed_fn
        rcnt_fn_list = []
        for suf in ['cnts2R.tsv','R_missing.err']:
            rcnt_fn = glob.glob(crunch_path + id + '*' + suf)
            rcnt_fn_list = rcnt_fn_list + rcnt_fn
        cts_beds.append(bed_fn_list)
        rcnts.append(rcnt_fn_list)
    return hts_beds,vcfs,trash,cts_beds,rcnts,subject_ids


def store_files(storage_path,hts_beds,vcfs,trash,cts_beds,rcnts):
    '''
    move selected output files to storage directory and remove selected trash
    files
    '''
    for set in hts_beds:
        for file in set:
            cmd = 'mv %s %s%s/' % (file,storage_path,'bed-files')
            subprocess.call(cmd, shell=True)
    for set in vcfs:
        for file in set:
            cmd = 'mv %s %s%s/' % (file,storage_path,'vcf-files')
            subprocess.call(cmd, shell=True)
    for set in trash:
        for file in set:
            cmd = 'rm %s' % (file)
            subprocess.call(cmd, shell=True)
    for set in cts_beds:
        for file in set:
            cmd = 'mv %s %s%s/' % (file,storage_path,'bed-files')
            subprocess.call(cmd, shell=True)
    for set in rcnts:
        for file in set:
            cmd = 'mv %s %s%s/' % (file,storage_path,'r-cnts')
            subprocess.call(cmd, shell=True)


def write_status(output_file,subject_ids):
    with open(output_file,'w+') as out_f:
        out_f.write('\n'.join(subject_ids))


def set_latest(original,revised):
    cmd = ' '.join(['cp',original,revised])
    subprocess.call(cmd, shell=True)


def main(args):
    input_file = args.input_file
    dt = args.dt
    set_num = args.set_num
    data_home = '/scratch/chd5n/aneuploidy/'
    seq_home = data_home + 'raw-data/sequencing/'
    crunch_path = seq_home + 'crunch/'
    storage_path = '/scratch/chd5n/aneuploidy/hetsites-data/'
    output_file = storage_path + '_'.join(['coad-read_set', set_num + 'v3', dt + '.txt'])
    output_latest = storage_path + 'latest_data.txt'
    hts_beds,vcfs,trash,cts_beds,rcnts,subject_ids = get_filenames(input_file,crunch_path)
    store_files(storage_path,hts_beds,vcfs,trash,cts_beds,rcnts)
    write_status(output_file,subject_ids)
    set_latest(output_file,output_latest)


parser = argparse.ArgumentParser(description='Move hetsites data to storage')
parser.add_argument('--in_file',help='input file with set information',action='store',dest='input_file')
parser.add_argument('--datetime',help='date and time of execution',action='store',dest='dt')
parser.add_argument('--set_num',help='sample set',action='store',dest='set_num')
args = parser.parse_args()
main(args)
