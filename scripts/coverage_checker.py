#!/usr/bin/env python3

## coverage_checker.py :: input raw bam + reference exome bed :: output coverage bed with observed counts at each position in ref exome

    # usage: coverage_checker.py --bam [normal/tumor].bam --refex_bed reference_exome.bed > subject_[normal/tumor]_coverage.bed

    # this script takes a bam file and a bed file and extracts from the bam
    # the exome coverage at refex_bed-specified intervals; finds coverage
    # across entire interval and then takes average based on interval length;
    # then writes avg coverage for each interval as well as average coverage
    # over entire reference exome to output bedfile


import pysam
import argparse


in_fields = ('chrom','start','stop')
out_fields = ('chrom','start','stop','length','total_coverage','mean_coverage')


def get_read_depth(bamfile, chrom, zero_based_start, zero_based_stop):
    '''
    return a list of depths at each position along a given interval

    the pileup call initiates an iterator over each position (i.e. column)
    in the set of positions (columns) composed of all positions in all reads
    that overlap each position in the interval
    '''
    depth_at_pos = []
    for pileupcolumn in bamfile.pileup(chrom, zero_based_start, zero_based_stop):
        if pileupcolumn.reference_pos >= zero_based_start and pileupcolumn.reference_pos <= zero_based_stop:
            depth_at_pos.append(len(pileupcolumn.pileups))
    return depth_at_pos


parser = argparse.ArgumentParser(description='Count reads at each position in exome intervals and take average')
parser.add_argument('--bam',help='target bam file',action='store',dest='bam_file')
parser.add_argument('--refex_bed',help='intervals of exome capture',action='store',dest='refex_bed')
args = parser.parse_args()


## open bamfile for getting depth at each position in exome
bamfile = pysam.AlignmentFile(args.bam_file,'rb')


## scan exome for depth of coverage
with open(args.refex_bed,'r') as in_f:
    print('\t'.join(out_fields))
    exome_reads = []
    exome_length = []
    for in_line in in_f:
        if in_line[0:5] == 'track' or in_line[0:5] == 'brows':
            continue
        in_line = in_line.strip('\n')
        ## get interval positions
        interval = dict(list(zip(in_fields,in_line.split('\t'))))
        ## get interval length
        interval_length = int(interval['stop']) - int(interval['start']) + 1
        ## get observed reads at each position in interval
        interval_reads = get_read_depth(bamfile,interval['chrom'],int(interval['start']),int(interval['stop']))
        ## get average coverage over each interval
        interval_mean = round(sum(interval_reads)/interval_length)
        ## tally all reads across exome
        exome_reads = exome_reads + interval_reads
        exome_length.append(interval_length)
        ## report depth
        out_str = '\t'.join([str(interval_length),str(sum(interval_reads)),str(interval_mean)])
        print('\t'.join([in_line,out_str]))
    print('whole exome coverage = %sx' % (round(sum(exome_reads)/sum(exome_length))))


## close bamfile once coverage calculated
bamfile.close()




### testing
# bam = '/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/TCGA-T9-A92H-10A-01D-A370-10_Illumina_gdc_realn.bam'
# refex_bed = '/scratch/chd5n/aneuploidy/exome-kits/SureSelect_CRE_V2_hg38/S30409818_Regions.bed'
# arg_str = ' '.join(['--bam',bam,'--refex_bed',refex_bed])
# args = parser.parse_args(arg_str.split())
# bamfile = pysam.AlignmentFile(args.bam_file,'rb')
#
# counter = 0
# with open(args.refex_bed,'r') as in_f:
#     print('\t'.join(out_fields))
#     exome_reads = []
#     exome_length = []
#     for in_line in in_f:
#         if in_line[0:5] == 'track' or in_line[0:5] == 'brows':
#             continue
#         in_line = in_line.strip('\n')
#         ## get interval positions
#         interval = dict(list(zip(in_fields,in_line.split('\t'))))
#         ## get interval length
#         interval_length = int(interval['stop']) - int(interval['start']) + 1
#         ## get observed reads at each position in interval
#         interval_reads = get_read_depth(bamfile,interval['chrom'],int(interval['start']),int(interval['stop']))
#         ## get average coverage over each interval
#         interval_mean = round(sum(interval_reads)/interval_length)
#         ## tally all reads across exome
#         exome_reads = exome_reads + interval_reads
#         exome_length.append(interval_length)
#         ## report depth
#         out_str = '\t'.join([str(interval_length),str(sum(interval_reads)),str(interval_mean)])
#         print('\t'.join([in_line,out_str]))
#         if counter > 10:
#             break
#         counter = counter + 1
#     print('whole exome coverage = %sx' % (round(sum(exome_reads)/sum(exome_length))))
# bamfile.close()
