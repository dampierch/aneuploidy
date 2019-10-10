#!/usr/bin/env python3

## count_hetalleles.py :: input raw bam + normal hetsites bed :: output TSV

    # usage: count_hetalleles.py --bam [normal/tumor].bam --hetsites_bed normal_hetsites.bed > subject_[normal/tumor]_hetcnts.tsv 2> subject_[normal/tumor].err_cnts

    # this script takes a bam file and a bed file with hetsites in the normal
    # exome sequences and extracts from the bam allele counts at each hetsite;
    # then writes allele counts at each hetsite to file

    # stderr reports locations with known heterozygotes < 80% of reads


import pysam
import sys
import argparse
import re


bases = ('A','C','G','T')
field_names = ('chrom','start','stop','type','score','strand','info')


def get_known_alleles(info_str):
    '''
    info_str: AF=0.500;AD::A=100,C=50
    '''
    known_alleles = []
    known_dict = {}
    for info in info_str.split(';'):
        if re.match('AD::',info):
            allele_cnts = info.split('::')[1].split(',')
            allele1 = allele_cnts[0].split('=')[0]
            allele2 = allele_cnts[1].split('=')[0]
            known_alleles.append(allele1)
            known_dict[allele1] = True
            known_alleles.append(allele2)
            known_dict[allele2] = True
    return known_alleles, known_dict


def get_alleles(bamfile, chrom, zero_based_pos):
    '''
    Return a list of alleles observed in the reads aligned at a locus of
    interest defined by chrom and one_based_pos
    '''
    alleles_at_pos = []
    for pileupcolumn in bamfile.pileup(chrom, zero_based_pos, zero_based_pos + 1):
        for pileupread in pileupcolumn.pileups:
            ## pysam includes all positions from all alignments that
            ## overlap the locus of interest, so you have to exclude
            ## loci that you don't care about
            if pileupcolumn.pos != zero_based_pos:
                continue
            ## query allele is None if is_del or is_refskip is set.
            if not pileupread.is_del and not pileupread.is_refskip:
                aligned_read_allele = pileupread.alignment.query_sequence[pileupread.query_position]
                alleles_at_pos.append(aligned_read_allele)
    return alleles_at_pos


def allele_counts(alleles, known_alleles):
    '''
    Count alleles where alleles is a list of A,C,G,T
    '''
    counts = {x:0 for x in bases}  ## reset counter
    known_total = total = 0
    for allele in alleles:
        if allele in counts:
            counts[allele] = counts[allele] + 1
            total = total + 1
        if allele in known_alleles:
            known_total = known_total + 1
    return (counts, known_total, total)


parser = argparse.ArgumentParser(description='Count alleles at heterozygouse sites')
parser.add_argument('--bam',help='target bam file',action='store',dest='bam_file')
parser.add_argument('--hetsites_bed',help='positions of hetsites',action='store',dest='hetsites_bed_file')
parser.add_argument('--min_count',help='minimum number of total reads at given position',action='store',dest='min_count',type=int,default=100)
parser.add_argument('--min_score',help='minimum score value for given position',action='store',dest='min_score',type=float,default=100.) ## wrp originally had 200 here but 100 in find_hetsites
args = parser.parse_args()


## open bamfile for reading bases at location
bamfile = pysam.AlignmentFile(args.bam_file,'rb')

## scan bed file of heterozygous sites
with open(args.hetsites_bed_file,'r') as in_f:
    for in_line in in_f:
        in_line = in_line.strip('\n')
        ## get het site positions
        data = dict(list(zip(field_names,in_line.split('\t'))))
        ## get known alleles, func returns a list with the major, then minor allele, but not the counts
        known_alleles, known_dict = get_known_alleles(data['info'])
        ## get observed alleles in bam, func returns a list of alleles at given chrom + position
        alleles = get_alleles(bamfile,data['chrom'],int(data['start']))
        ## get observed allele counts
        (counts, known_total, all_total) = allele_counts(alleles,known_dict)

        ## for x in known_alleles only returns counts for alleles in vcf reference file
        ## does not report when there are large numbers of other counts -- should probably check/report that

        if known_total/all_total > 0.8:
            allele_str = "SD::" + ','.join(["%s=%d"%(x,counts[x]) for x in known_alleles])
            print(";".join((in_line,allele_str)))
        else:
            allele_str = "BAD::" + ','.join(["%s=%d"%(x,counts[x]) for x in bases])
            sys.stderr.write('%s;%s ***\n'%(in_line,allele_str))

bamfile.close()


### testing
# bamfile = pysam.AlignmentFile('/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/TCGA-AF-3400-01A-01D-1989-10_gapfillers_Illumina_gdc_realn.bam','rb')
# bamfile = pysam.AlignmentFile('/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/TCGA-AF-3400-11A-01D-1554-10_Illumina_gdc_realn.bam','rb')
# in_line = [data[field] for field in out_fields]
# d1 = dict(list(zip(field_names,in_line)))
# known_alleles, known_dict = get_known_alleles(d1['info'])
# alleles = get_alleles(bamfile,d1['chrom'],int(d1['start']))
#
#
# limit = 30
# it = 0
# for pileupcolumn in bamfile.pileup(d1['chrom'], int(d1['start']), int(d1['start']) + 1):
#     for pileupread in pileupcolumn.pileups:
#         # while it < limit:
#         #     print(pileupcolumn.pos)
#         if pileupcolumn.pos == int(d1['start']):
#             it = it + 1
#
#
#
# f = pysam.AlignmentFile('/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/TCGA-AF-3400-01A-01D-1989-10_gapfillers_Illumina_gdc_realn.bam', 'rb')
# h = f.head(n=5, multiple_iterators=True)
# for read in h:
#     print(read)
# f.close()
