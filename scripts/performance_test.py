#!/usr/bin/env python3


## count_hetalleles.py test
    ## not much difference; bigger difference due to which called first


import pysam
import re
import time


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


def get_alleles_chd(bamfile, chrom, zero_based_pos):
    alleles_at_pos = []
    for pileupcolumn in bamfile.pileup(chrom, zero_based_pos, zero_based_pos + 1):
        if pileupcolumn.reference_pos != zero_based_pos:
            continue
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                aligned_read_allele = pileupread.alignment.query_sequence[pileupread.query_position]
                alleles_at_pos.append(aligned_read_allele)
    return alleles_at_pos


def get_alleles_wrp(bamfile, chrom, zero_based_pos):
    alleles_at_pos = []
    for pileupcolumn in bamfile.pileup(chrom, zero_based_pos, zero_based_pos+1):
        for pileupread in pileupcolumn.pileups:
            if pileupcolumn.pos != zero_based_pos:
                continue
            if not pileupread.is_del and not pileupread.is_refskip:
                aligned_read_allele = pileupread.alignment.query_sequence[pileupread.query_position]
                alleles_at_pos.append(aligned_read_allele)
    return alleles_at_pos


bamfile = pysam.AlignmentFile('/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/TCGA-T9-A92H-10A-01D-A370-10_Illumina_gdc_realn.bam','rb')


start = time.time()
with open('/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/TCGA-T9-A92H-10A-01D-A370-10_Illumina_gdc_realn_hetsites.bed') as in_f:
    counter = 0
    for in_line in in_f:
        counter = counter + 1
        in_line = in_line.strip('\n')
        ## get het site positions
        data = dict(list(zip(field_names,in_line.split('\t'))))
        ## get known alleles, func returns a list with the major, then minor allele, but not the counts
        known_alleles, known_dict = get_known_alleles(data['info'])  ## returns known_alleles, known_dict
        ## WRP
        alleles_wrp = get_alleles_wrp(bamfile,data['chrom'],int(data['start']))

end = time.time()
runtime = end - start
print('with get_alleles_wrp: %s sec' % (runtime))


start = time.time()
with open('/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/TCGA-T9-A92H-10A-01D-A370-10_Illumina_gdc_realn_hetsites.bed') as in_f:
    counter = 0
    for in_line in in_f:
        counter = counter + 1
        in_line = in_line.strip('\n')
        ## get het site positions
        data = dict(list(zip(field_names,in_line.split('\t'))))
        ## get known alleles, func returns a list with the major, then minor allele, but not the counts
        known_alleles, known_dict = get_known_alleles(data['info'])  ## returns known_alleles, known_dict
        ## CHD
        alleles_chd = get_alleles_chd(bamfile,data['chrom'],int(data['start']))

end = time.time()
runtime = end - start
print('with get_alleles_chd: %s sec' % (runtime))


bamfile.close()


print(len(alleles_chd))
print(len(alleles_wrp))
