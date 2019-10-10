#!/usr/bin/env python3

## find_hetsites.py :: input VCF :: output BED

    # usage: find_hetsites.py *.vcf > *_hetsites.bed

    # this script takes tab-delim output from gatk_haplo.bash (*.vcf_L) and
    # extracts positions of SNP alt alleles with AF=0.500 && sufficient read depth;
    # writes bed with location, AF, alleles

import fileinput
import re

qual_thresh = 100
dp_thresh = 100
out_fields = ('CHROM', 'POS0', 'POS', 'name', 'QUAL', 'strand', 'info_str')

def allele_info(data, genotype, alt_only):
    '''
    this function takes data line as dict; takes genotype data from last
    field in line as string, converts to dict; returns string with each
    allele and associated counts
    '''
    genotype_keys = data['FORMAT'].split(':')
    genotype_values = genotype.split(':')
    allele_data = dict(list(zip(genotype_keys, genotype_values)))
    allele_cnts = allele_data['AD'].split(',')
    if alt_only:
        out_str = "AD::%s=%s,%s=%s" % (data['ALT'].split(',')[0],allele_cnts[1],data['ALT'].split(',')[1],allele_cnts[2])
    else:
        out_str = "AD::%s=%s,%s=%s" % (data['REF'],allele_cnts[0],data['ALT'],allele_cnts[1])
    return out_str

with fileinput.input() as in_f:
    for in_line in in_f:
        if in_line[0]=='#':  ## use this test with continue to skip processing leading comments
            if in_line[0:6]=='#CHROM':
                field_names = in_line.strip('\n').split('\t')
                field_names[0] = 'CHROM'  ## this gets rid of '#'
                last_field = field_names[-1]  ## variable name with subject_id in first three '-' delim fields
            continue  ## pops out to next iteration of for loop to skip processing leading comments
        data = dict(list(zip(field_names,in_line.strip('\n').split('\t'))))
        data['POS0'] = str(int(data['POS'])-1)
        data['name'] = 'vcf'
        data['strand'] = '.'
        data['info_str'] = ''
        ## skip non-snps
        if len(data['REF']) > 1:
            if re.search(r',',data['REF']):
                if len(data['REF'].split(',')[0]) > 1 or len(data['REF'].split(',')[1]) > 1:
                    continue
            else:
                continue
        if len(data['ALT']) > 1:
            if re.search(r',',data['ALT']):
                if len(data['ALT'].split(',')[0]) > 1 or len(data['ALT'].split(',')[1]) > 1:
                    continue
            else:
                continue
        ## skip low quality snps
        if float(data['QUAL']) < qual_thresh:
            continue
        good_dp = False
        good_af = False
        alt_only = False
        for info in data['INFO'].split(';'):
            (key,val) = info.split('=')
            if key=='DP' and int(val) > dp_thresh:
                good_dp = True
            if key=='AF' and (re.search(r',',val) or float(val)==0.5):
                if re.search(r',',val):
                    alt_only = True
                good_af = True
                data['info_str'] = "%s;%s" % (info,allele_info(data,data[last_field],alt_only))
        if good_dp and good_af:
            print('\t'.join([data[field] for field in out_fields]))


## for testing
# with fileinput.input('/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/TCGA-AF-3400-11A-01D-1554-10_Illumina_gdc_realn.snp.indel.vcf_L') as in_f:
