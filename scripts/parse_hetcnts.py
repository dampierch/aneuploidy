## README
## aim:
##   -- calculate major allele fraction and then filter, sort, and intercalate
##      allele counts at heterozygous sites for tumor and normal samples
##   -- filter out sites with zero expected allele counts or that are absent from
##      one tissue type
## input:
##   -- two bed files with simple allele counts at sites of interest, one for
##      each tissue type (i.e. normal, tumor)
## output:
##   -- one tabular (i.e. tab separated values) file with filtered and sorted
##      counts for sites in both normal and tumor samples
## sketch:
##   -- parse_hetcnts.py < [norm/tumor]_hetcnts.bed > hetcnts.tsv
## usage:
##   -- intended to be called from within het_counter.sh
##   -- parse_hetcnts.py \
##          --subject_id ${subject_id} \
##          --normal_hetcnts_bed ${crunch_dir}${subject_id}_normal_hetcnts.bed \
##          --tumor_hetcnts_bed ${crunch_dir}${subject_id}_tumor_hetcnts.bed
##   -- currently a stand-alone rule that processes all samples at once
##   -- python parse_hetncts.py --info {input.info}


import gzip
import argparse


modname = 'parse_hetcnts'


def read_info(target):
    ''' read list of subject_ids from info file and make bed filenames '''
    target_path = ''.join(
        ['/'.join(target.split('/')[0:-1]), '/hetsites-data/']
    )
    sufx = '_hetcnts.bed.gz'
    l = []
    with open(target, 'r') as infile:
        for inline in infile:
            l.append(inline.strip('\n').split('\t')[0])
    l_id = sorted(set(l))
    d_norm = {i: ''.join([target_path, i, '_normal', sufx]) for i in l_id}
    d_tumor = {i: ''.join([target_path, i, '_tumor', sufx]) for i in l_id}
    d = {'subject_id': l_id, 'bed_norm': d_norm, 'bed_tumor': d_tumor}
    return d


def read_chroms():
    ''' read into dict absolute start position of each chromosome '''
    infields = ['chrom', 'length', 'pos0']
    ref_home = '/scratch/chd5n/reference-genome/assembly/tcga/'
    target = ''.join([ref_home, 'GRCh38.d1.vd1.chr1-XY.size.tsv'])
    d = {}
    with open(target, 'r') as infile:
        for inline in infile:
            data = dict(list(zip(infields, inline.strip('\n').split('\t'))))
            d[data['chrom']] = data['pos0']
    return d


def parse_infostr(info_str):
    '''
    parses the info string to return major and minor allele (simple) depth at
    given position
    '''
    af, ad, sd = info_str.split(';')
    a0, a1 = sd.split('::')[1].split(',')
    major = a0.split('=')
    minor = a1.split('=')
    return major, minor


def read_hetcnts(subject_id, target, tissue_type, d_chrom):
    in_fields = ['chrom', 'pos0', 'pos', 'an_type', 'score', 'strand',
        'info_str']
    d = {}
    with gzip.open(target, 'rt') as in_f:
        for in_line in in_f:
            data = dict(list(zip(in_fields, in_line.strip('\n').split('\t'))))
            data['t_type'] = tissue_type
            data['subject_id'] = subject_id
            data['pos_abs'] = int(d_chrom[data['chrom']]) + int(data['pos'])
            major, minor = parse_infostr(data['info_str'])
            data['maj'] = major[0]
            data['maj_cnt'] = major[1]
            data['min'] = minor[0]
            data['min_cnt'] = minor[1]
            dp = int(data['maj_cnt']) + int(data['min_cnt'])
            if dp == 0:
                continue
            else:
                data['maj_fract'] = '%.5f' % (int(data['maj_cnt'])/dp)
            d[data['pos_abs']] = data
    return d


def write_output(subject_id, d_norm, d_tumor):
    outfields = ['subject_id', 'chrom', 'pos', 't_type', 'maj', 'maj_cnt',
        'min', 'min_cnt', 'maj_fract']
    target_path = '/scratch/chd5n/aneuploidy/hetcnts/'
    target = ''.join([target_path, subject_id, '_hetcnts.tsv'])
    with open(target, 'w') as outfile:
        outfile.write(''.join(['\t'.join(outfields), '\n']))
        ## confirm that both sets of data have same coordinates
        for k in sorted(d_norm.keys()):
            s1 = ''.join(
                ['\t'.join([d_norm[k][f] for f in outfields]), '\n']
            )
            s2 = ''.join(
                ['\t'.join([d_tumor[k][f] for f in outfields]), '\n']
            )
            outfile.write(''.join([s1, s2]))
    print('file written to', target)


def main(args):
    d = read_info(args.info)
    d_chrom = read_chroms()
    for subject_id in d['subject_id']:
        bed_norm = d['bed_norm'][subject_id]
        bed_tumor = d['bed_tumor'][subject_id]
        ## prepare data (normal or tumor) as dict of dicts
        ## keys are absolute genomic position (for filter and sort)
        ## values are dicts of data
        d_norm = read_hetcnts(subject_id, bed_norm, 'norm', d_chrom)
        d_tumor = read_hetcnts(subject_id, bed_tumor, 'tumor', d_chrom)
        ## filter out unpaired sites
        d_norm = {k: d_norm[k] for k in d_norm if k in d_tumor}
        d_tumor = {k: d_tumor[k] for k in d_tumor if k in d_norm}
        write_output(subject_id, d_norm, d_tumor)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='filter and sort data')
    parser.add_argument(
        '--info', help='file with list of subject_ids to be analyzed',
        action='store', dest='info', default=False
    )
    args = parser.parse_args()
    main(args)
else:
    print('functions loaded for', modname)
