import gzip
from glob import glob
from pandas import read_csv, merge


modname = 'sum_sites'


def extract_value(l, count_type):
    if count_type == 'sites':
        v = len(l)
    elif count_type == 'chroms':
        v = len(set(l))
    else:
        v = -1
    return v


def count_data(count_type):
    '''
    parameters:
    -- count_type: either 'sites' or 'chroms'; sites counts total sites tested; chroms
    counts unique chromosomes for which sites are tested
    implementation:
    -- function reads 1761 files so takes 1-2 mins per count_type
    '''
    target_path = '/scratch/chd5n/aneuploidy/hetsites-data/'
    sufx = ['_cnts2R.tsv.gz', '_normal_hetcnts.bed.gz', '_tumor_hetcnts.bed.gz']
    l = glob(''.join([target_path, '*', sufx[0]]))
    l_subs = [i[len(target_path):len(i)-len(sufx[0])] for i in l]
    d = {}
    for i in l_subs:
        d[i] = {'bed_norm': -1, 'bed_tumor': -1, 'tsv_norm': -1, 'tsv_tumor': -1}
        for s in sufx:
            target = ''.join([target_path, i, s])
            if s == '_cnts2R.tsv.gz':
                header = True
                l_norm = []
                l_tumor = []
            else:
                header = False
                infields = ['chr', 'p0', 'p1', 'an_type', 'score', 'strd', 'info']
                l = []
            with gzip.open(target, 'rt') as infile:
                for inline in infile:
                    if header:
                        infields = inline.strip('\n').split('\t')
                        header = False
                        continue
                    else:
                        data = inline.strip('\n').split('\t')
                        data = dict(list(zip(infields, data)))
                    if 't_type' in data:
                        if data['t_type'] == 'norm':
                            l_norm.append(data['chrom'])
                        else:
                            l_tumor.append(data['chrom'])
                    else:
                        l.append(data['chr'])
            if s == '_cnts2R.tsv.gz':
                d[i]['tsv_norm'] = extract_value(l_norm, count_type)
                d[i]['tsv_tumor'] = extract_value(l_tumor, count_type)
            elif s == '_normal_hetcnts.bed.gz':
                d[i]['bed_norm'] = extract_value(l, count_type)
            else:
                d[i]['bed_tumor'] = extract_value(l, count_type)
    return d


def write_summary(target, d):
    with open(target, 'w') as outfile:
        s1 = '\t'.join([k for k in d[list(d.keys())[0]]])
        s2 = ''.join(['subject_id', '\t', s1, '\n'])
        outfile.write(s2)
        for k1 in d:
            s1 = '\t'.join([str(d[k1][k2]) for k2 in d[k1]])
            s2 = ''.join([k1, '\t', s1, '\n'])
            outfile.write(s2)
    print('file written to', target)


def update_ann(target_path):
    '''
    takes previously compiled phenotype annotations and adds site and chrom
    summaries created by this script
    '''
    targets = ['pheno', 'summary_sites', 'summary_chroms']
    d = {}
    for i in targets:
        target = ''.join([target_path, i, '.tsv'])
        d[i] = read_csv(target, sep='\t')
    d['file_info'] = read_csv(
        ''.join([target_path, 'coad-read.file_info']), sep='\t', header=None,
        names=['drop1', 't_type', 'drop2', 'file_id']
    )
    df1 = merge(d['file_info'], d['pheno'], how='left', on='file_id')
    df2 = merge(
        d['summary_sites'], d['summary_chroms'], how='inner',
        on='subject_id', suffixes=('_sites', '_chroms')
    )
    df = merge(df1, df2, how='left', on='subject_id')
    df = df.drop(labels=['drop1', 'drop2'], axis=1)
    target = ''.join([target_path, 'ann.tsv'])
    df.to_csv(target, sep='\t', index=False)
    print('file written to', target)
    return df


def main():
    target_path = '/scratch/chd5n/aneuploidy/'
    for count_type in ['sites', 'chroms']:
        d = count_data(count_type)
        target = ''.join([target_path, 'summary_', count_type, '.tsv'])
        write_summary(target, d)
    df = update_ann(target_path)


if __name__ == '__main__':
    main()
else:
    print('functions loaded for', modname)
