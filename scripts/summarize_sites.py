import gzip
from glob import glob


modname = 'summarize_sites'


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


def main():
    target_path = '/scratch/chd5n/aneuploidy/'
    for count_type in ['sites', 'chroms']:
        d = count_data(count_type)
        target = ''.join([target_path, 'summary_', count_type, '.tsv'])
        write_summary(target, d)


if __name__ == '__main__':
    main()
else:
    print('functions loaded for', modname)
