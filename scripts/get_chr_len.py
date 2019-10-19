##
## this script extracts the length in bases of canonical chromosomes
## from a reference fasta file if such lengths are included
##

## usage: python ~/projects/aneuploidy/scripts/get_chr_len.py >
##    /scratch/chd5n/Reference_genome/GRCh38.d1.vd1.chr1-XY.size.tsv


in_file_path = '/scratch/chd5n/Reference_genome/GRCh38.d1.vd1.fa'
in_fields = ('chrom','AC','gi','LN','rl','M5','AS')
data_list = []


def make_chrom_canon(start,stop,sex):
    chr_list=[]
    for i in range(start,stop):
        chr_list.append('chr'+str(i))
    if sex:
        chr_list.append('chrX')
        chr_list.append('chrY')
    return(chr_list)


with open(in_file_path) as in_f:
    for in_line in in_f:
        if 'chr' in in_line.strip('\n').split('\t')[0]:
            data_list.append(dict(list(zip(in_fields, in_line.strip('\n').split()))))


chr_canon = make_chrom_canon(1,23,True)
cum = 0
for each in data_list:
    if each['chrom'][1:] in chr_canon:
        print('\t'.join([each['chrom'][1:],each['LN'][3:],str(cum)]))
        cum = cum + int(each['LN'][3:])
