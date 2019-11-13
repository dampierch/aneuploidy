## chd coverage

import glob
import numpy as np

data_dir = '/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/'
suffix = '_coverage.bed'

file_set = glob.glob(data_dir+'*'+suffix)

set_cov = []
file_cov = {}
for file in file_set:
    sample = file.split('/')[-1].split('_coverage')[0]
    file_cov['sample'] = sample
    line_num = 0
    int_cov = []
    with open(file,'r') as in_f:
        for in_line in in_f:
            if line_num > 0 and in_line[0:5] != 'whole':
                int_cov.append(int(in_line.strip('\n').split('\t')[5]))
            line_num = line_num + 1
        file_cov['int_cov'] = int_cov
    set_cov.append(file_cov)

all_cov = []
for each in set_cov:
    all_cov = all_cov + each['int_cov']
data = np.array(all_cov)
data = data[data > 2]

import matplotlib.pyplot as plt
import pandas as pd
data = pd.Series(set_cov[0]['int_cov'])
data_kde = data.plot.kde()
plt.show()


data = np.array(set_cov[3]['int_cov'])
data = data[data > 2]

n, bins, patches = plt.hist(x=data, bins='auto', color='#0504aa',
                            alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel('Coverage (x)')
plt.ylabel('# targets with given coverage')
plt.title('Histogram of mean coverage per target')
maxfreq = n.max()
# Set a clean upper y-axis limit.
plt.ylim(ymax=maxfreq+100)
plt.xlim(xmax=300, xmin=0)
plt.show()


## notes
## broad on filtering
## https://software.broadinstitute.org/gatk/documentation/article.php?id=4721
## https://software.broadinstitute.org/gatk/documentation/article.php?id=6925
## broad recommends starting with QD<2 if using hard filter
##
## broad on coverage thresholds
## https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_diagnostics_diagnosetargets_DiagnoseTargets.php
## DiagnoseTargets calls low_coverage at 5x
## https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
## HaplotypeCaller --minReadsPerAlignmentStart default is 10 but this is only used in downsampling
##
## WUSTL VarScan for reference
## http://varscan.sourceforge.net/support-faq.html#usage-wgs-exome-rna
## defaults optimized for WXS expecting 10-20x coverage
## default --min-coverage set to 8
## http://varscan.sourceforge.net/somatic-calling.html
## TCGA use of VarScan
## https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#variant-call-command-line-parameters
## --min-coverage set to 8 for normal, 6 for tumor


## QD distrubution in my data
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT TCGA-T9-A92H-10A-01D-A370-10
# chr1 17385 . G A 176.77 . AC=1;AF=0.500;AN=2;BaseQRankSum=-2.255;ClippingRankSum=0.000;DP=39;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=39.04;MQRankSum=-1.732;QD=4.53;ReadPosRankSum=0.201;SOR=0.099 GT:AD:DP:GQ:PL 0/1:30,9:39:99:205,0,1033
#
# QUAL is field 5 ('QUAL')
# QD is field 7 ('INFO') subfield 13 ('QD')

from itertools import compress
import pandas as pd
import matplotlib.pyplot as plt

feat_0 = []
feat_1 = []
with open('/scratch/chd5n/aneuploidy/raw-data/sequencing/crunch/TCGA-AF-3400-11A-01D-1554-10_Illumina_gdc_realn.snp.indel.vcf_L','r') as in_f:
    for in_line in in_f:
        if in_line[0]=='#':  ## use this test with continue to skip processing leading comments
            if in_line[0:6]=='#CHROM':
                field_names = in_line.strip('\n').split('\t')
                field_names[0] = 'CHROM'  ## this gets rid of '#'
                last_field = field_names[-1]  ## variable name with subject_id in first three '-' delim fields
            continue  ## pops out to next iteration of for loop to skip processing leading comments
        data = dict(list(zip(field_names,in_line.strip('\n').split('\t'))))
        feat_0.append(float(data['QUAL']))
        idx_bool = ['QD=' in x for x in data['INFO'].split(';')]
        if sum(idx_bool) > 0:
            idx_int = list(compress(range(len(idx_bool)), idx_bool))[0]
            feat_1.append(float(data['INFO'].split(';')[idx_int].split('=')[1]))

dat_0 = pd.Series(feat_0)
dat_1 = pd.Series(feat_1)
dat_0.describe()
dat_1.describe()

dat_0_low = dat_0[dat_0 < 100]
len(dat_0_low) ## 3424
dat_1_low = dat_1[dat_1 < 2]
len(dat_1_low) ## 691 at 2, 2338 at 5

dat = dat_0
dat_kde = dat.plot.kde()
plt.show()

dat = dat_1
dat_kde = dat.plot.kde()
plt.show()
