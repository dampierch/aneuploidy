# TCGA COAD/READ exome capture kits

## exome capture kits used

| Kit | *n* |
| :--: | :--: |
| SeqCap EZ HGSC VCRome | 878 |
| SureSelect Human All Exon 38 Mb v2 | 186 |
| Gapfiller_7m | 94 |
| SeqCap EZ Human Exome Library v2.0 | 49 |
| VCRome V2.1 | 47 |
| VCRomeV2.1-PKv1 | 31 |
| Custom V2 Exome Bait, 48 RXN X 16 tubes | 2 |
| SeqCap EZ Exome V2.0 | 1 |
| NaN | 15 |

## dict for python
```
{
'SeqCap EZ HGSC VCRome': 'seqcap_vcrome',
'SureSelect Human All Exon 38 Mb v2': 'sureselect_v2',
'Gapfiller_7m': 'seqcap_vcrome',
'SeqCap EZ Human Exome Library v2.0': 'seqcap_v3',
'VCRome V2.1': 'seqcap_vcrome',
'VCRomeV2.1-PKv1': 'seqcap_vcrome',
'Custom V2 Exome Bait, 48 RXN X 16 tubes': 'agilent_bi',
'SeqCap EZ Exome V2.0': 'seqcap_v3',
'NaN': 'seqcap_vcrome'
}
```

## SeqCap EZ HGSC VCRome
* 45.1 Mb coverage
* [roche website](https://sequencing.roche.com/en/products-solutions/by-category/target-enrichment/hybridization/seqcap-ez-hgsc-vcrome.html)

## SureSelect Human All Exon 38 Mb v2
* 38 Mb coverage
* [biostars guidance](https://www.biostars.org/p/57675/)
* [agilent suredesign website](https://earray.chem.agilent.com/suredesign/index.htm?sessiontimeout=true)
  1. create login and workgroup
  2. find designs > SureSelect > agilent catalog
* earliest in online catalogue is SureSelect Human All Exon V4, which has 51 Mb coverage (hg19)
* emailed informatics_support@agilent.com for V2 on 9/24/2019
* Mandar Bedse responded with S0293689_Covered.bed for V2
* we want the "regions" bed file as the primary target file; the covered file is the bait
* for V4, regions and covered are the same
* Mandar Bedse confirmed that regions and covered are same for V2 as well 9/25/2019

### login info
* login: chd5n@virginia.edu
* password: 1surepass
* workgroup: tcga
* tech support password: 1Surepass!

### guidance for interpreting Agilent:
1. [biostars 1](https://www.biostars.org/p/5187/)
2. [biostars 2](https://www.biostars.org/p/220939/)

## Gapfiller_7m
* 7 Mb coverage
* Nimblegen, Roche technical support inquiry placed 9/24/2019
* [obscure lecture website](https://www.lanl.gov/conferences/finishfuture/pdfs/2012%20talks/sfaf12-muzny.pdf)
* Muzny talk from 2012 acknowledges TCGA exome gapfiller design with 7Mb genomic region targeted (10,200 genes)
* email sent to questions@hgsc.bcm.tmc.edu on 9/24/2019

## SeqCap EZ Human Exome Library v2.0
* can only find V3, 64 Mb coverage, hg19
* [website for v3.0](https://sequencing.roche.com/en/products-solutions/by-category/target-enrichment/hybridization/seqcap-ez-exome-v3-kit.html)

## Custom V2 Exome Bait, 48 RXN X 16 tubes
* email sent to genomics@broadinstitute.org
* Agilent design for Broad Institute
* Carrie Cibulskis <carrie@broadinstitute.org> responded with name of file: whole_exome_agilent_1.1_refseq_plus_3_boosters aka whole_exome_agilent_2.1
* file found on: https://bitbucket.org/cghub/cghub-capture-kit-info/src/master/BI/vendor/Agilent/
* the [CGHub bitbucket site](https://bitbucket.org/cghub/cghub-capture-kit-info/src/master/) seems like a great resource
* the [CGHub bitbucket wiki](https://bitbucket.org/cghub/cghub-capture-kit-info/wiki/Home) also seems useful
* targetIntervals is the file to use (as opposed to baitIntervals)

## Plan
* create hg38 exome references with liftOver
* assume all baseline references are hg19
* alternative: download current exome capture reference file with hg38
* [agilent suredesign website](https://earray.chem.agilent.com/suredesign/index.htm?sessiontimeout=true)
* login: chd5n@virginia.edu
* password: 1surepass
* workgroup: tcga
* tech support password: 1Surepass!
* SureSelect Human All Exon V7: 48.2 Mb
* SureSelect Clinical Research Exome V2: 67.3 Mb *i like this one*
* SureSelect Human All Exon V6+UTR r2: 91 Mb
* SureSelect Human All Exon V6 r2: 60 Mb
* SureSelect Human All Exon V6+COSMIC r2: 66 Mb


```
kit_path = '/scratch/chd5n/aneuploidy/exome-kits/'
kit_dirs = ['Agilent_BI/','SeqCapEZ_Exome_v3.0/','SeqCapEZ_VCRome_2.1/','SureSelect_V2/']
hg19_files = [
    'whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed',
    'SeqCap_EZ_Exome_v3_hg19_primary_targets.bed',
    'VCRome_2_1_hg19_primary_targets.bed',
    'SureSelect_V2_Covered.bed'
]



```

## code
kit_path = '/scratch/chd5n/aneuploidy/exome-kits/SeqCapEZ_Exome_v3.0/'

file = kit_path + 'SeqCap_EZ_Exome_v3_hg19_capture_targets.bed'
x = pd.read_table(file, skiprows=2, usecols=[0,1,2], header=None,
    names=['chr','start','stop'])
x['int'] = x.stop - x.start
x.int.describe()

file = kit_path + 'SeqCap_EZ_Exome_v3_hg19_primary_targets.bed'
y = pd.read_table(file, skiprows=2, usecols=[0,1,2], header=None,
    names=['chr','start','stop'])
y['int'] = y.stop - y.start
y.int.describe()
