# Aneuploidy Working Group Meetings

## 2021-01-15
* [ ] get cms classifications for tcga samples
* at least one group has generated them for this [manuscript](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6221003/)
* [ ] show rank by aneuploidy and frequency of gene mutations
* [ ] show kaplan meier curves for most vs least aneuploidy tumors
* exclude msi-h from analysis
* focus on chrom 18
* look at SMAD4 mutations and other mutations enriched on chrom 18

## 2020-12-17
* done with CRC classification wherein tumors scored by aneuploidy
* goal: find pattern among chromosome specific segments (chrom arms p and q)
  - look for recurrent chromosomal events in CRC
  - check for correlation with contamination
  - check for correlation with MSI
* expect chr5 (APC), chr17 (TP53) to show aneuploidy
* want to know what percentage of tumors have given pattern of aneuploidy
* want to know what percentage of CMS subtypes have given pattern of aneuploidy
* Her2 FISH could be used for something
* next meeting Jan 15 3pm

## 2020-03-18
* Goal: classify heterozygous (i.e. diploid), LOH, and aneuploid states correctly
* Resolution: chromosome arm
* Pankaj joins team
* Pearson update:
  1. done fitting normal distribution to central diploid peak and peripheral LOH peaks
  2. unsure what to do with aneuploid signal between extreme LOH peaks and diploid peak
  3. working on fitting beta distribution for aneuploid signal
    - beta useful because can change shape to be flatter than normal distribution
  4. goal is 3 state HMM in Python
* Future: Current Protocols in Bioinformatics
* TODO:
  1. add normal SD curve to aneuploid rank plot to help with interpretation
  2. use VAAST to identify genes of interest
  3. for follow-up analyses, include MSS + MSI-L together
  4. share AF density plots with Stukenberg, Pearson
  5. share AF counts at hetsites with Pankaj

## 2019-07-10
* introduction to method: identification of sites with heterozygous SNPs in WXS data
* review of TCGA breast cancer findings
* plan for simple internal cancer center grant submission
