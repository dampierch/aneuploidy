# Aneuploidy in CRC

An investigation into the drivers of CIN in CRC using TCGA WXS.

## Sample selection
1. From the [GDC data portal](https://portal.gdc.cancer.gov/), we select all colon, rectosigmoid, and rectum samples from TCGA-COAD and TCGA-READ without a filter for tumor type (i.e. we do not care whether the diagnosis was adenocarcinoma or cystic neoplasm).
2. We add all WXS files to our cart, and we find there are a total of 1,303 files from 608 unique cases (i.e. subjects); this data takes up 31.31 TB of disk space.
3. From our cart, we download the following annotation files:
  * biospecimen.cart.2019-07-13.json
  * biospecimen.cart.2019-07-13.tar.gz
  * clinical.cart.2019-07-13.json
  * clinical.cart.2019-07-13.tar.gz
  * gdc_sample_sheet.2019-07-13.tsv
  * metadata.cart.2019-07-13.json
4. Using gdc_sample_sheet.2019-07-13.tsv along with [samples.py](scripts/samples.py), we explore the basic characteristics of our cohort.

| Sample type | *n* |
| :--: | :--: |
| Primary Tumor | 628 |
| Blood Derived Normal | 560 |
| Solid Tissue Normal | 112 |
| Recurrent Tumor | 2 |
| Metastatic | 1 |

5. Among identifiers in the sample sheet, there are 1,259 unique 'Sample ID's, 608 unique 'Case ID's, 1,303 unique 'File Name's and 'File ID's.

### Normal tissue
1. We start by examining the normal tissue from which we will call SNVs.
2. There are 672 samples of normal tissue WXS.
3. There are 601 unique cases among those samples; 71 samples are non-unique.
4. There are 68 cases with multiple tissue samples, 65 with just two and three with a total of three samples each, accounting for all 71 non-unique samples (139 files).
5. We may have a soft preference for blood-derived normal WXS over adjacent normal tissue, but the choice at this point is arbitrary between duplicate blood-derived normals due to lack of quality parameters.
6. There are 533 cases with unique tissue samples, 477 blood-derived and 56 solid tissue (i.e. adjacent normal).

### Tumor tissue
1. We next examine the tumor tissue we will test for functional aneuploidy.
2. There are 628 samples of tumor tissue WXS.
3. There are 594 unique cases among those samples; 34 samples are non-unique.
4. There are 28 cases with multiple tissue samples, 22 with just two and six with a total of three samples each, accounting for all 34 non-unique samples (62 files).
5. At this point, the choice is arbitrary among the duplicates, as all are Primary Tumor, and we do not have quality parameters.
6. There are 566 cases with unique tissue samples, all primary tumor.

### Pilot
1. We create a pilot set of five cases from the set of unique normal samples and match them with their five tumor counterparts. We use [samples.py](scripts/samples.py) to generate a pilot manifest.
2. We use the [gdc-client](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) to download BAM files to a `/scratch` directory.
3. After checking library format with [libFormat.slurm](scripts/libFormat.slurm), we use [Biobambam2](https://www.sanger.ac.uk/science/tools/biobambam) to revert the BAM files to FASTQ and [pigz](https://zlib.net/pigz/) to compress them.

## Variants

1. From our normal WXS data, we want to call germline SNVs. We could also call somatic mutations, but the TCGA VCF files should suffice for that.
2. We start with our subset of five normal tissue samples that are unique and have a unique tumor sample partner. There are 510 of these unique pairs in total.
3. Our [download script](download.slurm) does not seem to work on the compute nodes, so we run it on a login node. We expect the download to be about 240 GB (average file size of [31 210 000 000 000 bytes / 1303 files] = 24 GB * 10 files). The actual file size is 218 GB for 10 BAM files.
4. We call variants with...
