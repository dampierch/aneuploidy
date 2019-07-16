#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6000
#SBATCH --time=2:00:00
#SBATCH --partition=standard
#SBATCH --account=cphg_caseylab
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chd5n@virginia.edu

# To get the GDC-Client:
# cd ~/Apps
# wget https://gdc.cancer.gov/system/files/authenticated%20user/0/gdc-client_v1.4.0_Ubuntu_x64.zip
# unzip gdc-client_v1.4.0_Ubuntu_x64.zip

home=/home/chd5n/
aneuploidy_home=/scratch/chd5n/aneuploidy/
aneuploidy_raw=${aneuploidy_home}raw-data/
d=${aneuploidy_raw}sequencing/
m=${aneuploidy_raw}annotations/dn-pilot-manifest_2019-07-13.tsv
t=${aneuploidy_raw}gdc-user-token.2019-07-14T23_26_20.271Z.txt
logfile=${aneuploidy_raw}sequencing/gdc-download_2_2019-07-14.out

${home}/Apps/gdc-client/gdc-client download -d ${d} --log-file ${logfile} -m ${m} -t ${t} --debug
