#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=32gb,walltime=168:00:00
# Local system locations
module remove R
module add R/3.6.2 bwa samtools java sambamba circos

install_dir=/home/wong.b/projects/onj_DURATION/gridss-purple-linx/gridss-purple-linx-master/gridss-purple-linx-v1.2.0
ref_data=/home/wong.b/projects/onj_DURATION/gridss-purple-linx/gridss-purple-linx-master/hg19
data_dir=/home/wong.b/onj_DURATION/processing/MESO/g-p-l/four_normal_gridss/15B26293_1B


GRIDSS_VERSION=2.9.3
AMBER_VERSION=3.3
COBALT_VERSION=1.8
PURPLE_VERSION=2.43
LINX_VERSION=1.9

export GRIDSS_JAR=/home/wong.b/projects/onj_DURATION/gridss-purple-linx/gridss-purple-linx-master/gridss-purple-linx-v1.2.0/gridss/gridss-2.9.3-gridss-jar-with-dependencies.jar
export AMBER_JAR=/home/wong.b/projects/onj_DURATION/gridss-purple-linx/gridss-purple-linx-master/gridss-purple-linx-v1.2.0/hmftools/amber-3.3.jar
export COBALT_JAR=/home/wong.b/projects/onj_DURATION/gridss-purple-linx/gridss-purple-linx-master/gridss-purple-linx-v1.2.0/hmftools/cobalt-1.8.jar 
export PURPLE_JAR=/home/wong.b/projects/onj_DURATION/gridss-purple-linx/gridss-purple-linx-master/gridss-purple-linx-v1.2.0/hmftools/purple-2.43.jar
export LINX_JAR=/home/wong.b/projects/onj_DURATION/gridss-purple-linx/gridss-purple-linx-master/gridss-purple-linx-v1.2.0/hmftools/sv-linx_1.9.jar

/home/wong.b/projects/onj_DURATION/gridss-purple-linx/gridss-purple-linx-master/gridss-purple-linx-v1.2.0/gridss-purple-linx/gridss-purple-linx.sh 	-o /home/wong.b/onj_DURATION/processing/MESO/g-p-l/four_normal_gridss/15B26293_1B 	-n /home/wong.b/onj_DURATION/processing/BGI_public/g-p-l/CNX0043488_process/CL1000074821_merged_mark.bam 	-t /home/wong.b/onj_DURATION/processing/MESO/g-p-l/15B26293_1B_process/15B26293_1B_sort_mark.bam  	--nosnvvcf 	-s 15B26293_1B 	--normal_sample CL1000074821 	--tumour_sample 15B26293_1B 	--ref_dir /home/wong.b/projects/onj_DURATION/gridss-purple-linx/gridss-purple-linx-master/hg19 	--install_dir /home/wong.b/projects/onj_DURATION/gridss-purple-linx/gridss-purple-linx-master/gridss-purple-linx-v1.2.0 

