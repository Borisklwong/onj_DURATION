#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=32gb,walltime=168:00:00
# Local system locations
module remove R
module add R/3.6.2 bwa samtools java sambamba circos

install_dir=~/projects/onj_DURATION/gridss-purple-linx/gridss-purple-linx-master/gridss-purple-linx-v1.2.0
ref_data=~/projects/onj_DURATION/gridss-purple-linx/gridss-purple-linx-master/hg19
data_dir=~/projects/onj_DURATION/gridss-purple-linx/gridss-purple-linx-master/smoke_test

echo rm -r $data_dir/amber $data_dir/cobalt $data_dir/gridss $data_dir/logs $data_dir/purple

GRIDSS_VERSION=2.9.2
AMBER_VERSION=3.3
COBALT_VERSION=1.8
PURPLE_VERSION=2.43
LINX_VERSION=1.9

export GRIDSS_JAR=$install_dir/gridss/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar
export AMBER_JAR=$install_dir/hmftools/amber-${AMBER_VERSION}.jar
export COBALT_JAR=$install_dir/hmftools/cobalt-${COBALT_VERSION}.jar 
export PURPLE_JAR=$install_dir/hmftools/purple-${PURPLE_VERSION}.jar
export LINX_JAR=$install_dir/hmftools/sv-linx_${LINX_VERSION}.jar

$install_dir/gridss-purple-linx/gridss-purple-linx.sh \
	-o $data_dir \
	-n $data_dir/CPCT12345678R.bam \
	-t $data_dir/CPCT12345678T.bam  \
	-v $data_dir/CPCT12345678T.somatic_caller_post_processed.vcf.gz \
	--snvvcf $data_dir/CPCT12345678T.somatic_caller_post_processed.vcf.gz \
	-s CPCT12345678 \
	--normal_sample CPCT12345678R \
	--tumour_sample CPCT12345678T \
	--ref_dir $ref_data \
	--install_dir $install_dir \
	
