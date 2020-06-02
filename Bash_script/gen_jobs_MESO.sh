#!/bin/bash

walltime="168:00:00"
gridss_pondir=/home/wong.b/onj_DURATION/gridss_files/pon3792v1/
#1 ref
#2 normal
#3 tumour
#4 job
create_jobs () {
	ref=$1
	normal=$2
	tumour=$3
	sample=$4
	echo "Processing $4 $1 $2 $3"
	wd=$PWD/gridss/$sample
	mkdir -p $wd
	gridss_dir=/home/wong.b/onj_DURATION/gridss_files/gridss-2.9.2/
	cat > $wd/job_gridss_$sample.sh << EOF
#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=32gb,walltime=$walltime
module remove R
module add R/3.6.2 bwa samtools java
/usr/bin/time -o $wd/time.log $gridss_dir/gridss.sh \
	--jar $gridss_dir/gridss-2.9.2-gridss-jar-with-dependencies.jar \
	--blacklist /home/wong.b/onj_DURATION/hg19_BW/wgEncodeDacMapabilityConsensusExcludable.bed \
	--reference $ref \
	--output $wd/$sample.gridss.vcf \
	--workingdir $wd \
	--assembly $wd/$4.assembly.bam \
	--threads 8 \
	$normal \
	$tumour \
	# --sanityCheck


Rscript $gridss_dir/gridss_somatic_filter.R \
	--pondir /home/wong.b/onj_DURATION/gridss_files/pon3792v1/ \
	--input $wd/$sample.gridss.vcf \
	--output $wd/$sample.gridss.somatic.vcf \
	--fulloutput $wd/$sample.gridss.somaticfull.vcf \
	--normalordinal 1 \
	--tumourordinal 2 \
	--ref BSgenome.Hsapiens.UCSC.hg19 \
	--scriptdir $gridss_dir


exit 0
	
EOF

}

# Truth set: BGI MESO WGS DATA
SAMPLE1=15B26293_1B
SAMPLE2=15T92457_B
SAMPLE3=17B_17424_1_4
SAMPLE4=MESO_122
SAMPLE5=MESO_455
SAMPLE6=MESO_474
SAMPLE7=MESO_514
SAMPLE8=MESO_1174
SAMPLE9=MESO_1451
SAMPLE10=MESO_1935
SAMPLE11=MESO_2493
SAMPLE12=MESO_2542
SAMPLE13=MESO_5242
SAMPLE14=MESO_5365
SAMPLE15=MESO_11815

create_jobs /home/wong.b/onj_DURATION/hg19_BW/hg19.fa /home/wong.b/onj_DURATION/processing/BGI_public/CL1000074821_merged.bam /home/wong.b/onj_DURATION/processing/MESO/${SAMPLE1}_process/${SAMPLE1}_sort_mark.bam $SAMPLE1
create_jobs /home/wong.b/onj_DURATION/hg19_BW/hg19.fa /home/wong.b/onj_DURATION/processing/BGI_public/CL1000074821_merged.bam /home/wong.b/onj_DURATION/processing/MESO/${SAMPLE2}_process/${SAMPLE2}_sort_mark.bam $SAMPLE2
create_jobs /home/wong.b/onj_DURATION/hg19_BW/hg19.fa /home/wong.b/onj_DURATION/processing/BGI_public/CL1000074821_merged.bam /home/wong.b/onj_DURATION/processing/MESO/${SAMPLE3}_process/${SAMPLE3}_sort_mark.bam $SAMPLE3
create_jobs /home/wong.b/onj_DURATION/hg19_BW/hg19.fa /home/wong.b/onj_DURATION/processing/BGI_public/CL1000074821_merged.bam /home/wong.b/onj_DURATION/processing/MESO/${SAMPLE4}_process/${SAMPLE4}_sort_mark.bam $SAMPLE4
create_jobs /home/wong.b/onj_DURATION/hg19_BW/hg19.fa /home/wong.b/onj_DURATION/processing/BGI_public/CL1000074821_merged.bam /home/wong.b/onj_DURATION/processing/MESO/${SAMPLE5}_process/${SAMPLE5}_sort_mark.bam $SAMPLE5
create_jobs /home/wong.b/onj_DURATION/hg19_BW/hg19.fa /home/wong.b/onj_DURATION/processing/BGI_public/CL1000074821_merged.bam /home/wong.b/onj_DURATION/processing/MESO/${SAMPLE6}_process/${SAMPLE6}_sort_mark.bam $SAMPLE6
create_jobs /home/wong.b/onj_DURATION/hg19_BW/hg19.fa /home/wong.b/onj_DURATION/processing/BGI_public/CL1000074821_merged.bam /home/wong.b/onj_DURATION/processing/MESO/${SAMPLE7}_process/${SAMPLE7}_sort_mark.bam $SAMPLE7
create_jobs /home/wong.b/onj_DURATION/hg19_BW/hg19.fa /home/wong.b/onj_DURATION/processing/BGI_public/CL1000074821_merged.bam /home/wong.b/onj_DURATION/processing/MESO/${SAMPLE8}_process/${SAMPLE8}_sort_mark.bam $SAMPLE8
create_jobs /home/wong.b/onj_DURATION/hg19_BW/hg19.fa /home/wong.b/onj_DURATION/processing/BGI_public/CL1000074821_merged.bam /home/wong.b/onj_DURATION/processing/MESO/${SAMPLE9}_process/${SAMPLE9}_sort_mark.bam $SAMPLE9
create_jobs /home/wong.b/onj_DURATION/hg19_BW/hg19.fa /home/wong.b/onj_DURATION/processing/BGI_public/CL1000074821_merged.bam /home/wong.b/onj_DURATION/processing/MESO/${SAMPLE10}_process/${SAMPLE10}_sort_mark.bam $SAMPLE10
create_jobs /home/wong.b/onj_DURATION/hg19_BW/hg19.fa /home/wong.b/onj_DURATION/processing/BGI_public/CL1000074821_merged.bam /home/wong.b/onj_DURATION/processing/MESO/${SAMPLE11}_process/${SAMPLE11}_sort_mark.bam $SAMPLE11
create_jobs /home/wong.b/onj_DURATION/hg19_BW/hg19.fa /home/wong.b/onj_DURATION/processing/BGI_public/CL1000074821_merged.bam /home/wong.b/onj_DURATION/processing/MESO/${SAMPLE12}_process/${SAMPLE12}_sort_mark.bam $SAMPLE12
create_jobs /home/wong.b/onj_DURATION/hg19_BW/hg19.fa /home/wong.b/onj_DURATION/processing/BGI_public/CL1000074821_merged.bam /home/wong.b/onj_DURATION/processing/MESO/${SAMPLE13}_process/${SAMPLE13}_sort_mark.bam $SAMPLE13
create_jobs /home/wong.b/onj_DURATION/hg19_BW/hg19.fa /home/wong.b/onj_DURATION/processing/BGI_public/CL1000074821_merged.bam /home/wong.b/onj_DURATION/processing/MESO/${SAMPLE14}_process/${SAMPLE14}_sort_mark.bam $SAMPLE14
create_jobs /home/wong.b/onj_DURATION/hg19_BW/hg19.fa /home/wong.b/onj_DURATION/processing/BGI_public/CL1000074821_merged.bam /home/wong.b/onj_DURATION/processing/MESO/${SAMPLE15}_process/${SAMPLE15}_sort_mark.bam $SAMPLE15
























