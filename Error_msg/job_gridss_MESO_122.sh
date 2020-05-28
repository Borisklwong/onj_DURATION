#!/bin/bash
#PBS -l nodes=1:ppn=8,mem=32gb,walltime=168:00:00
module remove R
module add R/3.6.2 bwa samtools java
/usr/bin/time -o /home/wong.b/onj_DURATION/processing/MESO/gridss/MESO_122/time.log /home/wong.b/onj_DURATION/gridss_files/gridss-2.9.2//gridss.sh 	--jar /home/wong.b/onj_DURATION/gridss_files/gridss-2.9.2//gridss-2.9.2-gridss-jar-with-dependencies.jar 	--blacklist /home/wong.b/onj_DURATION/hg19_BW/wgEncodeDacMapabilityConsensusExcludable.bed 	--reference /home/wong.b/onj_DURATION/hg19_BW/hg19.fa 	--output /home/wong.b/onj_DURATION/processing/MESO/gridss/MESO_122/MESO_122.gridss.vcf 	--workingdir /home/wong.b/onj_DURATION/processing/MESO/gridss/MESO_122 	--assembly /home/wong.b/onj_DURATION/processing/MESO/gridss/MESO_122/MESO_122.assembly.bam 	--threads 8 	/home/wong.b/onj_DURATION/processing/BGI_public/CL1000074821_merged.bam 	/home/wong.b/onj_DURATION/processing/MESO/MESO_122_process/MESO_122_sort_mark.bam 	# --sanityCheck


Rscript /home/wong.b/onj_DURATION/gridss_files/gridss-2.9.2//gridss_somatic_filter.R 	--pondir /home/wong.b/onj_DURATION/gridss_files/pon3792v1/ 	--input /home/wong.b/onj_DURATION/processing/MESO/gridss/MESO_122/MESO_122.gridss.vcf 	--output /home/wong.b/onj_DURATION/processing/MESO/gridss/MESO_122/MESO_122.gridss.somatic.vcf 	--fulloutput /home/wong.b/onj_DURATION/processing/MESO/gridss/MESO_122/MESO_122.gridss.somaticfull.vcf 	--pairoutput  /home/wong.b/onj_DURATION/processing/MESO/gridss/MESO_122/MESO_122.gridss.paired.sv.bedpe 	--pairoutput_details /home/wong.b/onj_DURATION/processing/MESO/gridss/MESO_122/MESO_122.gridss.paired.csv 	--ref BSgenome.Hsapiens.UCSC.hg19 	--scriptdir /home/wong.b/onj_DURATION/gridss_files/gridss-2.9.2/


exit 0
	
