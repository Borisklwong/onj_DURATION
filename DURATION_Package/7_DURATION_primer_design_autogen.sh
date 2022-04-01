#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=DURATION_gripss_to_primer
#SBATCH --time=48:00:00
#SBATCH --mem=16G
#Boris Wong
#
# Local system locations
date=23MAR2022_MESO
	file_dir=/home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing/${date}
	suffix=T.gripss.somatic.filtered.vcf.gz
	job_dir=${file_dir}/job
	primer3_dir=/home/users/allstaff/wong.b/onj_DURATION/tools/primer3-2.4.0
	R_script_dir=/home/users/allstaff/wong.b/onj_DURATION/R_scpript/DURATION_package_R
	blat_dir=/home/users/allstaff/wong.b/onj_DURATION/tools/blat
	mkdir -p ${job_dir}
for vcf in $(find $file_dir -name *$suffix); do
        sample_name=$(basename $vcf $suffix)
	vcf_dir=$(dirname $vcf)
	assembly_bam=$(find $(dirname $vcf_dir) -name *assembly.bam.sv.bam)
	job_file="${job_dir}/${sample_name}_$(basename $(dirname $vcf))_post_gridss.job"
	mkdir -p $vcf_dir
	echo  "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=${sample_name}_post_gripss
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --error=${sample_name}_slurm-%j.out
# THIS IS A GENERATED SCRIPT
# post_gripss

cd ${vcf_dir}
module remove R
module add R/4.0.0

Rscript ${R_script_dir}/1_DURATION.circos.R $sample_name $vcf_dir $vcf
echo "${sample_name} circos plot ready"
Rscript ${R_script_dir}/2_extract_seq_from_gripss_FFPE_only.R $sample_name $vcf_dir $vcf $assembly_bam
echo "${sample_name} GRIDSS assembly sequence ready"
# design primers for top 30 GRIDSS call
cd ${primer3_dir}/src
head -n 120 ${vcf_dir}/${sample_name}_vcf_to_p3.txt | ./primer3_core --strict_tags --p3_settings_file=${primer3_dir}/settings_files/biorad_settings_with_mispriming_library.txt --output=${vcf_dir}/${sample_name}_p3_out 
echo "${sample_name} primer designed"
head -n 120 ${vcf_dir}/${sample_name}_vcf_to_p3_overlap.txt | ./primer3_core --strict_tags --p3_settings_file=${primer3_dir}/settings_files/biorad_settings_with_mispriming_library.txt --output=${vcf_dir}/${sample_name}_p3_out_overlap 
echo "${sample_name} overlap primer designed"
Rscript ${R_script_dir}/3_p3_out_to_fasta.R ${sample_name} ${vcf_dir}
echo "${sample_name} fasta for top 10  SVs ready"
cd ${blat_dir}
./blat -fastMap -tileSize=8 -oneOff=1 -minMatch=1 -minScore=10 -repMatch=16000 -t=DNA \
-q=DNA ./hg19_blat/hg19.fa \
${vcf_dir}/${sample_name}_p3_out_full.fasta \
${vcf_dir}/${sample_name}_p3_out_full.psl
echo "${sample_name} primer blat result ready"
Rscript ${R_script_dir}/4_blat_filter.R ${sample_name} ${vcf_dir}
echo "${sample_name} blat result filtered"" > $job_file
	echo sbatch $job_file
	
done