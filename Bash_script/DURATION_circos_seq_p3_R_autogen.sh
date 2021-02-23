#!/bin/bash

#Gen job from post g-p-l to primer fasta

# Local system locations
	vcf_dir=/home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing
	primer3_dir=/home/users/allstaff/wong.b/onj_DURATION/tools/primer3-2.4.0
	job_dir=${vcf_dir}/.job
	R_script_dir=/home/users/allstaff/wong.b/onj_DURATION/R_scpript
	mkdir -p ${job_dir}
	cd ${vcf_dir}
for vcf in $vcf_dir/*/*/gripss/*somatic.filtered.vcf.gz ; do
        sample_name=$(basename $vcf T.gripss.somatic.filtered.vcf.gz)
	job_file="${job_dir}/${sample_name}.job"

	echo  "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=26
#SBATCH --job-name=${sample_name}_GRIPSS_circos_plot
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --error=${sample_name}_slurm-%j.out
# THIS IS A GENERATED SCRIPT
# MESO_GRIPSS_circos_plot

cd ${vcf_dir}
module remove R
module add R/4.0.0

Rscript ${R_script_dir}/DURATION.GRIPSS_circos_plot.R ${sample_name} ${vcf_dir}/${sample_name}_process/${sample_name}_g-p-l_out
echo "${sample_name} circos plot ready"
Rscript ${R_script_dir}/extract_seq_from_gripss.R ${sample_name} ${vcf_dir}/${sample_name}_process/${sample_name}_g-p-l_out
echo "${sample_name} GRIDSS assembly sequence ready"
#design primers for top 30 GRIDSS call
cd ${primer3_dir}/src
head -n 120 ${vcf_dir}/${sample_name}_process/${sample_name}_g-p-l_out/gripss/${sample_name}_vcf_to_p3.txt | ./primer3_core --strict_tags --p3_settings_file=${primer3_dir}/settings_files/biorad_settings_with_mispriming_library.txt --output=${vcf_dir}/${sample_name}_process/${sample_name}_g-p-l_out/gripss/${sample_name}_p3_out 
echo "${sample_name} primer designed"
head -n 120 ${vcf_dir}/${sample_name}_process/${sample_name}_g-p-l_out/gripss/${sample_name}_vcf_to_p3_overlap.txt | ./primer3_core --strict_tags --p3_settings_file=${primer3_dir}/settings_files/biorad_settings_with_mispriming_library.txt --output=${vcf_dir}/${sample_name}_process/${sample_name}_g-p-l_out/gripss/${sample_name}_p3_out_overlap 
echo "${sample_name} overlap primer designed"
Rscript ${R_script_dir}/p3_out_to_fasta.R ${sample_name} ${vcf_dir}/${sample_name}_process/${sample_name}_g-p-l_out
echo "${sample_name} fasta for top 5  SVs ready"" > $job_file
	echo sbatch  $job_file
	
done
