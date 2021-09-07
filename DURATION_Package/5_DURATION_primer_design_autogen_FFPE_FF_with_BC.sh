#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=DURATION_gripss_to_primer
#SBATCH --time=48:00:00
#SBATCH --mem=16G
#Boris Wong
#
#For post stand alone GRIDSS and GRIPSS

# Local system locations
	vcf_dir=/home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing
	job_dir=${vcf_dir}/job
	primer3_dir=/home/users/allstaff/wong.b/onj_DURATION/tools/primer3-2.4.0
	R_script_dir=/home/users/allstaff/wong.b/onj_DURATION/R_scpript/FFPE_vs_FF_R/FFPE_only_R
	mkdir -p ${job_dir}
for vcf in $vcf_dir/*_S_process/*_g-p-l_out/gripss/*_49_PON_gridss_raw_with_hartwig_PON_new/*.gripss.somatic.49.PON.filtered.vcf.gz ; do
        sample_name=$(basename $vcf T.gripss.somatic.49.PON.filtered.vcf.gz)
	out_dir=${vcf_dir}/post_gridss_FFPE_only_primers
	gripss_vcf_dir=${vcf_dir}/${sample_name}_process/${sample_name}_g-p-l_out/gripss/${sample_name}_49_PON_gridss_raw_with_hartwig_PON_new
	blat_dir=/home/users/allstaff/wong.b/onj_DURATION/tools/blat
	job_file="${job_dir}/${sample_name}_post_gridss_FFPE_only.job"
	mkdir -p $out_dir
	echo  "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=${sample_name}_GRIPSS_circos_plot
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --error=${sample_name}_slurm-%j.out
# THIS IS A GENERATED SCRIPT
# post_gridss_FFPE_only

cd ${out_dir}
module remove R
module add R/4.0.0

Rscript ${R_script_dir}/DURATION.circos_FFPE_only.R ${sample_name} ${vcf_dir}/${sample_name}_process/${sample_name}_g-p-l_out
echo "${sample_name} circos plot ready"
Rscript ${R_script_dir}/extract_seq_from_gripss_FFPE_only.R ${sample_name} ${vcf_dir}/${sample_name}_process/${sample_name}_g-p-l_out
echo "${sample_name} GRIDSS assembly sequence ready"
#design primers for top 30 GRIDSS call
cd ${primer3_dir}/src
head -n 120 ${gripss_vcf_dir}/${sample_name}_vcf_to_p3.txt | ./primer3_core --strict_tags --p3_settings_file=${primer3_dir}/settings_files/biorad_settings_with_mispriming_library.txt --output=${gripss_vcf_dir}/${sample_name}_p3_out 
echo "${sample_name} primer designed"
head -n 120 ${gripss_vcf_dir}/${sample_name}_vcf_to_p3_overlap.txt | ./primer3_core --strict_tags --p3_settings_file=${primer3_dir}/settings_files/biorad_settings_with_mispriming_library.txt --output=${gripss_vcf_dir}/${sample_name}_p3_out_overlap 
echo "${sample_name} overlap primer designed"
Rscript ${R_script_dir}/p3_out_to_fasta.R ${sample_name} ${vcf_dir}/${sample_name}_process/${sample_name}_g-p-l_out ${out_dir}
echo "${sample_name} fasta for top 10  SVs ready"
cd ${blat_dir}
./blat -fastMap -tileSize=8 -oneOff=1 -minMatch=1 -minScore=10 -repMatch=16000 -t=DNA \
-q=DNA ./hg19_blat/hg19.fa \
${vcf_dir}/${sample_name}_process/${sample_name}_g-p-l_out/gripss/${sample_name}_49_PON_gridss_raw_with_hartwig_PON_new/${sample_name}_p3_out_full.fasta \
${vcf_dir}/${sample_name}_process/${sample_name}_g-p-l_out/gripss/${sample_name}_49_PON_gridss_raw_with_hartwig_PON_new/${sample_name}_p3_out_full.psl
echo "${sample_name} primer blat result ready"
Rscript ${R_script_dir}/blat_filter.R ${sample_name} ${vcf_dir}/${sample_name}_process/${sample_name}_g-p-l_out
echo "${sample_name} blat result filtered"" > $job_file
	echo sbatch $job_file
	
done