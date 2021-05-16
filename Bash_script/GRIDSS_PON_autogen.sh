vcf_folder=/home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing
grdiss_dir=/home/users/allstaff/wong.b/onj_DURATION/tools/gridss/gridss-2.11.1-gridss-jar-with-dependencies.jar
ref_genome=/home/users/allstaff/wong.b/onj_DURATION/tools/ref_data/gpl_ref_data_hg37/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta
	
for vcf in $vcf_folder/PON_vcf/*.gripss.somatic.filtered.vcf.gz; do
        sample=$(basename $vcf .gripss.somatic.filtered.vcf.gz)
	PON_out_dir=${vcf_folder}/PON_out
	job_dir=${vcf_folder}/PON_job
        mkdir -p $job_dir
	mkdir -p $PON_out_dir
	cd ${vcf_folder}/PON_vcf
	job_file="${job_dir}/${sample}_PON.job"
        echo  "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=gripss
#SBATCH --time=48:00:00
#SBATCH --mem=12G
#Boris Wong
#
#Create PON with all other samples in the file
cd ${vcf_folder}/PON_vcf
java -Xmx8g \
-cp ${grdiss_dir} \
gridss.GeneratePonBedpe \
\$(ls -1 -I "*${sample}*vcf.gz" | awk '{print \"INPUT=\" \$0}') \
O=${PON_out_dir}/${sample}_pon_breakpoint.bedpe \
SBO=${PON_out_dir}/${sample}_pon_single_breakend.bed \
NORMAL_ORDINAL=0 \
REFERENCE_SEQUENCE=$ref_genome" > $job_file
echo sbatch  $job_file

done