date=23MAR2022_MESO

gripss_dir=/home/users/allstaff/wong.b/onj_DURATION/tools/gripss/gripss-1.11.jar
file_dir=/home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing/${date}
suffix=T.gridss.unfiltered.vcf.gz
job_dir=${file_dir}/job

for vcf in $(find $file_dir -name *$suffix); do
        sample=$(basename $vcf $suffix)
	tumour=${sample}T
	vcf_dir=$(dirname $vcf)
	vcf_out_dir=${file_dir}/${sample}_process/gripss_with_PON
	job_file="${job_dir}/$(basename $(dirname $(dirname $vcf)))_gripss_1st.job"
	mkdir -p $job_dir
	mkdir -p $vcf_dir
	mkdir -p $vcf_out_dir
        echo  "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=gripss
#SBATCH --time=48:00:00
#SBATCH --mem=12G
#Boris Wong
#
#Running gripss after GRIDSS for somatic filtering

java -Xms4G -Xmx16G -cp ${gripss_dir} com.hartwig.hmftools.gripss.GripssApplicationKt \
   -tumor ${tumour} \
   -ref_genome /home/users/allstaff/wong.b/onj_DURATION/tools/ref_data/gpl_ref_data_hg37/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta \
   -breakend_pon /home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing/PON_out_with_hartwig_PON+gridss_vcf/50_FFPE_pon_single_breakend.bed \
   -breakpoint_pon /home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing/PON_out_with_hartwig_PON+gridss_vcf/50_FFPE_pon_breakpoint.bedpe  \
   -breakpoint_hotspot /home/users/allstaff/wong.b/onj_DURATION/tools/gridss/gridss-2.12.1/known_fusions.37.bedpe \
   -input_vcf ${vcf} \
   -output_vcf ${vcf_out_dir}/${tumour}.gripss.somatic.vcf.gz" > $job_file
	echo sbatch  $job_file
	
done