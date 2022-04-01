gripss_dir=/home/users/allstaff/wong.b/onj_DURATION/tools/gripss/gripss-1.11.jar
date=23MAR2022_MESO
file_dir=/home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing/${date}
suffix=T.gridss.unfiltered.vcf.gz
job_dir=${file_dir}/job

for vcf in $(find $file_dir -name *$suffix); do
        sample=$(basename $vcf $suffix)
	tumour=${sample}T
	vcf_dir=$(dirname $vcf)
	vcf_out_dir=${file_dir}/${sample}_process/gripss_with_PON
	job_file="${job_dir}/${sample_name}_$(basename $(dirname $(dirname $vcf)))_gripss_2nd.job"
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

java -Xms4G -Xmx16G -cp ${gripss_dir} com.hartwig.hmftools.gripss.GripssHardFilterApplicationKt \
   -input_vcf ${vcf_out_dir}/${tumour}.gripss.somatic.vcf.gz \
   -output_vcf ${vcf_out_dir}/${tumour}.gripss.somatic.filtered.vcf.gz" > $job_file
	echo sbatch  $job_file
	
done