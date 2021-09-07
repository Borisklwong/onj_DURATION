#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=26
#SBATCH --job-name=fastp
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#Boris Wong
#
#GENERATING fastp
#align with 
#illumina Homo_sapiens.GRCh37


data_folder=/home/users/allstaff/wong.b/onj_DURATION/raw_data/BGI_DURATION/FFPE_FF_BC_120x/Raw_Fastq/WGS_2021_016

for fq1 in $data_folder/*/*_1.fq.gz ; do
        samplename=$(basename $(dirname $fq1))
        filename=$(basename $fq1 .fq.gz)
        fq2=${fq1/_1.fq/_2.fq}
	fq1_clean=${filename}_clean_1.fq.gz
	fq2_clean=${filename}_clean_2.fq.gz
	out_dir=/home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing/120x/${samplename}_process
	job_dir=${out_dir}/.job
        job_file=${job_dir}/${sample_name}_fastp.sh
		echo Generating ${job_file}
        mkdir -p ${out_dir}
	mkdir -p ${job_dir}
echo  "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=26
#SBATCH --job-name=fastp
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#Boris Wong
#
#GENERATING fastp
echo "running fastp" for ${samplename}
cd ${out_dir}
~/onj_DURATION/tools/fastp/fastp \
-i ${fq1} \
-I ${fq2} \
--adapter_sequence=AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
--adapter_sequence_r2=AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
-o ${fq1_clean} \
-O ${fq2_clean} \
-j ${filename}_fastp.json \
-h ${filename}_fastp.html" > $job_file
	sbatch $job_file
	
done