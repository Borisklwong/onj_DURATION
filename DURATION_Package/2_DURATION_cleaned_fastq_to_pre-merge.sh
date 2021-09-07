#!/bin/bash
#GENERATING DURATION_cleaned_fastq_to_pre-merge
#align with illumina Homo_sapiens.GRCh37
threads=8
ref=/home/users/allstaff/wong.b/onj_DURATION/tools/ref_data/gpl_ref_data_hg37/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta
process_data_folder=/home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing/120x

for fq1 in $process_data_folder/*/*_1.fq.gz ; do
        samplenames=$(basename $(dirname $fq1))
	samplename=${samplenames/_process/}
        filename=$(basename $fq1 .fq.gz)
        fq2=${fq1/_1.fq/_2.fq}
        scp_outdir=/home/users/allstaff/wong.b/onj_DURATION/scripts/BGI_DURATION/120x/${samplename}_script
        scriptname=$scp_outdir/${filename}_wgs_to_pre-merge.sh
	out_dir=/home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing/120x/${samplename}_process
        echo Generating $scriptname
        mkdir -p ${scp_outdir}
        mkdir -p ${out_dir}
        cat > $scriptname << EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=$threads
#SBATCH --job-name=DURATION_cleaned_fastq_to_pre-merge
#SBATCH --time=48:00:00
#SBATCH --mem=10G
# THIS IS A GENERATED SCRIPT
#DURATION_cleaned_fastq_to_pre-merge
#
echo ${samplename}
echo ${fq1}
echo ${fq2}
module add bwa/0.7.15 samtools/1.9 R/4.0.0 
mkdir -p ${out_dir}
#
cd ${out_dir}
# 
echo "Aligning with hg19 and runing fixamte"
bwa mem -t $threads -R '@RG\tID:${filename}\tSM:${samplename}' $ref $fq1 $fq2 | samtools fixmate -m -O BAM - ${filename}.bam
#
echo "Sorting bam"
samtools sort -@ $threads -O BAM -l 0 ${filename}.bam > ${filename}_sort.bam
#
EOF
done