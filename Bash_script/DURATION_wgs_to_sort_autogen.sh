#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=26
#SBATCH --job-name=DURATION_fastq_bwa-mem_sort
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#Boris Wong
#
#GENERATING DURATION_fastq_bwa-mem_sort
#align with 
#illumina Homo_sapiens.GRCh37

ref=/home/users/allstaff/wong.b/onj_DURATION/tools/ref_data/gpl_ref_data_hg37/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta

data_folder=/home/users/allstaff/wong.b/onj_DURATION/raw_data/BGI_DURATION/raw_data

for fq1 in $data_folder/*/*_1.fq.gz ; do
        samplename=$(basename $(dirname $fq1))
        filename=$(basename $fq1 .fq.gz)
        fq2=${fq1/_1.fq/_2.fq}
        outdir=/home/users/allstaff/wong.b/onj_DURATION/scripts/BGI_DURATION/${samplename}_script
        scriptname=$outdir/${filename}_fastq_bwa-mem_sort.sh
        echo Generating $scriptname
        mkdir -p $outdir
        cat > $scriptname << EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=26
#SBATCH --job-name=DURATION_fastq_bwa-mem_sort
#SBATCH --time=48:00:00
#SBATCH --mem=64G
# THIS IS A GENERATED SCRIPT
#URATION_fastq_bwa-mem_sort, for downstream megre and 
#
module add bwa/0.7.15 samtools/1.9
#
mkdir -p ~/onj_DURATION/processing/BGI_DURATION_processing/${samplename}_process
#
cd ~/onj_DURATION/processing/BGI_DURATION_processing/${samplename}_process
#
if [[ ! -f "${filename}.bam" ]]; then
	echo "Aligning with hg19 and runing fixamte"
	bwa mem -t 26 -R '@RG\tID:${filename}\tSM:${samplename}' $ref $fq1 $fq2 | samtools fixmate -m -O BAM - ${filename}.bam || exit 1 
fi 
if [[ ! -f "${filename}_sort.bam" ]]; then
	echo "Sorting bam"
	samtools sort -@ 26 -O BAM -l 0 ${filename}.bam > ${filename}_sort.bam || exit 1
fi
EOF
done