#!/bin/bash
#GENERATING DURATION_merge_to_pre-g-p-l
#align with 
#illumina Homo_sapiens.GRCh37

ref=/home/users/allstaff/wong.b/onj_DURATION/tools/ref_data/gpl_ref_data_hg37/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta

bam_dir=/home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing/120x

for folder in $bam_dir/* ; do
        samplenames=$(basename $folder)
	samplename=${samplenames/_process}
        scp_outdir=/home/users/allstaff/wong.b/onj_DURATION/scripts/BGI_DURATION/120x/${samplename}_script
        scriptname=$scp_outdir/${samplename}_merge_to_pre-g-p-l.sh
	out_dir=/home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing/120x/${samplename}_process
        echo Generating $scriptname
        mkdir -p ${scp_outdir}
	mkdir -p ${out_dir}
        cat > $scriptname << EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --job-name=DURATION_merge_to_pre-g-p-l
#SBATCH --time=48:00:00
#SBATCH --mem=32G
# THIS IS A GENERATED SCRIPT
#DURATION_merge_to_pre-g-p-l
#
module add bwa/0.7.15 samtools/1.9
#
mkdir -p ${out_dir}
#
cd ${out_dir}
#
echo "Merge all sorted bam"
samtools merge -@ 26 ${samplename}_merged.bam ./*_sort.bam
#
echo "Marking duplicates"
samtools markdup -@ 26 -O BAM ${samplename}_merged.bam  ${samplename}_merged_mark.bam
#
echo "Indexing bam"
samtools index -@26 ${samplename}_merged_mark.bam

EOF
done