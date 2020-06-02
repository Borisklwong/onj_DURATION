#!/bin/bash
#PBS -l nodes=1:ppn=26,mem=32gb,walltime=144:00:00
#Boris Wong
#
#GENERATING WGS raw data to Pre-GRIDSS script
#
ref=~/projects/reference_genomes/human/hg19.fa

data_folder=~/projects/onj_DURATION/raw_data/BGI_samples/BGI_MESO/Clean

for fq1 in $data_folder/*/*_1.fq.gz ; do
        name=$(basename $(dirname $fq1))
        fq2=$(dirname $fq1)/$(basename $fq1 _1.fq.gz)_2.fq.gz
        outdir=~/projects/onj_DURATION/scripts/MESO/${name}_script
        scriptname=$outdir/$(basename $fq1 _1.fq.gz)_wgs_to_pregridss.sh
        echo Generating $scriptname
        mkdir -p $outdir
        cat > $scriptname << EOF
#!/bin/bash
#PBS -l nodes=1:ppn=26,mem=32gb,walltime=144:00:00
# THIS IS A GENERATED SCRIPT
#WGS raw data to Pre-GRIDSS
#
module add bwa/0.7.15 samtools/1.9
#
mkdir -p ~/onj_DURATION/processing/MESO/${name}_process
#
cd ~/onj_DURATION/processing/MESO/${name}_process 
#
if [[ ! -f "${name}.bam" ]]; then
	echo "Aligning with hg19 and runing fixamte"
	bwa mem -t 26 $ref $fq1 $fq2 | samtools fixmate -m -O BAM - ${name}.bam || exit 1 
fi 
if [[ ! -f "${name}_sort.bam" ]]; then
	echo "Sorting bam"
	samtools sort -@ 26 -O BAM -l 0 ${name}.bam > ${name}_sort.bam || exit 1
fi
if [[ ! -f "${name}_sort_mark.bam" ]]; then
	echo "Marking duplicates"
	samtools markdup -@ 26 -O BAM ${name}_sort.bam ${name}_sort_mark.bam || exit 1
fi
if [[ ! -f "${name}*.bai" ]]; then
	echo "Indexing bam"
	samtools index -@26 ${name}_sort_mark.bam
fi


EOF
done

