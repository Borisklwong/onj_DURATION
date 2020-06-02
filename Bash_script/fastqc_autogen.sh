#!/bin/bash
#PBS -l nodes=1:ppn=26,mem=32gb,walltime=144:00:00
#Boris Wong
#
#GENERATING fastqc scrpit
#

data_folder=~/projects/onj_DURATION/raw_data/BGI_samples/BGI_MESO/Clean

for fq1 in $data_folder/*/*_1.fq.gz ; do
        name=$(basename $(dirname $fq1))
        fq2=$(dirname $fq1)/$(basename $fq1 _1.fq.gz)_2.fq.gz
        outdir=~/projects/onj_DURATION/scripts/MESO/${name}_script
        scriptname=$outdir/$(basename $fq1 _1.fq.gz)_fastqc.sh
        echo Generating $scriptname
        mkdir -p $outdir
        cat > $scriptname << EOF
#!/bin/bash
#PBS -l nodes=1:ppn=26,mem=32gb,walltime=144:00:00
# THIS IS A GENERATED SCRIPT
module add fastqc/0.11.8

mkdir -p ~/onj_DURATION/processing/MESO/${name}_process

fastqc -t 26 $fq1 $fq2 -o ~/onj_DURATION/processing/MESO/${name}_process

EOF
done

