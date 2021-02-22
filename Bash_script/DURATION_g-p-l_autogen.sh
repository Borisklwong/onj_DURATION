#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=26
#SBATCH --job-name=DURATION_g-p-l
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#Boris Wong
#
#GENERATING g-p-l 
#align with 
#illumina Homo_sapiens.GRCh37

# Local system locations
	bam_folder=/home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing
	install_dir=/home/users/allstaff/wong.b/onj_DURATION/tools/gridss-purple-linx/package
	ref_data=/home/users/allstaff/wong.b/onj_DURATION/tools/ref_data/gpl_ref_data_hg37
	
for bam in $bam_folder/*/*merged_mark.bam ; do
        sample=$(basename $bam _merged_mark.bam)
        outdir=/home/users/allstaff/wong.b/onj_DURATION/scripts/BGI_DURATION/${sample}_script
	data_dir=/home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing/${sample}_process
	run_dir=/home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing/${sample}_process/${sample}_g-p-l_out
        scriptname=$outdir/${sample}_g-p-l.sh
        echo Generating $scriptname
        mkdir -p $outdir
        cat > $scriptname << EOF
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=26
#SBATCH --job-name==DURATION_g-p-l
#SBATCH --time=48:00:00
#SBATCH --mem=64G
# THIS IS A GENERATED SCRIPT
# merged_sort.bam to g-p-l

cd /home/users/allstaff/wong.b/onj_DURATION/tools/gridss-purple-linx/
module remove R
module add R/4.0.0 bwa samtools java sambamba circos


rm -rf $run_dir/amber $run_dir/cobalt $run_dir/gridss $run_dir/gripss $run_dir/logs $run_dir/purple


export GRIDSS_JAR=$(find $install_dir -name '*gridss*jar')
export GRIPSS_JAR=$(find $install_dir -name '*gripss*jar')
export AMBER_JAR=$(find $install_dir -name '*amber*jar')
export COBALT_JAR=$(find $install_dir -name '*cobalt*jar')
export PURPLE_JAR=$(find $install_dir -name '*purple*jar')
export LINX_JAR=$(find $install_dir -name '*linx*jar')

bash -x gridss-purple-linx.sh \
	-o $run_dir \
	-t $data_dir/${sample}_merged_mark.bam \
	--nosnvvcf \
	-s $sample \
	--tumour_sample ${sample}T \
	--ref_dir $ref_data \
	--install_dir $install_dir \
	--cobalt_args "-tumor_only_diploid_bed DiploidRegions.hg19.bed.gz" \
		
EOF
done

