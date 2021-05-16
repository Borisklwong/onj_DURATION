vcf_folder=/home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing
gripss_dir=/home/users/allstaff/wong.b/onj_DURATION/tools/gripss/gripss-1.11.jar
	
for vcf in $vcf_folder/*_process/*_g-p-l_out/gridss/*.gridss.unfiltered.vcf.gz ; do
        sample=$(basename $vcf T.gridss.unfiltered.vcf.gz)
	tumour=${sample}
	vcf_out_dir=${vcf_folder}/${tumour}_process/${tumour}_g-p-l_out/gripss/${tumour}_49_PON
        scriptname=$outdir/DURATION_gripss_49_PON_${tumour}.sh
	job_dir=${vcf_folder}/49_PON_job
        mkdir -p $job_dir
	mkdir -p $vcf_out_dir
	job_file="${job_dir}/${tumour}_gripss_49_PON.job"
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
   -tumor ${tumour}T \
   -ref_genome /home/users/allstaff/wong.b/onj_DURATION/tools/ref_data/gpl_ref_data_hg37/refgenomes/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta \
   -breakend_pon /home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing/PON_out/${tumour}T_pon_single_breakend.bed \
   -breakpoint_pon /home/users/allstaff/wong.b/onj_DURATION/processing/BGI_DURATION_processing/PON_out/${tumour}T_pon_breakpoint.bedpe  \
   -breakpoint_hotspot /home/users/allstaff/wong.b/onj_DURATION/tools/gridss/known_fusions.37.bedpe \
   -input_vcf ${vcf_folder}/${tumour}_process/${tumour}_g-p-l_out/gridss/${tumour}T.gridss.unfiltered.vcf.gz \
   -output_vcf ${vcf_out_dir}/${tumour}T.gripss.somatic.49.PON.vcf.gz" > $job_file
	echo sbatch  $job_file
	
done