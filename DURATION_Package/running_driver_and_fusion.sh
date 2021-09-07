#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=${sample_name}_GRIPSS_circos_plot
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --error=${sample_name}_slurm-%j.out

module remove R
module add R/4.0.0

Rscript ~/onj_DURATION/R_scpript/120x/driver_and_fusion_list_from_linx_purple_NSCLC.R