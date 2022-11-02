#!/bin/bash
#SBATCH -J GM_Dend_Genos
#SBATCH --mem 100GB
#SBATCH -t 4:00:00

# Performing Marker QC on GigaMUGA samples
singularity run docker://sjwidmay/muga_qc:latest code/GigaMUGA_FounderDendGenos.R $1