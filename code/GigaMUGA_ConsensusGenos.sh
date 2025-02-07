#!/bin/bash
#SBATCH -J GM_BackgroundChecks
#SBATCH --mem 100GB
#SBATCH -t 4:00:00

# Performing Background QC on GigaMUGA samples
singularity run docker://sjwidmay/muga_qc:latest code/GigaMUGA_ConsensusGenos.R $1