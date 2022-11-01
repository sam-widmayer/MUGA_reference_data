#!/bin/bash
#SBATCH -J GM_QC_Test
#SBATCH --mem 100GB
#SBATCH -t 4:00:00

# Testing GigaMUGA QC Prep Steps
singularity run docker://sjwidmay/muga_qc:latest code/GigaMUGA_Reference_QC_FilePrep.R $1