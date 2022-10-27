#!/bin/bash
#SBATCH -J GM_QC_Test
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH --mem=100M

# Testing GigaMUGA QC Prep Steps
singularity run docker://sjwidmay/muga_qc:latest code/GigaMUGA_Reference_QC_FilePrep.R
