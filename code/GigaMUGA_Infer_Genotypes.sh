#!/bin/bash
#SBATCH -J GM_InferGenos
#SBATCH --mem 100GB
#SBATCH -t 1:00:00

# Performing Background QC on GigaMUGA samples
singularity run docker://sjwidmay/muga_qc:latest code/GigaMUGA_Infer_Genotypes.R $1
