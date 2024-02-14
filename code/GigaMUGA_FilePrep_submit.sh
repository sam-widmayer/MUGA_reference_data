#!/bin/bash
#SBATCH -J MUGA_REF_FILEPREP
#SBATCH --mem 100GB
#SBATCH -t 1:00:00

# submit from MUGA_reference_data directory using:
# config=/projects/compsci/vmp/USERS/widmas/MUGA_reference_data/data/GigaMUGA/chrs.txt
# chrs=$(wc -l < ${config})
# sbatch --array=1-${chrs} code/GigaMUGA_FilePrep_submit.sh

# chromosome config file
config=/projects/compsci/vmp/USERS/widmas/MUGA_reference_data/data/GigaMUGA/chrs.txt

# chromosome
# find the line in the config file
chr=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${config})
echo ${chr}
echo "Running Chromosome ${chr}"

# run GigaMUGA file prep
singularity exec docker://sjwidmay/muga_qc:latest code/GigaMUGA_Reference_QC_FilePrep.R ${chr}