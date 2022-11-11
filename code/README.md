# Code

## Overview

## ![](images/readme_diagram.png)

As shown above (and described [here](https://github.com/sam-widmayer/MUGA_reference_data/blob/main/README.md)), GigaMUGA reference genotypes had quality control analyses performed using a series of scripts deployed on the HPC infrastructure at JAX. Here we describe the function and usage of each script in the process of generating GigaMUGA reference genotypes.

### Dependencies

These QC steps were performed on an HPC using SLURM/sbatch routines. Recreating the analysis on other scheduling infrastructures may require adapting the following workflow.

We built a [Docker container](https://hub.docker.com/repository/docker/sjwidmay/muga_qc) with all required software to run all analyses. It is pulled within each shell script using singularity:

```{bash}
singularity run docker://sjwidmay/muga_qc:latest
```

All of the following commands were executed within the **MUGA_reference_data** directory.

### Preparing Chromosome-level Genotype Files

The script below runs `GigaMUGA_Reference_QC_FilePrep.R` for each chromosome in the input genotype file, generating chromosome-level genotypes stored as .fst files for more efficient reading/writing.

```{bash}
bash code/GigaMUGA_FilePrep_submit.sh
```

### Identifying Poor Markers and Performing Sex Checks

### Generating Chromosome-level Genotype Files for Dendrogram of Samples

### Aggregating Consensus Genotypes for CC/DO Founder Strains

### Summarizing Results and Generating Reference Files
