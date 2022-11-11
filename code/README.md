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

The expected output is the `data/GigaMUGA/GigaMUGA_reference_genotypes/` directory, containing chromosome-level genotype files for all reference samples in .fst format.

### Identifying Poor Markers and Performing Sex Checks

`GigaMUGA_sampleQC.R` identifies markers that fail to perform across many samples, flags samples with a high no-call rate across all makers, and infers sexes from sample-aggregated sex chromosome probe intensity values.

The script can be run in an interactive session using singularity:

```{bash}
singularity run docker://sjwidmay/muga_qc:latest code/GigaMUGA_sampleQC.R
```

If this exceeds the memory of the machine, the following can be used in a SLURM scheduling infrastructure and singularity.

```{bash}
bash code/GigaMUGA_sampleQC.sh
```

A successful execution of this script produces three output files:

-   `GigaMUGA_QC_Results.RData`: genotype frequencies for each marker across all samples, a list of samples and their no-call rate across all markers, and a list of samples from CC/DO founder strains and crosses between them.

-   `GigaMUGA_BadSamples_BadMarkers.RData`: a list of markers and samples with poor genotyping quality.

-   `GigaMUGA_SexCheck_Results.RData`: sex chromosome probe intensities for each sample, each sample's predicted sex based on sample ID, and their inferred sex from their sex chromosome intensities.

These three files are referenced in the next two processes, as well as the final summaries in `GigaMUGA_Reference_QC.Rmd`.

### Generating Chromosome-level Genotype Files for Dendrogram of Samples

The script below runs `GigaMUGA_FounderDendGenos.R` for each chromosome, which takes "long" chromosome-level genotype .fst files, filters them to high-quality markers and CC/DO founder strain samples, and makes them "wide", with each marker as a row and each column as a sample. These files feed into the dendrogram constructed in the `GigaMUGA_Reference_QC.Rmd`.

```{bash}
bash code/GigaMUGA_FounderDendGenos_submit.sh
```

The expected output from these scripts is the `data/GigaMUGA/GigaMUGA_founder_sample_dendrogram_genos/` directory, containing chromosome-level "wide" genotypes for constructing the sample dendrogram.

### Aggregating Consensus Genotypes for CC/DO Founder Strains

The script below runs `GigaMUGA_ConsensusGenos.R` for each chromosome. For each of these files, this script identifies samples derived from CC/DO founder strains and determines the consensus genotype for each founder at each high-quality marker. Sometimes the consensus genotype is fixed among all samples for a given strain, while other times there is a "N" call or misassigned genotype. An example of how these assignments are fixed is shown [here](https://sam-widmayer.github.io/MUGA_reference_data/MegaMUGA_Reference_QC.html#Validating_reference_sample_genetic_backgrounds) using MegaMUGA genotype data. Once consensus genotypes are called for each strain, the actual genotypes of each CC/DO founder sample, including F1 hybrids among them, are cross-validated against the expected genetic background according to the consensus genotypes. The concordance between the observed and expected genotypes is plotted in the final markdown.

```{bash}
bash code/GigaMUGA_ConsensusGenos_submit.sh
```

There are three expected outputs from this script:

-   `data/GigaMUGA/GigaMUGA_founder_consensus_genotypes/`: Chromosome-level CC/DO founder strain consensus genotypes for all high-quality markers.

-   `data/GigaMUGA/GigaMUGA_founder_sample_concordance/`: Chromosome-level measurements of the percentage of genotypes from each sample that match the expected consensus genotype.

-   `data/GigaMUGA/GigaMUGA_founder_sample_genotypes/`: Because some sample genotypes may have been re-inferred from probe intensities based on clustering with other strains, chromosome-level sample genotype files for each CC/DO founder sample are also generated.

### Summarizing Results and Generating Reference Files

`GigaMUGA_Reference_QC.Rmd` synthesizes the results contained in the `data/GigaMUGA/` directories and displays the results. It also generates all desired output files found on the [main workflowr page](https://sam-widmayer.github.io/MUGA_reference_data/index.html).
