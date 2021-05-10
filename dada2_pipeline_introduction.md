# **Tutorial for metabarcode analysis**
In this repository, you will find the tutotrial explaining how to process metabarcoding data with the [`dada2`](https://www.nature.com/articles/nmeth.3869) suite as implementated in [R](https://benjjneb.github.io/dada2/dada-installation.html) and generate tables of amplicon sequence variants (i.e., ASV). This tutorial is adapted from the `dada2` [tutorial](https://benjjneb.github.io/dada2/tutorial.html).

## **Overview of `dada2` pipeline** ###
The starting point of the `dada2` pipeline is a set of demultiplexed fastq files corresponding to the samples in our amplicon sequencing study. That is, `dada2` expects there to be an individual fastq file for each sample (or two fastq files, one forward and one reverse for each sample). Demultiplexing is often performed at the sequencing center, but if that has not been done there are a variety of tools to accomplish this. Once demultiplexed fastq files without non-biological nucleotides are in hand, the `dada2` pipeline proceeds in eight steps as follow:

* [Step 1: Inspect read quality profiles](https://github.com/benalric/Metab_pipeline_v2/blob/main/dada2_quality_plot.md)
* [Step 2: Filtering and trimming](https://github.com/benalric/Metab_pipeline_v2/blob/main/dada2_filtering_trimming.md)
* [Step 3: Learn error rates](https://github.com/benalric/Metab_pipeline_v2/blob/main/dada2_learn_error.md)
* [Step 4: Dereplicate](https://github.com/benalric/Metab_pipeline_v2/blob/main/dada2_denoising_merging.md)
* [Step 5: Infer sample composition](https://github.com/benalric/Metab_pipeline_v2/blob/main/dada2_denoising_merging.md)
* [Step 6: Merge paired reads](https://github.com/benalric/Metab_pipeline_v2/blob/main/dada2_denoising_merging.md)
* [Step 7: Make sequence table](https://github.com/benalric/Metab_pipeline_v2/blob/main/dada2_denoising_merging.md)
* [Step 8: Remove chimeras](https://github.com/benalric/Metab_pipeline_v2/blob/main/dada2_chimeras.md)
  
The output of pipeline is a table of ASVs with rows corresponding to samples and columns to amplicon. In this matrix, the value of each entry is the number of times ASV was observed in that sample. This table is analogous to the traditional OTU table, except at higher resolution, i.e. exact amplicon sequence variants rather than (usually 97%) clusters of sequencing reads.

| <span> |
| :------------------------------------------------------------------------------------------------------------ |
| **NOTE:** Each step of our analysis pipeline is composed in two parts, a first part (in bash) to launch the second part (in R) which is the core part of the analysis. Each step should be made separately. |
| <span> |
