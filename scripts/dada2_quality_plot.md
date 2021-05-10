## **Step 1: Inspect read quality profiles**
It is important to get a feel for the quality of the data that we are using. To do this, we will plot the quality of samples. These plots allow us to define where the quality of the sequences falls.

| <span> |
| :------------------------------------------------------------------------------------------------------------ |
| **NOTE:** Plots summarizing the quality of the reads are generated at sample scale. |
| <sapn> |
  
  *Description of quality plots:*
  >In gray-scale is a heat map of the frequency of each quality score at each base position. The median quality score at each position id shown by the green line, and the quartiles of the quality score distribution by the orange lines. Th ered line shown the scaled proportion of reads that extend to at least that position (this is more useful for other sequecing technologies, as Illumina reads are typically all the same length, hence the flat red line).

```bash
#!/bin/bash
#SBATCH --job-name=dada2_01_quality_plot.qsub                 # job name
#SBATCH --mail-type=BEGIN,END                                 # send a mail at the begining/end of job
#SBATCH --mail-user=balric@sb-roscoff.fr                      # where to send mail
#SBATCH --cpus-per-task 1                                     # number of CPUs required per task
#SBATCH --mem 4GB                                             # memory per processor
#SBATCH --workdir=/shared/projects/indigene/PHYTOPORT/18SV4   # set working directory
#SBATCH -p fast                                               # partition
#SBATCH -o dada2_01_quality_plot.out                          # output file
#SBATCH -e dada2_01_quality_plot.err                          # error file

module load r/4.0.3

srun Rscript /shared/projects/indigene/PHYTOPORT/18SV4/script/dada2_01_quality_plot.R
```

```r
# --------------------------------------------------------------- #
# This section corresponds to the Rscript dada2_01_quality_plot.R #
# --------------------------------------------------------------- #


# Some packages must be installed and loaded
library(dada2)
library(magrittr)
library(ggplot2)


# Set up pathway to forward and reverse fastq files
pathF <- "/shared/projects/indigene/PHYTOPORT/18SV4/FWD" 
pathR <- "/shared/projects/indigene/PHYTOPORT/18SV4/REV"

# Create folders for quality profiles
dir.create(paste0(pathF, "/qualityplot"))
dir.create(paste0(pathR, "/qualityplot"))

# Set file paths where quality plots will be stored
filtpathF <- file.path(pathF, "qualityplot") 
filtpathR <- file.path(pathR, "qualityplot")

# Get a list of all fastq files in the work directory and separate FWD and REV 
# Forward and reverse fastq files have format: SAMPLENAME_Cut1_trimmed.fastq.gz and SAMPLENAME_Cut2_trimmed.fastq.gz
fastqFs <- sort(list.files(pathF, pattern = "fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern = "fastq.gz"))

# Select file with a size above to 1000 bytes
fastqFs <- fastqFs[file.size(file.path(pathF, fastqFs)) > 1000]
fastqRs <- fastqRs[file.size(file.path(pathR, fastqRs)) > 1000]
if (length(fastqFs) != length(fastqRs)) {
  stop("Forward and reverse files do not match.")
}

# Identify the name of the run (here one run: PHYTOPORT)
runs <- sub("^PHP_[B, F]\\d+_\\d+_\\d+_([^_]+).+$", "\\1", fastqFs) %>% unique

# Quality plots to the sample level
pdf(file.path(filtpathF, "indiv_F_Qplots.pdf"))
for (i in file.path(pathF, fastqFs)) {
  print(plotQualityProfile(i))
}
dev.off()

pdf(file.path(filtpathR, "indiv_R_Qplots.pdf"))
for (i in file.path(pathR, fastqRs)) {
  print(plotQualityProfile(i))
}
dev.off()

# Version of packages used to build this document
sessionInfo()
```
**Visualization example of the quality plot**

![Figure 1. The example dataset: predictor variables and occurrence for
four species.](https://github.com/benalric/Metab_pipeline_v2/tree/main/figures)

<img src="figures/indiv_F_Qplots.png" />


*Description of quality plots:*
>In gray-scale is a heat map of the frequency of each quality score at each base position. The median quality score at each position id shown by the green line, and the quartiles of the quality score distribution by the orange lines. Th ered line shown the scaled proportion of reads that extend to at least that position (this is more useful for other sequecing technologies, as Illumina reads are typically all the same length, hence the flat red line).


[Step 2: Filtering and trimming](https://github.com/benalric/Metab_pipeline_v2/tree/main/src/scripts/dada2_filtering_trimming.md)
