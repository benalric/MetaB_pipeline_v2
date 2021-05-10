### **Step 2: Filtering and trimming**
In sequence data, low-quality sequences can contain unexpected and misleading errors, and Illumina sequencing quality tends to drop off at the end of reads.
Therefore, before chosing sequence variants, we trim reads where their quality scores begin to drop (the `truncLen` and `truncQ` values) and remove any low-quality reads that are left over after we have finished trimming (the `maxEE` value).

In `truncLen=c(xxx, yyy)` of the function `filterAndTrim()`, `xxx` refers to the forward read truncation length, whereas `yyy` refers to the reverse read truncation length. These parameters are determined from the inspection of quality profiles.

These parameters should be stored in a .csv file (*dada2_pipeline_parameters.csv*) located in your workspace.

| <span> |
| :------------------------------------------------------------------------------------------------------------ |
| **WARNING:** THESE PARAMETERS ARE NOT OPTIMAL FOR ALL DATASETS. Make sure you determine the trim and filtering parameters for your data. The settings used here are appropriates for MiSeq runs that are 2x250 bp. |
| <span> |

```bash
#!/bin/bash
#SBATCH --job-name=dad2_02_filter_trim_a.qsub                 # job name
#SBATCH --mail-type=BEGIN,END                                 # send a mail at the begining/end of job
#SBATCH --mail-user=balric@sb-roscoff.fr                      # where to send mail
#SBATCH --cpus-per-task 4                                     # number of CPUs required per task
#SBATCH --mem 10GB                                            # memory per processor
#SBATCH --workdir=/shared/projects/indigene/PHYTOPORT/18SV4   # set working directory
#SBATCH -p fast                                               # partition
#SBATCH -o dada2_02_filter_trim_a.out                         # output file
#SBATCH -e dada2_02_filter_trim_a.err                         # error file

module load r/4.0.3

srun Rscript /shared/projects/indigene/PHYTOPORT/18SV4/script/dada2_02_filter_trim_a.R
```

```r
# ---------------------------------------------------------------- #
# This section corresponds to the Rscript dada2_02_filter_trim_a.R #
# ---------------------------------------------------------------- #


# Some packages must be installed and loaded
library(dada2)
library(data.table)
library(magrittr)


# # Set up pathway to forward and reverse fastq files
path <- "/shared/projects/indigene/PHYTOPORT/18SV4"
pathF <- "/shared/projects/indigene/PHYTOPORT/18SV4/FWD" 
pathR <- "/shared/projects/indigene/PHYTOPORT/18SV4/REV"

# Create folders for dada2 results
dir.create((paste0(path, "/dada2/log")))

# Load the table of dada2 parameters
dada2_param <- read.csv2(paste0(path, "/dada2_pipeline_parameters.csv"), header = TRUE, stringsAsFactors = FALSE)

# Get a list of all fastq files in the work directory and separate FWD and REV 
# Forward and reverse fastq files have format: SAMPLENAME_Cut1_trimmed.fastq.gz and SAMPLENAME_Cut1_trimmed.fastq.gz
fastqFs <- sort(list.files(pathF, pattern = "fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern = "fastq.gz"))

# Select file with a size above to 1000 bytes
fastqFs <- fastqFs[file.size(file.path(pathF, fastqFs)) > 1000]
fastqRs <- fastqRs[file.size(file.path(pathR, fastqRs)) > 1000]
if (length(fastqFs) != length(fastqRs)) {
  stop("Forward and reverse files do not match.")
}

# Arguments for filter and trim
args <- dada2_param[dada2_param$step == "step01", ]
trunc_fwd <- args[args$vairable == "TRUNC_FWD", 3] %>% as.numeric
trunc_rev <- args[args$vairable == "TRUNC_REV", 3] %>% as.numeric
maxee <- args[args$vairable == "MAXEE", 3] %>% as.numeric

# Set file paths where files will be stored
pathF <- "/shared/projects/indigene/PHYTOPORT/18SV4/FWD" 
pathR <- "/shared/projects/indigene/PHYTOPORT/18SV4/REV"
filtpathF <- file.path(pathF, "filtered") 
filtpathR <- file.path(pathR, "filtered")

# Filtering and trimming
for (i in fastqFs) {
  dada2::filterAndTrim(fwd = file.path(pathF, i), filt = file.path(filtpathF, i),
                                     rev = file.path(pathR, i), filt.rev = file.path(filtpathR, i), 
                       truncLen = c(trunc_fwd, trunc_rev),
                       maxEE = maxee, maxN = 0, compress = TRUE, verbose = TRUE, multithread = FALSE) -> out
  
  # Save filtered and trimmed table in the folder of dada2 results
  data.table(out, keep.rownames = TRUE) %>% 
    setnames(names(data.table(out, keep.rownames = TRUE))[1], c("sample")) %>% 
    fwrite(paste0(path, "/dada2/log/filter_", sub("_trimmed.+$","",basename(i)), ".csv"), sep = ";", col.names = TRUE)
rm(out)
}

# Version of packages used to build this document
sessionInfo()
```
| <span> |
| :------------------------------------------------------------------------------------------------------------ |
| **WARNING:** The script dada2_02_filter_trim_a.R produces for each sample one file .csv. Here we merge all files in one table with rows as samples and three columns  correspoding to sample names, reads.in, reads.out. |
| <span> |

```bash
#!/bin/bash
#SBATCH --job-name=dada2_02_filter_trim_b.qsub                # job name
#SBATCH --mail-type=BEGIN,END                                 # send a mail at the begining/end of job
#SBATCH --mail-user=balric@sb-roscoff.fr                      # where to send mail
#SBATCH --cpus-per-task 1                                     # number of CPUs required per task
#SBATCH --mem 4GB                                             # memory per processor
#SBATCH --workdir=/shared/projects/indigene/PHYTOPORT/18SV4   # set working directory
#SBATCH -p fast                                               # partition
#SBATCH -o dada2_02_filter_trim_b.out                         # output file
#SBATCH -e dada2_02_filter_trim_b.err                         # error file

module load r/4.0.3

srun Rscript /shared/projects/indigene/PHYTOPORT/18SV4/script/dada2_02_filter_trim_b.R
```

```r
# ---------------------------------------------------------------- #
# This section corresponds to the Rscript dada2_02_filter_trim_b.R #
# ---------------------------------------------------------------- #


# Some packages must be installed and loaded
library(data.table)
library(magrittr)


# Merge all files in one table
path <- "/shared/projects/indigene/PHYTOPORT/18SV4/dada2/log"
myfiles <- list.files(path = path, pattern = "*.csv", full.names = TRUE)
filter <- lapply(myfiles, read.csv2)
filter <- do.call('rbind', filter)
filter %>% 
  data.table() %>%
  fwrite(paste0(path, "/filter.csv"), sep = ";", col.names = TRUE)

# Version of packages used to build this document
sessionInfo()
```
