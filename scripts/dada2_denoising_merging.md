## **Steps 4-7: Denoising, and merging**
After it understands the error rates, a dereplication step is required to condense the data by collapsing together all reads that encode the same sequence, which significantly reduces later computation times. Then, using the dereplicated data and error rates, `dada2` will infer the ASVs in our data.

| <span> |
| :------------------------------------------------------------------------------------------------------------ |
| **NOTE:** In `dada` algorithm, for each run, the samples are pooled together with the argument `pool` = TRUE, which can increase the sensitivity to rare variants. |
| <span> |

Finally, we will merge the corresponding forqard and reverse reads to create a list of the fully denoised sequences and create a sequence table for the result.

| <span> |
| :------------------------------------------------------------------------------------------------------------ |
| **NOTE:** Merging is performed by aligning the denoised forward reads with the reverse-compelment of the corresponding denoised reverse reads, and then constructing the merged "contig" sequences. By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region. |
| <sapn> |

```bash
#!/bin/bash
#SBATCH --job-name=dada2_04_denoising_merging.qsub            # job name
#SBATCH --mail-type=BEGIN,END                                 # send a mail at the begining/end of job
#SBATCH --mail-user=balric@sb-roscoff.fr                      # where to send mail
#SBATCH --cpus-per-task 6                                     # number of CPUs required per task
#SBATCH --mem 10GB                                            # memory per processor
#SBATCH --workdir=/shared/projects/indigene/PHYTOPORT/18SV4   # set working directory
#SBATCH -p fast                                               # partition
#SBATCH -o dada2_04_denoising_merging.out                     # output file
#SBATCH -e dada2_04_denoising_merging.err                     # error file

module load r/4.0.3

srun Rscript /shared/projects/indigene/PHYTOPORT/18SV4/script/dada2_04_denoising_merging.R
```

```r
# -------------------------------------------------------------------- #
# This section corresponds to the Rscript dada2_04_denoising_merging.R #
# -------------------------------------------------------------------- #


# Some packages must be installed and loaded
library(dada2)
library(data.table)
library(doParallel)
library(foreach)
library(magrittr)
registerDoParallel(cores = 6)


# Load the table of dada2 parameters
path <- "/shared/projects/indigene/PHYTOPORT/18SV4"
dada2_param <- read.csv2(paste0(path, "/dada2_pipeline_parameters.csv"), header = TRUE, stringsAsFactors = FALSE)

# Define the following path variable so that it points to the extracted directory
pathF <- "/shared/projects/indigene/PHYTOPORT/18SV4/FWD" 
pathR <- "/shared/projects/indigene/PHYTOPORT/18SV4/REV"
pathErr <- "/shared/projects/indigene/PHYTOPORT/18SV4/dada2/log"

# Arguments for denoising and merging
args <- dada2_param[dada2_param$step == "step02", ]
min_read_nb <- args[args$vairable == "MIN_READ_NUM", 3]

# Track reads through the pipeline
# As a final check of our progress, we'll look at the number of reads that made it through each step in the pipeline
getN <- function(x) sum(getUniques(x))
getA <- function(x) length(getUniques(x))

# selection based on the number of reads
tmp <- fread(paste0(path, "/dada2/log/filter.csv"))
tmpR <- tmpF <- tmp[reads.out >= min_read_nb, sample]

# Identify the name of the run (here one run: PHYTOPORT)
runs <- sub("^PHP_[B, F]\\d+_\\d+_\\d+_([^_]+).+$", "\\1", tmp$sample) %>% unique
runs <- c(paste(runs, "Cut1", sep = "_"), paste(runs, "Cut2", sep = "_"))

# File parsing
filtpathF <- file.path(pathF, "filtered") 
filtpathR <- file.path(pathR, "filtered")
filtFsall <- paste(filtpathF, tmpF, sep = "/")
filtRsall <- paste(filtpathR, tmpR, sep = "/")
sample.names.all <- sub("_trimmed.+$", "", basename(filtFsall))
sample.namesR.all <- sub("_trimmed.+$" ,"", basename(filtRsall))
if(!identical(sample.names.all, sample.namesR.all)) {
  stop("Forward and reverse files do not match.")
}
names(filtFsall) <- sample.names.all
names(filtRsall) <- sample.names.all
filtErr <- sort(list.files(pathErr, pattern = ".rda"))

# Set seed to ensure that randomized steps are replicatable
set.seed(100)

# Loop allowing to determinate learn error rates and infer sample composition
foreach(i = runs, .packages = c("data.table","dada2")) %dopar% {
  filtFs <- filtFsall[grep(i, names(filtFsall))]
  filtRs <- filtRsall[grep(i, names(filtRsall))]
  sample.names <- sub("_trimmed.+$","", basename(filtFs))
  # Select the FWD and REV error files
  err <- filtErr[grep(i, filtErr)]
  load(paste(pathErr, err[grep("errF", err)], sep = "/"))
  load(paste(pathErr, err[grep("errR", err)], sep = "/"))
  # Sample inference and merger of paired-end reads
  mergers <- vector("list", length(sample.names))
  names(mergers) <- sample.names
  track <- vector("list", length(sample.names))
  names(track) <- sample.names
  for(sam in sample.names) {
    cat("Processing:", sam, "\n")
    # Dereplication of forward reads+
    derepF <- derepFastq(filtFs[[sam]])
    # Infer sample composition of forward reads
    ddF <- dada(derepF, err = errF, multithread = FALSE, pool = TRUE)
    # Dereplication of revers reads
    derepR <- derepFastq(filtRs[[sam]])
    # Infer sample composition of reverse reads
    ddR <- dada(derepR, err = errR, multithread = FALSE, pool = TRUE)
    # Merge forward/reverse reads
    merger <- mergePairs(ddF, derepF, ddR, derepR, trimOverhang = TRUE, verbose = TRUE)
    mergers[[sam]] <- merger
    track[[sam]] <- data.table(sample = sam, 
                               denoisedF.read = getN(ddF), denoisedR.read = getN(ddR), merged.read = getN(merger),
                               denoisedF.seq = getA(ddF), denoisedR.seq = getA(ddR), merged.seq = getA(merger))
  }
  rm(derepF); rm(derepR)
  # Construct sequence table
  track <- rbindlist(track)
  seqtab <- makeSequenceTable(mergers)
  saveRDS(track, paste0(path, "/dada2/log/track_",i,".rds")) # CHANGE ME to where you want sequence table saved
  saveRDS(seqtab, paste0(path, "/dada2/seqtab_",i,".rds")) # CHANGE ME to where you want sequence table saved
  rm(track)
  rm(seqtab)
  for(i in 1:length(mergers)) {
    mergers[[i]] <- data.frame(samples = names(mergers)[i], mergers[[i]])
  }
  mergers <- rbindlist(mergers)
  saveRDS(mergers, paste0(path, "/dada2/log/mergers_",i,".rds")) # CHANGE ME to where you want sequence table saved
}

# Version of packages used to build this document
sessionInfo()
```

[**Next, Step 8: Remove chimeras**](https://github.com/benalric/MetaB_pipeline_v2/blob/main/scripts/dada2_chimeras.md)
