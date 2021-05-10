### **Step 8: Remove Chimeras**
Although `dada2` has searched for indel errors and substitutions, there may still be chimeric sequences in our dataset that are another important source of spurious sequences in amplicon sequencing. Chimeras are sequences that are derived form forward and reverse sequences from two different organisms becoming fused together during PCR and/or sequencing. To identify chimeras, we search for rare sequences variants that can be reconstructed by combining left-hand and right-hand segments from two more abudant "parent" sequences.

| <span> |
| :------------------------------------------------------------------------------------------------------------ |
| **NOTE:** The chimera detection method has to be `pooled` instead of `consensus`. |
| <span> |

```bash
#!/bin/bash
#SBATCH --job-name=dada2_05_chimeras.qsub                     # job name
#SBATCH --mail-type=BEGIN,END                                 # send a mail at the end of job
#SBATCH --mail-user=balric@sb-roscoff.fr                      # where to send mail
#SBATCH --cpus-per-task 4                                     # number of CPUs required per task
#SBATCH --mem 100GB                                           # memory per processor
#SBATCH --workdir=/shared/projects/indigene/PHYTOPORT/18SV4   # set working directory
#SBATCH -p fast                                               # partition
#SBATCH -o dada2_05_chimeras.out                              # output file
#SBATCH -e dada2_05_chimeras.err                              # error file

module load r/4.0.3

srun Rscript /shared/projects/indigene/PHYTOPORT/18SV4/script/dada2_05_chimeras.R
```

```r
# ----------------------------------------------------------- #
# This section corresponds to the Rscript dada2_05_chimeras.R #
# ---------------------------------------------- ------------ #

# Some packages must be installed and loaded
library(dada2); packageVersion("dada2")
library(data.table)
library(magrittr)
library(digest)


# Define the following path variable so that it points to the extracted directory
pathFilt <- "/shared/projects/indigene/PHYTOPORT/18SV4/dada2/log"

# Merge multiple runs (if necessary)
x <- grep("seqtab_.+\\.rds$", dir("/shared/projects/indigene/PHYTOPORT/18SV4/dada2/"), value = TRUE)
x <- paste0("/shared/projects/indigene/PHYTOPORT/18SV4/dada2/", x)
st.all <- mergeSequenceTables(tables = x)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, method = "pooled", multithread = TRUE)
seqtab.nochim2 <- seqtab.nochim %>% t %>% data.table

# load statistics
x <- grep("track_.+\\.rds$",dir("/shared/projects/indigene/PHYTOPORT/18SV4/dada2/log/"), value = TRUE)
x <- paste0("/shared/projects/indigene/PHYTOPORT/18SV4/dada2/log/", x)

stattab <- lapply(x, readRDS) %>%
  rbindlist

stattab[, nochim.read := sapply(sample, function(X){
  sum(seqtab.nochim2[, get(X)])
})]

stattab[, nochim.seq := sapply(sample, function(X){
  sum(seqtab.nochim2[, get(X)] != 0)
})]

stattab <- stattab[, list(sample, denoisedF.read, denoisedR.read, merged.read, 
                          nochim.read, denoisedF.seq, denoisedR.seq, merged.seq, nochim.seq)]

filter <- read.csv2(paste0(pathFilt, "/filter.csv"), header = TRUE, stringsAsFactors = FALSE)
filter$sample <- sub("_trimmed.+$","", filter$sample)
stattab <- merge(stattab, filter, by.x = "sample", by.y = "sample")
stattab <- stattab[, c("sample", "reads.in", "reads.out",
                       "denoisedF.read", "denoisedR.read", "merged.read", "nochim.read",
                       "denoisedF.seq", "denoisedR.seq", "merged.seq", "nochim.seq")]
fwrite(stattab, "/shared/projects/indigene/PHYTOPORT/18SV4/dada2/log/statistics.tsv", sep = "\t")

rm(stattab, seqtab.nochim2)

# Construct sequence table
tmp <- seqtab.nochim %>% t %>% data.table(keep.rownames = TRUE)
setnames(tmp, "rn", "sequence")
rm(seqtab.nochim)
tmp <- data.table(amplicon = sapply(tmp[, sequence], digest, algo = "sha1"), tmp)
fwrite(tmp, "/shared/projects/indigene/PHYTOPORT/18SV4/dada2/seqtab_all.tsv", sep = "\t")

# Version of packages used to build this document
sessionInfo()
```
At this stage, you have on ASV table which is the final product of the `dada2` pipeline. In this ASV table each row corresponds to a processes sample, and each column corresponds to an non-chimeric sample sequence (a more precise analogue to the common "OTU table"). But, the sequences variants are not yet annoted, i.e. assigning a taxonomy to the sequence variants. The `dada2` package provides a native implementation of the [naive Bayesian classifer method](https://aem.asm.org/content/73/16/5261.long) for this purpose.