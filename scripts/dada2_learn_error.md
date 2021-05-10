## **Step 3: Learn error rates**
Errors can be introduced by PCR amplification and sequencing. In this part of the pipeline `dada2` will learn to distinguish error from biological differences using a subset of our data as a training set. The `learnErrors()` method learns this error model from the data, by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution. he error parameters typically vary between sequencing runs and PCR protocols, so this method provides a way to estimate those parameters from the data itself.

| <span> |
| :------------------------------------------------------------------------------------------------------------ |
| **NOTE:** We want to make sure that the machine learning algorithm is learning the erroe rates properly. In plot of learning error rates, the red line represents what we should expect the learned error rates to look like for each of the 16 possible base transitions (A->A, A->C, A->G, etc.). Grey points are the obeserved error rates for each consensus quality score and black lines shown the error rates expected under the nominal definiton of the Q-score. The expected error rates should have a good fit with the observed rates, and the error rates must drop with increased quality. If black lines and red lines are very far off from each other, it may be a good idea to indrease the `nbases` parameter. This allows the machine leraning algorithm to train on a larger portion of your data and may help improve the fit. |
| <span> |

```bash
#!/bin/bash
#SBATCH --job-name=dada2_03_learn_error.qsub                  # job name
#SBATCH --mail-type=BEGIN,END                                 # send a mail at the begining/end of job
#SBATCH --mail-user=balric@sb-roscoff.fr                      # where to send mail
#SBATCH --cpus-per-task 6                                     # number of CPUs required per task
#SBATCH --mem 10GB                                            # memory per processor
#SBATCH --workdir=/shared/projects/indigene/PHYTOPORT/18SV4   # set working directory
#SBATCH -p fast                                               # partition
#SBATCH -o dada2_03_learn_error.out                           # output file
#SBATCH -e dada2_03_learn_error.err                           # error file

module load r/4.0.3

srun Rscript /shared/projects/indigene/PHYTOPORT/18SV4/script/dada2_03_learn_error.R
```

```r
# -------------------------------------------------------------- #
# This section corresponds to the Rscript dada2_03_learn_error.R #
# -------------------------------------------------------------- #


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
runs <- sub("^PHP_[B, F, H]\\d+_\\d+_\\d+_([^_]+).+$", "\\1", tmp$sample) %>% unique
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

# Set seed to ensure that randomized steps are replicatable
set.seed(100)

# Loop allowing to determinate learn error rates and infer sample composition
foreach(i = runs, .packages = c("data.table","dada2")) %dopar% {
  filtFs <- filtFsall[grep(i, names(filtFsall))]
  filtRs <- filtRsall[grep(i, names(filtRsall))]
  sample.names <- sub("_trimmed.+$","", basename(filtFs))
  # Learn forward error rates
  errF <- learnErrors(filtFs, nbases = 1e8, multithread = FALSE)
  pdf(paste0("dada2/log/errF_", i, ".pdf"))
  print(plotErrors(errF, nominalQ = TRUE))
  dev.off()
  # Learn reverse error rates
  errR <- learnErrors(filtRs, nbases = 1e8, multithread = FALSE)
  pdf(paste0("dada2/log/errR_", i, ".pdf"))
  print(plotErrors(errR, nominalQ = TRUE))
  dev.off()
  save(errF, file = paste0(path, "/dada2/log/errF_",i,".rda")) # CHANGE ME to where you want sequence table saved
  save(errR, file = paste0(path, "/dada2/log/errR_",i,".rda")) # CHANGE ME to where you want sequence table saved
}

# Version of packages used to build this document
sessionInfo()
```

**Visualization example of the estimates error rates**

![image](https://github.com/benalric/Metab_pipeline_v2/blob/main/figures/errF_PHYTOPORT_Cut1.png)

| <span> |
| :--- |
| **NOTE:** We want to make sure that the machine learning algorithm is learning the erroe rates properly. In plot of learning error rates, the red line represents what we should expect the learned error rates to look like for each of the 16 possible base transitions (A->A, A->C, A->G, etc.). Grey points are the obeserved error rates for each consensus quality score and black lines shown the error rates expected under the nominal definiton of the Q-score. The expected error rates should have a good fit with the observed rates, and the error rates must drop with increased quality. If black lines and red lines are very far off from each other, it may be a good idea to indrease the `nbases` parameter. This allows the machine leraning algorithm to train on a larger portion of your data and may help improve the fit. |
| <span> |

