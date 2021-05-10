### **Step 2_BIS: Filtering and trimming**
| <span> |
| :------------------------------------------------------------------------------------------------------------ |
| **NOTE:** The script dada2_02_filter_trim_a.R produces for each sample one file .csv. Here we merge all files in one table with rows as samples and three columns  correspoding to sample names, reads.in, reads.out. |
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