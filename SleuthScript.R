source("http://bioconductor.org/biocLite.R")
biocLite("remotes")
biocLite("devtools")
biocLite("pachterlab/sleuth")

library('sleuth')

help(package = 'sleuth')

base_dir <- "~/OneDrive - University of Florida/MastersThesisPub/Data/"
sample_id <- dir(file.path(base_dir,"KallistoQuantRun"))
sample_id

kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "KallistoQuantRun", id, "kallisto"))
kal_dirs

s2c <- read.table(file.path(base_dir,"hiseq_info.txt"), header=TRUE, stringsAsFactors = FALSE)
s2c <- dplyr::select(s2c, sample = sample, condition)
s2c

s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)

so <- sleuth_prep(s2c, ~ condition)
so <- sleuth_fit(so)
so <- sleuth_wt(so, 'conditionexposed')
models(so)

sleuth_live(so)
results_table <- sleuth_results(so, 'conditionexposed')
results_table