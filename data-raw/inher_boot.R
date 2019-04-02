#assign workspace
setwd('/home/crumplecup/work')

library(data.table)
library(magrittr)
library(muddier)
library(parallel)
options(mc.cores=detectCores())

# rank list of sites class 'DF' obs > 3
rl <- rank_list('DF')

rbit <- rl[[1]]
bit <- convo_rank(rbit)
also <- convo_lis(rbit, char_pmfs, years)


trunc_sum <- function(mat)  {
  index <- as.numeric(rownames(mat))
  vec <- rowSums(mat)
  trunc_0(vec, index)
}

# boot list of ranks with replacement
bl <- boot_ids(12, rl)
pmfl <- pmf_by_rank(rl)

# convolved pmfs of ranks in boot list
cbl <- mclapply(bl, convo_lis)

# convolved pmfs from list to data.table
cdt <- lapply(cbl, rack) %>% rack
colnames(cdt) <- lapply(bl, function(a) lapply(a, function(b) id_by_rank(b))) %>% unlist

# cis of convolved pmfs
xp <- cis_of_boot(cdt)



