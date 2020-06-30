#assign workspace
setwd('/home/crumplecup/work')

library(data.table)
library(magrittr)
library(muddier)
library(parallel)
options(mc.cores=detectCores())
set.seed(10101)

# summary data frame
df <- data.frame(age = sort(index), df = df_pmf, ff = ff_pmf, fg = fg_pmf)
write.csv(df, 'inherage.csv')

# debris flow inherited age
index <- char_pmfs %>% rownames %>% as.numeric
df_ranks <- rank_list('DF')
df_cl <- convo_list(df_ranks)
df_trc <- lapply(df_cl, function(a) trunc_list(a, sort(index)))
df_trc <- lapply(df_trc, rack)
df_pmf <- df_trc[[1]]

for (i in 2:7) {
  df_pmf <- cbind(df_pmf, df_trc[[i]])
}

df_pmf <- df_pmf %>% rowSums %>% normalize
df_cdf <- df_pmf %>% to_cdf
df_exc <- df_pmf %>% to_exceed

setwd('/home/crumplecup/work/muddier')
usethis::use_data(df_pmf, overwrite = T)
usethis::use_data(df_cdf, overwrite = T)
usethis::use_data(df_exc, overwrite = T)

# fluvial fines inherited age
index <- char_pmfs %>% rownames %>% as.numeric %>% sort
ff_ranks <- rank_list('FF')
ff_cl <- convo_list(ff_ranks)
ff_trc <- lapply(ff_cl, function(a) trunc_list(a, index))
ff_pmf <- to_mat(ff_trc) %>% rowSums %>% normalize
ff_cdf <- ff_pmf %>% to_cdf
ff_exc <- ff_pmf %>% to_exceed

setwd('/home/crumplecup/work/muddier')
usethis::use_data(ff_pmf, overwrite = T)
usethis::use_data(ff_cdf, overwrite = T)
usethis::use_data(ff_exc, overwrite = T)


# fluvial gravels inherited age
index <- char_pmfs %>% rownames %>% as.numeric %>% sort
fg_ranks <- rank_list('FG')
fg_cl <- convo_list(fg_ranks)
fg_trc <- lapply(fg_cl, function(a) trunc_list(a, index))
fg_pmf <- to_mat(fg_trc) %>% rowSums %>% normalize
fg_cdf <- fg_pmf %>% to_cdf
fg_exc <- fg_pmf %>% to_exceed

setwd('/home/crumplecup/work/muddier')
usethis::use_data(fg_pmf, overwrite = T)
usethis::use_data(fg_cdf, overwrite = T)
usethis::use_data(fg_exc, overwrite = T)

plot(ff_cdf)

# bootstrap cis on inherited age
df_bl <- boot_ids(10000, df_ranks)
ff_bl <- boot_ids(10000, ff_ranks)
fg_bl <- boot_ids(10000, fg_ranks)

usethis::use_data(df_bl)
usethis::use_data(ff_bl)
usethis::use_data(fg_bl)


# fluvial fines bootstrap

begin_og <- Sys.time()

index <- sort(as.numeric(rownames(char_pmfs)))
ff_c1 <- convo_list(ff_bl[[1]])
save(ff_c1, file = 'ff_c1.rda', overwrite = T)
ff_t1 <- lapply(ff_c1, function(a) trunc_list(a, index))
save(ff_t1, file = 'ff_t1.rda', overwrite = T)
rm(ff_c1)
ff_ta <- lapply(ff_t1, rack) %>% rack
save(ff_ta, file = 'ff_ta.rda', overwrite = T)
rm(ff_t1, ff_ta)

ff_c2 <- convo_list(ff_bl[[2]])
save(ff_c2, file = 'ff_c2.rda', overwrite = T)
ff_t2 <- lapply(ff_c2, function(a) trunc_list(a, index))
save(ff_t2, file = 'ff_t2.rda', overwrite = T)
rm(ff_c2)
ff_tb <- lapply(ff_t2, rack) %>% rack
save(ff_tb, file = 'ff_tb.rda', overwrite = T)
rm(ff_t2, ff_tb)

ff_c3 <- convo_list(ff_bl[[3]])
save(ff_c3, file = 'ff_c3.rda', overwrite = T)
ff_t3 <- lapply(ff_c3, function(a) trunc_list(a, index))
save(ff_t3, file = 'ff_t3.rda', overwrite = T)
rm(ff_c3)
ff_tc <- lapply(ff_t3, rack) %>% rack
save(ff_tc, file = 'ff_tc.rda', overwrite = T)
rm(ff_t3, ff_tc)

end_og <- Sys.time()

end_og - begin_og  # Time difference of 3.124383 hours

load('ff_t1.rda')

ff_t1 <- lapply(ff_t1, rackb)
ff_b1 <- ff_t1
save(ff_b1, file = 'ff_b1.rda')
ff_b1 <- ff_b1 %>% rack

load('ff_t2.rda')
ff_t2 <- lapply(ff_t2, rackb)
ff_b2 <- ff_t2
save(ff_b2, file = 'ff_b2.rda')
ff_b2 <- ff_b2 %>% rack

load('ff_t3.rda')
ff_t3 <- lapply(ff_t3, rackb)
ff_b3 <- ff_t3
save(ff_b3, file = 'ff_b3.rda')
ff_b3 <- ff_b3 %>% rack

ff_b <- cbind(ff_b1, cbind(ff_b2, ff_b3))
save(ff_b, file = 'ff_b.rda')

ff_b0 <- matrix(0, ncol = ncol(ff_b1), nrow = nrow(ff_b1))
for (i in 1:nrow(ff_b0))  {
  ff_b0[i, ] <- (ff_b1[i, ] + ff_b2[i, ] + ff_b3[i, ]) / 3
}

ff_pmfs_cis <- matrix(0, ncol = 3, nrow = nrow(ff_b))
for (i in 1:nrow(ff_b))  {
  print(i)
  ff_pmfs_cis[i, ] <- ff_b[i, ] %>% get_cis
}

ff_pmfs <- ff_b0
ff_cdfs <- apply(ff_pmfs, 2, to_cdf)
ff_cdfs_cis <- apply(ff_cdfs, 1, get_cis)
ff_excs <- 1 - ff_cdfs
ff_excs_cis <- apply(ff_excs, 1, get_cis)

# fluvial gravels bootstrap

begin <- Sys.time()

fg_c1 <- convo_list(fg_bl[[1]])
save(fg_c1, file = 'fg_c1.rda')
fg_t1 <- lapply(fg_c1, function(a) trunc_list(a, index))
save(fg_t1, file = 'fg_t1.rda', overwrite = T)
rm(fg_c1)
fg_ta <- lapply(fg_t1, rack) %>% rack
save(fg_ta, file = 'fg_ta.rda', overwrite = T)
rm(fg_t1, fg_ta)

fg_c2 <- convo_list(fg_bl[[2]])
save(fg_c2, file = 'fg_c2.rda')
fg_t2 <- lapply(fg_c2, function(a) trunc_list(a, index))
save(fg_t2, file = 'fg_t2.rda', overwrite = T)
rm(fg_c2)
fg_tb <- lapply(fg_t2, rack) %>% rack
save(fg_tb, file = 'fg_tb.rda', overwrite = T)
rm(fg_t2, fg_tb)

fg_c3 <- convo_list(fg_bl[[3]])
save(fg_c3, file = 'fg_c3.rda')
fg_t3 <- lapply(fg_c3, function(a) trunc_list(a, index))
save(fg_t3, file = 'fg_t3.rda', overwrite = T)
rm(fg_c3)
fg_tc <- lapply(fg_t3, rack) %>% rack
save(fg_tc, file = 'fg_tc.rda', overwrite = T)
rm(fg_t3, fg_tc)

fg_c4 <- convo_list(fg_bl[[4]])
save(fg_c4, file = 'fg_c4.rda')
fg_t4 <- lapply(fg_c4, function(a) trunc_list(a, index))
save(fg_t4, file = 'fg_t4.rda', overwrite = T)
rm(fg_c4)
fg_td <- lapply(fg_t4, rack) %>% rack
save(fg_td, file = 'fg_td.rda', overwrite = T)
rm(fg_t4, fg_td)

end <- Sys.time()

end - begin
end - begin_og


# assemble fluvial fines bootstraps into single matrix

load('ff_ta.rda')
load('ff_tb.rda')
load('ff_tc.rda')

ff_tab <- cbind(ff_ta, ff_tb)
save(ff_tab, file = 'ff_tab.rda')

load('ff_tab.rda')

ff_tab <- cbind(ff_tab, ff_tc)
save(ff_tab, file = 'ff_tabc.rda')

# transpose matrix
load('ff_tabc.rda')
ff_tab <- ff_tab %>% t
ff_t <- ff_tab
save(ff_t, file = 'ff_t.rda')

# convert from matrix to data.table
ff_t <- as.data.frame(ff_t)
ff_dt <- ff_t
save(ff_dt, file = 'ff_dt.rda')

# extract cis
ff_cis <- lapply(ff_dt, get_cis)
ff_cis <- ff_cis %>% rack
save(ff_cis, file = 'ff_cis.rda')

# convert to cdf
load('ff_tabc.rda')
ff_tab <- as.data.table(ff_tab)
ff_tab <- lapply(ff_tab, to_cdf)
ff_cdfs <- ff_tab
save(ff_cdfs, file = 'ff_cdfs.rda')

# convert cdfs from dt to matrix and transpose
load('ff_cdfs.rda')
ff_cdfs <- as.data.frame(ff_cdfs)

ff_cdfs_cis <- array(0, c(nrow(ff_cdfs), 3))
for (i in 1:nrow(ff_cdfs))  {
  print(i)
  ff_cdfs_cis[i, ] <- ff_cdfs[i, ] %>% unlist %>% get_cis
}

save(ff_cdfs_cis, file = 'ff_cdfs_cis.rda')


ff_cdfs <- as.matrix(ff_cdfs)
ff_cdfs <- 1 - ff_cdfs
ff_excs <- ff_cdfs

save(ff_excs, file = 'ff_excs')
load('ff_excs')

ff_excs_cis <- array(0, c(nrow(ff_excs), 3))
for (i in 1:nrow(ff_excs))  {
  print(i)
  ff_excs_cis[i, ] <- ff_excs[i, ] %>% get_cis
}

save(ff_excs_cis, file = 'ff_excs_cis.rda')

setwd('/home/crumplecup/work/muddier')
usethis::use_data(ff_pmfs_cis)
usethis::use_data(ff_cdfs_cis)
usethis::use_data(ff_excs_cis)

setwd('/home/crumplecup/work')
ff_b <- ff_b0
save(ff_b, file ='ff_b.rda')

# assemble fluvial gravel bootstraps into single matrix

load('fg_ta.rda')
load('fg_tb.rda')

fg_tab <- cbind(fg_ta, fg_tb)
save(fg_tab, file = 'fg_tab.rda')

load('fg_tab.rda')
load('fg_ta.rda')
load('fg_tb.rda')
load('fg_tc.rda')

fg_tab <- cbind(fg_tab, fg_tc)
fg_tabc <- fg_tab
save(fg_tabc, file = 'fg_tabc.rda')

load('fg_tabc.rda')
load('fg_td.rda')

fg_tabc <- cbind(fg_tabc, fg_td)
fg_pmfs <- fg_tabc
save(fg_pmfs, file = 'fg_pmfs.rda')

load('fg_pmfs.rda')
write.csv(fg_pmfs, 'fg_pmfs.csv')

fg_pmfs <- apply(fg_pmfs, 2, to_cdf)

fg_pmfs_cis <- array(0, c(nrow(fg_pmfs), 3))
for (i in 1:nrow(fg_pmfs))  {
  print(i)
  fg_pmfs_cis[i, ] <- fg_pmfs[i, ] %>% get_cis
}

save(fg_pmfs_cis, file = 'fg_pmfs_cis.rda')


load('fg_t1.rda')
load('fg_t2.rda')
load('fg_t3.rda')
load('fg_t4.rda')

fg_t1 <- lapply(fg_t1, rackb) %>% rack
fg_t2 <- lapply(fg_t2, rackb) %>% rack
fg_t3 <- lapply(fg_t3, rackb) %>% rack
fg_t4 <- lapply(fg_t4, rackb) %>% rack

fg_b <- matrix(0, ncol = 10000, nrow = nrow(fg_t1))

for (i in 1:nrow(fg_t1))  {
  fg_b[i,] <- (fg_t1[i,] + fg_t2[i,] + fg_t3[i,] + fg_t4[i,]) / 4
}

fg_pmfs_cis <- apply(fg_b, 1, get_cis)
fg_cdfs <- apply(fg_b, 2, to_cdf)
fg_cdfs_cis <- apply(fg_cdfs, 1, get_cis)
fg_excs <- 1 - fg_cdfs
fg_excs_cis <- apply(fg_excs, 1, get_cis)

setwd('/home/crumplecup/work/muddier')
usethis::use_data(fg_pmfs_cis)
usethis::use_data(fg_cdfs_cis)
usethis::use_data(fg_excs_cis)


# convert fluvial gravel pmfs to cdfs

load('fg_ta.rda')
fg_cdfs_a <- apply(fg_ta, 2, to_cdf)
save(fg_cdfs_a, file = 'fg_cdfs_a.rda')

load('fg_tb.rda')
fg_cdfs_b <- apply(fg_tb, 2, to_cdf)
save(fg_cdfs_b, file = 'fg_cdfs_b.rda')

load('fg_tc.rda')
fg_cdfs_c <- apply(fg_tc, 2, to_cdf)
save(fg_cdfs_c, file = 'fg_cdfs_c.rda')

load('fg_td.rda')
fg_cdfs_d <- apply(fg_td, 2, to_cdf)
save(fg_cdfs_d, file = 'fg_cdfs_d.rda')

# merge fluvial gravel cdfs into one matrix

load('fg_cdfs_a.rda')
load('fg_cdfs_b.rda')
fg_cdfs_a <- cbind(fg_cdfs_a, fg_cdfs_b)
fg_cdfs_ab <- fg_cdfs_a
save(fg_cdfs_ab, file = 'fg_cdfs_ab.rda')

load('fg_cdfs_c.rda')
load('fg_cdfs_d.rda')
fg_cdfs_c <- cbind(fg_cdfs_c, fg_cdfs_d)
fg_cdfs_cd <- fg_cdfs_c
save(fg_cdfs_cd, file = 'fg_cdfs_cd.rda')

load('fg_cdfs_ab.rda')
load('fg_cdfs_cd.rda')
fg_cdfs_ab <- cbind(fg_cdfs_ab, fg_cdfs_cd)
fg_cdfs <- fg_cdfs_ab
save(fg_cdfs, file = 'fg_cdfs.rda')


load('fg_cdfs.rda')
write.csv(fg_cdfs, 'fg_cdfs.csv')


fg_cdfs_cis <- array(0, c(nrow(fg_cdfs), 3))
for (i in 1:nrow(fg_cdfs))  {
  print(i)
  fg_cdfs_cis[i, ] <- fg_cdfs[i, ] %>% get_cis
}

save(fg_cdfs_cis, file = 'fg_cdfs_cis.rda')







data(fg_cdf, package = 'muddier')
load('fg_cdfs_cis.rda')

index <- sort(as.numeric(rownames(char_pmfs)))
setwd('/home/crumplecup/work')
png('fg_cdfs.png', height = 6, width = 6, units = 'in', res = 300)
plot(index, fg_cdf, xlim = c(0,10000), type = 'l', col = 'slateblue', lwd = 2,
     main = 'fluvial gravels inherited age cdf',
     xlab = 'inherited age', ylab = 'proportion age <= x')
lines(index, fg_cdfs_cis[1,], lty = 2)
lines(index, fg_cdfs_cis[3,], lty = 2)
lines(index, fg_cdfs_cis[2,], lty = 2, col = 'magenta3')
legend('bottomright', legend = c('observed', 'median', '95% CIs'),
       fill = c('slateblue', 'magenta3', 'black'))
dev.off()

png('fg_excs.png', height = 6, width = 6, units = 'in', res = 300)
plot(index, fg_exc, type='l', col='slateblue', lwd=2, log='xy', ylim=c(10^-2,10^0),
     main = 'fluvial gravels inherited age exceedance',
     xlab = 'inherited age', ylab = 'exceedance probability')
lines(index, fg_excs_cis[1,], lty = 2)
lines(index, fg_excs_cis[3,], lty = 2)
lines(index, fg_excs_cis[2,], lty = 2, col = 'magenta3')
legend('topright', legend = c('observed', 'median', '95% CIs'),
       fill = c('slateblue', 'magenta3', 'black'))
dev.off()


data(fg_pmf, package = 'muddier')
load('fg_pmfs_cis.rda')
plot(index, fg_pmfs_cis[,3], type = 'l', xlim = c(0,10000), lty = 2)
lines(index, fg_pmf)
lines(index, fg_pmfs_cis[,1], lty = 2)
lines(index, fg_pmfs_cis[,2], lty = 3, col = 'slateblue')

data(ff_pmf, package = 'muddier')
load('ff_cis.rda')
plot(index, ff_cis[3,], type = 'l', xlim = c(0,10000), lty = 2)
lines(index, ff_pmf)
lines(index, ff_cis[1,], lty = 2)
lines(index, ff_cis[2,], lty = 3, col = 'slateblue')

data(ff_cdf, package = 'muddier')
data('ff_cdfs_cis', package = 'muddier')
pal <- get_palette(c('coral', 'ocean', 'charcoal'), .7)
png('ff_cdfs_cis.png', height = 6, width = 6, units = 'in', res = 300)
plot(index, ff_cdfs_cis[3,], xlim = c(0,3500), lty=2, type = 'l', lwd = 2,
  main = 'fluvial fines inherited age cdf',
  xlab = 'inherited age', ylab = 'probability age <= x', col = pal[3])
lines(index, ff_cdfs_cis[1,], lty = 2, lwd = 2, col = pal[3])
lines(index, ff_cdf, lwd = 2, col = pal[1])
lines(index, ff_cdfs_cis[2,], lty = 2, col = pal[2], lwd = 2)
legend('bottomright', legend = c('observed', 'median', '95% CIs'),
       fill = pal)
dev.off()

data(ff_exc, package = 'muddier')
load('ff_excs_cis.rda')
png('ff_exc_cis.png', height = 4, width = 6, units = 'in', res = 300)
plot(index, ff_excs_cis[3,], log='xy', type='l', lty=2, ylim = c(10^-2,10^0),
     main = 'fluvial fines inherited age exceedance',
     xlab = 'inherited age', ylab = 'exceedance probability')
lines(index, ff_excs_cis[1,], lty = 2)
lines(index, ff_exc)
lines(index, ff_excs_cis[2,], lty = 2, col = 'slateblue')
dev.off()




# bootstrapping by facies but not site

df_ranks <- rank_list('DF')
dfr <- df_ranks %>% unlist
dfb <- boot_ids(100, list(dfr))
dfc <- mclapply(dfb, convo_list)
#save(dfc, file = 'dfc.rda')
index <- sort(as.numeric(rownames(char_pmfs)))
dfm <- lapply(dfc, rack)
dft <- lapply(dfm, function(a) trunc_list(a, index))
dfn <- dft[[1]] %>% rack
dfo <- apply(dfn, 2, function(a) trunc_norm(a, index))
dfp <- apply(dfo, 2, to_cdf)
dfq <- apply(dfp, 1, get_cis)

plot(index, dfq[3,])
lines(index, dfq[2,], lwd = 2, col = 'slateblue')
lines(index, dfq[1,], lwd = 2, lty = 2)



#  testing truncate and normalize functions

trunc_0 <- function(pmf, index)  {
  trunc <- pmf
  for (i in seq_along(index))  {
    if (index[i] < 0)  trunc[i] <- 0
  }
  trunc
}


trunc_1 <- function(pmf, index)  {
  trunc <- pmf
  trunc[index < 0] <- 0
  trunc
}


normalize <- function(probs)  {
  norm <- array(0, length(probs))
  for (i in seq_along(probs))  {
    norm[i] <- probs[i] / sum(probs)
  }
  norm
}


normalize1 <- function(probs)  {
  tot <- sum(probs)
  lapply(probs, function(a) a/tot) %>% unlist
}


normalize2 <- function(probs)  {
  tot <- sum(probs)
  for (i in seq_along(probs))  {
    probs[i] <- probs[i] / tot
  }
  probs
}




mat <- dfm[[1]] %>% rack
tet <- mat[1,] > 0

sub <- mat[,tet]

bit <- sub[,7] %>% trunc_1(index)
bb <- bit %>% normalize3


microbenchmark::microbenchmark(
  trunc0 = sub[,7] %>% trunc_0(index),
  trunc1 = sub[,7] %>% trunc_1(index)

)


microbenchmark::microbenchmark(
  norm = sub[,7] %>% trunc_1(index) %>% normalize,
  norm1 = sub[,7] %>% trunc_1(index) %>% normalize1,
  norm2 = sub[,7] %>% trunc_1(index) %>% normalize2
)


microbenchmark::microbenchmark(
  one = sub[,7] %>% trunc_norm(index)
)










read_row <- function(id, dt)  {
  nm <- paste0('row_',id,'.rda')
  vec <- dt[id, ]
  save(vec, file = nm)
}

read_rows <- function(dt, out_dir)  {
  og_dir <- getwd()
  setwd(out_dir)
  for (i in 1:nrow(dt))  {
    read_row(i, dt)
  }
  setwd(og_dir)
}



out_path <- '/home/crumplecup/work/ff_cdfs_rows/'
read_rows(ff_cdfs, out_path)

setwd(out_path)
files <- dir()
ff_cdfs_cis <- array(0, c(length(files), 3))


begin <- Sys.time()
for (i in 1:nrow(ff_cdfs_cis))  {
  load(files[i])
  ff_cdfs_cis[i,] <- get_cis(unlist(vec))
  rm(vec)
}
end <- Sys.time()
end - begin

plot(ff_cdfs_cis[,1])



ff_cdfs_cis <- array(0, c(nrow(ff_cdfs), 3))
for (i in 1:nrow(ff_cdfs))  {
  ff_cdfs_cis[i,] <- get_cis(ff_cdfs[i,])
}



ff_cdfs <- ff_cdfs %>% rack
ff_cdfs <- ff_cdfs %>% t
ff_cdfs_mat <- ff_cdfs
save(ff_cdfs_mat, file = 'ff_cdfs_mat.rda')

# convert cdfs from matrix to data.table
ff_cdfs_mat <- as.data.frame(ff_cdfs_mat)
ff_cdfs_dt <- ff_cdfs_mat
save(ff_cdfs_dt, file = 'ff_cdfs_dt.rda')

# extract cis for cdfs
ff_cdfs_cis <- lapply(ff_cdfs_dt, get_cis)
ff_cdfs_cis <- ff_cdfs_cis %>% rack
save(ff_cdfs_cis, file = 'ff_cdfs_cis.rda')


plot(ff_cis[3,])





data(df_bl, package = 'muddier')
sub <- list(df_bl[[1]], df_bl[[2]])


convolve_list1 <- function(ranks, pmfs = char_pmfs)  {
  index <- as.numeric(rownames(char_pmfs))
  rank_pmfs <- pmf_by_rank(ranks, pmfs)
  lapply(rank_pmfs, function(a)
    setDT(lapply(a, function(b, c) convo(unlist(b), unlist(c), index), c = a[,1])))
}

convolve_list2 <- function(ranks, pmfs = char_pmfs)  {
  index <- as.numeric(rownames(char_pmfs))
  rank_pmfs <- pmf_by_rank(ranks, pmfs)
  lapply(rank_pmfs, function(a)
    setDT(mclapply(a, function(b, c) convo(unlist(b), unlist(c), index), c = a[,1])))
}

convolve_list3 <- function(ranks, pmfs = char_pmfs)  {
  index <- as.numeric(rownames(char_pmfs))
  rank_pmfs <- pmf_by_rank(ranks, pmfs)
  mclapply(rank_pmfs, function(a)
    setDT(lapply(a, function(b, c) convo(unlist(b), unlist(c), index), c = a[,1])))
}

microbenchmark::microbenchmark(
  con1 = convolve_list1(df_bl[[1]]),
  con2 = convolve_list2(df_bl[[1]]),
  con3 = convolve_list2(df_bl[[1]]), times = 10L
)

# Unit: seconds
# expr      min       lq     mean   median       uq      max neval
# con1 12.37846 12.64029 12.69617 12.71025 12.80028 12.94136    10
# con2 20.27604 22.14754 22.44325 22.89713 23.09889 23.84933    10
# con3 20.74371 21.61118 22.59161 23.00127 23.44725 24.32304    10

boot_dir <- '/home/crumplecup/work/boot/'
setwd(boot_dir)

begin <- Sys.time()
df_c1 <- convo_list(df_bl[[1]])
end <- Sys.time()
end - begin  # Time difference of 43.02693 mins
save(df_c1, file = 'df_c1.rda')
rm(df_c1)
gc()

begin <- Sys.time()
df_c2 <- convo_list(df_bl[[2]])
end <- Sys.time()
end - begin  # Time difference of 49.74569 mins
save(df_c2, file = 'df_c2.rda')
rm(df_c2)
gc()

begin <- Sys.time()
df_c3 <- convo_list(df_bl[[3]])
end <- Sys.time()
end - begin  # Time difference of 37.76244 mins
save(df_c3, file = 'df_c3.rda')
rm(df_c3)
gc()

begin <- Sys.time()
df_c4 <- convo_list(df_bl[[4]])
end <- Sys.time()
end - begin  # Time difference of 44.56297 mins
save(df_c4, file = 'df_c4.rda')
rm(df_c4)
gc()

begin <- Sys.time()
df_c5 <- convo_list(df_bl[[5]])
end <- Sys.time()
end - begin  # Time difference of 34.94688 mins
save(df_c5, file = 'df_c5.rda')
rm(df_c5)
gc()

begin <- Sys.time()
df_c6 <- convo_list(df_bl[[6]])
end <- Sys.time()
end - begin  # Time difference of 44.63116 mins
save(df_c6, file = 'df_c6.rda')
rm(df_c6)
gc()

begin <- Sys.time()
df_c7 <- convo_list(df_bl[[7]])
end <- Sys.time()
end - begin  # Time difference of 43.02693 mins
save(df_c7, file = 'df_c7.rda')
rm(df_c7, df_bl)
gc()

load('df_c1.rda')
index <- sort(as.numeric(rownames(char_pmfs)))
begin <- Sys.time()
df_t1 <- lapply(df_c1, function(a) trunc_list(a, index))
end <- Sys.time()
end - begin  # Time difference of 27.79817 mins
save(df_t1, file = 'df_t1.rda')
rm(df_c1, df_t1)
gc()

load('df_c2.rda')
begin <- Sys.time()
df_t2 <- lapply(df_c2, function(a) trunc_list(a, index))
end <- Sys.time()
end - begin  # Time difference of 29.99942 mins
save(df_t2, file = 'df_t2.rda')
rm(df_c2, df_t2)
gc()

load('df_c3.rda')
begin <- Sys.time()
df_t3 <- lapply(df_c3, function(a) trunc_list(a, index))
end <- Sys.time()
end - begin  # Time difference of 18.27289 mins
save(df_t3, file = 'df_t3.rda')
rm(df_c3, df_t3)
gc()

load('df_c4.rda')
begin <- Sys.time()
df_t4 <- lapply(df_c4, function(a) trunc_list(a, index))
end <- Sys.time()
end - begin  # Time difference of 31.00294 mins
save(df_t4, file = 'df_t4.rda')
rm(df_c4, df_t4)
gc()

load('df_c5.rda')
begin <- Sys.time()
df_t5 <- lapply(df_c5, function(a) trunc_list(a, index))
end <- Sys.time()
end - begin  # Time difference of 27.82506 mins
save(df_t5, file = 'df_t5.rda')
rm(df_c5, df_t5)
gc()

load('df_c6.rda')
begin <- Sys.time()
df_t6 <- lapply(df_c6, function(a) trunc_list(a, index))
end <- Sys.time()
end - begin  # Time difference of 31.02878 mins
save(df_t6, file = 'df_t6.rda')
rm(df_c6, df_t6)
gc()

load('df_c7.rda')
begin <- Sys.time()
df_t7 <- lapply(df_c7, function(a) trunc_list(a, index))
end <- Sys.time()
end - begin  # Time difference of 25.12493 mins
save(df_t7, file = 'df_t7.rda')
rm(df_c7, df_t7)
gc()


setwd('/home/crumplecup/work/boot/')
library(muddier)
load('df_t1.rda')
df_t1 <- rackl(df_t1)
df_t1_cdf <- apply(df_t1, 2, to_cdf)
save(df_t1_cdf, file = 'df_t1_cdf.rda')
rm(df_t1_cdf)
gc()

load('df_t2.rda')
df_t2 <- rackl(df_t2)
df_t2_cdf <- apply(df_t2, 2, to_cdf)
save(df_t2_cdf, file = 'df_t2_cdf.rda')
rm(df_t2_cdf)
gc()

df_t1 <- cbind(df_t1, df_t2)
df_m1 <- df_t1
save(df_m1, file = 'df_m1.rda')
rm(df_m1, df_t1, df_t2)
gc()

setwd('/home/crumplecup/work/boot/')
library(muddier)
load('df_t1_cdf.rda')
load('df_t2_cdf.rda')
df_t1_cdf <- cbind(df_t1_cdf, df_t2_cdf)
df_m1_cdf <- df_t1_cdf
save(df_m1_cdf, file = 'df_m1_cdf.rda')
rm(df_t1_cdf, df_t2_cdf, df_m1_cdf)
gc()


setwd('/home/crumplecup/work/boot')
library(muddier)
load('df_t3.rda')
df_t3 <- rackl(df_t3)
df_t3_cdf <- apply(df_t3, 2, to_cdf)
save(df_t3_cdf, file = 'df_t3_cdf.rda')
rm(df_t3_cdf)
gc()

load('df_t4.rda')
df_t4 <- rackl(df_t4)
df_t4_cdf <- apply(df_t4, 2, to_cdf)
save(df_t4_cdf, file = 'df_t4_cdf.rda')
rm(df_t4_cdf)
gc()

df_t3 <- cbind(df_t3, df_t4)
df_m2 <- df_t3
save(df_m2, file = 'df_m2.rda')
rm(df_t3, df_t4, df_m2)
gc()

setwd('/home/crumplecup/work/boot/')
library(muddier)
load('df_t3_cdf.rda')
load('df_t4_cdf.rda')
df_t3_cdf <- cbind(df_t3_cdf, df_t4_cdf)
df_m2_cdf <- df_t3_cdf
save(df_m2_cdf, file = 'df_m2_cdf.rda')
rm(df_t3_cdf, df_t4_cdf, df_m2_cdf)
gc()

setwd(boot_dir)
library(muddier)
load('df_m1.rda')
load('df_m2.rda')

df_m1 <- cbind(df_m1, df_m2)
save(df_m1, file = 'df_m1.rda')
rm(df_m1, df_m2)
gc()

load('df_m1_cdf.rda')
load('df_m2_cdf.rda')
df_m1_cdf <- cbind(df_m1_cdf, df_m2_cdf)
save(df_m1_cdf, file = 'df_m1_cdf.rda')
rm(df_m1_cdf, df_m2_cdf)
gc()

setwd(boot_dir)
library(muddier)
load('df_t5.rda')
df_t5 <- rackl(df_t5)
df_t5_cdf <- apply(df_t5, 2, to_cdf)
save(df_t5_cdf, file = 'df_t5_cdf.rda')
rm(df_t5_cdf)
gc()

load('df_t6.rda')
df_t6 <- rackl(df_t6)
df_t6_cdf <- apply(df_t6, 2, to_cdf)
save(df_t6_cdf, file = 'df_t6_cdf.rda')
rm(df_t6_cdf)
gc()

df_t5 <- cbind(df_t5, df_t6)
df_m2 <- df_t5
save(df_m2, file = 'df_m2.rda')
rm(df_t5, df_t6, df_m2)

setwd('/home/crumplecup/work/boot/')
library(muddier)
load('df_t5_cdf.rda')
load('df_t6_cdf.rda')
df_t5_cdf <- cbind(df_t5_cdf, df_t6_cdf)
df_m2_cdf <- df_t5_cdf
save(df_m2_cdf, file = 'df_m2_cdf.rda')
rm(df_t5_cdf, df_t6_cdf, df_m2_cdf)
gc()


load('df_t7.rda')
df_t7 <- rackl(df_t7)
df_t7_cdf <- apply(df_t7, 2, to_cdf)
save(df_t7_cdf, file = 'df_t7_cdf.rda')
rm(df_t7_cdf)
gc()

load('df_m2.rda')
df_m2 <- cbind(df_m2, df_t7)
save(df_m2, file = 'df_m2.rda')
rm(df_m2, df_t7)
gc()

setwd('/home/crumplecup/work/boot/')
library(muddier)
load('df_m2_cdf.rda')
load('df_t7_cdf.rda')
df_m2_cdf <- cbind(df_m2_cdf, df_t7_cdf)
save(df_m2_cdf, file = 'df_m2_cdf.rda')
rm(df_m2_cdf, df_t7_cdf)
gc()

setwd('/home/crumplecup/work/boot/')
library(muddier)
load('df_m1.rda')
df_m1a <- df_m1[1:2000, ]
save(df_m1a, file = 'df_m1a.rda')
df_m1a <- df_m1[2001:4430,]
df_m1b <- df_m1a
save(df_m1b, file = 'df_m1b.rda')
rm(df_m1, df_m1a, df_m1b)
gc()

load('df_m1_cdf.rda')
df_m1a_cdf <- df_m1_cdf[1:2000, ]
save(df_m1a_cdf, file = 'df_m1a_cdf.rda')
df_m1a_cdf <- df_m1_cdf[2001:4430,]
df_m1b_cdf <- df_m1a_cdf
save(df_m1b_cdf, file = 'df_m1b_cdf.rda')
rm(df_m1_cdf, df_m1a_cdf, df_m1b_cdf)
gc()


load('df_m2.rda')
df_m2a <- df_m2[1:2000, ]
save(df_m2a, file = 'df_m2a.rda')
df_m2a <- df_m2[2001:4430,]
df_m2b <- df_m2a
save(df_m2b, file = 'df_m2b.rda')
rm(df_m2, df_m2a, df_m2b)
gc()

load('df_m2_cdf.rda')
df_m2a_cdf <- df_m2_cdf[1:2000, ]
save(df_m2a_cdf, file = 'df_m2a_cdf.rda')
df_m2a_cdf <- df_m2_cdf[2001:4430,]
df_m2b_cdf <- df_m2a_cdf
save(df_m2b_cdf, file = 'df_m2b_cdf.rda')
rm(df_m2_cdf, df_m2a_cdf, df_m2b_cdf)
gc()


load('df_m1a.rda')
load('df_m2a.rda')
df_m1a <- cbind(df_m1a, df_m2a)
df_ma <- df_m1a
save(df_ma, file = 'df_ma.rda')

df_ma_cis <- array(0, c(nrow(df_ma), 3))
begin <- Sys.time()
for (i in 1:nrow(df_ma))  {
  df_ma_cis[i,] <- get_cis(df_ma[i,])
}
end <- Sys.time()
end - begin  # Time difference of 2.191562 mins
save(df_ma_cis, file = 'df_ma_cis.rda')
rm(df_m1a, df_m2a, df_ma, df_ma_cis)
gc()


setwd('/home/crumplecup/work/boot/')
library(muddier)
load('df_m1b.rda')
load('df_m2b.rda')
df_m1b <- cbind(df_m1b, df_m2b)
df_mb <- df_m1b
save(df_mb, file = 'df_mb.rda')

df_mb_cis <- array(0, c(nrow(df_mb), 3))
begin <- Sys.time()
for (i in 1:nrow(df_mb))  {
  df_mb_cis[i,] <- get_cis(df_mb[i,])
}
end <- Sys.time()
end - begin  # Time difference of 1.563305 mins
save(df_mb_cis, file = 'df_mb_cis.rda')
rm(df_m1b, df_m2b, df_mb)
gc()

load('df_ma_cis.rda')
df_ma_cis <- rbind(df_ma_cis, df_mb_cis)
df_pmfs_cis <- df_ma_cis
save(df_pmfs_cis, file = 'df_pmfs_cis.rda')
rm(df_ma_cis, df_mb_cis)
gc()

index <- sort(as.numeric(rownames(char_pmfs)))
plot(index, df_pmfs_cis[,3], type = 'l', lwd = 3, col = 'slateblue', log = 'x',
     ylim = c(0, .05))
lines(index, df_pmfs_cis[,1], lwd = 3, lty = 2)
lines(index, df_pmfs_cis[,2], lwd = 3)
lines(index, df_pmf, lwd = 3, col = 'forestgreen')
sum(df_pmfs_cis[,2])

sum(df_pmf)
lines(index,df_pmf)

setwd('/home/crumplecup/work/boot')
library(muddier)
load('df_m1a_cdf.rda')
load('df_m2a_cdf.rda')
df_m1a_cdf <- cbind(df_m1a_cdf, df_m2a_cdf)
df_ma_cdf <- df_m1a_cdf
save(df_ma_cdf, file = 'df_ma_cdf.rda')
rm(df_m1a_cdf, df_m2a_cdf)
gc()

df_ma_cdf_cis <- array(0, c(nrow(df_ma_cdf), 3))
for (i in 1:nrow(df_ma_cdf))  {
  df_ma_cdf_cis[i,] <- get_cis(df_ma_cdf[i,])
}
save(df_ma_cdf_cis, file = 'df_ma_cdf_cis.rda')
rm(df_ma_cdf, df_ma_cdf_cis)
gc()


setwd('/home/crumplecup/work/boot/')
library(muddier)
load('df_m1b_cdf.rda')
load('df_m2b_cdf.rda')
df_m1b_cdf <- cbind(df_m1b_cdf, df_m2b_cdf)
df_mb_cdf <- df_m1b_cdf
save(df_mb_cdf, file = 'df_mb_cdf.rda')
rm(df_m1b_cdf, df_m2b_cdf)
gc()

df_mb_cdf_cis <- array(0, c(nrow(df_mb_cdf), 3))
begin <- Sys.time()
for (i in 1:nrow(df_mb_cdf))  {
  df_mb_cdf_cis[i,] <- get_cis(df_mb_cdf[i,])
}
save(df_mb_cdf_cis, file = 'df_mb_cdf_cis.rda')
end <- Sys.time()
end - begin
rm(df_mb_cdf)
gc()


setwd('/home/crumplecup/work/boot/')
library(muddier)
load('df_ma_cdf_cis.rda')
load('df_mb_cdf_cis.rda')
df_ma_cdf_cis <- rbind(df_ma_cdf_cis, df_mb_cdf_cis)
df_cdfs_cis <- df_ma_cdf_cis
save(df_cdfs_cis, file = 'df_cdfs_cis.rda')

library(magrittr)
index <- sort(as.numeric(rownames(char_pmfs)))
lwr <- df_cdfs_cis[,1]
med <- df_cdfs_cis[,2]
upr <- df_cdfs_cis[,3]
emp <- df_cdf

pal <- get_palette(c('ocean', 'charcoal', 'coral'), .8)
png('df_cdf_cis_logx.png', height = 14, width = 17, unit = 'cm', res = 300)
plot(index, lwr, type = 'l', lwd = 3, lty = 3, col = pal[2], #log = 'x',
     main = 'Debris Flow Inherited Age CDF',
     xlab = 'Age Post Year 2000', ylab = 'Proportion Age <= X')
lines(index, upr, lwd = 3, lty = 3, col = pal[2])
lines(index, emp, lwd = 3, col = pal[1])
lines(index, med, lwd = 3, lty = 3, col = pal[3])
legend('topleft', legend = c('ecdf', 'median', '95% CIs'),
       fill = pal[c(1,3,2)])
dev.off()


##
setwd('/home/crumplecup/work/boot/')
library(magrittr)
load('df_t1.rda')
df_t1 <- lapply(df_t1, rackb) %>% rack
load('df_t2.rda')
df_t2 <- lapply(df_t2, rackb) %>% rack
load('df_t3.rda')
df_t3 <- lapply(df_t3, rackb) %>% rack
load('df_t4.rda')
df_t4 <- lapply(df_t4, rackb) %>% rack
load('df_t5.rda')
df_t5 <- lapply(df_t5, rackb) %>% rack
load('df_t6.rda')
df_t6 <- lapply(df_t6, rackb) %>% rack
load('df_t7.rda')
df_t7 <- lapply(df_t7, rackb) %>% rack

df_pmfs <- matrix(0, nrow = nrow(df_t1), ncol = ncol(df_t1))

for (i in 1:nrow(df_pmfs))  {
  df_pmfs[i,] <- (df_t1[i,] + df_t2[i,] + df_t3[i,] + df_t4[i,] +
                    df_t5[i,] + df_t6[i,] + df_t7[i,]) / 7
}

df_pmfs_cis <- apply(df_pmfs, 1, get_cis)
df_cdfs <- apply(df_pmfs, 2, to_cdf)
df_cdfs_cis <- apply(df_cdfs, 1, get_cis)
df_excs <- 1 - df_cdfs
df_excs_cis <- apply(df_excs, 1, get_cis)

setwd('/home/crumplecup/work/muddier/')
usethis::use_data(df_pmfs, overwrite = T)
usethis::use_data(df_pmfs_cis, overwrite = T)
usethis::use_data(df_cdfs, overwrite = T)
usethis::use_data(df_cdfs_cis, overwrite = T)
usethis::use_data(df_excs, overwrite = T)
usethis::use_data(df_excs_cis, overwrite = T)

#plot results

data(df_pmf, package = 'muddier')
data(ff_pmf, package = 'muddier')
data(fg_pmf, package = 'muddier')

data(df_cdf, package = 'muddier')
data(ff_cdf, package = 'muddier')
data(fg_cdf, package = 'muddier')

load('ff_cis.rda')
load('fg_pmfs_cis.rda')
load('df_pmfs_cis.rda')

load('ff_cdfs_cis.rda')
load('fg_cdfs_cis.rda')
load('df_cdfs_cis.rda')

index <- sort(as.numeric(rownames(char_pmfs)))

png('df_cdf.png', width = 6, height = 4, units = 'in', res = 300)
plot(index, df_cdf, type = 'l', col = 'slateblue', log = 'x',
     ylim = c(0,1), main = 'Debris Flow Inherited Age CDF',
     xlab = 'Inherited Age (years)', ylab = 'Proportion Age <= X')
lines(index, df_cdfs_cis[1,], lty = 2)
lines(index, df_cdfs_cis[3,], lty = 2)
lines(index, df_cdfs_cis[2,], lty = 2, col = 'magenta3')
legend('bottomright', legend = c('observed', 'median', 'cis'),
       fill = c('slateblue', 'magenta3', 'black'))
dev.off()

png('df_exc.png', width = 6, height = 4, units = 'in', res = 300)
plot(index, df_exc, type = 'l', col = 'slateblue', ylim = c(10^-2,10^0), log='xy',
     main = 'Debris Flow Inherited Age Exceedance',
     xlab = 'Inherited Age (years)', ylab = 'Exceedance Probability')
lines(index, df_excs_cis[1,], lty = 2)
lines(index, df_excs_cis[3,], lty = 2)
lines(index, df_excs_cis[2,], lty = 2, col = 'magenta3')
legend('bottomleft', legend = c('observed', 'median', 'cis'),
       fill = c('slateblue', 'magenta3', 'black'))
dev.off()


png('ff_cdf.png', width = 6, height = 4, units = 'in', res = 300)
plot(index, ff_cdf, type = 'l', lwd = 3, col = 'slateblue',
     xlim = c(0,10000), ylim = c(0,1), main = 'Fluvial Fines Inherited Age CDF',
     xlab = 'Inherited Age (years)', ylab = 'Proportion of Sample at or Below Age X')
lines(index, ff_cdfs_cis[,1], lwd = 3, col = 'coral', lty = 2)
lines(index, ff_cdfs_cis[,3], lwd = 3, col = 'coral', lty = 2)
lines(index, ff_cdfs_cis[,2], lwd = 3, lty = 2)
legend('bottomright', legend = c('observed', 'median', 'cis'),
       fill = c('slateblue', 'black', 'coral'))
dev.off()

png('fg_cdf.png', width = 6, height = 4, units = 'in', res = 300)
plot(index, fg_cdf, type = 'l', lwd = 3, col = 'slateblue',
     xlim = c(0,11000), ylim = c(0,1), main = 'Fluvial Gravels Inherited Age CDF',
     xlab = 'Inherited Age (years)', ylab = 'Proportion of Sample at or Below Age X')
lines(index, fg_cdfs_cis[,1], lwd = 3, col = 'coral', lty = 2)
lines(index, fg_cdfs_cis[,3], lwd = 3, col = 'coral', lty = 2)
lines(index, fg_cdfs_cis[,2], lwd = 3, lty = 2)
legend('bottomright', legend = c('observed', 'median', 'cis'),
       fill = c('slateblue', 'black', 'coral'))
dev.off()

png('df_pmf.png', width = 6, height = 4, units = 'in', res = 300)
plot(index, df_pmf, type = 'l', lwd = 3, col = 'black',
     xlim = c(0,10000), main = 'Debris Flow Inherited Age PMF',
     xlab = 'Inherited Age (years)', ylab = 'Probability of Inherited Age X')
lines(index, df_pmfs_cis[,2], col = 'slateblue', lwd = 3)
#lines(index, df_pmfs_cis[,1], lwd = 3, col = 'coral')
lines(index, df_pmfs_cis[,3], lwd = 3, col = 'coral')
legend('topright', legend = c('observed', 'median', 'cis'),
       fill = c('black', 'slateblue', 'coral'))
dev.off()

png('ff_pmf.png', width = 6, height = 4, units = 'in', res = 300)
plot(index, ff_cis[3,], type = 'l', lwd = 3, col = 'coral',
     xlim = c(0,10000), main = 'Fluvial Fines Inherited Age PMF',
     xlab = 'Inherited Age (years)', ylab = 'Probability of Inherited Age X')
lines(index, ff_cis[2,], col = 'slateblue', lwd = 3)
#lines(index, ff_cis[1,], lwd = 3, col = 'coral')
lines(index, ff_pmf, lwd = 3, col = 'black')
legend('topright', legend = c('observed', 'median', 'cis'),
       fill = c('black', 'slateblue', 'coral'))
dev.off()

png('fg_pmf.png', width = 6, height = 4, units = 'in', res = 300)
plot(index, fg_pmfs_cis[,3], type = 'l', lwd = 3, col = 'coral',
     xlim = c(0,10000), main = 'Fluvial Gravels Inherited Age PMF',
     xlab = 'Inherited Age (years)', ylab = 'Probability of Inherited Age X')
lines(index, fg_pmfs_cis[,2], col = 'slateblue', lwd = 3)
#lines(index, fg_pmfs_cis[,1], lwd = 3, col = 'coral')
lines(index, fg_pmf, lwd = 3, col = 'black')
legend('topright', legend = c('observed', 'median', 'cis'),
       fill = c('black', 'slateblue', 'coral'))
dev.off()












begin <- Sys.time()
df_t1a <- list()
for (i in 1:2000) {
  df_t1a[[i]] <- df_t1[[i]]
}
df_t1a <- mclapply(df_t1a, rack) %>% rack
end <- Sys.time()
end - begin


begin <- Sys.time()
df_t1 <- mclapply(df_t1, rack) %>% rack
end <- Sys.time()
end - begin


dfb_pmf <- lapply(df_t, to_mat) %>% rack
dfb_pmf_cis <- apply(dfb_pmf, 1, get_cis)
dfb_cdf <- apply(dfb_pmf, 2, to_cdf)
dfb_cdf_cis <- apply(dfb_cdf, 1, get_cis)
dfb_exc <- apply(dfb_pmf, 2, to_exceed)
dfb_exc_cis <- apply(dfb_exc, 1, get_cis)

# plot results
plot(sort(index), dfb_pmf_cis[3,], type = 'l', lwd = 3, col = 'brown', xlim = c(0,2000), main = 'DF Inherited Age', xlab = 'inherited age', ylab = 'pmf')
lines(sort(index), dfb_pmf_cis[2,], lwd = 3, col = 'black')
lines(sort(index), dfb_pmf_cis[1,], lwd = 3, col = 'slateblue')
lines(sort(index), df_dt, lwd = 3, col = 'goldenrod')
legend('topright', legend = c('obs','upr','med','lwr'), fill = c('goldenrod', 'brown', 'black', 'slateblue'))


plot(sort(index), dfb_cdf_cis[3,], type = 'l', lwd = 3, col = 'brown', xlim = c(0,7000), main = 'DF Inherited Age', xlab = 'inherited age', ylab = 'cdf')
lines(sort(index), dfb_cdf_cis[2,], lwd = 3, col = 'black')
lines(sort(index), dfb_cdf_cis[1,], lwd = 3, col = 'slateblue')
lines(sort(index), df_cdf, lwd = 3, col = 'goldenrod')
legend('bottomright', legend = c('obs','upr','med','lwr'), fill = c('goldenrod', 'brown', 'black', 'slateblue'))

plot(sort(index), dfb_exc_cis[3,], type = 'l', lwd = 3, col = 'brown', xlim = c(0,11000), main = 'DF Inherited Age', xlab = 'inherited age', ylab = 'log10 exc')
lines(sort(index), dfb_exc_cis[2,], lwd = 3, col = 'black')
lines(sort(index), dfb_exc_cis[1,], lwd = 3, col = 'slateblue')
lines(sort(index), df_exc, lwd = 3, col = 'goldenrod')
legend('topright', legend = c('obs','upr','med','lwr'), fill = c('goldenrod', 'brown', 'black', 'slateblue'))



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



