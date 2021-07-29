library(magrittr)
library(parallel)
options(mc.cores=parallel::detectCores())
setwd('/media/erik/catacomb/research')
setwd('/home/erik/output')

debris_flows <- charcoal$mn[charcoal$facies == 'DF'] + 50
gravels <- charcoal$mn[charcoal$facies == 'FG'] + 50
fines <- charcoal$mn[charcoal$facies == 'FF'] + 50

# debris-flow pmf and index
storage_index <- char_pmfs %>% rownames %>% as.numeric %>% rev + 50
dfc <- emp_cdf(charcoal$mn[charcoal$facies == 'DF'])
dfi <- dfc[ , 1]
dfi <- dfi + 50
dfp <- to_pmf(dfc[ , 2])

# gravels pmf and index
grc <- emp_cdf(gravels)
gri <- grc[ , 1]
grp <- to_pmf(grc[ , 2])

muddier::fines_synth(.12, .12, 208, source_index = dfi, source_prob = dfp, storage_index = storage_index,
                     gravel_index = gri, gravel_prob = grp)

muddier::test_fit(muddier::gof(
  muddier::gravel_synth(.12, .12, 208, source_index = dfi, source_prob = dfp, storage_index = storage_index),
  gravels
))

muddier::gravel_fit_n(10, .12, .12, 208, source_index = dfi, source_prob = dfp, storage_index = storage_index, gravels = gravels)
muddier::fines_fit_n(10, .12, .12, 208, source_index = dfi, source_prob = dfp, storage_index = storage_index,
                     gravel_index = gri, gravel_prob = grp, fines = fines)

rec <- muddier::fluvial_fit(10, 200, 0, 1, 0, 1, 50, 10000, storage_index, dfi, dfp, gri, grp, fines, fines = T)

dur <- as.difftime(50, units = 'hours')
begin <- Sys.time()
end <- begin + dur
while (Sys.time() < end) {
  rec <- rbind(rec, muddier::fluvial_fit(
    batch = 10,
    n = 200,
    min_cap = 0,
    max_cap = 1,
    min_stor = 0,
    max_stor = 1,
    min_turn = 50,
    max_turn = 10000,
    storage_index = storage_index,
    source_index = dfi,
    source_prob = dfp,
    gravel_index = gri,
    gravel_prob = grp,
    obs = fines,
    fines = F))
  save(rec, file = 'fl_200f.rds')
}


rec[rec$ad_a == max(rec$ad_a), ]
plot3D::points3D(rec$capture, rec$storage, rec$ad_a, ticktype = 'detailed', pch = 20,
                 phi = 40, theta = 30,
                 xlab = 'capture rate', ylab = 'storage rate', zlab = 'ad alpha')

rec[rec$ch_a == max(rec$ch_a), ]
plot3D::points3D(rec$capture, rec$storage, rec$ch_a, ticktype = 'detailed', pch = 20,
                 phi = 20, theta = 290,
                 xlab = 'capture rate', ylab = 'storage rate')

rec[rec$kp_a == max(rec$kp_a), ]
plot3D::points3D(rec$capture, rec$storage, rec$kp_a, ticktype = 'detailed', pch = 20,
                 phi = 30, theta = 310,
                 xlab = 'capture rate', ylab = 'storage rate')

rec[rec$ks_a == max(rec$ks_a), ]
plot3D::points3D(rec$capture, rec$storage, rec$ks_a, ticktype = 'detailed', pch = 20,
                 phi = 40, theta = 310,
                 xlab = 'capture rate', ylab = 'storage rate')

magicaxis::magplot(rec$turnover, rec$ks_a, pch = 20, col = get_palette('ocean'),
                   xlab = 'turnover period [y]', ylab = 'k-s fit')
abline(v = 191, col = get_palette('violet', .9), lwd = 2)
abline(v = 208, col = get_palette('crimson', .9), lwd = 2)
abline(v = 293, col = get_palette('gold', .9), lwd = 2)
abline(v = 318, col = get_palette('ocean', .9), lwd = 2)
