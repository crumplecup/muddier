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

load('fl_200f.rds')
recf <- rec
load('fl_200a.rds')
recg <- rec
load('fl_200b.rds')
recg <- rbind(recg, rec)
rm(rec)


rec <- recg[recg$turnover < 350, ]
rec <- recf[recf$turnover < 350, ]



rec[rec$ad_a == max(rec$ad_a), ]
plot3D::points3D(rec$capture, rec$storage, rec$ad_a, ticktype = 'detailed', pch = 20,
                 phi = 40, theta = 310,
                 xlab = 'capture rate', ylab = 'storage rate', zlab = 'ad alpha')

rec[rec$ch_a == max(rec$ch_a), ]
rec[rec$ch_b == max(rec$ch_b), ]
plot3D::points3D(rec$capture, rec$storage, rec$ch_a, ticktype = 'detailed', pch = 20,
                 phi = 20, theta = 310,
                 xlab = 'capture rate', ylab = 'storage rate')

rec[rec$kp_a == max(rec$kp_a), ]
plot3D::points3D(rec$capture, rec$storage, rec$kp_a, ticktype = 'detailed', pch = 20,
                 phi = 30, theta = 310,
                 xlab = 'capture rate', ylab = 'storage rate')

rec[rec$ks_a == max(rec$ks_a), ]
rec[rec$ks_b == max(rec$ks_b), ]
plot3D::points3D(rec$capture, rec$storage, rec$ks_b, ticktype = 'detailed', pch = 20,
                 phi = 40, theta = 310,
                 xlab = 'capture rate', ylab = 'storage rate')

magicaxis::magplot(rec$turnover, rec$ks_a, pch = 20, col = get_palette('ocean'),
                   xlab = 'turnover period [y]', ylab = 'k-s fit')
abline(v = 191, col = get_palette('violet', .9), lwd = 2)
abline(v = 208, col = get_palette('crimson', .9), lwd = 2)
abline(v = 293, col = get_palette('gold', .9), lwd = 2)
abline(v = 318, col = get_palette('ocean', .9), lwd = 2)


gr_ch <- data.table::fread('/home/erik/reservoirs/data/debris_flow_transits_ch.csv')



gr_ks <- data.table::fread('/home/erik/output/stereotype_gravels_ks2.csv')
gr_kp <- data.table::fread('/home/erik/output/stereotype_gravels_kp_1000.csv')
gr_tr <- data.table::fread('/home/erik/output/fines_stereo_kp_1000.csv')

gr_ad <- data.table::fread('/home/erik/reservoirs/data/debris_flow_transits_ad.csv')
gr_ch <- data.table::fread('/home/erik/reservoirs/data/debris_flow_transits_ch.csv')
gr_kp <- data.table::fread('/home/erik/reservoirs/data/debris_flow_transits_kp.csv')
gr_ks <- data.table::fread('/home/erik/reservoirs/data/debris_flow_transits_ks.csv')
gr_ad <- data.table::fread('/home/erik/output/stereotype_gravels_ad_1007.csv')
gr_ks <- data.table::fread('/home/erik/output/stereotype_gravels_ks_1002.csv')
gr_ch <- data.table::fread('/home/erik/output/stereotype_gravels_ch_1006.csv')
gr_kp <- data.table::fread('/home/erik/output/stereotype_gravels_kp_1017.csv')


png('reservoirs_gravels_log.png', width = 21, height = 17, units = 'cm', res = 300)
magicaxis::magplot(emp_cdf(charcoal$mn[charcoal$facies == 'FF'] + 50), pch = 20, col = get_palette('coral'),
                   xlab = 'Charcoal Age [y]', ylab = 'CDF of samples', xlim = c(1, 30000), log = 'x')
points(emp_cdf(debris_flows), pch = 20, col = get_palette('hardwood'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FG'] + 50), pch = 20, col = get_palette('charcoal'))
lines(emp_cdf(gr_ch$V1), lwd = 2.5, col = get_palette('gold', .7))
lines(emp_cdf(gr_kp$V1), lwd = 2.5, col = get_palette('violet', .7))
lines(emp_cdf(gr_ks$V1), lwd = 2.5, col = get_palette('crimson', .7))
lines(emp_cdf(gr_ad$V1), lwd = 2.5, col = get_palette('ocean', .7))
legend('topleft', legend = c('Debris-Flows', 'Gravels', 'Fines', 'Anderson-Darling', 'Chi-squared', 'Kolmogorov-Smirnov', 'Kuiper'),
       lty = c(NA, NA, NA, 1, 1, 1, 1), lwd = c(NA, NA, NA, 2, 2, 2, 2),
       pch = c(20, 20, 20, NA, NA, NA, NA), col = get_palette(c('hardwood', 'charcoal', 'coral', 'ocean', 'gold', 'crimson', 'violet'), .9))
dev.off()

png('reservoirs_gravels.png', width = 21, height = 17, units = 'cm', res = 300)
magicaxis::magplot(emp_cdf(charcoal$mn[charcoal$facies == 'FF'] + 50), pch = 20, col = get_palette('coral'),
                   xlab = 'Charcoal Age [y]', ylab = 'CDF of samples', xlim = c(0, 25000))
points(emp_cdf(debris_flows), pch = 20, col = get_palette('hardwood'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FG'] + 50), pch = 20, col = get_palette('charcoal'))
lines(emp_cdf(gr_ad$V1), lwd = 2.5, col = get_palette('ocean', .7))
lines(emp_cdf(gr_kp$V1), lwd = 2.5, col = get_palette('violet', .7))
lines(emp_cdf(gr_ks$V1), lwd = 2.5, col = get_palette('crimson', .7))
lines(emp_cdf(gr_ch$V1), lwd = 2.5, col = get_palette('gold', .7))
legend('bottomright', legend = c('Debris-Flows', 'Gravels', 'Fines', 'Anderson-Darling', 'Chi-squared', 'Kolmogorov-Smirnov', 'Kuiper'),
       lty = c(NA, NA, NA, 1, 1, 1, 1), lwd = c(NA, NA, NA, 2, 2, 2, 2),
       pch = c(20, 20, 20, NA, NA, NA, NA), col = get_palette(c('hardwood', 'charcoal', 'coral', 'ocean', 'gold', 'crimson', 'violet'), .9))


dev.off()

gr_ad <- data.table::fread('/home/erik/output/gravels_cdf_ad_1000.csv') %>% unlist
gr_ch <- data.table::fread('/home/erik/output/gravels_cdf_ch_1000.csv') %>% unlist
gr_kp <- data.table::fread('/home/erik/output/gravels_cdf_kp_1000.csv') %>% unlist
gr_ks <- data.table::fread('/home/erik/output/gravels_cdf_ks_1000.csv') %>% unlist


png('reservoirs_gravels.png', width = 21, height = 17, units = 'cm', res = 300)
magicaxis::magplot(emp_cdf(charcoal$mn[charcoal$facies == 'FF'] + 50), pch = 20, col = get_palette('coral'),
                   xlab = 'Charcoal Age [y]', ylab = 'CDF of samples', xlim = c(0, 25000))
points(emp_cdf(debris_flows), pch = 20, col = get_palette('hardwood'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FG'] + 50), pch = 20, col = get_palette('charcoal'))
lines(cumsum(gr_ad), lty = 2, lwd = 1.5, col = get_palette('ocean', .7))
lines(cumsum(gr_kp), lty = 2, lwd = 1.5, col = get_palette('violet', .7))
lines(cumsum(gr_ks), lty = 2, lwd = 1.5, col = get_palette('crimson', .7))
lines(cumsum(gr_ch), lty = 2, lwd = 1.5, col = get_palette('gold', .7))

lines(emp_cdf(gr_ks$V1), lwd = 2, col = get_palette('crimson', .7))
lines(emp_cdf(gr_kp$V1), lwd = 2, col = get_palette('violet', .7))
lines(emp_cdf(gr_ad$V1), lwd = 2, col = get_palette('ocean', .7))
lines(emp_cdf(gr_ch$V1), lwd = 2, col = get_palette('gold', .7))
legend('bottomright', legend = c('Debris-Flows', 'Gravels', 'Fines', 'Anderson-Darling', 'Chi-squared', 'Kolmogorov-Smirnov', 'Kuiper', 'Best', 'Mean'),
       lty = c(NA, NA, NA, 1, 1, 1, 1, 1, 2), lwd = c(NA, NA, NA, 2, 2, 2, 2, 1, 1),
       pch = c(20, 20, 20, NA, NA, NA, NA, NA, NA), col = get_palette(c('hardwood', 'charcoal', 'coral', 'ocean', 'gold', 'crimson', 'violet', 'charcoal', 'charcoal'), .9))
dev.off()

gr_src <- (gr_ad + gr_ch + gr_ks) / 3
gr_src[gr_src < 0] <- 0
lines(cumsum(gr_src))
write.csv(gr_src, file = 'gravel_cdf_src.csv', row.names = F)

gr_ad[gr_ad < 0] <- 0
gr_ch[gr_ch < 0] <- 0
gr_kp[gr_kp < 0] <- 0
gr_ks[gr_ks < 0] <- 0

gr_ch <- data.table::fread('/home/erik/output/stereotype_fines_ch_1011.csv')
gr_ad <- data.table::fread('/home/erik/output/stereotype_fines_ad_1004.csv')
gr_kp <- data.table::fread('/home/erik/output/stereotype_fines_kp_1002.csv')
gr_ks <- data.table::fread('/home/erik/output/stereotype_fines_ks_1011.csv')

gr_ad <- data.table::fread('/home/erik/output/fines_cdf_ad_1000.csv') %>% unlist
gr_ch <- data.table::fread('/home/erik/output/fines_cdf_ch_1000.csv') %>% unlist
gr_kp <- data.table::fread('/home/erik/output/fines_cdf_kp_1000.csv') %>% unlist
gr_ks <- data.table::fread('/home/erik/output/fines_cdf_ks_1000.csv') %>% unlist

png('reservoirs_fines.png', height = 17, width = 21, units = 'cm', res = 300)
magicaxis::magplot(emp_cdf(charcoal$mn[charcoal$facies == 'FF'] + 50), pch = 20, col = get_palette('coral'),
                   xlab = 'Charcoal Age [y]', ylab = 'CDF of samples', xlim = c(0, 20000))
points(emp_cdf(debris_flows), pch = 20, col = get_palette('hardwood'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FG'] + 50), pch = 20, col = get_palette('charcoal'))
lines(emp_cdf(gr_ks$V1), lwd = 2, col = get_palette('crimson', .7))
lines(emp_cdf(gr_kp$V1), lwd = 2, col = get_palette('violet', .7))
lines(emp_cdf(gr_ad$V1), lwd = 2, col = get_palette('ocean', .7))
lines(emp_cdf(gr_ch$V1), lwd = 2, col = get_palette('gold', .7))

lines(cumsum(gr_ad), lty = 2, lwd = 1.5, col = get_palette('ocean', .7))
lines(cumsum(gr_kp), lty = 2, lwd = 1.5, col = get_palette('violet', .7))
lines(cumsum(gr_ks), lty = 2, lwd = 1.5, col = get_palette('crimson', .7))
lines(cumsum(gr_ch), lty = 2, lwd = 1.5, col = get_palette('gold', .7))
legend('bottomright', legend = c('Debris-Flows', 'Gravels', 'Fines', 'Anderson-Darling', 'Chi-squared', 'Kolmogorov-Smirnov', 'Kuiper', 'Best', 'Mean'),
       lty = c(NA, NA, NA, 1, 1, 1, 1, 1, 2), lwd = c(NA, NA, NA, 2, 2, 2, 2, 1, 1),
       pch = c(20, 20, 20, NA, NA, NA, NA, NA, NA), col = get_palette(c('hardwood', 'charcoal', 'coral', 'ocean', 'gold', 'crimson', 'violet', 'charcoal', 'charcoal'), .9))
dev.off()

lines(cumsum(unlist(gr_ad)), lwd = 2.5, col = get_palette('ocean', .7))
lines(cumsum(unlist(gr_ch)), lwd = 2.5, col = get_palette('gold', .7))
lines(cumsum(unlist(gr_kp)), lwd = 2.5, col = get_palette('violet', .7))
lines(cumsum(unlist(gr_ks)), lwd = 2.5, col = get_palette('crimson', .7))

lines(emp_cdf(gr_ks$V1), lwd = 2.5, col = get_palette('crimson', .7))
lines(emp_cdf(gr_kp$V1), lwd = 2.5, col = get_palette('violet', .7))

write.csv(gr_ad, file = 'gravels_cdf_ad.csv', row.names = F, col.names = F)
write.csv(gr_ch, file = 'gravels_cdf_ch.csv', row.names = F, col.names = F)
write.csv(gr_kp, file = 'gravels_cdf_kp.csv', row.names = F, col.names = F)
gr_ks1 <- gr_ks %>% unlist %>% as.numeric
write.csv(gr_ks1, file = 'gravels_cdf_ks.csv', row.names = F, col.names = F)

df_ad <- data.table::fread('/home/erik/reservoirs/data/debris_flow_transits_ad.csv') %>% unlist
df_ch <- data.table::fread('/home/erik/reservoirs/data/debris_flow_transits_ch.csv') %>% unlist
df_kp <- data.table::fread('/home/erik/reservoirs/data/debris_flow_transits_kp.csv') %>% unlist
df_ks <- data.table::fread('/home/erik/reservoirs/data/debris_flow_transits_ks.csv') %>% unlist
lines(cumsum(df_ad), lwd = 2.5, col = get_palette('ocean', .7))

png('reservoirs_debris_flows.png', width = 21, height = 17, units = 'cm', res = 300)
magicaxis::magplot(emp_cdf(charcoal$mn[charcoal$facies == 'FF'] + 50), pch = 20, col = get_palette('coral'),
                   xlab = 'Charcoal Age [y]', ylab = 'CDF of samples', xlim = c(0, 20000))
points(emp_cdf(debris_flows), pch = 20, col = get_palette('hardwood'))
points(emp_cdf(charcoal$mn[charcoal$facies == 'FG'] + 50), pch = 20, col = get_palette('charcoal'))
lines(cumsum(df_ad), lwd = 1.5, col = get_palette('ocean', .7))
lines(cumsum(df_kp), lwd = 1.5, col = get_palette('violet', .7))
lines(cumsum(df_ks), lwd = 1.5, col = get_palette('crimson', .7))
lines(cumsum(df_ch), lwd = 1.5, col = get_palette('gold', .7))

legend('bottomright', legend = c('Debris-Flows', 'Gravels', 'Fines', 'Anderson-Darling', 'Chi-squared', 'Kolmogorov-Smirnov', 'Kuiper', 'Observed', 'Modeled Mean'),
       lty = c(NA, NA, NA, 1, 1, 1, 1, NA, 1), lwd = c(NA, NA, NA, 2, 2, 2, 2, NA, 2),
       pch = c(20, 20, 20, NA, NA, NA, NA, 20, NA), col = get_palette(c('hardwood', 'charcoal', 'coral', 'ocean', 'gold', 'crimson', 'violet', 'slate', 'slate'), .9))
dev.off()


fn_src <- convo_add(unlist(gr_ad), unlist(df_ad), 0:20000)
lines(0:20000, cumsum(fn_src))

fn_src_ad <- unlist((df_ad + gr_ad))
fn_src_ch <- unlist((df_ch + gr_ch))
fn_src_kp <- unlist((df_kp + gr_kp))
fn_src_ks <- unlist((df_ks + gr_ks))
fn_src_ks[fn_src_ks < 0] <- 0

write.csv(fn_src_ad, file = 'fines_source_ad.csv', row.names = F, col.names = F)
write.csv(fn_src_ch, file = 'fines_source_ch.csv', row.names = F, col.names = F)
write.csv(fn_src_kp, file = 'fines_source_kp.csv', row.names = F, col.names = F)
write.csv(fn_src_ks, file = 'fines_source_ks.csv', row.names = F, col.names = F)

?write.csv

load('fn_200.rds')
recf2 <- rec
load('fn_1000.rds')
recf2 <- rbind(rec, recf2)
rec <- recf2


load('gr_1000.rds')



magicaxis::magplot(rec$capture, rec$ks, pch = 20, col = get_palette('ocean'),
                    xlab = 'capture rate', ylab = 'k-s fit')
magicaxis::magplot(rec$storage, rec$ks, pch = 20, col = get_palette('ocean'),
                    xlab = 'storage rate', ylab = 'k-s fit')

gr_st <- bin_dual_stat(rec$capture, rec$storage, rec$ks, 11)
gr_st <- bin_dual_stat(rec$capture_rate_gravels, rec$storage_rate_gravels, rec$ks2, 35)
gr_st <- bin_dual_stat(rec$capture_rate_gravels, rec$storage_rate_gravels, rec$ks1, 35)
plot3D::points3D(gr_st$rng1, gr_st$rng2, gr_st$mns, ticktype = 'detailed', pch = 20,
                 phi = 40, theta = 310)
mean(gr_st$ns[!is.na(gr_st$ns)])

gr_min <- gr_st[!is.nan(gr_st$mns) & !is.na(gr_st$mns), ]
(gr_min <- gr_min[gr_min$mns == min(gr_min$mns), ])
gr_cap <- gr_st$rng1[gr_st$lwr <= gr_min$upr]
range(gr_cap[!is.na(gr_cap)])
gr_str <- gr_st$rng2[gr_st$lwr <= gr_min$upr]
range(gr_str[!is.na(gr_str)])

bins <- 60
gr_st_ad <- bin_dual_stat(rec$capture_rate_gravels, rec$storage_rate_gravels, rec$ad1, bins)
gr_st_ch <- bin_dual_stat(rec$capture_rate_gravels, rec$storage_rate_gravels, rec$ch, bins)
gr_st_kp <- bin_dual_stat(rec$capture_rate_gravels, rec$storage_rate_gravels, rec$kp, bins)
gr_st_ks <- bin_dual_stat(rec$capture_rate_gravels, rec$storage_rate_gravels, rec$ks1, bins)
gr_st <- gr_st_ad
gr_st <- gr_st_ch
gr_st <- gr_st_kp
gr_st <- gr_st_ks
plot3D::points3D(gr_st$rng1, gr_st$rng2, gr_st$mns, ticktype = 'detailed', pch = 20,
                 phi = 40, theta = 310)
mean(gr_st$ns[!is.na(gr_st$ns)])

gr_max <- gr_st[!is.nan(gr_st$mns) & !is.na(gr_st$mns), ]
(gr_max <- gr_max[gr_max$mns == max(gr_max$mns), ][1,])
gr_inc <- gr_st[gr_st$upr >= gr_max$lwr & gr_st$upr >= gr_max$lwr, ]
gr_cap <- gr_st$rng1[gr_st$upr >= gr_max$lwr]
range(gr_cap[!is.na(gr_cap)])
gr_str <- gr_st$rng2[gr_st$upr >= gr_max$lwr]
range(gr_str[!is.na(gr_str)])

gr_inc_ad <- gr_inc
gr_inc_ch <- gr_inc
gr_inc_kp <- gr_inc
gr_inc_ks <- gr_inc

gr_inc_rng1 <- range(c(gr_inc_ad$rng1, gr_inc_ch$rng1, gr_inc_kp$rng1, gr_inc_ks$rng1))
gr_inc_rng2 <- range(c(gr_inc_ad$rng2, gr_inc_ch$rng2, gr_inc_kp$rng2, gr_inc_ks$rng2))

magicaxis::magplot(gr_inc$rng1, gr_inc$rng2, col = get_palette('slate', .001),
                   xlim = gr_inc_rng1, ylim = gr_inc_rng2,
                   xlab = 'Capture Rate', ylab = 'Storage Rate')
points(gr_inc_ad$rng1, gr_inc_ad$rng2, pch = 20, col = get_palette('ocean'))
points(gr_inc_ch$rng1, gr_inc_ch$rng2, pch = 20, col = get_palette('gold'))
points(gr_inc_kp$rng1, gr_inc_kp$rng2, pch = 20, col = get_palette('violet'))
points(gr_inc_ks$rng1, gr_inc_ks$rng2, pch = 20, col = get_palette('crimson'))

rec <- data.table::fread('/home/erik/output/gravels_200_1000.csv')
rec <- data.table::fread('/home/erik/output/gravel_hits_200_1000.csv')
rec <- rbind(rec, data.table::fread('/home/erik/output/gravel_hits_ks_200_1001.csv'))
rec <- rbind(rec, data.table::fread('/home/erik/output/gravel_hits_ks_200_1002.csv'))
rec <- rbind(rec, data.table::fread('/home/erik/output/gravel_hits_ks_200_1003.csv'))
rec <- rbind(rec, data.table::fread('/home/erik/output/gravel_hits_ks_200_1004.csv'))
rec <- data.table::fread('/home/erik/output/gravel_hits_kp_1000.csv')
rec <- data.table::fread('/home/erik/output/fines_hits_ad_1000.csv')
rec <- data.table::fread('/home/erik/output/fines_hits_200_1001.csv')
rec <- rbind(rec, data.table::fread('/home/erik/output/fines_hits_200_1002.csv'))
rec <- rbind(rec, data.table::fread('/home/erik/output/fines_hits_200_1003.csv'))
rec <- data.table::fread('/home/erik/output/fines_hits_ks_1000.csv')



rec[rec$ad1 == max(rec$ad1), ]
plot3D::points3D(rec$capture_rate_gravels, rec$storage_rate_gravels, rec$ad1, ticktype = 'detailed', pch = 20,
                 phi = 20, theta = 310,
                 xlab = 'capture rate', ylab = 'storage rate', zlab = 'ad alpha')

rec[rec$ch == max(rec$ch), ]

plot3D::points3D(rec$capture_rate_gravels, rec$storage_rate_gravels, rec$ch, ticktype = 'detailed', pch = 20,
                 phi = 20, theta = 340,
                 xlab = 'capture rate', ylab = 'storage rate')

rec[rec$kp == max(rec$kp), ]
plot3D::points3D(rec$capture_rate_gravels, rec$storage_rate_gravels, rec$kp, ticktype = 'detailed', pch = 20,
                 phi = 30, theta = 310,
                 xlab = 'capture rate', ylab = 'storage rate')

rec[rec$ks1 == max(rec$ks1), ]
plot3D::points3D(rec$capture_rate_gravels, rec$storage_rate_gravels, rec$ks1, ticktype = 'detailed', pch = 20,
                 phi = 25, theta = 340,
                 xlab = 'capture rate', ylab = 'storage rate')


# cedar creek

cd <- creeks[creeks$creek_name == 'cedar', ]
names(cd)
plot(cd$DebrisFlow, cd$xsec_area)
