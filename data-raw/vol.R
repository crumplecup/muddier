library(data.table)
library(magrittr)
library(parallel)
library(plot3D)
options(mc.cores=detectCores())
set.seed(10101)

# figure out which charcoal sample corresponds to which node in the creeks data
crk <- creeks[creeks$creek_name %in% c('bear', 'knowles'), ]
crk <- crk[crk$NODE_ID %in% creeks_radio$node_ids, ]

setwd('/home/crumplecup/work/muddier')
# sa_nodes <- crk
# usethis::use_data(sa_nodes)
# setwd('/home/crumplecup/work/')


crks <- creeks[creeks$creek_name %in% c('bear', 'knowles'), ]
wt <- lapply(crks$DebrisFlow, function(x) weight_by_delprob(x)) %>% unlist
crks$dp <- crks$DebrisFlow * wt

crks$lvl <- crks$xsec_area / crks$valley_width

# global input/output rates
gi <- 0.11
go <- 0.08

# localized weighted rate
# nrow(crk) is the number of sampled nodes in the study areas
# not the total nodes in the study area crks

ri_ave <- gi / nrow(crk)
ri_ttl <- ri_ave * nrow(crks)
ri_dp <- ri_ave * crks$dp
crks$ri <- ri_dp * (ri_ttl / sum(ri_dp))

bit <- (log(crks$dp) - mean(log(crks$dp)))
bit <- (1/bit)


# relative optimized delivery probability
crks$rdp <- crks$dp / max(crks$dp)

# output rate is dependent on streampower
# relative streampower based on contr area and slope

crks$pwr <- log(crks$contr_area+1) * log(crks$slope + 1)

pwr_df <- data.frame(CA = log(crks$contr_area),
                     slope = log(crks$slope),
                     creek = crks$creek_name)
pwr_mod <- lm(CA ~ slope + creek, data = pwr_df)
summary(pwr_mod)

pwr_prd <- predict(
  pwr_mod, newdata = data.frame(slope = log(crks$slope), creek = crks$creek_name)
)

creek_bool <- rep(0, nrow(crks))
creek_bool[crks$creek_name == 'knowles'] <- 1

lca_obs <- log(crks$slope) * pwr_mod$coefficients[2] +
  creek_bool * pwr_mod$coefficients[3]

lca_dif <- log(crks$contr_area) - pwr_prd
difs <- (exp(lca_dif))
# values of difs from zero to one had negative logged difs
# no need to normalize to use as a weight, as negatives are deflated, positive inflated
# streampower coefficient (difs) is normalized to use as a weight

# output weighted by inverse dp and streampower coefficient
ro_ave <- go / nrow(crk)
ro_ttl <- ro_ave * nrow(crks)
ro_wt <- (1 - crks$dp / max(crks$dp)) * difs
crks$ro <- ro_wt * (ro_ttl / sum(ro_wt))
crks$ro[crks$ro <= 0] <- 0.00000001

crks$ro <- difs * (ro_ttl / sum(difs))
bat <- (1 - crks$dp / max(crks$dp))
crks$ro <- bat * (ro_ttl / sum(bat))

ro_dif <- crks$ri - ro_ave
length(ro_dif[ro_dif > 0]) / length(ro_dif)
plot(ro_dif)

crks$ri <- ri_ave
crks$ro <- ro_ave

# weighted by streampower
crks$ro <- (crks$ro - crks$ri) * difs + crks$ri
crks$ro <- crks$ro * ro_ttl / sum(crks$ro)
# weighted by inverse delivery probability
crks$ro <- (crks$ro - crks$ri) * bat + crks$ri
crks$ro <- crks$ro * ro_ttl / sum(crks$ro)

ldp <- log(crks$dp)
ladp <- mean(ldp)
ldp_difs <- ldp - ladp

# set up input/output rates

# both streampower and inverse delivery prob
ro_scl <- 1
crks$ri <- ri_ave
crks$ro <- ro_ave
crks$ro <- (crks$ro - crks$ri) * difs * bat * ro_scl + crks$ri
crks$ro[crks$ro <= 0] <- min(crks$ro[crks$ro > 0])
crks$ro <- crks$ro * ro_ttl / sum(crks$ro)


ri_scl <- 1
crks$ri <- ri_ave
ri_dp <- crks$ro - (crks$ro - crks$ri) * ldp_difs * ri_scl
ri_dp[ri_dp <= 0] <- min(ri_dp[ri_dp > 0])
crks$ri <- ri_dp * (ri_ttl / sum(ri_dp))


length(crks$ri[crks$ri > crks$ro]) / length(crks$ri)
plot(crks$ri - crks$ro, pch = 20, col = get_palette('charcoal'))
abline(h = 0, lty = 2)

# simulate accumulation record for creek nodes

record <- search_topology(
  crks, #[crks$creek_name == 'knowles', ],
  c(0,1), c(0,300), bt = 100, by = 'vol',
  w_scl = 4, backfill = F)
record[record[, 3] < .10, ]
scatter3D(record[,1], record[,2], record[,3],
          ticktype = 'detailed', pch = 20, theta = 150, phi = 50)

rec <- record


png('vol_9-18.png', height = 20, width = 20, units = 'cm', res = 300)
scatter3D(record[,1], record[,2], record[,3],
          ticktype = 'detailed', pch = 20, theta = 150, phi = 50)
dev.off()


best_so_far <- rec
rec <- best_so_far
linked <- rec

save(best_so_far, file = 'bookmark_9-18.rds')
load('bookmark_9-18.rds')

best_fit <- record[record[,3] < .10, ]

pred <- fit_volumes(crks, 10000, best_fit[2,1], best_fit[2,2], 20)
png('vol_fit_9-18.png', height = 15, width = 17, units = 'cm', res = 300)
plot(crks$xsec_area %>% emp_cdf, pch = 20, col = get_palette('ocean', .2),
     xlab = 'unit cross-sectional area m2', ylab = 'CDF')
lines(pred$vol %>% emp_cdf, lwd = 3, col = get_palette('charcoal', .6))
legend('bottomright', legend = c('observed', 'fit'),
       fill = get_palette(c('ocean', 'charcoal'), .7))
dev.off()

plot(crks$lvl %>% emp_cdf)
lines(pred$lvl %>% emp_cdf, lwd = 3, col = get_palette('ocean'))


# linked-bucket model

record <- search_topology(
  crks[crks$creek_name == 'knowles', ], c(0,1), c(0,100),
  bt = 100, it = 20, by = 'vol', w_scl = 3, linked = T)
record[record[, 5] < .20, ]
rec <- record
# png('vol_fit_9-19b.png', height = 15, width = 17, units = 'cm', res = 300)
scatter3D(rec[,1], rec[,2], rec[,3],
          ticktype = 'detailed', pch = 20, theta = 270, phi = 30)
# dev.off()

mins <- rec[rec[,3] < .12, ]
scatter3D(mins[,1], mins[,2], mins[,7]*10,
          ticktype = 'detailed', pch = 20, theta = 320, phi = 50)


best_fit <- record[record[,3] == min(record[,3]), ]

pred <- fit_volumes1(crks[crks$creek_name == 'knowles', ], 10000, best_fit[1], best_fit[2], 10)
pred <- pred[[1]]
#png('vol_fit_9-18.png', height = 15, width = 17, units = 'cm', res = 300)
plot(crks$xsec_area %>% emp_cdf, pch = 20, col = get_palette('ocean', .2),
     xlab = 'unit cross-sectional area m2', ylab = 'CDF')
lines(pred$vol %>% emp_cdf, lwd = 3, col = get_palette('charcoal', .6))
legend('bottomright', legend = c('observed', 'fit'),
       fill = get_palette(c('ocean', 'charcoal'), .7))
#dev.off()

plot(crks$lvl[crks$creek_name == 'knowles'] %>% emp_cdf, pch = 20, col = get_palette('ocean', .2),
     xlab = 'average valley bottom depth (m)', ylab = 'CDF')
lines(pred$lvl %>% emp_cdf, lwd = 3, col = get_palette('charcoal', .7))
legend('bottomright', legend = c('observed', 'fit'),
       fill = get_palette(c('ocean', 'charcoal'), .7))

save(rec, file = 'linked_9-20.rds')
load('linked_9-20.rds')


# make dataframe of study area node coords, slope
crk_coords <- coordinates(crks)
cslp <- crks$slope
bslp <- crks$slope[crks$creek_name == 'bear']
bkm <- crks$ToMouth_km[crks$creek_name == 'bear']
kslp <- crks$slope[crks$creek_name == 'knowles']
kkm <- crks$ToMouth_km[crks$creek_name == 'knowles']


belev <- slope_to_elev(bslp, bkm)
kelev <- slope_to_elev(kslp, kkm)
elev <- c(belev,kelev)
plot(elev)

creeks_xyz <- data.frame(
  x_data = crk_coords[,1],
  y_data = crk_coords[,2],
  z_data = elev,
  slope = crks$slope
)

write.csv(creeks_xyz, file = 'creeks_xyz.csv')

crks$elev <- elev

options(scipen = 5)

# backfill model
record <- search_topology(
  crks[crks$creek_name == 'knowles', ], c(0,1), c(0,100),
  bt = 100, it = 20, by = 'vol', w_scl = 3)
record[record[, 3] < .20, ]
rec <- record
png('vol_fit_9-21.png', height = 15, width = 17, units = 'cm', res = 300)
scatter3D(rec[,1], rec[,2], rec[,5],
          ticktype = 'detailed', pch = 20, theta = 320, phi = 30)
dev.off()




best_fit <- record[record[,3] == min(record[,3]), ]

pred <- backfill(crks[crks$creek_name == 'knowles', ], 10000, best_fit[1], best_fit[2], 10)
pred <- pred[[1]]
#png('vol_fit_9-21a.png', height = 15, width = 17, units = 'cm', res = 300)
plot(crks$xsec_area %>% emp_cdf, pch = 20, col = get_palette('ocean', .2),
     xlab = 'unit cross-sectional area m2', ylab = 'CDF')
lines(pred$vol %>% emp_cdf, lwd = 3, col = get_palette('charcoal', .6))
legend('bottomright', legend = c('observed', 'fit'),
       fill = get_palette(c('ocean', 'charcoal'), .7))
#dev.off()

png('vol_fit_9-21c.png', height = 15, width = 17, units = 'cm', res = 300)
plot(crks$lvl[crks$creek_name == 'knowles'] %>% emp_cdf, pch = 20, col = get_palette('ocean', .2),
     xlab = 'average valley depth m2', ylab = 'CDF')
lines(pred$lvl %>% emp_cdf, lwd = 3, col = get_palette('charcoal', .7))
dev.off()

setwd('/home/crumplecup/work/')
save(rec, file = 'backfill_9-21.rds')
load('backfill_9-21.rds')



mins <- rec[rec[,3] < .12, ]
scatter3D(mins[,1], mins[,2], mins[,7]*10,
          ticktype = 'detailed', pch = 20, theta = 320, phi = 50)

kn <- creeks[creeks$NODE_ID %in% crks$NODE_ID, ]
kn <- kn[kn$creek_name == 'knowles', ]
kn <- creeks[creeks$creek_name == 'knowles', ]


begin <- Sys.time()
diff <- as.difftime('01:00:00', '%H:%M:%S', units = 'hour')
end <- begin + diff
i <- 0
while (Sys.time() < end) {
  res <- n_dimensional_search(
    kn,
    c(.1, .15), c(.05, .13), c(.05, .4), c(.05, .6),
    ti_scl = c(0.01, 3),
    to_scl = c(0.01, 3),
    ks_scl = c(0.01, 3),
    dp_scl = c(0.01, 2),
    m_scl = c(0.01, 1),
    b_scl = c(0.01, 1),
    it = 30,
    batch = 100
  )
  rec <- rbind(rec, res)
  save(rec, file = "run_20200929.rds")
}

# sb <- rec[rec[,9] < .2 &
#       rec[,1] > .1 & rec[,1] < .3, ]
# scatter3D(sb[,3], sb[,4], (sb[,9] + sb[,11])/2,
#           ticktype = 'detailed', pch = 20, theta = 320, phi = 40,
#           xlab = 'input rate', ylab = 'output rate', zlab = 'test value')

(good <- rec[rec[,10] + rec[,12] == min(rec[,10] + rec[,12]), ])
scatter3D(rec[,3], rec[,4], (rec[,10] + rec[,12])/2,
          ticktype = 'detailed', pch = 20, theta = 320, phi = 40,
          xlab = 'input rate', ylab = 'output rate', zlab = 'test value')

scatter3D(rec[,3], rec[,4], rec[,9],
          ticktype = 'detailed', pch = 20, theta = 320, phi = 40,
          xlab = 'input rate', ylab = 'output rate', zlab = 'test value')

# rec <- rec[1,]
# og_rec <- rec
rec <- res
# rec <- rbind(rec, og_rec)
# rownames(rec) <- NULL
(good <- rec[rec[,9] < .20 & rec[, 11] < .20, ])
(good <- rec[rec[,10] == min(rec[,10]), ])
(good <- rec[rec[,12] == min(rec[,12]), ])
(good <- rec[rec[,10] + rec[,12] == min(rec[,10] + rec[,12]), ])
good <- good[1,]
good <- res[3,]

best_fit <- as.numeric(good[, 1:15])
best_so <- rater(kn, best_fit[1], best_fit[2], best_fit[5], best_fit[6], best_fit[7], best_fit[8], best_fit[9], good$type)
pred <- backfill(best_so, 10000, best_fit[3], best_fit[4], .1, 3)
# pred <- fit_volumes1(best_so, 10000, best_fit[3], best_fit[4], 10)
pred[[2]]
pred <- pred[[1]]
# png('vol_fit_10-01.png', height = 15, width = 17, units = 'cm', res = 300)
plot(kn$ToMouth_km, kn$xsec_area, pch = 20, col = get_palette('forest'),
     xlab = 'distance to outlet (km)', ylab = 'unit cross-sectional area (km2)')
points(kn$ToMouth_km, pred$vol, pch = 20, col = get_palette('charcoal'))
legend('topright', legend = c('observed', 'fit'),
       fill = get_palette(c('forest', 'charcoal'), .7))
# dev.off()

# png('xsec_fit_10-01.png', height = 15, width = 17, units = 'cm', res = 300)
plot(best_so$xsec_area %>% emp_cdf, pch = 20, col = get_palette('ocean', .2),
     xlab = 'unit cross-sectional area m2', ylab = 'CDF')
lines(pred$vol %>% emp_cdf, lwd = 3, col = get_palette('charcoal', .6))
legend('bottomright', legend = c('observed', 'fit'),
       fill = get_palette(c('ocean', 'charcoal'), .7))
# dev.off()

# png('lvl_fit_10-01.png', height = 15, width = 17, units = 'cm', res = 300)
plot(pred$lvl %>% emp_cdf, type = 'l', lwd = 3, col = get_palette('charcoal', .7),
     xlab = 'average valley depth m2', ylab = 'CDF')
points(best_so$lvl %>% emp_cdf, pch = 20, col = get_palette('ocean', .2))
legend('bottomright', legend = c('observed', 'fit'),
       fill = get_palette(c('ocean', 'charcoal'), .7))
# dev.off()

# png('elev_fit_10-01.png', height = 15, width = 17, units = 'cm', res = 300)
plot(pred$elev %>% emp_cdf, type = 'l', lwd = 3, col = get_palette('charcoal', .7),
     xlab = 'elevation m', ylab = 'CDF')
points(best_so$elev %>% emp_cdf, pch = 20, col = get_palette('ocean', .2))
legend('bottomright', legend = c('observed', 'fit'),
       fill = get_palette(c('ocean', 'charcoal'), .7))
# dev.off()

setwd('/home/crumplecup/work/')
save(rec, file = 'n_dim_9-29a.rds')
load('n_dim_9-29a.rds')

#png('fit_9-22.png', height = 15, width = 17, units = 'cm', res = 300)
scatter3D(rec[,3], rec[,4], (rec[,9] + rec[,11])/2,
          ticktype = 'detailed', pch = 20, theta = 320, phi = 40,
          xlab = 'input rate', ylab = 'output rate', zlab = 'test value')
#dev.off()

# png('fit_9-29.png', height = 15, width = 17, units = 'cm', res = 300)
sub <- rec[rec[,3] > .1 & rec[,4] > .1, ]
scatter3D(sub[,3], sub[,4], (sub[,10] + sub[,12])/2,
          ticktype = 'detailed', pch = 20, theta = 320, phi = 40,
          xlab = 'input rate', ylab = 'output rate', zlab = 'test value')
# dev.off()

# png('fit_9-28a.png', height = 15, width = 17, units = 'cm', res = 300)
sub <- rec[rec[,9] < .2 & rec[,11] < .15, ]
scatter3D(sub[,3], sub[,4], (sub[,9] + sub[,11])/2,
          ticktype = 'detailed', pch = 20, theta = 320, phi = 40,
          xlab = 'input rate', ylab = 'output rate', zlab = 'test value')
# dev.off()


br_tribs <- c(386219, 386239, 386243, 386258, 386269, 386278, 386299, 386314,
              386331, 386348, 386359, 386378, 386389, 386409, 386440, 386452)
kn_tribs <- landmarks$nodeid
flags <- array(0, nrow(crks))
flags[crks$NODE_ID %in% c(br_tribs, kn_tribs)] <- 1
idx <- 1:nrow(crks) * flags
idx <- idx[idx > 0]

plot(crks$dp[crks$NODE_ID %in% c(br_tribs, kn_tribs,
                                 br_tribs + 1, kn_tribs + 1,
                                 br_tribs + 2, kn_tribs + 2,
                                 br_tribs + 3, kn_tribs + 3,
                                 br_tribs + 4, kn_tribs + 4)],
     ylab = 'vals')

plot(crks$dp, pch = 20, col = get_palette('charcoal'))
points(idx, crks$dp[crks$NODE_ID %in% c(br_tribs, kn_tribs)], pch = 20, col = get_palette('crimson'))

plot(kn_dif)
plot(log(best_so$contr_area) + log(best_so$slope))
plot(1/best_so$slope)

slp_dif <- 0
for (i in 2:nrow(best_so)) {
  slp_dif[i] <- best_so$slope[i] - best_so$slope[i-1]
}
slp_dif <- -slp_dif
slp_dif <- slp_dif - min(slp_dif)

kn_dif <- lca_dif[crks$creek_name == 'knowles']
kn_dif <- kn_dif - min(kn_dif)
plot(slp_dif)
plot(kn_dif)
plot(log(slp_dif))
slp_dif <- log(slp_dif+.01) - min(log(slp_dif+.01))
wt_dif <- kn_dif * slp_dif
plot(wt_dif)

plot(wt_mth)


bit <- voluminous2(best_so, 10000, best_fit[3], best_fit[4])


# inspect weight candidates for input
library(data.table)
kn <- creeks[creeks$creek_name == 'knowles', ]
best_so <- rater(kn, best_fit[1], best_fit[2], best_fit[5], best_fit[6], best_fit[7], best_fit[8], best_fit[9], good$type)

plot(best_so$lvl, type = 'l', lwd = 3, col = get_palette('hardwood', .7))

plot(best_so$xsec_area, type = 'l', lwd = 3, col = get_palette('hardwood', .7))
lines(best_so$lvl * 50, lwd = 3, col = get_palette('forest', .7))

plot(log(best_so$dp) - min(log(best_so$dp)))
lines(best_so$lvl, lwd = 3, col = get_palette('forest', .7))

lines(best_so$contr_area)
lines(best_so$xsec_area / 50, lwd = 3, col = get_palette('hardwood', .7))

pwr_df <- data.frame(CA = log(best_so$contr_area),
                     slope = log(best_so$slope))
pwr_mod <- lm(CA ~ slope, data = pwr_df)
pwr_prd <- predict(pwr_mod, newdata = data.frame(slope = log(best_so$slope)))
pwr_k <- log(best_so$contr_area) - pwr_prd

ti_ave <- best_fit[1] / nrow(sa_nodes)
to_ave <- best_fit[2] / nrow(sa_nodes)

ca_df <- data.frame(vol = best_so$xsec_area,
                    area = best_so$contr_area,
                    ca2 = best_so$contr_area^2,
                    dp = best_so$dp,
                    pwr = pwr_k - min(pwr_k))
ca_mod <- glm(vol ~ area + ca2 + dp + pwr, data = ca_df)
summary(ca_mod)
ca_prd <- predict(ca_mod, newdata = data.frame(ca_df))
plot(ca_prd)

plot(log(1 / (log(best_so$dp) - min(log(best_so$dp))) ))

iwt <- 1 / exp(pwr_k) *
  (log(best_so$dp) - min(log(best_so$dp))) *
  best_so$contr_area / max(best_so$contr_area)

plot(iwt, type = 'l', lwd = 2, col = get_palette('charcoal', .7))
lines(best_so$xsec_area / 40, lwd = 3, col = get_palette('hardwood', .7))


plot(1/exp(pwr_k), ylim = c(0, 4))
lines(best_so$xsec_area / 150, lwd = 3, col = get_palette('hardwood', .7))

lines(pwr_k - min(pwr_k)*10)

lelev <- log(best_so$elev+1)
llvl <- log(best_so$elev+1 + best_so$lvl)
ldist <- log(best_so$ToMouth_km)
plot(ldist, lelev)
celev <- lelev[-c(1:80)]
cdist <- ldist[-c(1:80)]
plot(cdist, celev)
lm(celev ~ cdist) %>% summary
plot(best_so$elev, type = 'l')
lines(best_so$elev + best_so$lvl)


br <- creeks[creeks$creek_name == 'bear', ]
best_so <- rater(br, best_fit[1], best_fit[2], best_fit[5], best_fit[6], best_fit[7], best_fit[8], best_fit[9], good$type)
pred <- backfill1(best_so, 10000, best_fit[3], best_fit[4], .1, .5)

lelev <- log(best_so$elev+1)
llvl <- log(best_so$elev+1 + best_so$lvl)
ldist <- log(best_so$ToMouth_km)
plot(ldist, lelev)
celev <- lelev[-c(1:80)]
cdist <- ldist[-c(1:80)]
plot(cdist, celev, type = 'l')
lm(celev ~ cdist) %>% summary
plot(best_so$elev)

plot(best_so$elev, type = 'l')
lines(best_so$elev + best_so$lvl)

begin <- Sys.time()
diff <- as.difftime('24:00:00', '%H:%M:%S', units = 'hour')
end <- begin + diff
i <- 0
while (Sys.time() < end) {
  res <- ndim_search(
    kn,
    c(.0674, .0674), c(.0658, .0658), c(.05, .2), c(.05, .5),
    tp = c(0.25, 0.25),
    lp = c(0.5, 0.5),
    it = 100,
    batch = 50
  )
  rec <- rbind(rec, res)
  save(rec, file = "link_20201114.rds")
}

# load("vol_20201008.rds")
load("link_20201114.rds")
# rec <- res
nrow(rec)
(good <- rec[rec$ks_lvl + rec$ks_vol == min(rec$ks_lvl + rec$ks_vol), ])
best_fit <- as.numeric(good)
best_so <- rater1(kn, best_fit[1], best_fit[2])


pred <- backfill1(best_so, 10000, best_fit[3], best_fit[4], best_fit[5], best_fit[6])
# pred <- fit_volumes1(best_so, 10000, best_fit[3], best_fit[4], 10)
pred[[2]]
pred <- pred[[1]]

# setwd('/home/crumplecup/work')
png('vol_fit_11-03.png', height = 15, width = 17, units = 'cm', res = 300)
plot(kn$ToMouth_km, kn$xsec_area, type = 'l', lwd = 2.5, col = get_palette('forest', .7),
     xlab = 'distance to outlet (km)', ylab = 'unit cross-sectional area (km2)')
lines(kn$ToMouth_km, rev(pred$vol), lwd = 2, col = get_palette('charcoal', .7))
legend('topright', legend = c('observed', 'fit'),
       fill = get_palette(c('forest', 'charcoal'), .7))
# dev.off()

# png('xsec_fit_11-03.png', height = 15, width = 17, units = 'cm', res = 300)
plot(best_so$xsec_area %>% emp_cdf, pch = 20, col = get_palette('ocean', .2),
     xlab = 'unit cross-sectional area m2', ylab = 'CDF')
lines(pred$vol %>% emp_cdf, lwd = 3, col = get_palette('charcoal', .6))
legend('bottomright', legend = c('observed', 'fit'),
       fill = get_palette(c('ocean', 'charcoal'), .7))
# dev.off()

# png('lvl_11-03.png', height = 15, width = 17, units = 'cm', res = 300)
plot(pred$lvl %>% emp_cdf, type = 'l', lwd = 3, col = get_palette('charcoal', .7),
     xlab = 'average valley depth m2', ylab = 'CDF', xlim = c(0, max(best_so$lvl)))
points(best_so$lvl %>% emp_cdf, pch = 20, col = get_palette('ocean', .2))
legend('bottomright', legend = c('observed', 'fit'),
       fill = get_palette(c('ocean', 'charcoal'), .7))
# dev.off()

# png('fit_11-14.png', height = 15, width = 17, units = 'cm', res = 300)
scatter3D(rec$vi, rec$vo, (rec$ks_vol + rec$ks_lvl)/2,
          ticktype = 'detailed', pch = 20, theta = 335, phi = 30,
          xlab = 'input rate', ylab = 'output rate', zlab = 'test value')
# dev.off()


pal <- get_palette(c('ocean', 'forest'), .7)
plot(creeks$valley_width, pch = 20, col = pal[1])
points(creeks$VALWIDTH_5, pch = 20, col = pal[2])

setwd('/home/crumplecup/work/')
png('creeks_valley_width.png', height = 24, width = 24, units = 'cm', res = 300)
par(mfrow = c(2, 2))
br <- creeks[creeks$creek_name == 'bear', ]
plot(br$ToMouth_km, br$valley_width, pch = 20, col = pal[1],
     xlab = 'distance from mouth (km)', ylab = 'valley width (m)',
     ylim = c(0, max(br$VALWIDTH_5)), main = 'bear')
points(br$ToMouth_km, br$VALWIDTH_5, pch = 20, col = pal[2])
legend('topleft', legend = c('survey', '10-m DEM'), fill = pal)

ce <- creeks[creeks$creek_name == 'cedar', ]
plot(ce$ToMouth_km, ce$valley_width, pch = 20, col = pal[1],
     xlab = 'distance from mouth (km)', ylab = 'valley width (m)',
     ylim = c(0, max(ce$VALWIDTH_5)), main = 'cedar')
points(ce$ToMouth_km, ce$VALWIDTH_5, pch = 20, col = pal[2])
legend('topright', legend = c('survey', '10-m DEM'), fill = pal)

ho <- creeks[creeks$creek_name == 'hoffman', ]
plot(ho$ToMouth_km, ho$valley_width, pch = 20, col = pal[1],
     xlab = 'distance from mouth (km)', ylab = 'valley width (m)',
     main = 'hoffman')
points(ho$ToMouth_km, ho$VALWIDTH_5, pch = 20, col = pal[2])
legend('bottomleft', legend = c('survey', '10-m DEM'), fill = pal)

kn <- creeks[creeks$creek_name == 'knowles', ]
plot(kn$ToMouth_km, kn$valley_width, pch = 20, col = pal[1],
     xlab = 'distance from mouth (km)', ylab = 'valley width (m)',
     ylim = c(0, max(kn$VALWIDTH_5)), main = 'knowles')
points(kn$ToMouth_km, kn$VALWIDTH_5, pch = 20, col = pal[2])
legend('topright', legend = c('survey', '10-m DEM'), fill = pal)
dev.off()



png('creeks_slope.png', height = 24, width = 24, units = 'cm', res = 300)
par(mfrow = c(2, 2))
plot(br$ToMouth_km, br$slope, pch = 20, col = pal[1],
     xlab = 'distance from mouth (km)', ylab = 'slope',
     #ylim = c(0, max(br$GRADIENT)),
     main = 'bear')
points(br$ToMouth_km, br$GRADIENT, pch = 20, col = pal[2])
legend('topleft', legend = c('survey', '10-m DEM'), fill = pal)

plot(ce$ToMouth_km, ce$slope, pch = 20, col = pal[1],
     xlab = 'distance from mouth (km)', ylab = 'slope',
     ylim = c(min(ce$slope), max(ce$GRADIENT)), main = 'cedar')
points(ce$ToMouth_km, ce$GRADIENT, pch = 20, col = pal[2])
legend('topleft', legend = c('survey', '10-m DEM'), fill = pal)

plot(ho$ToMouth_km, ho$slope, pch = 20, col = pal[1],
     xlab = 'distance from mouth (km)', ylab = 'slope',
     ylim = c(min(ho$GRADIENT), max(ho$GRADIENT)),
     main = 'hoffman')
points(ho$ToMouth_km, ho$GRADIENT, pch = 20, col = pal[2])
legend('topleft', legend = c('survey', '10-m DEM'), fill = pal)

plot(kn$ToMouth_km, kn$slope, pch = 20, col = pal[1],
     xlab = 'distance from mouth (km)', ylab = 'slope',
     ylim = c(0, max(kn$GRADIENT)),
     main = 'knowles')
points(kn$ToMouth_km, kn$GRADIENT, pch = 20, col = pal[2])
legend('topleft', legend = c('survey', '10-m DEM'), fill = pal)
dev.off()

cols <- creeks$creek_name
cols[cols == 'bear'] <- get_palette('coral')
cols[cols == 'cedar'] <- get_palette('gold')
cols[cols == 'hoffman'] <- get_palette('forest')
cols[cols == 'knowles'] <- get_palette('ocean')

creeks$slope[creeks$slope <= 0] <- .001
pwr_df <- data.frame(CA = log(creeks$contr_area),
                     slope = log(creeks$slope))
pwr_mod <- lm(CA ~ slope, data = pwr_df)
pwr_prd <- predict(pwr_mod, newdata = data.frame(slope = log(creeks$slope)))
pwr_k <- log(so$contr_area) - pwr_prd

so <- rater1(creeks, .11, .08)
wtp <- log(so$dp) - min(log(so$dp))
ca_wt <- so$contr_area / max(so$contr_area)
wt_i <- wtp * ca_wt / exp(pwr_k)

png('creeks_weights.png', height = 36, width = 24, units = 'cm', res = 300)
par(mfrow = c(4, 1))
plot(creeks$xsec_area, pch = 20, col = cols,
     ylab = 'cross-sectional area (m2)')
legend('topleft', legend = c('bear', 'cedar', 'hoffman', 'knowles'),
       fill = get_palette(c('coral', 'gold', 'forest', 'ocean'), .7))

plot(exp(pwr_k), pch = 20, col = cols, ylab = 'transport capacity')
plot(so$dp, pch = 20, col = cols, ylab = 'optimized delivery probability')
plot(wtp, pch = 20, col = cols, ylab = 'wt_p')
dev.off()

png('creeks_weights1.png', height = 36, width = 24, units = 'cm', res = 300)
par(mfrow = c(3, 1))
plot(so$lvl, pch = 20, col = cols,
     ylab = 'average valley depth (m)')
legend('topleft', legend = c('bear', 'cedar', 'hoffman', 'knowles'),
       fill = get_palette(c('coral', 'gold', 'forest', 'ocean'), .7))
plot(creeks$xsec_area, pch = 20, col = cols,
     ylab = 'cross-sectional area (m2)')
legend('topleft', legend = c('bear', 'cedar', 'hoffman', 'knowles'),
       fill = get_palette(c('coral', 'gold', 'forest', 'ocean'), .7))
plot(wt_i, ylim = c(0, 8), pch = 20, col = cols, ylab = 'wt_i')
legend('topleft', legend = c('bear', 'cedar', 'hoffman', 'knowles'),
       fill = get_palette(c('coral', 'gold', 'forest', 'ocean'), .7))
dev.off()

