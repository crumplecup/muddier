
library(magrittr)
library(sp)
library(muddier)

# here i scrub NAs and bad values that need to be checked / redone
sub <- creeks[!is.na(creeks$slope) & !(creeks$NODE_ID %in% cedar$NODE_ID[105:156]),]

# subset by attribute ranges
sub <- sub[sub$slope < .2 &
             sub$contr_area < 6 &
             sub$slope > 0 &
             sub$valley_width < 49 &
             sub$valley_width > 0 &
             sub$xsec_area > 0, ]

names(sub)
nrow(sub)


dp <- sub$DebrisFlow  # debris-flow delivery
xs <- sub$xsec_area  # cross-sectional area m^2
vw <- sub$valley_width  # valley width
ca <- sub$contr_area  # contributing area
sl <- sub$slope  # slope

lslope <- (sub$slope + .026) %>% log # log slope with offset for negative grads
lcontr <- sub$contr_area %>% log # log contr area
as_mod <- lm(lslope ~ lcontr)  # fit AS to power law
as_theta <- as_mod$coefficients[2]  # pull theta value

# estimate individual ks values using pulled theta
log_ks <- lslope - as_theta * lcontr
ks <- exp(log_ks)

# accumulation index options
ndp <- dp / max(dp) # normalized p
nks <- ks / max(ks)  # normalized k_s
acc_dif <- ndp - nks  # p - k_s accumulation difference
acc_rat <- ndp / nks  # p / k_s accumulation ratio


# gather variables into data.frame for linear modeling
dat <- data.frame(xsec_area = xs,  # 1
                  log_xsec = log(xs+.01),  # 2
                  val_width = vw, # 3
                  log_val_width = log(vw),  # 4
                  del_prob = dp,  # 5
                  slope = sl, #  6
                  contr_area = ca,  #  7
                  strm_pwr = ks,  #  8
                  acc_dif = acc_dif,  #  9
                  acc_rat = acc_rat)  #  10


modlist <- list(dat[ , c(1, 8)], dat[ , c(1, 5, 8)],
                dat[ , c(1, 9)], dat[ , c(1, 10)])

(tab <- stat_tab(modlist))

modlist <- list(dat[ , c(1, 3)],
                dat[ , c(1, 3, 5)],
                dat[ , c(1, 3, 6)],
                dat[ , c(1, 3, 7)],
                dat[ , c(1, 3, 8)],
                dat[ , c(1, 3, 9)],
                dat[ , c(1, 3, 10)],
                dat[ , c(2, 4)],
                dat[ , c(2, 4, 5)],
                dat[ , c(2, 4, 6)],
                dat[ , c(2, 4, 7)],
                dat[ , c(2, 4, 8)],
                dat[ , c(2, 4, 9)],
                dat[ , c(2, 4, 10)])

(tab <- stat_tab(modlist))

modlist <- list(dat[ , c(2, 4, 7)],
                dat[ , c(2, 4, 7, 5)],
                dat[ , c(2, 4, 7, 6)],
                dat[ , c(2, 4, 7, 5, 6)],
                dat[ , c(2, 4, 5, 6)])

(tab <- stat_tab(modlist))
mod_stat(dat[ , c(2, 4, 7, 5, 6)])

mod <- lm(log_xsec ~ log_val_width + contr_area + del_prob + slope, data = dat)

pred <- exp(predict(mod, dat))
sub$pred <- pred

min_dat <- dat
min_dat$del_prob <- min(dat$del_prob)
max_dat <- dat
max_dat$del_prob <- max(dat$del_prob)

min_pred <- exp(predict(mod, min_dat))
max_pred <- exp(predict(mod, max_dat))
difs <- max_pred - min_pred
quantile(difs)

wt <- weight_by(dat$xsec_area, pred, pred)
sub <- sub[order(pred), ]
sub$wt <- wt
#op <- (wt * sort(dat$del_prob)) / sum(wt * sort(dat$del_prob))
#op <- op / max(op)
#np <- sort(dat$del_prob) / max(dat$del_prob)

site <- creeks_radio$node_ids
counts <- 0
for (i in seq_along(site)) {
  counts[i] <- length(site[site == site[i]])
}
mults <- site[counts > 1] %>% factor %>% levels

crks <- creeks_radio[!creeks_radio$node_ids %in% mults, ]

for (i in seq_along(mults)) {
  sub <- creeks_radio[creeks_radio$node_ids == mults[i], ]
  crks <- rbind(crks, sub[1,])
}


crks$dp <- 0
for (i in seq_along(crks$dp)) {
  crks$dp[i] <- creeks$DebrisFlow[creeks$NODE_ID == crks$node_ids[i]]
}

df_crks <- crks[crks$Stratigraphy == 'DF', ]

pct_cdf_df <- percent_cdf(df_crks$dp, ln = length(crks$dp))
pct_cdf_crks <- percent_cdf(df_crks$dp, ln = nrow(df_crks))

cdf_df <- to_cdf(df_crks$dp)
cdf_crks <- to_cdf(crks$dp)


rho_df <- nrow(df_crks) / nrow(crks)

lg_cdf_df <- log(cdf_df)
lg_cdf_crks <- log(pct_cdf_crks)

lmod <- lm(lg_cdf_df ~ lg_cdf_crks)
lpred <- predict(lmod, data.frame(lg_cdf_crks = log(cdf_crks)))
epred <- exp(lpred)

edif <- cdf_dif(epred)
ddif <- cdf_dif(cdf_crks)
fx <- edif / ddif * (1/rho_df)

dp_vals <- sort(crks$dp)

df <- data.frame(df = lg_cdf_df, crk = lg_cdf_crks)

mod1 <- lm(df ~ crk, data = df[1:18,])
mod2 <- lm(df ~ crk, data = df[19:35,])
mod3 <- lm(df ~ crk, data = df[36:46,])

pred1 <- exp(predict(mod1, newdata = data.frame(crk = log(cdf_crks))))
pred2 <- exp(predict(mod2, newdata = data.frame(crk = log(cdf_crks))))
pred3 <- exp(predict(mod3, newdata = data.frame(crk = log(cdf_crks))))

# cdf_ids <- ragged_bin(1:length(cdf_crks), 3)
pred <- c(pred2[1:94], pred3[95:98])

edif <- cdf_dif(pred)
ddif <- cdf_dif(cdf_crks)
fx <- edif / ddif * (1/rho_df)








#older code

creeks_radio$dp <- 0
for (i in seq_along(creeks_radio$dp)) {
  creeks_radio$dp[i] <- creeks$DebrisFlow[creeks$NODE_ID == creeks_radio$node_ids[i]]
}

df_crks <- creeks_radio[creeks_radio$Stratigraphy == 'DF', ]

pct_cdf_df <- percent_cdf(df_crks$dp)
pct_cdf_crks <- percent_cdf(creeks_radio$dp)

cdf_df <- to_cdf(df_crks$dp)
cdf_crks <- to_cdf(creeks_radio$dp)

plot(sort(creeks_radio$dp), cdf_crks, type = 'l', lwd=2, col = 'slateblue')
lines(sort(df_crks$dp), cdf_df, lwd = 2, col = 'forestgreen')

plot(pct_cdf_crks, pct_cdf_df)
rho_df <- nrow(df_crks) / nrow(creeks_radio)

lg_cdf_df <- log(pct_cdf_df)
lg_cdf_crks <- log(pct_cdf_crks)

lmod <- lm(lg_cdf_df ~ lg_cdf_crks)
lpred <- predict(lmod, data.frame(lg_cdf_crks))
epred <- exp(lpred)



df_p <- creeks$DebrisFlow[creeks$NODE_ID %in% df_crks$node_ids]
crks_p <- creeks$DebrisFlow[creeks$NODE_ID %in% creeks_radio$node_ids]
df_cdf <- to_cdf(normalize(df_p))
crks_cdf <- to_cdf(normalize(crks_p))


creeks_radio$df <- 0
creeks_radio$df[creeks_radio$Stratigraphy == 'DF'] <- 1


creeks_radio$vw <- 0
creeks_radio$xs <- 0
creeks_radio$sl <- 0

for (i in 1:nrow(creeks_radio)) {
  creeks_radio$vw[i] <- creeks$valley_width[creeks$NODE_ID == creeks_radio$node_ids[i]]
  creeks_radio$xs[i] <- creeks$xsec_area[creeks$NODE_ID == creeks_radio$node_ids[i]]
  creeks_radio$sl[i] <- creeks$slope[creeks$NODE_ID == creeks_radio$node_ids[i]]
}
creeks_radio$lvw <- log(creeks_radio$vw)
creeks_radio$lxs <- log(creeks_radio$xs)

d <- data.frame(df = creeks_radio$df, dp = creeks_radio$dp,
                vw = creeks_radio$vw, lvw = creeks_radio$lvw,
                xs = creeks_radio$xs, lxs = creeks_radio$lxs,
                sl = creeks_radio$sl)
d <- d[order(d$dp), ]
d$ln <- 1:nrow(d)

dfcdf <- percent_cdf(d$df[d$df == 1], d$dp[d$df == 1])
ldfcdf <- log(dfcdf)
dpcdf <- percent_cdf(d$dp)
ldpcdf <- log(dpcdf)

lmod <- lm(ldfcdf ~ ldpcdf)
lpred <- predict(lmod, data.frame(ldpcdf))
epred <- exp(lpred)


edif <- cdf_dif(epred)
ddif <- cdf_dif(pct_cdf_crks)
fx <- edif / ddif * (1/rho_df)

dp_vals <- val_by_cdf(creeks_radio)


vec <- seq(.01, 1, .01)
plot(vec, vec, type = 'l', lwd = 2, lty = 2,
     xlab = 'proportion total area', ylab = 'proportion debris flow area')
lines(percent_cdf(crks_p), percent_cdf(df_p), type = 'l', lwd = 2, col = 'blue')
lines(percent_cdf(crks_p), percent_cdf(df_df$pred), type = 'l', lwd = 2, col = 'blue')

plot(percent_cdf(df_p), percent_cdf(df_p)/ percent_cdf(crks_p),
     type = 'l', lwd = 2, col = 'blue')


ragged_bin(length(df_p), 10)
ragged_bin(100)

(f_x <- bin_cdfs(log(df_p+1), log(crks_p+1), log(crks_p+1), rho_df, 10))
#(f_x <- bin_cdfs(df_p, crks_p, crks_p, rho_df, 4))
rng <- f_x[[4]]
# rng <- range(crks_p)
rng <- exp(rng) - 1
ddif <- f_x[[1]]


cdf_dif(cumsum(f_x[[2]])/sum(f_x[[2]])) / cdf_dif(cumsum(f_x[[3]])/sum(f_x[[3]]))

rho_df * f_x[[1]] * f_x[[3]]
f_x[[2]]



plot(rng[-1], ddif, xlim = c(0, max(rng)),
     xlab = 'delivery probability')

plot(cumsum(f_x[[3]]) / sum(f_x[[3]]), cumsum(f_x[[2]]) / sum(f_x[[2]]),
     xlab = 'proportion creek area', ylab = 'proportion debris flows',
     type = 'l', lwd = 2, col = get_palette('ocean'))
lines(vec, vec, lwd = 2, lty = 2)
points(cumsum(f_x[[3]]) / sum(f_x[[3]]), cumsum(f_x[[2]]) / sum(f_x[[2]]))

plot(vec, vec, type = 'l', lwd = 2, lty = 2,
     xlab = 'proportion total area', ylab = 'proportion debris flow area')
lines(percent_cdf(crks_n), percent_cdf(df_n, crks_n), type = 'l', lwd = 2, col = 'blue')

wt <- bin_cdfs(df_p, crks_p, crks_p, rho = rho_df)

plot(percent_cdf(df_p), (1 / rho_df) * (percent_cdf(df_p) / percent_cdf(crks_p)),
     type = 'l', lwd = 2, col = 'blue')

df_dat <- data.frame(df = percent_cdf(df_p), area = percent_cdf(crks_p))
mod <- lm(df ~ area, data = df_dat)
summary(mod)
pred_df <- predict(mod, newdata = df_dat)


logmod <- glm(df ~ dp + sl, family = binomial, data = d)
summary(logmod)
dfpred <- predict(logmod, newdata = d, type = 'response')
dfbin <- dfpred
dfbin[dfpred <= .5] <- 0
dfbin[dfpred > .5] <- 1
d$prd <- dfbin

d$ndp <- d$dp / max(d$dp)
ct <- seq(.1, 1, .1)

logit_stat(d[c(1, 2, 7)])

pal <- get_palette(c('ocean', 'forest', 'coral'))
plot(to_cdf(d$dp), to_cdf(d$dp), type = 'l', lwd = 2, col = pal[1])
plot(percent_cdf(d$dp),
      percent_cdf(d$df[d$df == 1], d$dp[d$df == 1]),
      lwd = 2, col = pal[2], type = 'l',
     xlab = 'proportion creek area', ylab = 'proportion debris flows')
lines(percent_cdf(d$dp),
      percent_cdf(d$df[d$df == 1], d$dp[d$df == 1]),
      lwd = 2, col = pal[3])


vdp <- val_by_cdf(d$dp)

descdist(df_dat[,1] / df_dat[,2], boot = 100)
df_gam <- fitdistrplus::fitdist(as.numeric(df_dat[,1] / df_dat[,2]), 'gamma')
df_nm <- fitdistrplus::fitdist(as.numeric(df_dat[,1] / df_dat[,2]), 'norm')
df_lnm <- fitdistrplus::fitdist(as.numeric(df_dat[,1] / df_dat[,2]), 'lnorm')
df_wb <- fitdistrplus::fitdist(as.numeric(df_dat[,1] / df_dat[,2]), 'weibull')

descdist(dfcdf / dpcdf, boot = 100)
df_gam <- fitdistrplus::fitdist(as.numeric(dfcdf / dpcdf), 'gamma')
df_nm <- fitdistrplus::fitdist(as.numeric(dfcdf / dpcdf), 'norm')
df_lnm <- fitdistrplus::fitdist(as.numeric(dfcdf / dpcdf), 'lnorm')
df_wb <- fitdistrplus::fitdist(as.numeric(dfcdf / dpcdf), 'weibull')


leg <- c('gamma', 'normal', 'lognormal', 'weibull')
denscomp(list(df_gam, df_nm, df_lnm), legendtext = leg)
qqcomp(list(df_gam, df_nm, df_lnm), legendtext = leg)
cdfcomp(list(df_gam, df_nm, df_lnm), legendtext = leg)
ppcomp(list(df_gam, df_nm, df_lnm), legendtext = leg)

plot(dfit)

library(fitdistrplus)

#subset debris flow deposit ages
df <- as.matrix(char_pmfs)
df <- df[ , charcoal$facies == 'DF']

#determine which sites have multiple samples
site <- charcoal$family[charcoal$facies == 'DF']
counts <- 0
for (i in seq_along(site)) {
  counts[i] <- length(site[site == site[i]])
}
mults <- site[counts > 1] %>% factor %>% levels

#select youngest age from each set of multiples
dfmn <- charcoal$mn[charcoal$facies == 'DF']
mins <- vector(length(site), mode = 'numeric')
for (i in seq_along(mults)) {
  mins[dfmn == min(dfmn[site == mults[i]])] <- 1
}

#screen out all but youngest from multiples
df <- df[ , counts == 1 | (counts > 1 & mins > 0)]

#order from youngest to oldest
dfmn <- dfmn[counts == 1 | (counts > 1 & mins > 0)]
df <- df[ ,order(dfmn)]
dfmns <- dfmn + 50

bin_mns <- discrete_bin(dfmn+50, 500)

# subset bins with sufficient counts to fit an exp distribution
# sm_bins <- list(bin_mns[[1]], bin_mns[[2]], bin_mns[[3]],
#                bin_mns[[4]], bin_mns[[5]], bin_mns[[6]])
sm_bins <- list(bin_mns[[1]], bin_mns[[2]], bin_mns[[3]])
df_exp <- fit_bins(df, dfmns, sm_bins)
df_gam <- fit_bins(df, dfmns, sm_bins, 'gamma')
df_wbl <- fit_bins(df, dfmns, sm_bins, 'weibull')

exp_est <- t(df_exp[[3]])

# xlab <- c('0-200', '201-400', '401-600', '601-800', '801-1000', '1001-1200')
xlab <- c('0-500', '501-1000', '1001-1500')


(gof <- gof_tab(list(df_exp[[2]], df_gam[[2]], df_wbl[[2]]),
                c('exp', 'gamma', 'weibull')))


# fit exponential distribution to interarrival times
descdist(as.numeric(df_exp[[1]][[3]]))
plot(sort(df_exp[[1]][[3]]))

labs <- c('exp', 'gamma', 'weibull')
denscomp(ft = list(df_exp[[2]][[1]], df_gam[[2]][[1]], df_wbl[[2]][[1]]), legendtext = labs)
cdfcomp(ft = list(df_exp[[2]][[1]], df_gam[[2]][[1]], df_wbl[[2]][[1]]), legendtext = labs)
# qqcomp(ft = list(df_exp, df_gam, df_wbl), legendtext = labs)
# ppcomp(ft = list(df_exp, df_gam, df_wbl), legendtext = labs)


library(sp)
library(rgdal)
library(raster)
library(data.table)
# Lancaster survey data
setwd('/home/crumplecup/work')
kn_vol <- fread('kn_vol.csv')
br_vol <- fread('bear.csv')
cd_vol <- fread('cedar.csv')
hf_vol <- fread('hoffman_sedvol.csv')
gr_vol <- read.csv('grc_vol.csv', skipNul = TRUE)

# contributing area per study area
br_ca <- (max(br_vol$corr_contr_area_m2) - min(br_vol$corr_contr_area_m2)) / 1e6
cd_ca <- (max(cd_vol$corr_contr_area_m2) - min(cd_vol$corr_contr_area_m2)) / 1e6
gr_ca <- max(gr_vol$contr_area_km2) - min(gr_vol$contr_area_km2)
hf_ca <- (max(hf_vol$corr_contr_area_m2) - min(hf_vol$corr_contr_area_m2)) / 1e6
kn_ca <- max(kn_vol$contr_area_km2) - min(kn_vol$contr_area_km2)
cas <- c(br_ca, cd_ca, gr_ca, hf_ca, kn_ca)

# stream length per study area
br_ln <- max(br_vol$outlet_dist_m)
cd_ln <- max(cd_vol$corr_hipchain_m)
gr_ln <- max(gr_vol$hipchain_m) - min(gr_vol$hipchain_m)
hf_ln <- max(hf_vol$corr_hipchain_m)
kn_ln <- max(kn_vol$corrected_hipchain) - min(kn_vol$corrected_hipchain)
lns <- c(br_ln, cd_ln, gr_ln, hf_ln, kn_ln)

# nodes per study area
nds <- c(nrow(bear), nrow(cedar), nrow(hoffman), nrow(transects))

nd_ln <- lns[-3] / nds
nd_ln <- c(nd_ln[1:2], mean(nd_ln), nd_ln[3:4])  # use mean for grc

# total contributing area represented by nodes with charcoal sampled debris flows
nd_ca <- sum(cas / lns * nd_ln * ncol(df))

# arrival rate per square kilometer of contributing area
ar_ca <- df_exp[[3]] / nd_ca

# develop color ramp to map probabilities to colors for plotting
ramp <- colorRampPalette(c(get_palette(c('slate', .1)),
                           get_palette('rose', .5)), alpha = T)
rose_ramp <- colorRampPalette(c(get_palette('white', .01),
                                get_palette('rose')), alpha = T)
rose_bush <- colorRampPalette(c(get_palette('hardwood', .01),
                                get_palette('rose')), alpha = T)
forest_ramp <- colorRampPalette(c(get_palette('white', .01),
                                  get_palette('forest')), alpha = T)
tree_ramp <- colorRampPalette(c(get_palette('hardwood', .01),
                                  get_palette('forest')), alpha = T)
blood_ramp <- colorRampPalette(c(get_palette('rose', .1),
                                get_palette('crimson', .5)), alpha = T)
col_ramp <- ramp(255)
rose_ramp <- rose_ramp(255)
rose_bush <- rose_bush(255)
forest_ramp <- forest_ramp(255)
tree_ramp <- tree_ramp(255)
blood_ramp <- blood_ramp(255)

sub <- creeks
sub <- sub[order(sub$DebrisFlow), ]
sub$p_cdf <- (1:sub_n) / sub_n
sub_n <- nrow(sub)

sub$wt <- 0
for(i in seq_along(dp_crk)) {
  if (i == 1) {
    sub$wt[sub$DebrisFlow <= dp_crk[i]] <- fx[i]
  }
  if (i > 1 & i < length(dp_crk)) {
    sub$wt[sub$DebrisFlow > dp_crk[i-1] &
             sub$DebrisFlow <= dp_crk[i]] <- fx[i]
  }
  if (i == length(dp_crk)) {
    sub$wt[sub$DebrisFlow > dp_crk[i-1]] <- fx[i]
  }
}

scl <- 255
sub$dp_cols <- rose_bush[round((sub$DebrisFlow / max(sub$DebrisFlow)) * 255)]
sub$wt_cols <- 0
sub$wt_cols[sub$wt > 1] <- forest_ramp[round((sub$wt[sub$wt > 1] /
                                                max(sub$wt[sub$wt > 1])) * scl)]
wt_inv <- round(((1 / sub$wt) / max(1 / sub$wt)) * scl)
for (i in seq_along(wt_inv)) {
  wt_inv[i] <- blood_ramp[wt_inv[i]]
}
sub$wt_cols[sub$wt <= 1] <- wt_inv[sub$wt <= 1]









sub$op <- (sub$DebrisFlow * sub$wt) / sum(sub$DebrisFlow * sub$wt)
sub$np <- sub$DebrisFlow / sum(sub$DebrisFlow)

sub$op_cols <- rose_bush[round((sub$op / max(sub$op)) * 255)]

np_cols <- round(np * 255)
op_cols <- round(op * 255)
sub <- sub[order(sub$DebrisFlow),]
sub$np_cols <- col_ramp[np_cols]
sub$op_cols <- col_ramp[op_cols]
sub$op <- op
sub$dif <- op - np
sub$dif_neg <- NA
sub$dif_neg[sub$dif < 0] <- (abs(sub$dif[sub$dif < 0]) / abs(min(sub$dif))) * 255
sub$dif_pos <- NA
sub$dif_pos[sub$dif > 0] <- (sub$dif[sub$dif > 0] / max(sub$dif)) * 255
sub$dif_neg[sub$dif < 0] <- rose_ramp[round(sub$dif_neg[sub$dif < 0])]
sub$dif_pos[sub$dif > 0] <- forest_ramp[round(sub$dif_pos[sub$dif > 0])]

# Miller/Burnett model output
nodes <- readOGR('nodes_debrisflow.shp')

frame <- extent(transects)

pad_ext <- function(ext, scale = c(.15, .1), h = TRUE, v = TRUE) {
  if (length(scale) == 1) scale <- c(scale, scale)
  if (h) {
    adj <- (ext[2] - ext[1]) * scale[1]
    ext[1] <- round(ext[1] - adj)
    ext[2] <- round(ext[2] + adj)
  }
  if (v) {
    adj <- (ext[4] - ext[3]) * scale[2]
    ext[3] <- round(ext[3] - adj)
    ext[4] <- round(ext[4] + adj)
  }
  ext
}





# plots for export
setwd('/home/crumplecup/work')
png('cor_dp_xs.png',
    width = 7, height = 5, units = 'in', res = 300)
cor_creeks(xs, dp,
           creeks = sub$creek_name,
           creek_labs = levels(factor(sub$creek_name)),
           labs = c('cross-sectional area km2', 'delivery probability'),
           cols = c('rose', 'forest', 'ocean', 'violet'),
           leg_pos = 'topright', a = .25)
dev.off()


png('cor_dp_vw.png',
    width = 7, height = 5, units = 'in', res = 300)
cor_creeks(vw, dp,
           creeks = sub$creek_name,
           creek_labs = levels(factor(sub$creek_name)),
           labs = c('valley width m', 'delivery probability'),
           cols = c('rose', 'forest', 'ocean', 'violet'),
           leg_pos = 'topright', a = .25)
dev.off()

png('cor_dp_ca.png',
    width = 7, height = 5, units = 'in', res = 300)
cor_creeks(ca, dp,
           creeks = sub$creek_name,
           creek_labs = levels(factor(sub$creek_name)),
           labs = c('contributing area km2', 'delivery probability'),
           cols = c('rose', 'forest', 'ocean', 'violet'),
           leg_pos = 'topright', a = .25)
dev.off()

png('cor_dp_sl.png',
    width = 7, height = 5, units = 'in', res = 300)
cor_creeks(sl, dp,
           creeks = sub$creek_name,
           creek_labs = levels(factor(sub$creek_name)),
           labs = c('slope', 'delivery probability'),
           cols = c('rose', 'forest', 'ocean', 'violet'),
           leg_pos = 'topright', a = .25)
dev.off()

png('cor_xs_ca.png',
    width = 7, height = 5, units = 'in', res = 300)
cor_creeks(ca, xs,
           creeks = sub$creek_name,
           creek_labs = levels(factor(sub$creek_name)),
           labs = c('contributing area km2',
                    'cross-sectional area m2'),
           cols = c('rose', 'forest', 'ocean', 'violet'),
           leg_pos = 'topright', a = .25)
dev.off()

png('cor_xs_sl.png',
    width = 7, height = 5, units = 'in', res = 300)
cor_creeks(sl, xs,
           creeks = sub$creek_name,
           creek_labs = levels(factor(sub$creek_name)),
           labs = c('slope', 'cross-sectional area m2'),
           cols = c('rose', 'forest', 'ocean', 'violet'),
           leg_pos = 'topright', a = .25)
dev.off()

png('cor_xs_vw.png',
    width = 7, height = 5, units = 'in', res = 300)
cor_creeks(vw, xs,
           creeks = sub$creek_name,
           creek_labs = levels(factor(sub$creek_name)),
           labs = c('valley width', 'cross-sectional area m2'),
           cols = c('rose', 'forest', 'ocean', 'violet'),
           leg_pos = 'topright', a = .25)
dev.off()

png('cor_xs_ks.png',
    width = 7, height = 5, units = 'in', res = 300)
cor_creeks(ks, xs,
           creeks = sub$creek_name,
           creek_labs = levels(factor(sub$creek_name)),
           labs = c('streampower coefficient',
                    'cross-sectional area m2'),
           cols = c('rose', 'forest', 'ocean', 'violet'),
           leg_pos = 'topright', a = .25)
dev.off()

png('cor_xs_acc_dif.png',
    width = 7, height = 5, units = 'in', res = 300)
cor_creeks(acc_dif, xs,
           creeks = sub$creek_name,
           creek_labs = levels(factor(sub$creek_name)),
           labs = c('p - k_s',
                    'cross-sectional area m2'),
           cols = c('rose', 'forest', 'ocean', 'violet'),
           leg_pos = 'topright', a = .25)
dev.off()

png('cor_xs_acc_rat.png',
    width = 7, height = 5, units = 'in', res = 300)
cor_creeks(acc_rat, xs,
           creeks = sub$creek_name,
           creek_labs = levels(factor(sub$creek_name)),
           labs = c('p / k_s',
                    'cross-sectional area m2'),
           cols = c('rose', 'forest', 'ocean', 'violet'),
           leg_pos = 'topright', a = .25)
dev.off()

png('pred_xs.png',
    width = 7, height = 5, units = 'in', res = 300)
plot(max_pred[order(pred)],
     type = 'l', lwd = 1.5, col = 'slateblue',
     ylab = 'cross-sectional area m2')
lines(abs(dat$xsec_area[order(pred)] - pred[order(pred)]),
      lwd = 1.5, col = 'gray')
lines(min_pred[order(pred)], lwd = 1.5, col = 'slateblue')
lines(pred[order(pred)], col = 'goldenrod', lwd = 2)
lines(difs[order(pred)], lwd = 1.5, col = 'coral3')
legend('topleft',
       legend = c('predicted', 'min & max del_prob',
                  'max - min del_prob', '|observed - predicted|'),
       fill = c('goldenrod', 'slateblue', 'coral3', 'gray'))
dev.off()

png('wt_xs.png',
    width = 7, height = 5, units = 'in', res = 300)
plot_weight(num = dat$xsec_area,
            den = pred,
            by = dat$del_prob,
            bin_no = 10,
            labs = c('delivery probability', 'weight'))
dev.off()

chg <- (op[order(np)] - np[order(np)]) / sort(np)

png('op_chg.png',
    width = 7, height = 5, units = 'in', res = 300)
plot(np, op - np, pch = 20, col = get_palette('hardwood', .2),
     xlab = 'original delivery probability',
     ylab = 'optimal - original')
abline(h = 0, lwd = 2, lty = 2)
dev.off()

png('op_np.png',
    width = 7, height = 5, units = 'in', res = 300)
plot(np, pch = 20, col = get_palette('hardwood', .2),
     ylab = 'delivery probability')
points(op, pch = 20, col = get_palette('rose', .2))
legend('topleft', legend = c('original', 'optimized'),
       fill = c(get_palette('hardwood'), get_palette('crimson')))
dev.off()

png('op_cdf.png',
    width = 7, height = 5, units = 'in', res = 300)
plot(sort(sub$op), to_cdf(sort(sub$op)), type = 'l', lwd = 3, col = get_palette('rose'),
     xlab = 'delivery probability', ylab = 'CDF F(x)')
lines(sort(sub$np), to_cdf(sort(sub$np)),
      lwd = 3, col = get_palette('hardwood'))
legend('bottomright', legend = c('original', 'optimized'),
       fill = c(get_palette('hardwood'), get_palette('rose')))
dev.off()

png('delivery_probs_adjusted.png',
    width = 7, height = 5, units = 'in', res = 300)
plot(ord$p_ks, ord$pO, pch = 20, log = 'x',
     col = rgb(red=.267, green=.671, blue=.471, alpha=.33),
     ylab = 'delivery probability', xlab = 'p / k_s',
     ylim = c(min(ord$pO - ord$DebrisFlow/sum(ord$DebrisFlow)), max(ord$pO)))
points(ord$p_ks, ord$DebrisFlow/sum(ord$DebrisFlow),
       pch = 20, col = rgb(red=.306, green=.569, blue=.831, alpha=.33))
points(ord$p_ks, ord$pO - ord$DebrisFlow/sum(ord$DebrisFlow),
       pch = 20, col = rgb(red=1, green=.592, blue=.478, alpha=.33))
legend('bottomright', legend = c('p', 'pO', 'pO-p'),
       fill = c(rgb(red=.306, green=.569, blue=.831, alpha=.33),
                rgb(red=.267, green=.671, blue=.471, alpha=.33),
                rgb(red=1, green=.592, blue=.478, alpha=.33)))
dev.off()


png('interarrival_rates_df.png',
    width = 7, height = 5, units = 'in', res = 300)
plot(ar_ca[,1], type = 'l', lwd = 2, col = get_palette('ocean', .5),
     xlab = 'age range of local sample window', ylab = 'interarrival rate',
     ylim = c(0, max(ar_ca[,4])), axes=F)
lines(ar_ca[,2], lwd = 2, col = get_palette('forest', .8), lty = 2)
lines(ar_ca[,3], lwd = 2, col = get_palette('rose', .8), lty = 2)
lines(ar_ca[,4], lwd = 2, col = get_palette('rose', .8), lty = 2)
axis(1, at=1:length(xlab), labels=xlab)
axis(2)
legend('topleft', legend = c('Mean', 'Median', '95% CI'),
       fill = c(get_palette('ocean', .5),
                get_palette('forest', .8),
                get_palette('rose', .8)))
dev.off()

sub$dp_cols[sub$DebrisFlow <= 0.009] <- get_palette('crimson', .075)
sub$dp_cols[sub$DebrisFlow > 0.009 & sub$DebrisFlow <= 0.013] <- get_palette('rose', .075)
sub$dp_cols[sub$DebrisFlow > 0.013 & sub$DebrisFlow <= 0.022] <- get_palette('gold', .075)
sub$dp_cols[sub$DebrisFlow > 0.022 & sub$DebrisFlow <= 0.031] <- get_palette('leaf', .075)
sub$dp_cols[sub$DebrisFlow > 0.031] <- get_palette('forest', .075)


png('kn_p.png',
    width = 7, height = 5, units = 'in', res = 300)
frame <- pad_ext(extent(transects), c(1.25,.15))
pic <- crop(nodes, frame)
pts <- crop(sub, frame)
plot(pic, pch = 20, col = get_palette('ocean', .05))
points(pts, col = pts$dp_cols, pch = 20)
legend('bottomleft', legend = c('p < .009', '.01 to .013', '.013 to .022', '.022 to .031', '.031+'),
       fill = c(get_palette(c('crimson', 'rose', 'gold', 'leaf', 'forest'), .7)))
dev.off()


sub$wt_cols[sub$wt <= 0.1] <- get_palette('crimson', .075)
sub$wt_cols[sub$wt > 0.1 & sub$wt <= 1] <- get_palette('rose', .075)
sub$wt_cols[sub$wt > 1 & sub$wt <= 2] <- get_palette('gold', .075)
sub$wt_cols[sub$wt > 2 & sub$wt <= 3] <- get_palette('leaf', .075)
sub$wt_cols[sub$wt > 3] <- get_palette('forest', .075)

png('kn_wt.png',
    width = 7, height = 5, units = 'in', res = 300)
frame <- pad_ext(extent(transects), c(1.25,.15))
pic <- crop(nodes, frame)
pts <- crop(sub[sub$creek_name == 'knowles',], frame)
plot(pic, pch = 20, col = get_palette('ocean', .05))
points(pts, col = pts$wt_cols, pch = 20)
legend('bottomleft', legend = c('wt < .01', '.01 to 1', '1 to 2', '2 to 3', '3+'),
       fill = c(get_palette(c('crimson', 'rose', 'gold', 'leaf', 'forest'), .7)))
dev.off()

png('kn_dif.png',
    width = 7, height = 5, units = 'in', res = 300)
frame <- pad_ext(extent(transects), c(1.25,.15))
pic <- crop(nodes, frame)
pts <- crop(sub, frame)
plot(pic, pch = 20, col = get_palette('ocean', .05))
points(pts, col = pts$dif_pos, pch = 20)
points(pts, col = pts$dif_neg, pch = 20)
dev.off()

png('br_p.png',
    width = 7, height = 5, units = 'in', res = 300)
frame <- pad_ext(extent(bear), c(1.25,.15))
pic <- crop(nodes, frame)
pts <- crop(sub, frame)
plot(pic, pch = 20, col = get_palette('ocean', .05))
points(pts, col = pts$dp_cols, pch = 20)
legend('bottomleft', legend = c('hi p', 'low p'),
       fill = c(get_palette('rose', .8), get_palette('hardwood', .8)))
dev.off()

png('br_wt.png',
    width = 7, height = 5, units = 'in', res = 300)
frame <- pad_ext(extent(bear), c(1.25,.15))
pic <- crop(nodes, frame)
pts <- crop(sub, frame)
plot(pic, pch = 20, col = get_palette('ocean', .05))
points(pts, col = pts$wt_cols, pch = 20)
legend('bottomleft', legend = c('wt < .01', '.01 to 1', '1 to 2', '2 to 3', '3+'),
       fill = c(get_palette(c('crimson', 'rose', 'gold', 'leaf', 'forest'), .7)))
dev.off()

png('br_dif.png',
    width = 7, height = 5, units = 'in', res = 300)
frame <- pad_ext(extent(bear), c(1.25,.15))
pic <- crop(nodes, frame)
pts <- crop(sub, frame)
plot(pic, pch = 20, col = get_palette('ocean', .05))
points(pts, col = pts$dif_pos, pch = 20)
points(pts, col = pts$dif_neg, pch = 20)
dev.off()

png('cd_p.png',
    width = 7, height = 5, units = 'in', res = 300)
frame <- pad_ext(extent(cedar), c(1.25,.15))
pic <- crop(nodes, frame)
pts <- crop(sub, frame)
plot(pic, pch = 20, col = get_palette('ocean', .05))
points(pts, col = pts$dp_cols, pch = 20)
legend('bottomleft', legend = c('hi p', 'low p'),
       fill = c(get_palette('rose', .8), get_palette('hardwood', .8)))
dev.off()

png('cd_wt.png',
    width = 7, height = 5, units = 'in', res = 300)
frame <- pad_ext(extent(cedar), c(1.25,.15))
pic <- crop(nodes, frame)
pts <- crop(sub, frame)
plot(pic, pch = 20, col = get_palette('ocean', .05))
points(pts, col = pts$wt_cols, pch = 20)
legend('bottomleft', legend = c('wt < .01', '.01 to 1', '1 to 2', '2 to 3', '3+'),
       fill = c(get_palette(c('crimson', 'rose', 'gold', 'leaf', 'forest'), .7)))
dev.off()

png('cd_dif.png',
    width = 7, height = 5, units = 'in', res = 300)
frame <- pad_ext(extent(cedar), c(1.25,.15))
pic <- crop(nodes, frame)
pts <- crop(sub, frame)
plot(pic, pch = 20, col = get_palette('ocean', .05))
points(pts, col = pts$dif_pos, pch = 20)
points(pts, col = pts$dif_neg, pch = 20)
dev.off()

png('hf_p.png',
    width = 7, height = 5, units = 'in', res = 300)
frame <- pad_ext(extent(hoffman), c(1.25,.15))
pic <- crop(nodes, frame)
pts <- crop(sub, frame)
plot(pic, pch = 20, col = get_palette('ocean', .05))
points(pts, col = pts$dp_cols, pch = 20)
legend('bottomleft', legend = c('hi p', 'low p'),
       fill = c(get_palette('rose', .8), get_palette('hardwood', .8)))
dev.off()

png('hf_wt.png',
    width = 7, height = 5, units = 'in', res = 300)
frame <- pad_ext(extent(hoffman), c(1.25,.15))
pic <- crop(nodes, frame)
pts <- crop(sub, frame)
plot(pic, pch = 20, col = get_palette('ocean', .05))
points(pts, col = pts$wt_cols, pch = 20)
legend('bottomleft', legend = c('wt < .01', '.01 to 1', '1 to 2', '2 to 3', '3+'),
       fill = c(get_palette(c('crimson', 'rose', 'gold', 'leaf', 'forest'), .7)))
dev.off()

png('hf_dif.png',
    width = 7, height = 5, units = 'in', res = 300)
frame <- pad_ext(extent(hoffman), c(.6,.15))
pic <- crop(nodes, frame)
pts <- crop(sub, frame)
plot(pic, pch = 20, col = get_palette('ocean', .05))
points(pts, col = pts$dif_pos, pch = 20)
points(pts, col = pts$dif_neg, pch = 20)
dev.off()

png('pred_wt.png',
    width = 7, height = 5, units = 'in', res = 300)
plot(sub$pred, sub$wt, col = get_palette('white', .001), pch = 20,
     xlab = 'cross-sectional area m2',
     ylab = 'weighting function w')
abline(h = 1, lwd = 2, lty = 2, col = get_palette('charcoal', .9))
points(sub$xsec_area, sub$wt, col = get_palette('charcoal', .15), pch = 20)
points(sub$pred, sub$wt, col = get_palette('ocean', .15), pch = 20)
legend('topright', legend = c('predicted', 'observed'),
       fill = c(get_palette('ocean'), get_palette('charcoal')))
dev.off()


png('wt_dp.png',
    width = 7, height = 5, units = 'in', res = 300)
pal <- get_palette(c('ocean', 'charcoal'))
dp_crk <- val_by_cdf(crks$dp, ln = nrow(crks))
plot(dp_crk, fx, pch = 20, col = pal[2], lwd = 2, type = 'l',
     xlab = 'delivery probability', ylab = 'weighting function')
points(dp_crk, fx, col = pal[1])
legend('bottomright', legend = c('fit', 'bins'), fill = pal[2:1])
dev.off()

png('wt_fit.png',
    width = 7, height = 5, units = 'in', res = 300)
vec <- seq(.01, 1, .01)
plot(vec, vec, type = 'l', lwd = 2, lty = 2, col = get_palette('slate'),
     xlab = 'proportion total samples', ylab = 'proportion debris flows')
lines(cdf_crks, pred, lwd = 2, col = get_palette('charcoal'))
points(pct_cdf_crks, cdf_df, pch = 20, col = get_palette('ocean'))
legend('topleft', legend = c('observed', 'fit'),
       fill = get_palette(c('ocean', 'charcoal')[c(2,1)]))
dev.off()

png('df_a_cdfs.png',
    width = 7, height = 5, units = 'in', res = 300)
pal <- get_palette(c('ocean', 'forest', 'coral'), .5)
dp_vals = val_by_cdf(crks$dp, ln = nrow(df_crks))
plot(sort(dp_vals), pct_cdf_crks, type = 'l', col = pal[1], lwd = 2, ylim = c(0,1),
     xlab = 'delivery probability p', ylab = 'F(p)')
lines(sort(dp_vals), cdf_df, col = pal[2], lwd = 2)
legend('bottomright', legend = c('debris flows', 'total samples'),
       fill = pal[2:1])
dev.off()














