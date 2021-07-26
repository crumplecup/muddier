
library(magrittr)
library(sp)
library(muddier)

# screen out older samples from sites with multiples
# node ids of sites with charcoal samples
site <- creeks_radio$node_ids
# count number of samples at each site
counts <- 0
for (i in seq_along(site)) {
  counts[i] <- length(site[site == site[i]])
}
# site ids with multiple samples
mults <- site[counts > 1] %>% factor %>% levels
# subset creek nodes with single samples only
crks <- creeks_radio[!creeks_radio$node_ids %in% mults, ]
# loop through sites with multiples and add the youngest from each
for (i in seq_along(mults)) {
  sub <- creeks_radio[creeks_radio$node_ids == mults[i], ]
  crks <- rbind(crks, sub[1,])
}

# table of multiple sites
site_ct <- charcoal[, .N, by = family]
char <- charcoal
char$ct <- 0
for (i in 1:nrow(char)) {
  char$ct[i] <- site_ct$N[site_ct$family == char$family[i]]
}
site_ct[site_ct$N > 2, ]
png('site_count.png', height = 17, width = 21, units = 'cm', res = 300)
plot(emp_cdf(char$ct), xlim = c(1, max(char$ct)), ylim = c(0, 1), type = 'l',
     col = get_palette('ocean', .7), lwd = 2,
     xlab = 'samples per site', ylab = 'CDF of samples')
points(emp_cdf(char$ct), pch = 20, col = get_palette('charcoal', .5))
dev.off()

site_ct$N[site_ct$N > 2] %>% sum
site_ct$N[site_ct$N > 1] %>% sum

index <- char_pmfs %>% rownames %>% as.numeric %>% rev
png('inherited_age.png', height = 17, width = 21, units = 'cm', res = 300)
plot(index, df_cdf, log = 'x', type = 'l', lwd = 2, col = get_palette('charcoal', .7),
     xlab = 'inherited age [years]', ylab = 'CDF of samples', ylim = c(0.1,1))
lines(index, fg_cdf, lwd = 2, col = get_palette('ocean', .7))
lines(index, ff_cdf, lwd = 2, col = get_palette('crimson', .7))
legend('topleft', legend = c('debris flows', 'fines', 'gravels'),
       fill = get_palette(c('charcoal', 'crimson', 'ocean'), .9))
dev.off()

# remove samples from the cedar fan
crks <- crks[crks$creek_name != 'cedar', ]


# add optimized MB delivery probability
# add MB delivery probability
crks_so <- rater1(creeks, 0.7, 0.7)
crks$dp <- 0
crks$odp <- 0

for (i in seq_along(crks$dp)) {
  crks$dp[i] <- creeks$DebrisFlow[creeks$NODE_ID == crks$node_ids[i]]
  crks$odp[i] <- crks_so$dp[crks_so$NODE_ID == crks$node_ids[i]]
}

# subset debris-flow deposits
df_crks <- crks[crks$stratigraphy == 'DF', ]
# cdf_crks <- to_cdf(crks$dp)

# debris-flow deposit density
rho_df <- nrow(df_crks) / nrow(crks)

# takes pmf, returns cdf binned
pct_cdf <- function(vec, by = vec, ln = 100) {
  vec <- obs_to_cdf(normalize(vec), by = by)
  n <- length(vec)
  cdf <- vector(length = ln, mode = 'numeric')
  for (i in 1:ln) {
    cdf[i] <- length(vec[vec <= i/ln]) / n
  }
  return(cdf)
}

# takes pmf, returns cdf
obs_to_cdf <- function(pmf, by = pmf)  {
  cdf <- array(0,length(pmf))
  if (sum(pmf) != 1) pmf <- normalize(pmf)
  pmf <- pmf[order(by)]
  for (i in seq_along(pmf)) {
    if (i == 1) cdf[i] <- pmf[i]
    if (i > 1)  cdf[i] <- pmf[i] + cdf[i-1]
  }
  cdf
}

# cdf of debris-flow deposits
dfcdf <- obs_to_cdf(df_crks$dp)
# cdf of total deposits
crkcdf <- obs_to_cdf(crks$dp)
# match cdf length of total deposits to cdf of debris-flows
pct_crkcdf <- pct_cdf(crks$dp, ln = length(dfcdf))

# log both cdfs to fit to model
lg_cdf_df <- log(dfcdf)
lg_cdf_crks <- log(pct_crkcdf)

lmod <- lm(log(dfcdf) ~ log(pct_crkcdf))
lpred <- predict(lmod, data.frame(pct_crkcdf = crkcdf))
epred <- exp(lpred)

edif <- cdf_dif(epred)
ddif <- cdf_dif(crkcdf)
fx <- edif / ddif * (1/rho_df)

dp_vals <- sort(crks$dp)
plot(dp_vals, fx)

# same for optimized delivery probability
# cdf of debris-flow deposits
odfcdf <- obs_to_cdf(df_crks$odp)

# match cdf length of total deposits to cdf of debris-flows
opct_crkcdf <- pct_cdf(crks$odp, ln = length(odfcdf))

olmod <- lm(log(odfcdf) ~ log(opct_crkcdf))
olpred <- predict(olmod, data.frame(opct_crkcdf = crkcdf))
oepred <- exp(olpred)

oedif <- cdf_dif(oepred)
ofx <- oedif / ddif * (1/rho_df)

odp_vals <- sort(crks$odp)
plot(odp_vals, ofx)


df <- data.frame(df = lg_cdf_df, crk = lg_cdf_crks)

mod1 <- lm(df ~ crk, data = df[1:35,])
mod2 <- lm(df ~ crk, data = df[30:39,])
mod3 <- lm(df ~ crk, data = df[36:45,])

pred1 <- exp(predict(mod1, newdata = data.frame(crk = log(cdf_crks))))
pred2 <- exp(predict(mod2, newdata = data.frame(crk = log(cdf_crks))))
pred3 <- exp(predict(mod3, newdata = data.frame(crk = log(cdf_crks))))

# cdf_ids <- ragged_bin(1:length(cdf_crks), 3)
pred <- c(pred1[1:78], pred2[79:85], pred3[86:89])

# save data
setwd('/home/crumplecup/work/muddier/')

df_wt <- data.frame(wt = fx, dp = dp_vals)
usethis::use_data(df_wt, overwrite = T)

dfcdf_prd <- data.frame(dfcdf_prd = pred, cdf_crks = cdf_crks)
usethis::use_data(dfcdf_prd, overwrite = T)


setwd('/home/crumplecup/work/')
png('wt_dp1.png', height = 17, width = 23, units = 'cm', res = 300)
plot(emp_cdf(crks$dp)[,1], fx, type = 'l', lwd = 3, col = get_palette('ocean', .7),
     xlab = 'MB Delivery Probability Index', ylab = 'Relative Density of Debris Flows',
     log = 'x')
points(emp_cdf(crks$dp)[,1], fx, col = get_palette('charcoal', .7))
dev.off()

s <- seq(0, 3.5, 0.1)
lines(s, s * (0.035 / 3.5))
lines(df_wt)

s <- seq(0, 1, 0.1)
setwd('/media/erik/catacomb/research')
png('delivery_wgt.png', height = 20, width = 24, units = 'cm', res = 300)
plot(s, s, type = 'l', lwd = 2, lty = 1, col = get_palette('charcoal'),
     xlab = 'CDF of total deposits',
     ylab = 'CDF of debris-flow deposits')
lines(pct_crkcdf, dfcdf, lwd = 1, lty = 2, col = get_palette('ocean'))
points(pct_crkcdf, dfcdf, pch = 20, col = get_palette('ocean'))
lines(opct_crkcdf, odfcdf, lwd = 1, lty = 2, col = get_palette('crimson'))
points(opct_crkcdf, odfcdf, pch = 20, col = get_palette('crimson'))
legend('topleft', legend = c('Unweighted', 'MB predictor', 'Optimized predictor'),
       fill = get_palette(c('charcoal', 'ocean', 'crimson'), .8))
dev.off()

png('weighting_function_df.png', height = 20, width = 24, units = 'cm', res = 300)
plot(df_wt$dp, df_wt$wt, type = 'l', lwd = 3, col = get_palette('charcoal'),
     xlab = 'delivery probability',
     ylab = 'weighting function', log = 'x')
points(df_wt$dp, df_wt$wt, pch = 20, col = get_palette('ocean'))
legend('topleft', legend = c('observed', 'fit'),
       fill = get_palette(c('ocean', 'charcoal'), .8))
dev.off()


# Miller/Burnett model output
library(rgdal)
nodes <- readOGR('nodes_debrisflow.shp')

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


# develop color ramp to map probabilities to colors for plotting
ramp <- colorRampPalette(c(get_palette('ocean', .75),
                           get_palette('crimson', .75)), alpha = T)


col_ramp <- ramp(255)

sub <- creeks
sub <- sub[order(sub$DebrisFlow), ]
sub_n <- nrow(sub)
sub$p_cdf <- (1:sub_n) / sub_n

dp_crk <- val_by_cdf(crks$dp, ln = nrow(crks))


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

sub$wt_cols <- 0
sub$wt_cols[sub$wt <= 0.1] <- get_palette('ocean', .075)
sub$wt_cols[sub$wt > 0.1 & sub$wt <= 1] <- get_palette('sky', .075)
sub$wt_cols[sub$wt > 1 & sub$wt <= 2] <- get_palette('gold', .075)
sub$wt_cols[sub$wt > 2 & sub$wt <= 3] <- get_palette('rose', .075)
sub$wt_cols[sub$wt > 3] <- get_palette('crimson', .075)

rgdal::writeOGR(sub, 'creeks.shp', 'creeks', driver = 'ESRI Shapefile')

library(raster)
png('kn_wt.png',
    width = 7, height = 5, units = 'in', res = 300)
frame <- pad_ext(extent(transects), c(1.25,.15))
pic <- crop(nodes, frame)
pts <- crop(sub[sub$creek_name == 'knowles',], frame)
plot(pic, pch = 20, col = get_palette('slate', .05))
points(pts, col = pts$wt_cols, pch = 20)
legend('bottomleft', legend = c('wt < .01', '.01 to 1', '1 to 2', '2 to 3', '3+'),
       fill = c(get_palette(c('ocean', 'sky', 'gold', 'rose', 'crimson'), .7)))
dev.off()

png('br_wt.png',
    width = 7, height = 5, units = 'in', res = 300)
frame <- pad_ext(extent(bear), c(1.25,.15))
pic <- crop(nodes, frame)
pts <- crop(sub, frame)
plot(pic, pch = 20, col = get_palette('slate', .05))
points(pts, col = pts$wt_cols, pch = 20)
legend('bottomleft', legend = c('wt < .01', '.01 to 1', '1 to 2', '2 to 3', '3+'),
       fill = c(get_palette(c('ocean', 'sky', 'gold', 'rose', 'crimson'), .7)))
dev.off()


hi <- .018
mid <- .0102
low <- .006
df_low <- df_crks[df_crks$dp < low, ]
crk_low <- crks[crks$dp < low, ]
df_mid <- df_crks[df_crks$dp >= low & df_crks$dp < mid, ]
crk_mid <- crks[crks$dp >= low & crks$dp < mid, ]
df_hi <- df_crks[df_crks$dp >= mid & df_crks$dp < hi, ]
crk_hi <- crks[crks$dp >= mid & crks$dp < hi, ]
tot_df <- nrow(df_low) + nrow(df_mid) + nrow(df_hi)
tot_crk <- nrow(crk_low) + nrow(crk_mid) + nrow(crk_hi)

png('weighting_eg.png', height = 12, width = 18, units = 'cm', res = 300)
barplot(c(nrow(df_low) / tot_df, nrow(crk_low) / tot_crk,
          nrow(df_mid) / tot_df, nrow(crk_mid) / tot_crk,
          nrow(df_hi) / tot_df, nrow(crk_hi) / tot_crk),
        col = pal[rep(2:1,3)],
        space = c(.5,.1),
        ylab = 'proportion of samples',
        names.arg = c('        p =', 'low       ',
                      '        p =', 'mid        ',
                      '        p =', 'hi        '))
legend('topleft', legend = c('DF', 'ALL'), fill = pal[2:1])
dev.off()

dfn <- c(nrow(df_low) / tot_df, nrow(df_mid) / tot_df, nrow(df_hi) / tot_df)
crkn <- c(nrow(crk_low) / tot_crk, nrow(crk_mid) / tot_crk, nrow(crk_hi) / tot_crk)

df_ave <- (nrow(df_hi) + nrow(df_mid) + nrow(df_hi)) / (nrow(crk_low) + nrow(crk_mid) + nrow(crk_hi))
df_inv <- df_ave^-1

dfns <- c(nrow(df_low), nrow(df_mid), nrow(df_hi))
crkns <- c(nrow(crk_low), nrow(crk_mid), nrow(crk_hi))
difs <- dfns / crkns
wt_p <- df_inv * difs
prd <- wt_p * df_ave * c(27, 32, 34)
prd

# ia_sv <- ia_rn
# ia_sv2 <- ia_rn
ia_rn <- rexp(73, 73/8095)
ia_df <- dfmns[-1] - dfmns[-73]
setwd('/home/crumplecup/work/')
png('ia_eg1.png', height = 16, width = 19, units = 'cm', res = 300)
plot(emp_cdf(ia_df), col = get_palette('charcoal', .7),
     xlab = 'interarrival time of deposit',
     ylab = 'CDF')
lines(emp_cdf(ia_df), col = get_palette('charcoal', .7))
lines(emp_cdf(ia_rn), col = get_palette('crimson', .7), lwd = 2)
lines(emp_cdf(ia_sv), col = get_palette('gold', .7), lwd = 2)
legend('topleft', legend = c('observed', 'random', 'cherry-picked'),
       fill = get_palette(c('charcoal', 'crimson', 'gold'), .7))
dev.off()

ia_syn <- cumsum(ia_rn)
ia_syn2 <- cumsum(ia_sv)

png('ia_eg2.png', height = 16, width = 19, units = 'cm', res = 300)
plot(emp_cdf(dfmns), col = get_palette('charcoal', .7),
     xlab = 'sample age', ylab = 'CDF')
lines(emp_cdf(dfmns), col = get_palette('charcoal', .7))
lines(emp_cdf(ia_syn), col = get_palette('crimson', .7), lwd = 3)
lines(emp_cdf(ia_syn2), col = get_palette('gold', .7), lwd = 3)
legend('bottomright', legend = c('observed', 'random', 'cherry-picked'),
       fill = get_palette(c('charcoal', 'crimson', 'gold'), .7))
dev.off()

png('ia_eg3.png', height = 16, width = 19, units = 'cm', res = 300)
plot(dfmns[-1], ia_df, pch = 20, col = get_palette('charcoal'),
     xlab = 'sample age', ylab = 'interarrival time')
points(cumsum(ia_sv), ia_sv, pch = 20, col = get_palette('gold'))
points(cumsum(ia_rn), ia_rn, pch = 20, col = get_palette('crimson'))
legend('topleft', legend = c('observed', 'random', 'cherry-picked'),
       fill = get_palette(c('charcoal', 'crimson', 'gold'), .7))
dev.off()

rt <- 73/8095

png('ia_eg4.png', height = 16, width = 19, units = 'cm', res = 300)
plot(emp_cdf(dfmns), col = get_palette('charcoal', .7),
     xlab = 'sample age', ylab = 'CDF')
lines(emp_cdf(dfmns), col = get_palette('charcoal', .7))
lines(emp_cdf(ia_syn), col = get_palette('crimson', .7), lwd = 2)
for (i in 1:5) {
  lines(emp_cdf(boot_rec(100, 73, (i+1)*rt, i*rt, 0)), col = get_palette('crimson', .7), lwd = 2)
}
text(7500, .81, 'in:out = 1:0')
text(5500, .86, '2:1')
text(4200, .86, '3:2')
text(3000, .85, '4:3')
text(2500, .85, '5:4')
text(2000, .85, '6:5')
legend('bottomright', legend = c('observed', 'synthetic'),
       fill = get_palette(c('charcoal', 'crimson'), .7))
dev.off()

boot_rec <- function(n, x, rin, rout, ia) {
  ar <- lapply(1:n, function(m) recorder(x, rin, rout, ia)) %>% rack
  apply(ar, 1, median)
}


png('ia_eg5.png', height = 16, width = 19, units = 'cm', res = 300)
plot(emp_cdf(ia_ar), col = get_palette('gold', .7),
     xlab = 't', ylab = 'CDF', log = 'x')
points(emp_cdf(dfmns), col = get_palette('charcoal', .7))
lines(emp_cdf(dfmns), col = get_palette('charcoal', .7))
lines(emp_cdf(ia_ar), col = get_palette('gold', .7))
legend('topleft', legend = c('deposit age', 'inherited age'),
       fill = get_palette(c('charcoal', 'gold'), .7))
dev.off()

plot(ksmin$ri, ksmin$ro, xlim = c(min(ksmin$ri), max(kpmin$ri)),
     ylim = c(min(ksmin$ro), max(kpmin$ro)),
     col = get_palette('ocean'), pch = 20)
points(kpmin$ri, kpmin$ro, col = get_palette('crimson'), pch = 20)


boot_restime <- function(n, tab, x, ia = 0) {
  res <- vector(length = nrow(tab), mode = 'numeric')
  for (i in 1:nrow(tab)) {
    res[i] <- boot_rec(n, x, tab[i, 1], tab[i, 2], ia) %>% mean
  }
  res
}

res_kp <- boot_restime(100, kpmin, 73)
res_ks <- boot_restime(100, ksmin, 73)


png('restime.png', height = 16, width = 19, units = 'cm', res = 300)
plot(emp_cdf(res_ks), pch = 20, col = get_palette('ocean'),
     xlim = c(min(res_kp), max(res_ks)),
     xlab = 'mean residence time (yrs)', ylab = 'CDF')
lines(emp_cdf(res_ks), lwd = 2, col = get_palette('ocean'))
points(emp_cdf(res_kp), pch = 20, col = get_palette('crimson'))
lines(emp_cdf(res_kp), lwd = 2, col = get_palette('crimson'))
legend('bottomright', legend = c('K-S fit', 'Kuiper fit'),
       fill = get_palette(c('ocean', 'crimson'), .7))
dev.off()




