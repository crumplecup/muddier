library(magrittr)

rec <- data.table::fread('/home/erik/output/df20k_1000.csv')
rec1 <- data.table::fread('/home/erik/output/df20k_1001.csv')
recf <- data.table::fread('/home/erik/output/ff20k_1000.csv')
recf1 <- data.table::fread('/home/erik/output/ff20k_1001.csv')
recf2 <- data.table::fread('/home/erik/output/ff20k_1002.csv')
recg <- data.table::fread('/home/erik/output/fg20k_1000.csv')
recg1 <- data.table::fread('/home/erik/output/fg20k_1001.csv')
rec <- rbind(rec, rec1)
recf <- rbind(recf, rbind(recf1, recf2))
recg <- rbind(recg, recg1)
rm(rec1, recf1, recf2, recg1)

trans <- data.table::fread('/home/erik/output/transit_fits.csv')

# total sediment storage in the study creeks
sedvol <- data.table::fread('/home/erik/data/sedvol.csv')
sfill <- sum(sedvol$`Volume increment (m3)`[sedvol$Reach == 'Cedar']) +
  sum(sedvol$`Volume increment (m3)`[sedvol$Reach == 'Hoffman']) +
  sum(sedvol$`Volume increment (m3)`[sedvol$Reach %in% c('LB', 'UB')]) +
  sum(sedvol$`Volume increment (m3)`[sedvol$Reach %in% c('DFK', 'LK', 'UK')])
# contributing area for study creeks
ein <- max(creeks$contr_area[creeks$creek_name == 'bear']) +
  max(creeks$contr_area[creeks$creek_name == 'cedar']) +
  max(creeks$contr_area[creeks$creek_name == 'hoffman']) +
  max(creeks$contr_area[creeks$creek_name == 'knowles'])

# reneau & dietrich sediment density averages 1.75 g cm-3
# total reservoir volume in Mg
sfill * 1.75

# frueh and lancaster 2014

# flux in kg / yr
fl_br_lwr_kg <- 6.07e4 # lower bear
fl_br_lwr_ca <- 2.23 # km2
fl_br_upr_kg <- 1.1e5 # upper bear
fl_br_upr_ca <- 1.35 # km2
fl_cd_kg <- 3.3e4 # cedar fan
fl_cd_ca <- 0.14 # km2
fl_gr_kg <- 3300 # golden ridge
fl_gr_ca <- 1.5 # km2

# estimated unit flux volume in CM
# frueh & lancaster use a density of 1.26 g / cm3
fl_br_lwr_kg * 1e3 / # to g
  1.26 / # divide by density 1.75 g cm-3 to get volume in cm3
  # divide by contr area to get average depth per unit area
  (fl_br_lwr_ca * 1e10) # contr area from km2 to cm2
fl_br_upr_kg * 1e3 / # to g
  1.26 / # divide by density 1.75 g cm-3 to get volume in cm3
  # divide by contr area to get average depth per unit area
  (fl_br_upr_ca * 1e10) # contr area from km2 to cm2
fl_cd_kg * 1e3 / # to g
  1.26 / # divide by density 1.75 g cm-3 to get volume in cm3
  # divide by contr area to get average depth per unit area
  (fl_cd_ca * 1e10) # contr area from km2 to cm2
fl_gr_kg * 1e3 / # to g
  1.26 / # divide by density 1.75 g cm-3 to get volume in cm3
  # divide by contr area to get average depth per unit area
  (fl_gr_ca * 1e10) # contr area from km2 to cm2

# kolmogorov-smirnov
# debris-flow input/output rate
drt <- rec$input[rec$ks == min(rec$ks)[1]]
# fines input/output rate
frt <- recf$input[recf$ks == min(recf$ks)[1]][1]
# gravels input/output rate
grt <- recg$input[recg$ks == min(recg$ks)[1]][1]

# debris-flow events remaining
der <- rec$n[rec$ks == min(rec$ks)[1]][1]
# fines events remaining
fer <- recf$n[recf$ks == min(recf$ks)[1]][1]
# gravel events remaining
ger <- recg$n[recg$ks == min(recg$ks)[1]][1]

# turnover period for debris-flow inputs
tp <- 205.27607176168172



event_stats <- function(drt, frt, grt, der, fer, ger, tp, fill, so) {
  # debris-flow events per turnover period
  det <- tp * drt
  # fluvial events per turnover period
  fet <- tp * frt
  # gravel events per turnover period
  get <- tp * grt
  # total events per turnover
  te <- det + fet + get
  # total events remaining
  tr <- der + fer + ger
  # fill per remaining event
  fe <- fill / tr
  # average volume per year flux
  af <- fe * te / tp
  # denudation rate
  dr <- af / (ein * 1e6)

  # create debris-flow weights
  # log delivery probability p
  ldp <- log(so$dp)
  # median log p
  md_ldp <- median(ldp)
  # distance from median p in log space
  ad_ldp <- ldp - md_ldp
  # range of variation in log space
  rng_ldp <- range(ad_ldp)[2]-range(ad_ldp)[1]
  # relative distance, scaled by range
  wt_dp <- lapply(ad_ldp, function(x) x / rng_ldp) %>% unlist
  # negated to indicate increasing significance as p goes down
  wt_dp <- -wt_dp

  # create fluvial weights
  rng_k <- range(so$k)[2] - range(so$k)[1]
  md_k <- median(so$k)
  wt_k <- lapply(so$k, function(x) (x - md_k) / rng_k) %>% unlist
  wt_k <- -wt_k

  df_wt <- sum(wt_dp[so$type == 'DF'] + 1)
  ff_wt <- sum(wt_k[so$type == 'FF'] + 1)
  fg_wt <- sum(wt_k[so$type == 'FG'] + 1)

  # weighted total events remaining
  wtr <- df_wt + ff_wt + fg_wt
  # scale sampled events to total expected events
  wt_sc <- tr / nrow(so) * wtr
  # fill per weighted remaining event
  wfe <- sfill / wt_sc
  # weighted average volume per year flux
  waf <- wfe * te / tp
  # weighted denudation rate
  wdr <- waf / (ein * 1e6)
  df <- data.frame(fe = fe,
                   af = af,
                   dr = dr,
                   wfe = wfe,
                   waf = waf,
                   wdr = wdr)
  df
}


sum(wt_dp[so$type == 'DF'] + 1) / nrow(so[so$type == 'DF'])
sum(wt_k[so$type == 'FF'] + 1) / nrow(so[so$type == 'FF'])
sum(wt_k[so$type == 'FG'] + 1) / nrow(so[so$type == 'FG'])


length(wt_dp[so$type == 'DF'])
(der - df_wt) / der
(fer - ff_wt) / fer
(ger - fg_wt) / ger
(df$dr - df$wdr) / df$dr


# look up optimized delivery probability and streampower coefficient k for charcoal samples
crks <- creeks_radio
crks_so <- rater1(creeks, .73, .73)
# convert lengths from to-mouth to hipchain from study area
crks_so$hip <- crks_so$ToMouth_km
crks_so$hip[crks_so$creek_name == 'knowles'] <- crks_so$hip[crks_so$creek_name == 'knowles'] - min(crks_so$hip[crks_so$creek_name == 'knowles'])

# subset bear then knowles
br <- crks_so[crks_so$creek_name == 'bear', ]
kn <- crks_so[crks_so$creek_name == 'knowles', ]

# from km to m
# switch orientation from mouth to initiation point
br$hp <- (max(br$hip) - br$hip) * 1000
kn$hp <- (max(kn$hip) - kn$hip) * 1000
# cumulative lengths
br$cm <- br$hp %>% rev %>% cumsum %>% rev
kn$cm <- kn$hp %>% rev %>% cumsum %>% rev
# upstream count
br$ct <- nrow(br):1
kn$ct <- nrow(kn):1
# average length traveled per segment
br$alt <- br$cm / br$ct
kn$alt <- kn$cm / kn$ct

# mean reach length [m]
mn_ln <- mean(c(br$hp, kn$hp))
# mean sediment distance per year [m/y]
mn_d <- mn_ln / tp

# fluvial inheritance
br$fi <- br$alt / mn_d
kn$fi <- kn$alt / mn_d


# match site ids between radiocarbon and transect datasets
# note which sites are named differently
crks$sid[!crks$sid %in% charcoal$site_id]
charcoal$site_id[grep('BC', charcoal$site_id)]
crks$sid[grep('BC-18B', crks$sid)] <- 'BC-18'
crks$sid[grep('BC-20B', crks$sid)] <- 'BC-20'
charcoal$site_id[grep('LK', charcoal$site_id)] %>% sort
crks$sid[grep('LK_118R', crks$sid)] <- 'LK_118r'
crks$sid[grep('LK_57b', crks$sid)] <- 'LK_57'
crks$sid[grep('LK_0.42Q', crks$sid)] <- 'LK_0.42q'
crks$sid[grep('LK_140', crks$sid)] <- 'LK140'
crks$sid[grep('LK_36O', crks$sid)] <- 'LK_36o'
crks$sid[grep('LK_0.16L', crks$sid)] <- 'LK_0.16l'
crks$sid[grep('LK_0.16M', crks$sid)] <- 'LK_0.16m'
crks$sid[grep('LK_151N', crks$sid)] <- 'LK_151n'
crks$sid[grep('LK_187L', crks$sid)] <- 'LK_187l'
crks$sid[grep('LK_145F', crks$sid)] <- 'LK_145f'
crks$sid[grep('LK_145M', crks$sid)] <- 'LK_145m'
crks$sid[grep('LK_145O', crks$sid)] <- 'LK_145o'
crks$sid[grep('LK_145P', crks$sid)] <- 'LK_145p'
charcoal$site_id[grep('DFK', charcoal$site_id)] %>% sort
crks$sid[grep('DFK_42', crks$sid)] <- 'DFK_142'
crks$sid[grep('DFK_110B', crks$sid)] <- 'DFK_110b'
crks$sid[grep('DFK_140L', crks$sid)] <- 'DFK_140l'
crks$sid[grep('DFK_140M', crks$sid)] <- 'DFK_140m'
crks$sid[grep('DFK_140N', crks$sid)] <- 'DFK_140n'
crks$sid[grep('DFK_140P', crks$sid)] <- 'DFK_140p'
crks$sid[grep('DFK_140Q', crks$sid)] <- 'DFK_140q'
crks$sid[grep('DFK_140R', crks$sid)] <- 'DFK_140r'
crks$sid[grep('DFK_140S', crks$sid)] <- 'DFK_140s'
crks$sid[grep('DFK_184B', crks$sid)] <- 'DFK_184b'
crks$sid[grep('DFK_280B', crks$sid)] <- 'DFK_280b'

crks$dp <- 0
crks$type <- 0
crks$hip <- 0
crks$mn <- 0
crks$tm <- 0
crks$id <- 0
for (i in 1:nrow(crks)) {
  crks$dp[i] <- crks_so$dp[crks_so$NODE_ID == crks$node_ids[i]]
  crks$k[i] <- crks_so$k[crks_so$NODE_ID == crks$node_ids[i]]
  crks$hip[i] <- crks_so$hip[crks_so$NODE_ID == crks$node_ids[i]]
  # crks$tm[i] <- crks_so$tm[crks_so$NODE_ID == crks$node_ids[i]]
  crks$type[i] <- charcoal$facies[charcoal$site_id == crks$sid[i]]
  crks$mn[i] <- charcoal$mn[charcoal$site_id == crks$sid[i]]
}


ks_stats <- event_stats(drt, frt, grt, der, fer, ger, tp, sfill, crks)

# debris-flow events per turnover period
det <- tp * drt
# fluvial events per turnover period
fet <- tp * frt
# gravel events per turnover period
get <- tp * grt

# total events per turnover
te <- det + fet + get
# total events remaining
tr <- der + fer + ger

# odds of flux selection
pf <- te / (te + tr)

# poisson distribution
# given a rate, time and number of events
# returns p(k) distribution
fish <- function(rt, t, k) {
  (rt*t)^k*exp(-rt*t) / factorial(k)
}


# megagram to grams 1:1e6
# grams per cm3 granite 2.64
# cm3 to m3: 1e6:1
# Mg km-2 yr-1 to m3 km-2 yr-1 == Mg * 1e6 / density(volume) / 1e6
# m to mm 1:1000
# m3 km-2 yr-1 to mm yr-1 = m3 / 1000

# needle branch min ssc in megagrams km-2 yr-1
mg_lwr <- 31
# flynn creek max ssc in megagrams km-2 yr-1
mg_upr <- 313

options(scipen = 7)
c(mg_lwr, mg_upr) / 2.64 / 1000

# from Hatten etal 2018
# contributing area in km2
fc_ca <- 2.19 # flynn creek
nb_ca <- 0.94 # needle branch
dc_ca <- 3.15 # deer creek

fcg_t <- c(2007:2015, 1959:1973)
# suspended sediment yield in Mg km-2 yr-1
fcg_mg <- c(92, 55, 64, 60, NA, 313, 108, 70, 102,
            30, 22, 114, 46, 38, 76, 428, 98, 44, 23, 48, 41, 64, 372, 20)
nblg_mg <- c(102, 53, 57, 33, NA, NA, 51, 31, 52,
             17, 12, 54, 41, 34, 53, 124, 106, 262, 142, 149, 67, 120, 150, 38)
dcg_mg <- c(127, 69, 127, NA, NA, NA, 92, NA, NA,
            31, 31, 115, 40, 55, 72, 361, 100, 251, 73, 29, 54, 50, 475, 44)

# suspended sediment concentrations mg L-1
fcg_mgl <- c(47.8, 130.1, 48.2, 78.2, 39.9, 7.6, 142, 32.2, 72.9, 65.5, 68.2)
nblg_mgl <- c(32.8, 127.7, 53.8, 65.1, 578.2, 6.3, 43.1, 29.1, 41.5)
dcg_mgl <- c(43.2, 79.6, 32.3, 41.9, 22.1, 1.7, 2.1, 28.7, 7.3, 51.4, 11.4)
all_mgl <- c(fcg_mgl, nblg_mgl, dcg_mgl)
mn_mgl <- mean(all_mgl)

all_ca <- c(rep(fc_ca, length(fcg_mg)),
            rep(nb_ca, length(nblg_mg)),
            rep(dc_ca, length(dcg_mg)))
all_mg = c(fcg_mg, nblg_mg, dcg_mg)
mn_mg = mean(all_mg[!is.na(all_mg)])

# convert Hatton numbers for plotting
# hatton includes larson and sidle numbers from 1959-1973
# obs from 2007-2015
ht_fc_t <- c(92, 55, 64, 60, 313, 108, 70, 102)
ht_fc_ca <- rep(fc_ca, length(ht_fc_t))
ht_nb_t <- c(102, 53, 57, 33, 51, 31, 52)
ht_nb_ca <- rep(nb_ca, length(ht_nb_t))
ht_dc_t <- c(127, 69, 127, 92)
ht_dc_ca <- rep(dc_ca, length(ht_dc_t))

ls_fc_t <- c(30, 22, 114, 46, 38, 76, 428, 98, 44, 23, 48, 41, 64, 372, 20)
ls_fc_ca <- rep(fc_ca, length(ls_fc_t))
ls_nb_t <- c(17, 12, 54, 41, 34, 53, 124, 106, 262, 142, 149, 67, 120, 150, 38)
ls_nb_ca <- rep(nb_ca, length(ls_nb_t))
ls_dc_t <- c(31, 31, 115, 40, 55, 72, 361, 100, 251, 73, 29, 54, 50, 475, 44)
ls_dc_ca <- rep(dc_ca, length(ls_dc_t))

fc_cas <- rep(fc_ca, 23)
nb_cas <- rep(nb_ca, 22)
dc_cas <- rep(dc_ca, 19)
fc_mg <- fcg_mg[!is.na(fcg_mg)]
nb_mg <- nblg_mg[!is.na(nblg_mg)]
dc_mg <- dcg_mg[!is.na(dcg_mg)]




# from O'Connor et al. 2014
# reach name
oc_nm <- c('Chetco', 'Hunter', 'Lower Applegate', 'Illinois', 'Lobster Creek',
           'Broadbent', 'Wilson', 'Miami', 'Tillamook', 'Trask', 'Kilchis', 'Nehalem',
           'Upper Applegate', 'Grants Pass', 'Merlin', 'Galice', 'Powers', 'Bridge',
           'Gravelford', 'Days Creek', 'Roseburg', 'Garden Valley',
           'Coast Range', 'North Umpqua')
# reach type (alluvial, mixed, or bedrock)
oc_type <- c(rep('alluvial', 12), rep('mixed', 10), rep('bedrock', 2))
# contributing area in km2
oc_ca <- c(900, 110, 1990, 1560, 13310, 640, 494, 86, 110, 424, 167, 1840,
           1370, 6470, 8890, 10290, 490, 800, 750, 1960, 4660, 8930, 10490, 3520)
# bed-material flux (including attrition) t yr-1
oc_tyr <- c(111226, 8552, 120292, 272843, 444404, 21462, 25352, 9949, 3396, 15146, 23426, 34692,
            100425, 1919, 132007, 225962, 19677, 5061, 1917, 113007, 165172, 145594, 63035, 68965)
# flux in t km-2 yr-1
oc_tkmy <- oc_tyr / oc_ca


# from Reneau & Dietrich 1991
# reach name
rd_nm <- c('Needle Branch', 'Needle Branch', 'Flynn Creek', 'Deer Creek', 'Deer Creek',
           'Cooper Creek', 'Sutherlin Creek', 'Olalla Creek',
           'Yaquina River', 'Alsea River', 'Siuslaw River')
rd_ca <- c(nb_ca, nb_ca, fc_ca, dc_ca, dc_ca, 11.4, 23.3, 158, 655, 865, 1523)
rd_t <- c(53, 146, 98, 97, 157, 179, 174, 95, 129, 187, 125)



# from Larson & Sidle 1980
# Alsea study area
# contributing area originally in mi^2, conver to km^2
tokm <- 2.589
ls_fc_ca <- 0.78 * tokm
ls_dc_ca <- 1.17 * tokm
ls_nb_ca <- 0.27 * tokm
# alsea numbers are reproduced by Hatton etal 2018
# oak creek
ls_oc_nm <- c('1978', '1979', '1980', 'Mean')
# oak creek contributing area, from mi2 to km2
ls_oc_ca <- rep(2.8 * tokm, length(ls_oc_nm))
# from tons per mi2 to tonnes per km2
toMg <- 0.907185
ls_oc_t <- c(51, 23, 24, 33) / tokm * toMg

# andrews data
ls_an <- data.table::fread('larson_sidle_andrews.csv')
ls_an_t <- ls_an$ss[!is.na(ls_an$ss) & ls_an$year != 'Mean'] / tokm * toMg
ls_an_ca <- ls_an$ca[!is.na(ls_an$ss) & ls_an$year != 'Mean'] * tokm
plot(ls_an_ca, ls_an_t, log = 'xy')

ls_an_ca %>% factor %>% levels
# percent bedload of total
ls_an_bd <- ls_an$bed[!is.na(ls_an$total)] / tokm * toMg
ls_an_tl <- ls_an$total[!is.na(ls_an$total)] / tokm * toMg
ls_an_bp <- ls_an_bd / ls_an_tl
ls_an_bca <- ls_an$ca[!is.na(ls_an$total)] * tokm

plot(ls_an_bca, ls_an_bp, log = 'x')


# coyote creek (near Roseburg)
ls_cy <- data.table::fread('larson_sidle_coyote.csv')
ls_cy_t <- ls_cy$ss[!is.na(ls_cy$ss) & ls_cy$year != 'Mean'] / tokm * toMg
ls_cy_ca <- ls_cy$ca[!is.na(ls_cy$ss) & ls_cy$year != 'Mean'] * tokm
plot(ls_cy_ca, ls_cy_t, log = 'xy')

# percent bedload of total
ls_cy_bd <- ls_cy$bed[!is.na(ls_cy$total)] / tokm * toMg
ls_cy_tl <- ls_cy$total[!is.na(ls_cy$total)] / tokm * toMg
ls_cy_bp <- ls_cy_bd / ls_cy_tl
ls_cy_bca <- ls_cy$ca[!is.na(ls_cy$total)] * tokm

png('bedload_fraction.png', height = 17, width = 21, units = 'cm', res = 300)
plot(ls_an_tl, ls_an_bp, log = 'x',
     xlab = 'flux [Mg km^-2]',
     ylab = 'bedload fraction',
     pch = 20, col = get_palette('ocean'))
points(ls_cy_tl, ls_cy_bp, pch = 20, col = get_palette('crimson', .5))
abline(h = mean(ls_cy_bp), lty = 2, col = get_palette('crimson', .9))
abline(h = mean(ls_an_bp), lty = 2, col = get_palette('ocean', .5))
abline(h = 1 - .7226, lty = 2, col = get_palette('gold', .9))
legend('bottomright', legend = c('Andrews EF', 'Coyote Creek', 'study area (predicted)', 'mean'),
       pch = c(20, 20, 20, NA), lty = c(NA, NA, NA, 2),
       col = get_palette(c('ocean', 'crimson', 'gold', 'charcoal'), .75))
text(8500, .21, labels = 'from Larson & Sidle 1980')
dev.off()



# plot transit times
plot(trans$rate, trans$mean)
plot(rec$input, rec$n)
trans_n <- 0
for (i in 1:nrow(trans)) {
  trans_n[i] <- mean(rec$n[rec$input > trans$rate[i]*0.9 & rec$input < trans$rate[i]*1.1])
}
png('event_transit_rate.png', height = 17, width = 21, units = 'cm', res = 300)
par(mar = c(5, 5, 3, 5))
plot(trans$rate, trans_n, pch = 20, col = get_palette('ocean'),
     xlab = 'interarrival rate',
     ylab = 'events remaining')
par(new = T)
plot(trans$rate, trans$mean, pch = 20, col = get_palette('crimson'),
     xaxt = 'n', yaxt = 'n',
     xlab = '', ylab = '')
axis(side = 4)
mtext('mean transit time [years]', side = 4, line = 3)
legend('top', legend = c('events', 'transit time'),
       fill = get_palette(c('ocean', 'crimson'), .7))
dev.off()

par(mar = c(5, 5, 3, 3))
png('event_transit_time.png', height = 17, width = 21, units = 'cm', res = 300)
plot(trans$mean, trans_n, pch = 20, col = get_palette('ocean'),
     xlab = 'mean transit time',
     ylab = 'events remaining')
dev.off()




png('bedload_fraction.png', height = 17, width = 21, units = 'cm', res = 300)
plot(ls_cy_bca, ls_cy_bp, pch = 20, col = get_palette('ocean', .5),
     xlab = 'contributing area [km2]',
     ylab = 'bedload fraction',
     xlim = c(0.1, 1e4), ylim = c(0,1), log = 'x')
points(ls_an_bca, ls_an_bp, pch = 20, col = get_palette('crimson', .5))
legend('topleft', legend = c('Coyote Creek', 'Andrews Forest'),
       fill = get_palette(c('ocean', 'crimson'), .7))
text(0.4, 0.87, label = 'from Larson & Sidle 2018')
dev.off()

setwd('/media/erik/catacomb/research/')
options(scipen = 7)
png('flux_Mg_per_km2.png', height = 17, width = 21, units = 'cm', res = 300)
pal <- get_palette(c('crimson', 'gold', 'hardwood', 'ocean', 'charcoal'), .5)
plot(oc_ca, oc_tkmy, pch = 20, col = pal[4],
     xlab = 'contributing area (km2)',
     ylab = 'flux (Mg km-2)', log = 'xy', xlim = c(0.1, 1e7), ylim = c(0.1, max(ls_an_t)))
points(rd_ca[-c(1:5)], rd_t[-c(1:5)], pch = 20, col = pal[1])
points(rd_ca[c(1:2)], rd_t[c(1:2)], pch = 1, col = pal[1])
points(rd_ca[c(3)], rd_t[c(3)], pch = 0, col = pal[1])
points(rd_ca[c(4:5)], rd_t[c(4:5)], pch = 2, col = pal[1])
points(ls_fc_ca, ls_fc_t, pch = 0, col = pal[3])
points(ls_nb_ca, ls_nb_t, pch = 1, col = pal[3])
points(ls_dc_ca, ls_dc_t, pch = 2, col = pal[3])
points(ls_oc_ca, ls_oc_t, pch = 3, col = pal[3])
points(ls_an_ca, ls_an_t, pch = 4, col = pal[3])
points(ls_cy_ca, ls_cy_t, pch = 5, col = pal[3])
points(ht_fc_ca, ht_fc_t, pch = 0, col = pal[2])
points(ht_nb_ca, ht_nb_t, pch = 1, col = pal[2])
points(ht_dc_ca, ht_dc_t, pch = 2, col = pal[2])
abline(h = 122.5, lty = 2, col = pal[1])
abline(h = mn_mg, lty = 2, col = pal[2])
legend('topright',
       legend = c("Hatton etal 2018", "O'Connor etal 2014", "Reneau & Dietrich 1991", "Larson & Sidle 1980",
                  "Flynn Creek", "Needle Branch", "Deer Creek", "Oak Creek",
                  "Andrews Forest", "Coyote Creek", "Mean"),
       col = c(pal[c(2, 4, 1, 3, rep(5, 8))]), pch = c(20, 20, 20, 20, 0:5, NA),
       lty = c(rep(NA, 10), 2)
       )
dev.off()

# turnover time based on Reneau & Dietrich sed transport rate
(sfill * 1.75) / (122.5 * ein) # 538.0983
# turnover time based on Hatton etal sed transport rate
(sfill * 1.75) / (mn_mg * ein) # 672.6228


palette <- get_palette(c('ocean', 'crimson', 'gold', 'charcoal'), .5)
png('hatten_creeks_Mg.png', height = 17, width = 21, units = 'cm', res = 300)
plot(emp_cdf(dcg_mg), pch = 20, col = palette[2],
     xlab = 'Mg km-2 yr-1', xlim = c(0, max(dcg_mg[!is.na(dcg_mg)])),
     ylab = 'CDF', ylim = c(0,1))
abline(v = mn_mg, lty = 2, col = palette[4])
lines(emp_cdf(fcg_mg), lwd = 2, col = palette[1])
points(emp_cdf(fcg_mg), pch = 20, col = palette[1])
lines(emp_cdf(nblg_mg), lwd = 2, col = palette[2])
points(emp_cdf(nblg_mg), pch = 20, col = palette[2])
lines(emp_cdf(dcg_mg), lwd = 2, col = palette[3])
legend('bottomright', legend = c('Flynn Creek', 'Needle Branch', 'Deer Creek', 'mean'),
       fill = palette)
text(425, .22, labels = c('from Hatten et al. 2018'))
dev.off()

png('hatten_creeks_mgl.png', height = 17, width = 21, units = 'cm', res = 300)
plot(emp_cdf(dcg_mgl), pch = 20, col = palette[3],
     xlab = 'mg L-1', xlim = c(0, max(all_mgl)),
     ylab = 'CDF', ylim = c(0,1))
abline(v = mn_mgl, lty = 2, col = palette[4])
lines(emp_cdf(fcg_mgl), lwd = 2, col = palette[1])
points(emp_cdf(fcg_mgl), pch = 20, col = palette[1])
lines(emp_cdf(nblg_mgl), lwd = 2, col = palette[2])
points(emp_cdf(nblg_mgl), pch = 20, col = palette[2])
lines(emp_cdf(dcg_mgl), lwd = 2, col = palette[3])
legend('bottomright', legend = c('Flynn Creek', 'Needle Branch', 'Deer Creek', 'mean'),
       fill = palette)
text(520, .22, labels = c('from Hatten et al. 2018'))
dev.off()


# function to split observations of a predictor into bins
# returns the mean sd and count of fit for each bin
bin_stat <- function(pred, fit, bins = 10) {
  # vector of means, sds and counts
  mns <- 0
  sds <- 0
  ns <- 0
  # divide range of pred into bins number of steps

  rng <- seq(min(pred), max(pred), (max(pred) - min(pred)) / bins)
  for (i in 1:bins) {
    # select fits where pred in is range
    if (i == 1) {
      bin <- fit[pred <= rng[i]]
    } else {
      bin <- fit[pred > rng[i-1] & pred <= rng[i]]
    }
    mns[i] <- mean(bin)
    sds[i] <- sd(bin)
    ns[i] <- length(bin)
  }
  data.frame(mns = mns, sds = sds, ns = ns, rng = rng[1:length(mns)])
}

ks_mns <- bin_stat(rec$input, rec$ks, 200)
ks_lwr <- ks_mns$mns - 1.96 * (ks_mns$sds / sqrt(ks_mns$ns))
ks_upr <- ks_mns$mns + 1.96 * (ks_mns$sds / sqrt(ks_mns$ns))

kp_mns <- bin_stat(rec$input, rec$kp, 200)
kp_lwr <- kp_mns$mns - 1.96 * (kp_mns$sds / sqrt(kp_mns$ns))
kp_upr <- kp_mns$mns + 1.96 * (kp_mns$sds / sqrt(kp_mns$ns))

ad_mns <- bin_stat(rec$input, rec$ad, 200)
ad_lwr <- ad_mns$mns - 1.96 * (ad_mns$sds / sqrt(ad_mns$ns))
ad_upr <- ad_mns$mns + 1.96 * (ad_mns$sds / sqrt(ad_mns$ns))

ch_mns <- bin_stat(rec$input, rec$ch, 200)
ch_lwr <- ch_mns$mns - 1.96 * (ch_mns$sds / sqrt(ch_mns$ns))
ch_upr <- ch_mns$mns + 1.96 * (ch_mns$sds / sqrt(ch_mns$ns))


png('ks_kp_fit.png', height = 17, width = 21, units = 'cm', res = 300)
par(mar = c(5, 5, 3, 5))
plot(ks_mns$rng, ks_mns$mns, type = 'l', lwd = 2, col = get_palette('crimson', .7),
    ylim = c(0.1805,0.184), xlim = c(0.4, 1.35),
    xlab = 'interarrival rate', ylab = 'Kolmogorov-Smirnov test')
lines(ks_mns$rng, ks_lwr, lwd = 1, lty = 2, col = get_palette('crimson', .7))
lines(ks_mns$rng, ks_upr, lwd = 1, lty = 2, col = get_palette('crimson', .7))
par(new = T)
plot(kp_mns$rng, kp_mns$mns, type = 'l', lwd = 2, col = get_palette('violet', .7),
     xlim = c(0.4, 1.35), ylim = c(0.2292, 0.231),
     xaxt = 'n', yaxt = 'n',
     xlab = '', ylab = '')
lines(kp_mns$rng, kp_lwr, lwd = 1, lty = 2, col = get_palette('violet', .7))
lines(kp_mns$rng, kp_upr, lwd = 1, lty = 2, col = get_palette('violet', .7))
axis(side = 4)
mtext('Kuiper test', side = 4, line = 3)
legend('bottomleft', legend = c('K-S test', 'Kuiper test', '95% CI'),
       lwd = 2, lty = c(1,1,2), col = get_palette(c('crimson', 'violet', 'charcoal'), .7))
dev.off()

png('ad_ch_fit.png', height = 17, width = 21, units = 'cm', res = 300)
par(mar = c(5, 5, 3, 5))
plot(ad_mns$rng, ad_mns$mns, type = 'l', lwd = 2, col = get_palette('ocean', .7),
     ylim = c(1.98,2.1), xlim = c(0.15, 0.7),
     xlab = 'interarrival rate', ylab = 'Anderson-Darling test')
lines(ad_mns$rng, ad_lwr, lwd = 1, lty = 2, col = get_palette('ocean', .7))
lines(ad_mns$rng, ad_upr, lwd = 1, lty = 2, col = get_palette('ocean', .7))
par(new = T)
plot(ch_mns$rng, ch_mns$mns, type = 'l', lwd = 2, col = get_palette('gold', .7),
     xlim = c(0.15, 0.7), ylim = c(5, 5.3),
     xaxt = 'n', yaxt = 'n',
     xlab = '', ylab = '')
lines(ch_mns$rng, ch_lwr, lwd = 1, lty = 2, col = get_palette('gold', .7))
lines(ch_mns$rng, ch_upr, lwd = 1, lty = 2, col = get_palette('gold', .7))
axis(side = 4)
mtext('Chi-squared test', side = 4, line = 3)
legend('bottomright', legend = c('A-D test', 'Chi-squared test', '95% CI'),
       lwd = 2, lty = c(1,1,2), col = get_palette(c('ocean', 'gold', 'charcoal'), .7))
dev.off()

# confidence intervals for statistical tests
# debris-flows
(ad_mn <- ad_mns$rng[ad_mns$mns == min(ad_mns$mns)])
ad_in <- ad_mns$rng[ad_lwr < ad_upr[ad_mns$rng == ad_mn]]
(ad_lw <- min(ad_in[!is.na(ad_in)]))
(ad_up <- max(ad_in[!is.na(ad_in)]))

(ch_mn <- ch_mns$rng[ch_mns$mns == min(ch_mns$mns)])
ch_in <- ch_mns$rng[ch_lwr < ch_upr[ch_mns$rng == ch_mn]]
(ch_lw <- min(ch_in[!is.na(ch_in)]))
(ch_up <- max(ch_in[!is.na(ch_in)]))

(kp_mn <- kp_mns$rng[kp_mns$mns == min(kp_mns$mns)])
kp_in <- kp_mns$rng[kp_lwr < kp_upr[kp_mns$rng == kp_mn]]
(kp_lw <- min(kp_in[!is.na(kp_in)]))
(kp_up <- max(kp_in[!is.na(kp_in)]))

(ks_mn <- ks_mns$rng[ks_mns$mns == min(ks_mns$mns)])
ks_in <- ks_mns$rng[ks_lwr < ks_upr[ks_mns$rng == ks_mn]]
(ks_lw <- min(ks_in[!is.na(ks_in)]))
(ks_up <- max(ks_in[!is.na(ks_in)]))

# confidence intervals for statistical tests
# fines
ksf_mns <- bin_stat(recf$input, recf$ks, 200)
ksf_lwr <- ksf_mns$mns - 1.96 * (ksf_mns$sds / sqrt(ksf_mns$ns))
ksf_upr <- ksf_mns$mns + 1.96 * (ksf_mns$sds / sqrt(ksf_mns$ns))

kpf_mns <- bin_stat(recf$input, recf$kp, 200)
kpf_lwr <- kpf_mns$mns - 1.96 * (kpf_mns$sds / sqrt(kpf_mns$ns))
kpf_upr <- kpf_mns$mns + 1.96 * (kpf_mns$sds / sqrt(kpf_mns$ns))

adf_mns <- bin_stat(recf$input, recf$ad, 200)
adf_lwr <- adf_mns$mns - 1.96 * (adf_mns$sds / sqrt(adf_mns$ns))
adf_upr <- adf_mns$mns + 1.96 * (adf_mns$sds / sqrt(adf_mns$ns))

chf_mns <- bin_stat(recf$input, recf$ch, 200)
chf_lwr <- chf_mns$mns - 1.96 * (chf_mns$sds / sqrt(chf_mns$ns))
chf_upr <- chf_mns$mns + 1.96 * (chf_mns$sds / sqrt(chf_mns$ns))

(adf_mn <- adf_mns$rng[adf_mns$mns == min(adf_mns$mns)])
adf_in <- adf_mns$rng[adf_lwr < adf_upr[adf_mns$rng == adf_mn]]
(adf_lw <- min(adf_in[!is.na(adf_in)]))
(adf_up <- max(adf_in[!is.na(adf_in)]))

(chf_mn <- chf_mns$rng[chf_mns$mns == min(chf_mns$mns)])
chf_in <- chf_mns$rng[chf_lwr < chf_upr[chf_mns$rng == chf_mn]]
(chf_lw <- min(chf_in[!is.na(chf_in)]))
(chf_up <- max(chf_in[!is.na(chf_in)]))

(kpf_mn <- kpf_mns$rng[kpf_mns$mns == min(kpf_mns$mns)])
kpf_in <- kpf_mns$rng[kpf_lwr < kpf_upr[kpf_mns$rng == kpf_mn]]
(kpf_lw <- min(kpf_in[!is.na(kpf_in)]))
(kpf_up <- max(kpf_in[!is.na(kpf_in)]))

(ksf_mn <- ksf_mns$rng[ksf_mns$mns == min(ksf_mns$mns)])
ksf_in <- ksf_mns$rng[ksf_lwr < ksf_upr[ksf_mns$rng == ksf_mn]]
(ksf_lw <- min(ksf_in[!is.na(ksf_in)]))
(ksf_up <- max(ksf_in[!is.na(ksf_in)]))


# confidence intervals for statistical tests
# gravels
ksg_mns <- bin_stat(recg$input, recg$ks, 200)
ksg_lwr <- ksg_mns$mns - 1.96 * (ksg_mns$sds / sqrt(ksg_mns$ns))
ksg_upr <- ksg_mns$mns + 1.96 * (ksg_mns$sds / sqrt(ksg_mns$ns))

kpg_mns <- bin_stat(recg$input, recg$kp, 200)
kpg_lwr <- kpg_mns$mns - 1.96 * (kpg_mns$sds / sqrt(kpg_mns$ns))
kpg_upr <- kpg_mns$mns + 1.96 * (kpg_mns$sds / sqrt(kpg_mns$ns))

adg_mns <- bin_stat(recg$input, recg$ad, 200)
adg_lwr <- adg_mns$mns - 1.96 * (adg_mns$sds / sqrt(adg_mns$ns))
adg_upr <- adg_mns$mns + 1.96 * (adg_mns$sds / sqrt(adg_mns$ns))

chg_mns <- bin_stat(recg$input, recg$ch, 200)
chg_lwr <- chg_mns$mns - 1.96 * (chg_mns$sds / sqrt(chg_mns$ns))
chg_upr <- chg_mns$mns + 1.96 * (chg_mns$sds / sqrt(chg_mns$ns))

(adg_mn <- adg_mns$rng[adg_mns$mns == min(adg_mns$mns)])
adg_in <- adg_mns$rng[adg_lwr < adg_upr[adg_mns$rng == adg_mn]]
(adg_lw <- min(adg_in[!is.na(adg_in)]))
(adg_up <- max(adg_in[!is.na(adg_in)]))

(chg_mn <- chg_mns$rng[chg_mns$mns == min(chg_mns$mns)])
chg_in <- chg_mns$rng[chg_lwr < chg_upr[chg_mns$rng == chg_mn]]
(chg_lw <- min(chg_in[!is.na(chg_in)]))
(chg_up <- max(chg_in[!is.na(chg_in)]))

(kpg_mn <- kpg_mns$rng[kpg_mns$mns == min(kpg_mns$mns)])
kpg_in <- kpg_mns$rng[kpg_lwr < kpg_upr[kpg_mns$rng == kpg_mn]]
(kpg_lw <- min(kpg_in[!is.na(kpg_in)]))
(kpg_up <- max(kpg_in[!is.na(kpg_in)]))

(ksg_mn <- ksg_mns$rng[ksg_mns$mns == min(ksg_mns$mns)])
ksg_in <- ksg_mns$rng[ksg_lwr < ksg_upr[ksg_mns$rng == ksg_mn]]
(ksg_lw <- min(ksg_in[!is.na(ksg_in)]))
(ksg_up <- max(ksg_in[!is.na(ksg_in)]))


# knowles creek pebble counts
peb <- data.table::fread('/media/erik/catacomb/research/knowles_pebble_count.csv')
pebs <- data.table::fread('/media/erik/catacomb/research/knowles_subsurface_pebbles.csv')

png('knowles_pebble_count.png', height = 17, width = 21, units = 'cm', res = 300)
plot(emp_cdf(peb$base_1a), type = 'l', lwd = 2, col = get_palette('ocean'),
     xlab = 'b-axis width [mm]', ylab = 'CDF of samples', log = 'x')
lines(emp_cdf(pebs$b), lwd = 2, col = get_palette('crimson'))
points(emp_cdf(pebs$b), pch = 20, col = get_palette('crimson'))
lines(exp(seq(2.71, 4.32, .01)), exp(-off + 1.683 * seq(2.71, 4.32, .01)))

points(emp_cdf(peb$base_1a), pch = 20, col = get_palette('ocean'))
lines(emp_cdf(peb$site_1b), lwd = 2, col = get_palette('crimson'))
points(emp_cdf(peb$site_1b), pch = 20, col = get_palette('crimson'))
lines(emp_cdf(peb$site_1a), lwd = 2, col = get_palette('crimson'))
points(emp_cdf(peb$site_1a), pch = 20, col = get_palette('crimson'))
lines(emp_cdf(peb$base_1b), lwd = 2, col = get_palette('ocean'))
points(emp_cdf(peb$base_1b), pch = 20, col = get_palette('ocean'))
legend('bottomright', legend = c('K38 +1260', 'K38 - 724'),
       fill = get_palette(c('crimson', 'ocean'), .7))
dev.off()

(bm <- lm(c ~ log(b) + exp(b), data = b5_9d))
summary(bm)
surf <- c(peb$site_1a, peb$site_1b, peb$base_1a, peb$base_1b)
# rarity weight based on median b-axis width
surf_w <- log(surf) / median(log(surf))
surf_w[surf_w > 1] <- 1 / surf_w[surf_w > 1]
# normalize weight into pmf
surf_n <- surf_w / sum(surf_w)
plot(sort(surf), cumsum(sort(surf_n)), type = 'l', lwd = 3,
     xlab = 'b-axis width [mm]',
     ylab = 'CDF of samples',
     col = get_palette('crimson'))
points(emp_cdf(surf), pch = 20, col = get_palette('charcoal'))
prd <- log(1:1000) / median(log(surf))
prd[prd > 1] <- 1 / prd[prd > 1]
prd <- prd / sum(prd)
prd <- cumsum(prd)
lines(1:1000, prd)
# intercept at surface d50
off <- median(log(surf)) * b5_9m$coefficients[2] - log(0.5)
exp(off / b5_9m$coefficients[2]) # max pred
exp((-exp(1) + off) / b5_9m$coefficients[2]) # min pred
# predicted size when half as rare
d75 <- exp((off - log(0.75)) / b5_9m$coefficients[2])
exp(log(20) * b5_9m$coefficients[2] - off)

rarity <- function(d, d50, mod, win = .95, r = 1) {
  # given a b-axis width, the d50, and a model fit
  # returns the relative rarity of the diameter

  # intercept at surface d50
  off <- log(d50) * mod$coefficients[2] - log(0.5)
  # pred width when less common
  upr <- exp((off + log(0.5 + 0.5 * win)) / mod$coefficients[2])
  # pred width when less rare
  lwr <- exp((off + log(0.5 - 0.5 * win)) / mod$coefficients[2])
  # lower and higher paths
  if (d < d50) {
    if (d >= lwr) {
      # cdf at d
      c <- exp(log(d) * mod$coefficients[2] - off)
      # relative rarity to d50
      r <- r * 0.5 / c
    }
    if (d < lwr) {
      r <- rarity(d, lwr, mod, r = r * 0.5 / (0.5 - 0.5 * win))
    }
  }
  if (d > d50) {
    if (d <= upr) {
      # cdf at d
      c <- exp(log(d) * mod$coefficients[2] - off)
      # relative rarity to d50
      r <- r *  0.5 / c
    }
    if (d > upr) {
      r <- rarity(d, upr, mod, r = r * 0.5 / (0.5 + 0.5 * win))
    }
  }
  r
}

rarer <- function(d, d50, mod, win = .95, r = 1) {
  # given a b-axis width, the d50, and a model fit
  # returns the relative rarity of the diameter

  # intercept at surface d50
  off <- log(d50) * mod$coefficients[2] - log(0.5)
  # pred width when less common
  upr <- exp((off + log(0.5 + 0.5 * win)) / mod$coefficients[2])
  # pred width when less rare
  lwr <- exp((off + log(0.5 - 0.5 * win)) / mod$coefficients[2])

  # test if d is in range, else recurse up or down
  if (d < lwr) {
    r <- rarer(d, lwr, mod, r = r * 0.5 / (0.5 - 0.5 * win))
  }

  if (d > upr) {
    r <- rarer(d, upr, mod, r = r * 0.5 / (0.5 + 0.5 * win))
  }

  # cdf at d
  c <- exp(log(d) * mod$coefficients[2] - off)
  # relative rarity to d50
  r <- r * 0.5 / c
  r
}

cdf_hook <- function(obs) {
  # given observations in random order
  # returns a vector of the cdf value of each ob
  vals <- rep(0, length(obs))
  cdf <- emp_cdf(obs)
  x <- cdf[ , 1]
  y <- cdf[ , 2]
  for (i in 1:length(x)) {
    if (i == 1) {
      vals[obs <= x[i]] <- y[i]
    }
    if (i > 1) {
      vals[obs > x[i-1] & obs <= x[i]] <- y[i]
    }
  }
  vals
}
b_cdf <- emp_cdf(pebs$b)
b5_9 <- pebs[pebs$b >= 5 & pebs$b <= 9, ]
b5_9c <- cdf_hook(log(b5_9$b))
b5_9d <- data.frame(b = b5_9$b, lb = log(b5_9$b), c = b5_9c)
b5_9m <- lm(c ~ lb, data = b5_9d)
summary(b5_9m)

b4 <- pebs[pebs$b >= 4.5 & pebs$b <= 7.5, ]
b4c <- cdf_hook(log(b4$b))
b4d <- data.frame(b = b4$b, lb = log(b4$b), c = b4c)
b4m <- lm(c ~ lb, data = b4d)
summary(b4m)




xs <- 1:500 # sequence of b-axis widths
rs <- 0 # rarity over xs
rw <- 0 # rarity weights over xs
rss <- 0 # rarity subsurface
rws <- 0 # rarity weights subsurface
for (i in seq_along(xs)) {
  rs[i] <- rarity(xs[i], median(surf), b4m)
  rw[i] <- log(xs[i]) / log(median(surf))
  rss[i] <- rarity(xs[i], median(pebs$b), b4m)
  rws[i] <- log(xs[i]) / log(median(pebs$b))
}
rr <- rs / sum(rs) # relative rarity surface
rrs <- rss / sum(rss) # relative rarity subsurface
rw[rw > 1] <- 1 / rw[rw > 1]
rws1 <- rws
rws1[rws < 1] <- 1
rws1[rws > 1] <- 1 / rws1[rws > 1]
rws <- rws1^3

mw <- log(xs)^3.6
mws <- log(xs)^4
wr <- (rr * rw^3 * mw) / sum(rr * rw^3 * mw) # weighted rarity surface
wrs <- (rrs * rws * mws) / sum(rrs * rws * mws) # weighted rarity subsurface


# estimate counts for fines fractions
# predict mass by width
gmod <- lm(log(wgt) ~ log(b), data = pebs[pebs$wgt > 0, ])
summary(gmod)

# 70.39 - 1.56 fines 2-4mm
# 51.55 - 1.56 fines < 2mm
d3d <- 1.3
d1d <- 0.33
d3 <- (70.39 - 1.56) / exp(gmod$coefficients[1] + gmod$coefficients[2] * d3d)
d1 <- (51.55 - 1.56) / exp(gmod$coefficients[1] + gmod$coefficients[2] * d1d)
d2 <- (70.39 + 51.55 + 1.56 * 2) / exp(gmod$coefficients[1] + gmod$coefficients[2])

subs <- c(pebs$b, rep(d3d, ceiling(d3)), rep(d1d, ceiling(d1)))
# subs <- c(pebs$b, rep(1, ceiling(d2)))

d <- data.frame(b = surf, c = cdf_hook(surf))
ds <- data.frame(b = pebs$b, c= cdf_hook(pebs$b))
d$r <- vapply(surf, function(x) rarity(x, median(surf), b4m), 0)
ds$r <- vapply(ds$b, function(x) rarity(x, median(pebs$b), b4m), 0)
d$r <- vapply(surf, function(x) rarer(x, median(surf), b4m), 0)
ds$r <- vapply(ds$b, function(x) rarer(x, median(pebs$b), b4m), 0)

plot(ds$b, ds$r, log = 'xy')
plot(d$b, d$r, log = 'xy')

d$rw <- log(d$b) / log(median(d$b))
d$rw[d$rw > 1] <- 1 / d$rw[d$rw > 1]
# d$m <- exp(predict(gmod, newdata = d))
# d$rw <- d$rw * (1 - d$m / max(d$m))
ds$rw <- log(ds$b) / log(median(pebs$b))

# umbrella distribution
ds$rw <- ds$b / median(pebs$b) # width relative to median
d$rw <- d$b / median(d$b) # width relative to median
# ds$rw[ds$rw < 1] <- median(pebs$b) / ds$b[ds$rw < 1]
# above median, increasing width decreases selection probability
# subtract 1 from rarity weight and max because median == 1 to get relative probability
# subtract from one to indicate decreasing selection likelihood with increasing width
ds$rw[ds$rw < 1] <- 1
ds$rw[ds$rw > 1] <- 1 - (ds$rw[ds$rw > 1] - 1) / (max(ds$rw) - 1)
d$rw[d$rw > 1] <- 1 - (d$rw[d$rw > 1] - 1) / (max(d$rw) - 1)
# below median, decreasing width decreases selection probability
# invert the width ratio
# ds$rw[ds$b < median(ds$b)] <- 1 / ds$rw[ds$b < median(ds$b)]
d$rw[d$b < median(d$b)] <- 1 / d$rw[d$b < median(d$b)]
# convert to relative probability using same transform as above median
# ds$rw[ds$rw > 1] <- 1 - (ds$rw[ds$rw > 1] - 1) / (max(ds$rw) - 1)
d$rw[d$rw > 1] <- 1 - (d$rw[d$rw > 1] - 1) / (max(d$rw) - 1)

# ds$rw[ds$rw > 1] <- 1 / log(ds$rw[ds$rw > 1])
# ds$m <- exp(predict(gmod, newdata = ds))
# ds$lm <- log(ds$m) - min(log(ds$m))
# ds$rw <- ds$rw * (1 - ds$lm / max(ds$lm))
plot(ds$b, ds$rw, log = 'x')

png('scale_weight.png', height = 17, width = 21, units = 'cm', res = 300)
plot(d$b, d$rw, log = 'x', pch = 20, col = get_palette('ocean'),
     xlab = 'b-axis widtth [mm]',
     ylab = 'scale weight')
points(ds$b, ds$rw, pch = 20, col = get_palette('crimson'))
legend('topright', legend = c('surface', 'subsurface'),
       fill = get_palette(c('ocean', 'crimson'), .77))
dev.off()

dx <- data.frame(b = 1:500)
dx$r <- vapply(dx$b, function(x) rarity(x, median(d$b), b4m), 0)
dx$rw <- log(dx$b) / log(median(d$b))
dx$rw[dx$rw > 1] <- 1 / dx$rw[dx$rw > 1]


m <- lm(c ~ log(r) + rw, data = d[d$r != 0, ])
summary(m)
ms <- lm(c ~ log(r) + rw, data = ds)
summary(ms)
m <- lm(c ~ log(b) + rw, data = d)
summary(m)
ms <- lm(c ~ log(b) + rw, data = ds)
summary(ms)

# preferred model
m <- lm(c ~ log(b) + rw + b, data = d)
summary(m)
ms <- lm(c ~ log(b) + rw + b, data = ds)
summary(ms)

m <- lm(c ~ log(b) + b, data = d)
summary(m)
ms <- lm(c ~ log(b) + b, data = ds)
summary(ms)

m <- lm(c ~ log(b), data = d)
summary(m)
ms <- lm(c ~ log(b), data = ds)
summary(ms)


evac <- log(ds$b) * ms$coefficients[4]
evac <- evac - min(evac)
evac <- 1 - evac / max(evac)

plot(ds$b, evac, log = 'x',
     xlab = 'b-axis width [mm]',
     ylab = 'relative evacuation probability')

msp1 <- predict(m, newdata = dx)
points(dx$b, msp1)

dp <- predict(m, newdata = d)
dp[dp < 0] <- 0
dp[dp > 1] <- 1

msp <- predict(ms, newdata = ds)
msp[msp < 0] <- 0
msp[msp > 1] <- 1

setwd('/media/erik/catacomb/research')
png('gravel_fit1.png', height = 17, width = 21, units = 'cm', res = 300)
plot(sort(d$b), sort(dp), log = 'x',
     xlab = 'b-axis width [mm]', xlim = c(2, 458),
     ylab = 'CDF of samples',
     type = 'l', lwd = 3, col = get_palette('crimson', .7))
points(emp_cdf(d$b), pch = 20, col = get_palette('charcoal', .2))
points(emp_cdf(ds$b), pch = 20, col = get_palette('hardwood', .2))
lines(sort(ds$b), sort(msp), lwd = 3, col = get_palette('crimson', .7))
legend('bottomright', legend = c('surface', 'subsurface', 'modeled'),
       fill = get_palette(c('charcoal', 'hardwood', 'crimson'), .8))
dev.off()

summary(b4m)
tp <- 205
# sternberg mass loss coefficient over time
# estimating mass from b-axis width
# turnover period from k-s test of debris-flow charcoal ages
st <- log(predict(gmod, newdata = data.frame(b = median(d$b)))) / tp
#
# look up optimized delivery probability and streampower coefficient k for charcoal samples
crks <- creeks_radio
crks_so <- rater1(creeks, .73, .73)
# convert lengths from to-mouth to hipchain from study area
crks_so$hip <- crks_so$ToMouth_km
crks_so$hip[crks_so$creek_name == 'knowles'] <- crks_so$hip[crks_so$creek_name == 'knowles'] - min(crks_so$hip[crks_so$creek_name == 'knowles'])

# subset bear then knowles
br <- crks_so[bk$creek_name == 'bear', ]
kn <- crks_so[bk$creek_name == 'knowles', ]

# from km to m
# switch orientation from mouth to initiation point
br$hp <- (max(br$hip) - br$hip) * 1000
kn$hp <- (max(kn$hip) - kn$hip) * 1000
library(magrittr)
# cumulative lengths
br$cm <- br$hp %>% rev %>% cumsum %>% rev
kn$cm <- kn$hp %>% rev %>% cumsum %>% rev
# upstream count
br$ct <- nrow(br):1
kn$ct <- nrow(kn):1
# average length traveled per segment
br$alt <- br$cm / br$ct
kn$alt <- kn$cm / kn$ct
# fluvial inheritance [y m^-1]
mn_d <- mean(c(br$hip, kn$hip)) # mean travel distance [km]
br$fi <- br$alt / (mn_d * 1000 / tp)
kn$fi <- kn$alt / (mn_d * 1000 / tp)

# sternberg mass loss coefficient over distance
# estimating mean travel distance from bear and knowles
# starting weight is d50 of surface gravel at knowles
stm <- log(predict(gmod, newdata = data.frame(b = median(d$b)))) / mn_d
# 0.814 comparable to O'Conner 2014

# mass remaining after turnover period
mr <- exp(log(sum(pebs$wgt)) - st * tp) # 477.58 g remain

# fines produced by measured gravel in turnover period [g] (50g observed)
fn <- sum(pebs$wgt) - mr # 1243.835 g produced
# ratio of fines to gravel
fn / sum(pebs$wgt)

plot(ds$b, log(ds$b), log = 'x')
log(4.5)
log(7.5)
pebs$c <- cdf_hook(pebs$b)
d4 <- pebs[log(pebs$b) > 1.5 & log(pebs$b) <= 1.75, ]
d6 <- pebs[log(pebs$b) > 1.75 & log(pebs$b) <= 2, ]
d9 <- pebs[log(pebs$b) > 2 & log(pebs$b) <= 2.25, ]

(sum(d6$wgt) / (max(d6$c) - min(d6$c))) / (sum(d4$wgt) / (max(d4$c) - min(d4$c)))
(sum(d9$wgt) / (max(d6$c) - min(d6$c))) / (sum(d4$wgt) / (max(d4$c) - min(d4$c)))

plot(pebs$wgt, cdf_hook(pebs$wgt), log ='x')
plot(pebs$b, pebs$c, log = 'x')
plot(ds$b, ds$r, log = 'x')
plot(d$b, d$r, log = 'x')
plot(emp_cdf(d$b)[, 1], emp_cdf(dp)[,2])
points(emp_cdf(d$b)[, 1], emp_cdf(dp)[,2])

dx <- data.frame(b = floor(min(d$b)):ceiling(max(d$b)))
dx$r <- vapply(dx$b, function(x) rarity(x, median(d$b), b4m), 0)
# umbrella distribution
dx$rw <- dx$b / median(d$b) # width relative to median
dx$rw[dx$rw > 1] <- 1 - (dx$rw[dx$rw > 1] - 1) / (max(dx$rw) - 1)
dx$rw[dx$b < median(d$b)] <- 1 / dx$rw[dx$b < median(d$b)]
dx$rw[dx$rw > 1] <- 1 - (dx$rw[dx$rw > 1] - 1) / (max(dx$rw) - 1)

dx$rw <- dx$rw * log(dx$b)^2 / max(log(dx$b)^2)
dx$rw <- dx$rw / max(dx$rw)
plot(dx$b, dx$rw, log = 'x')


dx$p <- predict(m, newdata = dx) # predict cdf
dx$pd <- log(dx$p + 0.5)
dx$pd[dx$pd < 1] <- dx$pd[dx$pd < 1] * (log(1.5) / abs(min(dx$pd)))
dx$pd[dx$pd > 1] <- dx$pd[dx$pd > 1] * (log(1.5) / abs(max(dx$pd)))
dx$pd <- dx$pd - min(dx$pd)
dx$pd <- dx$pd / max(dx$pd)
dx$pc <- emp_cdf(dx$p)[ , 2] # normalize
dx$pm <- to_pmf(dx$pc) # to pmf
plot(dx$b, dx$pd, log = 'x')
lines(emp_cdf(d$b))

synth <- 0
for (i in 1:nrow(d)) {
  synth[i] = min(dx$b[dx$pc >= runif(1)])
}
plot(emp_cdf(synth), log = 'x')
lines(emp_cdf(d$b))

plot(ds$b, ds$r, log = 'xy')
max(ds$r)
sum(pebs$wgt[pebs$b > 6 & pebs$b < 7.5])
sum(pebs$wgt[pebs$b > 4.5 & pebs$b < 6])
plot(d$b, d$r, log = 'xy')

plot(xs, cumsum(rr), log = 'x', ylim = c(0,1))
lines(xs, cumsum(rrs))
lines(xs, cumsum(wr))
lines(emp_cdf(pebs$b))
lines(emp_cdf(surf))
lines(xs, cumsum(wrs))


plot(xs, rw)

plot(emp_cdf(surf_w), log = 'x')
lines(emp_cdf(pebs$b), lwd = 2, col = get_palette('crimson'))
points(emp_cdf(pebs$b), pch = 20, col = get_palette('crimson'))

d50 <- 50
bm$coefficients[2] * log(d50) + bm$coefficients[2] * exp(d50)

pebs$rw <- pebs$wgt / sum(pebs$wgt)
pebs$cw <- cdf_hook(pebs$wgt)
summary(lm(cw ~ log(rw), data = pebs[pebs$wgt != 0, ]))

plot(fish(.25, 1:100, 1) %>% cumsum)

xs <- seq(2.71, 4.3, .01)

lines(exp(seq(2.71, 4.32, .01)), exp(-off + 1.87821 * seq(2.71, 4.32, .01)))


mass_mod <- lm(log(wgt) ~ log(b), data = pebs[pebs$wgt != 0, ])
mass_pred <- predict(mass_mod, newdata = pebs)
summary(mass_mod)
plot(exp(mass_pred), pebs$wgt, log = 'xy')



png('pebble_mass_baxis.png', height = 17, width = 21, units = 'cm', res = 300)
plot(pebs$b, pebs$wgt, log = 'xy',
     pch = 20, col = get_palette('ocean'),
     xlab = 'b-axis width [mm]', ylab = 'mass [g]')
lines(pebs$b, exp(mass_pred), lwd = 2.5, col = get_palette('gold', .5))
text(50, 400, labels = 'R2 = 0.92')
dev.off()

plot(emp_cdf(pebs$wgt), log = 'x',
     xlab = 'Mass [g]', ylab = 'CDF of samples')



