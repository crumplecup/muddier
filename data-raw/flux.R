#assign workspace
setwd('/home/crumplecup/work')

library(data.table)
library(magrittr)
library(rgdal)
library(raster)

crs_ref <- '+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'

uk <- raster('uk.tif')
ukh <- raster('ukh.tif')

#IMPORT VALLEY TRIBUTARY JUNCTION PTS

kn_vol <- read.csv('kn_vol.csv')
flowline <- read.csv('flowline.csv')[,2:3] %>% as.matrix
knowles_samples <- read.csv('knowles_samples.csv')
radio <- read.csv('radio.csv')
gpsing <- read.csv('gpsing.csv')
nodes <- readOGR('nodes_debrisflow.shp')
pdel <- raster('LSDel_sius.tif')


#crop map to upper knowles creek only, it is too big
#extent(uk)
frame_uk <- c(439316.7,442400,4865500,4870600)
crop_uk <- uk %>% crop(frame_uk)
crop_ukh <- ukh %>% crop(frame_uk)
crop_pdel <- pdel %>% crop(frame_uk)

nodes <- nodes %>% crop(frame_uk)
pdel <- pdel %>% crop(frame_uk)

# kn_vol <- fread('kn_vol.csv')
# knowles_samples <- fread('knowles_samples.csv')
# setwd('/home/crumplecup/work/muddier')
# usethis::use_data(kn_vol, overwrite = T)
# usethis::use_data(knowles_samples, overwrite = T)

#geotagged location of site
bear_spot <- c(439290,4869773)	#
grc_spot <- c(429747,4840441)
ccc_spot <- c(426574,4868903)	#
ccf_spot <- c(426555,4868905)
lk_spot <- c(440403,4865453)	#
grc_spot <- c(432159,4838575)
lk2_spot <- c(441479,4868850)	#lk

spots <- c(bear_spot, grc_spot, ccc_spot, ccf_spot,
           lk_spot, grc_spot, lk2_spot)

spotar <- array(spots,dim=c(2,7)) %>% t


labs <- gpsing[,3]%>%as.character
labs[c(15,length(labs))] <- c('32a','32b')


gpsing <- gpsing[,c(3,2,1)]
gpsing$hip <- 0

gpsing[15,4] <- -724
gpsing[16,4] <- -426
gpsing[19,4] <- -340
gpsing[21,4] <- 0

#from lancaster 2006a
lanc <- c(4867939.89,441376.51,793)
lanc <- lanc %>% rbind(c(4866229.08,440854.59,2907))

#from lancaster 2006b
h631 <- c(440662,4867502)
lancs <- lanc[,c(2,1)] %>% rbind(h631)


# extract knowles creek from node network of watershed
nodepath <- nodes
nodes <- nodes[nodes$NODE_ID>300000,]
nodes <- nodes[nodes$NODE_ID<384365,]
knodes <- nodes
nodes <- nodepath

knowles <- knodes
usethis::use_data(knowles, overwrite = T)

# Georef for Heron Creek (LB major trib)  #
knodes@coords[knodes$NODE_ID==384144,] %>% birdseye

# Georef for p1448 RB trib  #
knodes@coords[knodes$NODE_ID==384162,] %>% birdseye

# Georef for Boomer Creek (LB last major trib)  #
knodes@coords[knodes$NODE_ID==384190,] %>% birdseye

#  placing lancaster points  #
lancs[2,1:2] %>% birdseye(100)
text(knodes,knodes$NODE_ID)  #384291

lancs[1,1:2] %>% birdseye(150)
text(knodes,knodes$NODE_ID)  #384108


#build landmarks matrix 'gps'
gps <- gpsing[c(15:16,19,21),]

gps <- gps %>% rbind(c(knodes@coords[knodes$NODE_ID==384108,],'ref793',793))
gps <- gps %>% rbind(c(knodes@coords[knodes$NODE_ID==384144,],'heron',1239))
gps <- gps %>% rbind(c(knodes@coords[knodes$NODE_ID==384162,],'trib1448',1448))
gps <- gps %>% rbind(c(knodes@coords[knodes$NODE_ID==384190,],'boomer',1759))
gps <- gps %>% rbind(c(knodes@coords[knodes$NODE_ID==384291,],'ref2907',2907))
gps <- gps %>% rbind(c(knodes@coords[knodes$NODE_ID==384364,],'init',4000))

gps[,1] <- gps[,1] %>% as.numeric
gps[,2] <- gps[,2] %>% as.numeric

gps$nodeid <- 0
gps$nodeid[1] <- 383982
gps$nodeid[2] <- 384004
gps$nodeid[3] <- 384009
gps$nodeid[4] <- 384041
gps$nodeid[5] <- 384108
gps$nodeid[6] <- 384144
gps$nodeid[7] <- 384162
gps$nodeid[8] <- 384190
gps$nodeid[9] <- 384291
gps$nodeid[10] <- 384364



#pull kilometer to outlet values from linked node network at landmark points

gps$kmtout <- 0

for (i in 1:nrow(gps))	{
  gps$kmtout[i] <- knodes$ToMouth_km[knodes$NODE_ID==gps$nodeid[i]]
}
gps$outhip <- 0
gps$hip <- gps$hip %>% as.numeric
for (i in 2:nrow(gps))	{
  gps$outhip[i] <- (gps$kmtout[i]-gps$kmtout[i-1])*1000/
    (gps$hip[i]-gps$hip[i-1])
}
gps$outhip[10] <- gps$outhip[9]


# id tags
gps$id <- c('m724','m426','m340','ok38','ref793','heron','trib1448',
            'boomer','ref2907','init')


#landmarks <- gps
#usethis::use_data(landmarks, overwrite = T)

#calculate meters to outlet for transects

#based upon the 'outlet to hipchain' ratio
#estimate the distance from nearest landmark point for each transect

tout <- 0
up <- 0
down <- 0
flag <- 0

j <- 2
for (i in 1:length(kn_vol$corrected_hipchain))	{

  while(kn_vol$corrected_hipchain[i] > gps$hip[j]) j <- j+1
  up <- gps$kmtout[j] - ((gps$hip[j] - kn_vol$corrected_hipchain[i]) * gps$outhip[j]/1000)
  down <- gps$kmtout[j-1] + ((kn_vol$corrected_hipchain[i] - gps$hip[j-1]) * gps$outhip[j]/1000)
  flag <- 0
  if((gps$kmtout[j]-up)<(down-gps$kmtout[j-1])) flag <- 1
  if(flag) tout[i] <- up
  if(!flag) tout[i] <- down

}

#snap transects to nearest node
ntout <- 0
for (i in 1:length(tout))	{
  dif <- knodes$ToMouth_km - tout[i]
  ntout[i] <- knodes$NODE_ID[which.min(abs(dif))]
}

#subsect transect nodes
tranodes <- knodes[knodes$NODE_ID%in%ntout,]

#subset nodes within transect length
trans <- knodes[knodes$NODE_ID>=min(tranodes$NODE_ID) & knodes$NODE_ID<=max(tranodes$NODE_ID),]


tranodes$area <- 0  # cross-sectional area in meters
tranodes$contr <- 0  # contributing area
tranodes$top <- 0  # top valley width in meters
tranodes$basal <- 0  # basal valley width in meters

for (i in 1:nrow(tranodes))  {
  tranodes$area[i] <- kn_vol$total_area_m2[i]
  tranodes$contr[i] <- kn_vol$contr_area_km2[i]
  tranodes$top[i] <- kn_vol$top_width_m[i]
  tranodes$basal[i] <- kn_vol$basal_width_m[i]
}

#interpolate between transects
trans$area <- 0  # cross-sectional area in meters
trans$contr <- 0  # contributing area
trans$top <- 0  # top valley width in meters
trans$basal <- 0  # basal valley width in meters

seg <- 0

for (i in 1:(nrow(tranodes)-1))	{
  trans$area[trans$NODE_ID==tranodes$NODE_ID[i]]	<- tranodes$area[tranodes$NODE_ID==tranodes$NODE_ID[i]]
  trans$area[trans$NODE_ID==tranodes$NODE_ID[i+1]]	<- tranodes$area[tranodes$NODE_ID==tranodes$NODE_ID[i+1]]
  trans$contr[trans$NODE_ID==tranodes$NODE_ID[i]]	<- tranodes$contr[tranodes$NODE_ID==tranodes$NODE_ID[i]]
  trans$contr[trans$NODE_ID==tranodes$NODE_ID[i+1]]	<- tranodes$contr[tranodes$NODE_ID==tranodes$NODE_ID[i+1]]
  trans$top[trans$NODE_ID==tranodes$NODE_ID[i]]	<- tranodes$top[tranodes$NODE_ID==tranodes$NODE_ID[i]]
  trans$top[trans$NODE_ID==tranodes$NODE_ID[i+1]]	<- tranodes$top[tranodes$NODE_ID==tranodes$NODE_ID[i+1]]
  trans$basal[trans$NODE_ID==tranodes$NODE_ID[i]]	<- tranodes$basal[tranodes$NODE_ID==tranodes$NODE_ID[i]]
  trans$basal[trans$NODE_ID==tranodes$NODE_ID[i+1]]	<- tranodes$basal[tranodes$NODE_ID==tranodes$NODE_ID[i+1]]

  seg <- tranodes$NODE_ID[i+1] - tranodes$NODE_ID[i]
  for (j in 1:(seg-1))	{
    trans$area[trans$NODE_ID==(tranodes$NODE_ID[i]+j)] <- (
      #percent distance to node
      (trans$ToMouth_km[trans$NODE_ID==(tranodes$NODE_ID[i]+j)] -
         trans$ToMouth_km[trans$NODE_ID==tranodes$NODE_ID[i]]) /
        (trans$ToMouth_km[trans$NODE_ID==tranodes$NODE_ID[i+1]] -
           trans$ToMouth_km[trans$NODE_ID==tranodes$NODE_ID[i]])	) *
      #multiply by change in x-sec area
      (trans$area[trans$NODE_ID==tranodes$NODE_ID[i+1]] -
         trans$area[trans$NODE_ID==tranodes$NODE_ID[i]]) +
      #add starting area
      trans$area[trans$NODE_ID==tranodes$NODE_ID[i]]

    trans$contr[trans$NODE_ID==(tranodes$NODE_ID[i]+j)] <- (
      #percent distance to node
      (trans$ToMouth_km[trans$NODE_ID==(tranodes$NODE_ID[i]+j)] -
         trans$ToMouth_km[trans$NODE_ID==tranodes$NODE_ID[i]]) /
        (trans$ToMouth_km[trans$NODE_ID==tranodes$NODE_ID[i+1]] -
           trans$ToMouth_km[trans$NODE_ID==tranodes$NODE_ID[i]])	) *
      #multiply by change in contr area
      (trans$contr[trans$NODE_ID==tranodes$NODE_ID[i+1]] -
         trans$contr[trans$NODE_ID==tranodes$NODE_ID[i]]) +
      #add starting area
      trans$contr[trans$NODE_ID==tranodes$NODE_ID[i]]

    trans$top[trans$NODE_ID==(tranodes$NODE_ID[i]+j)] <- (
      #percent distance to node
      (trans$ToMouth_km[trans$NODE_ID==(tranodes$NODE_ID[i]+j)] -
         trans$ToMouth_km[trans$NODE_ID==tranodes$NODE_ID[i]]) /
        (trans$ToMouth_km[trans$NODE_ID==tranodes$NODE_ID[i+1]] -
           trans$ToMouth_km[trans$NODE_ID==tranodes$NODE_ID[i]])	) *
      #multiply by change in val
      (trans$top[trans$NODE_ID==tranodes$NODE_ID[i+1]] -
         trans$top[trans$NODE_ID==tranodes$NODE_ID[i]]) +
      #add starting val
      trans$top[trans$NODE_ID==tranodes$NODE_ID[i]]

    trans$basal[trans$NODE_ID==(tranodes$NODE_ID[i]+j)] <- (
      #percent distance to node
      (trans$ToMouth_km[trans$NODE_ID==(tranodes$NODE_ID[i]+j)] -
         trans$ToMouth_km[trans$NODE_ID==tranodes$NODE_ID[i]]) /
        (trans$ToMouth_km[trans$NODE_ID==tranodes$NODE_ID[i+1]] -
           trans$ToMouth_km[trans$NODE_ID==tranodes$NODE_ID[i]])	) *
      #multiply by change in val
      (trans$basal[trans$NODE_ID==tranodes$NODE_ID[i+1]] -
         trans$basal[trans$NODE_ID==tranodes$NODE_ID[i]]) +
      #add starting val
      trans$basal[trans$NODE_ID==tranodes$NODE_ID[i]]

  }
}


# Locate tribs and add 'debris fan effect' to delivery prob

#variable for new debris flow probs
trans$dfp <- 0

trans@coords[11,] %>% birdseye(200)
text(nodes@coords,labels = nodes$NODE_ID)

trans$dfp[11] <- trans$DebrisFlow[11] +
  nodes$DebrisFlow[nodes$NODE_ID==385005] +
  nodes$DebrisFlow[nodes$NODE_ID==385006] +
  nodes$DebrisFlow[nodes$NODE_ID==385007]

trans@coords[15,] %>% birdseye(200)
text(nodes@coords,labels = nodes$NODE_ID)

trans$dfp[15] <- trans$DebrisFlow[15] +
  nodes$DebrisFlow[nodes$NODE_ID==384965] +
  nodes$DebrisFlow[nodes$NODE_ID==384966] +
  nodes$DebrisFlow[nodes$NODE_ID==384967]

trans@coords[30,] %>% birdseye(250)
text(nodes@coords,labels = nodes$NODE_ID)

trans$dfp[30] <- trans$DebrisFlow[30] +
  nodes$DebrisFlow[nodes$NODE_ID==384948] +
  nodes$DebrisFlow[nodes$NODE_ID==384949] +
  nodes$DebrisFlow[nodes$NODE_ID==384950]

trans@coords[48,] %>% birdseye(250)
text(nodes@coords,labels = nodes$NODE_ID)

trans$dfp[48] <- trans$DebrisFlow[48] +
  nodes$DebrisFlow[nodes$NODE_ID==384807] +
  nodes$DebrisFlow[nodes$NODE_ID==384808] +
  nodes$DebrisFlow[nodes$NODE_ID==384809]

trans@coords[117,] %>% birdseye(250)
text(nodes@coords,labels = nodes$NODE_ID)

trans$dfp[117] <- trans$DebrisFlow[117] +
  nodes$DebrisFlow[nodes$NODE_ID==384759] +
  nodes$DebrisFlow[nodes$NODE_ID==384760] +
  nodes$DebrisFlow[nodes$NODE_ID==384761] +
  nodes$DebrisFlow[nodes$NODE_ID==384783] +
  nodes$DebrisFlow[nodes$NODE_ID==384784] +
  nodes$DebrisFlow[nodes$NODE_ID==384785]

trans@coords[151,] %>% birdseye(250)
text(nodes@coords,labels = nodes$NODE_ID)

trans$dfp[151] <- trans$DebrisFlow[151] +
  nodes$DebrisFlow[nodes$NODE_ID==384537] +
  nodes$DebrisFlow[nodes$NODE_ID==384538] +
  nodes$DebrisFlow[nodes$NODE_ID==384539]

trans@coords[168,] %>% birdseye(250)
text(nodes@coords,labels = nodes$NODE_ID)

trans$dfp[168] <- trans$DebrisFlow[168] +
  nodes$DebrisFlow[nodes$NODE_ID==384528] +
  nodes$DebrisFlow[nodes$NODE_ID==384529] +
  nodes$DebrisFlow[nodes$NODE_ID==384530]

trans@coords[195,] %>% birdseye(250)
text(nodes@coords,labels = nodes$NODE_ID)

trans$dfp[195] <- trans$DebrisFlow[195] +
  nodes$DebrisFlow[nodes$NODE_ID==384406] +
  nodes$DebrisFlow[nodes$NODE_ID==384407] +
  nodes$DebrisFlow[nodes$NODE_ID==384408]

trans@coords[220,] %>% birdseye(250)
text(nodes@coords,labels = nodes$NODE_ID)

trans$dfp[220] <- trans$DebrisFlow[220] +
  nodes$DebrisFlow[nodes$NODE_ID==384396] +
  nodes$DebrisFlow[nodes$NODE_ID==384397] +
  nodes$DebrisFlow[nodes$NODE_ID==384398]

trans@coords[233,] %>% birdseye(250)
text(nodes@coords,labels = nodes$NODE_ID)

trans$dfp[233] <- trans$DebrisFlow[233] +
  nodes$DebrisFlow[nodes$NODE_ID==384392] +
  nodes$DebrisFlow[nodes$NODE_ID==384393] +
  nodes$DebrisFlow[nodes$NODE_ID==384394]

# new variable overwriting trib points with new probs, same as M/B otherwise
trans$dfprob <- 0
trans$dfprob[trans$dfp>0] <- trans$dfp[trans$dfp>0]
trans$dfprob[trans$dfp==0] <- trans$DebrisFlow[trans$dfp==0]


setwd('/home/crumplecup/work/muddier/')
transects <- trans
usethis::use_data(transects, overwrite = T)



# weight delivery probs by transect deposit volumes

sum_area <- sum(trans$area)
area_wt <- trans$area / sum_area
area_p <- trans$dfprob * area_wt
rel_p <- area_p / sum(area_p)


# estimate valley sediment volume #

t_dist <- 0  # distance between transects
seg_vol <- 0  # segment volume
for (i in 2:nrow(kn_vol))	{
  t_dist[i] <- kn_vol$corrected_hipchain[i] - kn_vol$corrected_hipchain[i-1]
  seg_vol[i] <- kn_vol$total_area_m2[(i-1):i] %>% sum / 2 * t_dist[i]
}

tot_vol <- seg_vol %>% sum

# debris flow pct by bank height

bank_ht <- (knowles_samples$coarse_fluvial_m +
              knowles_samples$fine_alluvium_m +
              knowles_samples$df_deposit_m) %>% sum

df_ht <- knowles_samples$df_deposit_m %>% sum

# total debris flow volume (volume x pct debris flow)
df_vol <- tot_vol * df_ht / bank_ht


index <- char_pmfs %>% rownames %>% as.numeric %>% sort
cor_mns <- vector(ncol(cor_pmfs), mode = 'numeric')
for (i in 1:ncol(cor_pmfs))  {
  cor_mns[i] <- weighted.mean(index, cor_pmfs[,i])
}

mn_difs <- charcoal$mn - cor_mns

df_difs <- vector(length(mn_difs), mode = 'numeric')
ff_difs <- vector(length(mn_difs), mode = 'numeric')
fg_difs <- vector(length(mn_difs), mode = 'numeric')

df_difs[charcoal$facies == 'DF'] <- mn_difs[charcoal$facies == 'DF']
ff_difs[charcoal$facies == 'FF'] <- mn_difs[charcoal$facies == 'FF']
fg_difs[charcoal$facies == 'FG'] <- mn_difs[charcoal$facies == 'FG']

plot(fg_difs, pch = 20, col = 'coral3')
points(ff_difs, pch = 20, col = 'slateblue')
points(df_difs, pch = 20, col = 'forestgreen')

df_mn_dep <- weighted.mean(index, df_dep)
ff_mn_dep <- weighted.mean(index, ff_dep)
fg_mn_dep <- weighted.mean(index, fg_dep)

# total debris flow flux over mean residence time
df_ann <- df_vol / df_mn_dep

# expected debris flow flux per year
df_flux <- rel_p * df_ann

ave_df1 <- mean(df_flux)
ave_df2 <- df_ann / length(trans)

png('flux_by_delprob.png', width = 11, height = 8.5, res = 300, units = 'in')
plot(trans$DebrisFlow, df_flux, pch = 20, col = 'forestgreen',
     xlab = 'M/B Delivery Probability',
     ylab = 'Sediment Flux [m3/yr]')
abline(lm(df_flux ~ trans$DebrisFlow))
abline(h = ave_df1, lty = 2, col = 'slateblue')
legend('topright', legend = c('Sample Sites', 'Linear Fit', 'Average per Node'), fill = c('forestgreen','black','blue'))
dev.off()


# weighting function for delivery probs

# delivery prob index
p_index <- sort(transects$dfprob)

# cdf of delivery probs
p_cdf <- to_cdf(p_index/sum(p_index))

# change in delivery prob cdf per change in index
delta_p <- vector(length(p_index), mode = 'numeric')
for (i in 2:length(delta_p))  {
  delta_p[i] <- p_cdf[i] - p_cdf[i-1]
}

# cdf of volume by prob index
v_cdf <- vector(length(p_index), mode = 'numeric')
for (i in 2:length(v_cdf))  {
  v_cdf[i] <- sum(transects$area[transects$dfprob <= p_index[i]]) /
    sum(transects$area)
}

# change in volume cdf per change in index
delta_v <- vector(length(p_index), mode = 'numeric')
for (i in 2:length(delta_v))  {
  delta_v[i] <- v_cdf[i] - v_cdf[i-1]
}

# change in volume cdf per change in del prob cdf
p_wt <- delta_v / delta_p

# find ks values in S = ks(A^theta)
lslope <- transects$GRADIENT %>% log  # logged slope
lcontr <- transects$contr %>% log # logged contr area

as_mod <- lm(lslope ~ lcontr)  # fit AS to power law
as_theta <- as_mod$coefficients[2]  # pull theta value

# estimate individual ks values using pulled theta
log_ks <- lslope + as_theta * lcontr
ks <- exp(log_ks)
transects$p_ks <- transects$dfprob / ks

# cdf of ks / pdel ratio over pdel
p_ks_cdf <- vector(length(p_index), mode = 'numeric')
for (i in 2:length(p_ks_cdf))  {
  p_ks_cdf[i] <- sum(transects$p_ks[transects$dfprob <= p_index[i]]) /
    sum(transects$p_ks)
}

# change in p_ks cdf per change in index
delta_p_ks <- vector(length(p_index), mode = 'numeric')
for (i in 2:length(delta_p_ks))  {
  delta_p_ks[i] <- p_ks_cdf[i] - p_ks_cdf[i-1]
}

# change in volume cdf per change in del prob cdf
p_ks_wt <- delta_p_ks / delta_p

p_wt_cdf <- to_cdf(p_wt[-1] / sum(p_wt[-1]))
p_ks_wt_cdf <- to_cdf(p_ks_wt[-1] / sum(p_ks_wt[-1]))


# cdf of top width by prob index
top_cdf <- vector(length(p_index), mode = 'numeric')
for (i in 2:length(top_cdf))  {
  top_cdf[i] <- sum(transects$top[transects$dfprob <= p_index[i]]) /
    sum(transects$top)
}

# change in top width cdf per change in index
delta_top <- vector(length(p_index), mode = 'numeric')
for (i in 2:length(delta_top))  {
  delta_top[i] <- top_cdf[i] - top_cdf[i-1]
}

# change in top width cdf per change in del prob cdf
top_wt <- delta_top / delta_p
top_cdf <- to_cdf(top_wt[-1] / sum(top_wt[-1]))


# cdf of basal width by prob index
basal_cdf <- vector(length(p_index), mode = 'numeric')
for (i in 2:length(basal_cdf))  {
  basal_cdf[i] <- sum(transects$basal[transects$dfprob <= p_index[i]]) /
    sum(transects$basal)
}

# change in basal width cdf per change in index
delta_basal <- vector(length(p_index), mode = 'numeric')
for (i in 2:length(delta_basal))  {
  delta_basal[i] <- basal_cdf[i] - basal_cdf[i-1]
}

# change in basal width cdf per change in del prob cdf
basal_wt <- delta_basal / delta_p
basal_cdf <- to_cdf(basal_wt[-1] / sum(basal_wt[-1]))


plot(p_index[-1], p_ks_wt_cdf, type = 'l', lwd = 2, col = 'forestgreen',
     xlab = 'delivery probability', ylab = 'cdf')
lines(p_index[-1], p_wt_cdf, lwd = 2, col = 'slateblue')
lines(p_index[-1], top_cdf, lwd = 1, col = 'coral3')
lines(p_index[-1], basal_cdf, lwd = 1, col = 'goldenrod')
legend('bottomright', legend = c('volume', 'p/ks', 'top width', 'basal width'),
       fill = c('forestgreen', 'slateblue', 'coral3', 'goldenrod'))





