#assign workspace
setwd('/home/crumplecup/work')

library(data.table)
library(magrittr)
library(rgdal)
library(raster)

crs_ref <- '+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'

# Lancaster survey data
kn_vol <- fread('kn_vol.csv')
br_vol <- fread('bear.csv')
cd_vol <- fread('cedar.csv')
hf_vol <- fread('hoffman_sedvol.csv')
radio <- fread('radiocarbon.csv')
cd_lng <- fread('cedar_long.csv')
gr_vol <- read.csv('grc_vol.csv', skipNul = TRUE)

# Miller/Burnett model output
nodes <- readOGR('nodes_debrisflow.shp')

# identify trib junctions by change in contr area
delta_contr_kn <- contr_delta(kn_vol$contr_area_km2)
delta_contr_br <- contr_delta(br_vol$corr_contr_area_m2)
delta_contr_cd <- contr_delta(cd_vol$corr_contr_area_m2)

#geotagged location of site
bear_spot <- c(439290,4869773)	#
grc_spot <- c(429747,4840441)
ccc_spot <- c(426574,4868903)	#
ccf_spot <- c(426555,4868905)
lk_spot <- c(440403,4865453)	#
grc2_spot <- c(432159,4838575)
lk2_spot <- c(441479,4868850)	#lk

bear_ref <- snap(bear_spot, nodes@coords)
ccc_ref <- snap(ccc_spot, nodes@coords)
ccf_ref <- snap(ccf_spot, nodes@coords)
grc_ref <- snap(grc_spot, nodes@coords)
grc2_ref <- snap(grc2_spot, nodes@coords)
kn_ref <- snap(lk_spot, nodes@coords)
kn2_ref <- snap(lk2_spot, nodes@coords)

node_coords <- nodes@coords

bear_id <- nodes$NODE_ID[node_coords[,1] == bear_ref[1] &
                           node_coords[,2] == bear_ref[2] ]
ccc_id <- nodes$NODE_ID[node_coords[,1] == ccc_ref[1] &
                           node_coords[,2] == ccc_ref[2] ]
ccf_id <- nodes$NODE_ID[node_coords[,1] == ccf_ref[1] &
                           node_coords[,2] == ccf_ref[2] ]
grc_id <- nodes$NODE_ID[node_coords[,1] == grc_ref[1] &
                           node_coords[,2] == grc_ref[2] ]
grc2_id <- nodes$NODE_ID[node_coords[,1] == grc2_ref[1] &
                          node_coords[,2] == grc2_ref[2] ]
kn_id <- nodes$NODE_ID[node_coords[,1] == kn_ref[1] &
                          node_coords[,2] == kn_ref[2] ]
kn2_id <- nodes$NODE_ID[node_coords[,1] == kn2_ref[1] &
                          node_coords[,2] == kn2_ref[2] ]



plot(nodes[nodes$NODE_ID >= 386220 & nodes$NODE_ID <= 386800,])
plot(delta_contr_kn)
plot(delta_contr_br[delta_contr_br > 40000])



spot <- nodes@coords[nodes$NODE_ID == 386219]  # bear outlet to knowles
spot <- nodes@coords[nodes$NODE_ID == 386239]  # trib B1L lanc 224
spot <- nodes@coords[nodes$NODE_ID == 386243]  # trib B2R lanc 270
spot <- nodes@coords[nodes$NODE_ID == 386258]  # trib B3L lanc 445
spot <- nodes@coords[nodes$NODE_ID == 386269]  # trib B4L lanc 631
spot <- nodes@coords[nodes$NODE_ID == 386278]  # trib B5R lanc 682
spot <- nodes@coords[nodes$NODE_ID == 386299]  # trib B6L lanc 942
spot <- nodes@coords[nodes$NODE_ID == 386314]  # trib B7R lanc 1191
spot <- nodes@coords[nodes$NODE_ID == 386331]  # trib B8L lanc 1377
spot <- nodes@coords[nodes$NODE_ID == 386348]  # trib B9R lanc 1572
spot <- nodes@coords[nodes$NODE_ID == 386359]  # trib B10L lanc 1694
spot <- nodes@coords[nodes$NODE_ID == 386378]  # trib B11L lanc 1902
spot <- nodes@coords[nodes$NODE_ID == 386389]  # trib B12R lanc 2061
spot <- nodes@coords[nodes$NODE_ID == 386409]  # trib B13R lanc 2284
spot <- nodes@coords[nodes$NODE_ID == 386440]  # trib B14R lanc 2544
spot <- nodes@coords[nodes$NODE_ID == 386452]  # init lanc 2658
spot %>% plot_pt(200)

bear_nodes <- nodes[nodes$NODE_ID %in% 386219:386452,]

trib_ids <- c(386219, 386239, 386243, 386258, 386269, 386278, 386299, 386314,
              386331, 386348, 386359, 386378, 386389, 386409, 386440, 386452)
bear_lanc_hip <- c(0, 224, 270, 445, 631, 682, 942, 1191, 1377, 1572, 1694,
               1902, 2061, 2284, 2544, 2658)

bear_tribs <- nodes[nodes$NODE_ID %in% trib_ids, ]
bear_tribs$lanc_hip <- bear_lanc_hip

bear_tribs$out_hip <- 0

for (i in 2:nrow(bear_tribs))	{
  bear_tribs$out_hip[i] <- ((bear_tribs$ToMouth_km[i] -
      bear_tribs$ToMouth_km[i-1]) * 1000) / (
      bear_tribs$lanc_hip[i] - bear_tribs$lanc_hip[i-1])
}


#based upon the 'outlet to hipchain' ratio
#estimate the distance from nearest landmark point for each transect
#snap transects to nearest node


bear_out <- hip_to_out(br_vol$outlet_dist_m, bear_tribs$lanc_hip,
                       bear_tribs$ToMouth_km * 1000, bear_tribs$out_hip)

bear_ids <- snap_to_node(bear_out, bear_nodes$ToMouth_km * 1000, bear_nodes$NODE_ID)
bear_ids[20] <- bear_ids[20] + 1

bear_trib_out <- nodes$ToMouth_km[nodes$NODE_ID %in% bear_ids]


bear_xsec_raw <- unlist(br_vol[,10])
bear_xsec_raw[is.na(bear_xsec_raw)] <- 0

bear_xsec_area <- interp_by_node(vals = bear_xsec_raw,
                                 trib_out = bear_trib_out,
                                 trib_id = bear_ids,
                                 chan_out = bear_nodes$ToMouth_km,
                                 chan_id = bear_nodes$NODE_ID)

bear_contr_area <- interp_by_node(vals = br_vol$corr_contr_area_m2,
                                 trib_out = bear_trib_out,
                                 trib_id = bear_ids,
                                 chan_out = bear_nodes$ToMouth_km,
                                 chan_id = bear_nodes$NODE_ID)

bear_valley_width_raw <- br_vol$valley_width_m
bear_valley_width_raw[c(49:50, 52:53, 57:58, 60)] <- 0
bear_valley_width_raw <- bear_valley_width_raw %>% as.numeric

bear_valley_width <- interp_by_node(vals = bear_valley_width_raw,
                                 trib_out = bear_trib_out,
                                 trib_id = bear_ids,
                                 chan_out = bear_nodes$ToMouth_km,
                                 chan_id = bear_nodes$NODE_ID)

bear_slope <- interp_by_node(vals = br_vol$local_avg_bed_slope,
                                    trib_out = bear_trib_out,
                                    trib_id = bear_ids,
                                    chan_out = bear_nodes$ToMouth_km,
                                    chan_id = bear_nodes$NODE_ID)

bear <- bear_nodes
bear$xsec_area <- bear_xsec_area
bear$contr_area <- bear_contr_area
bear$valley_width <- bear_valley_width
bear$slope <- bear_slope

setwd('/home/crumplecup/work/muddier')
usethis::use_data(bear, overwrite = T)


# Cedar Creek
# C1.5R 418461
# C1.5R mouth 419085
# RB trib 419112 Lanc 275
# RB trib 419131 Lanc 479
# LB trib 419154 Lanc 726
# LB trib 419164 Lanc 852
# LB trib 419187 Lanc 1111
# init 419240 (lanc 1703 +13 dif at lanc 1111)

ced_nodes <- nodes[nodes$NODE_ID %in% 419085:419240,]

ced_trib_ids <- c(419085, 419112, 419131, 419154, 419164, 419187, 419240)
ced_lanc_hip <- c(0, 275, 479, 726, 852, 1111, 1717)

ced_tribs <- nodes[nodes$NODE_ID %in% ced_trib_ids, ]
ced_tribs$lanc_hip <- ced_lanc_hip

ced_tribs$out_hip <- 0

for (i in 2:nrow(ced_tribs))	{
  ced_tribs$out_hip[i] <- ((ced_tribs$ToMouth_km[i] -
                               ced_tribs$ToMouth_km[i-1]) * 1000) / (
                                 ced_tribs$lanc_hip[i] - ced_tribs$lanc_hip[i-1])
}


ced_out <- hip_to_out(unlist(cd_lng[,2]), ced_tribs$lanc_hip,
                       ced_tribs$ToMouth_km * 1000, ced_tribs$out_hip)

ced_ids <- snap_to_node(ced_out, ced_nodes$ToMouth_km * 1000, ced_nodes$NODE_ID)
# screen out oversampling within node sections
# 4-5, 10-11, 17-18, 24-25, 39-40, 43-44, 45-46, 47-48, 56-57, 58-59
# 62-63, 67-68, 74-75, 79-80, 82-83, 84-85, 86-88, 92-93, 102-103,
# 104-105, 106-107, 111-112
# compare outlet distances to find transect closest to node (drop other)
over_list <- list(c(4:5), c(10:11), c(17:18), c(24:25), c(39:40), c(43:44),
                  c(45:46), c(47:48), c(56:57), c(58:59), c(62:63), c(67:68),
                  c(74:75), c(79:80), c(82:83), c(84:85), c(86:88), c(92:93),
                  c(102:103), c(104:105), c(106:107), c(111:112))

difs <- 0
for (i in seq_along(ced_ids)) {
  difs[i] <- (nodes$ToMouth_km[nodes$NODE_ID == ced_ids[i]]*1000) - ced_out[i]
}

min_ids <- c(4, 10, 17, 25, 39, 44, 46, 48, 56, 59, 62, 68,
             74, 80, 83, 85, 88, 92, 103, 104, 106, 112)

# remove ids closest to node sites from oversample list
rem_ids <- 0
for (i in seq_along(min_ids)) {
  rem <- over_list[[i]]
  rem_ids <- c(rem_ids, rem[!rem %in% min_ids[i]])
}
rem_ids <- rem_ids[-1]

cd_trc <- cd_lng[-rem_ids]

ced_out_lng <- hip_to_out(unlist(cd_trc[,2]), ced_tribs$lanc_hip,
                      ced_tribs$ToMouth_km * 1000, ced_tribs$out_hip)

ced_ids_lng <- snap_to_node(ced_out_lng, ced_nodes$ToMouth_km * 1000, ced_nodes$NODE_ID)

# local average bed slope m
ced_bed <- cd_trc$chan_bed_elev_m
ced_hip <- cd_trc$corr_hipchain_m
ced_n <- length(ced_bed)
ced_slp <- vector(length = ced_n, mode = 'numeric')
for (i in seq_along(ced_bed)) {
  if (i == 1) {
    ced_slp[i] <- (ced_bed[i + 1] - ced_bed[i]) / (ced_hip[i + 1] - ced_hip[i])
  }
  if (i > 1 & i < ced_n) {
    slp_a <- (ced_bed[i] - ced_bed[i - 1]) / (ced_hip[i] - ced_hip[i - 1])
    slp_b <- (ced_bed[i + 1] - ced_bed[i]) / (ced_hip[i + 1] - ced_hip[i])
    ced_slp[i] <- (slp_a + slp_b) / 2
  }
  if (i == ced_n) {
    ced_slp[i] <- (ced_bed[i] - ced_bed[i - 1]) / (ced_hip[i] - ced_hip[i - 1])
  }
}

ced_out <- hip_to_out(unlist(cd_vol[,2]), ced_tribs$lanc_hip,
                      ced_tribs$ToMouth_km * 1000, ced_tribs$out_hip)

ced_ids <- snap_to_node(ced_out, ced_nodes$ToMouth_km * 1000, ced_nodes$NODE_ID)
ced_ids <- c(ced_ids, 419240)

ced_trib_out <- nodes$ToMouth_km[nodes$NODE_ID %in% ced_ids]

ced_xsec_raw <- c(unlist(cd_vol$xsec_area_m2), 0)
ced_xsec_area <- interp_by_node(vals = ced_xsec_raw,
                                 trib_out = ced_trib_out,
                                 trib_id = ced_ids,
                                 chan_out = ced_nodes$ToMouth_km,
                                 chan_id = ced_nodes$NODE_ID)

ced_contr_raw <- c(unlist(cd_vol$corr_contr_area_m2), 10000)
ced_contr_area <- interp_by_node(vals = ced_contr_raw,
                                  trib_out = ced_trib_out,
                                  trib_id = ced_ids,
                                  chan_out = ced_nodes$ToMouth_km,
                                  chan_id = ced_nodes$NODE_ID)

ced_valley_raw <- c(unlist(cd_vol[,8]), 1)
ced_valley_width <- interp_by_node(vals = ced_valley_raw,
                                    trib_out = ced_trib_out,
                                    trib_id = ced_ids,
                                    chan_out = ced_nodes$ToMouth_km,
                                    chan_id = ced_nodes$NODE_ID)

ced_slope <- interp_by_node(vals = ced_slp,
                            trib_out = ced_trib_out,
                            trib_id = ced_ids,
                            chan_out = ced_nodes$ToMouth_km,
                            chan_id = ced_nodes$NODE_ID)

cedar <- ced_nodes
cedar$xsec_area <- ced_xsec_area
cedar$contr_area <- ced_contr_area
cedar$valley_width <- ced_valley_width
cedar$slope <- ced_slope

setwd('/home/crumplecup/work/muddier')
usethis::use_data(cedar, overwrite = T)


# Hoffman Creek

# H1.7L 429008
# H1.7L mouth 431077
# RB trib 431104 Lanc 342
# RB trib 431116 Lanc 462
# LB trib 431146 Lanc 856
# RB trib 431155 Lanc 1013
# LB trib 431174 Lanc 1203
# LB trib 431179 Lanc 1302
# LB trib 431194 Lanc 1468
# LB trib 431219 Lanc 1736
# LB trib 431228 Lanc 1857
# RB trib 431240 Lanc 1971 "bloody finger trib"
# RB trib 431252 Lanc 2103
# LB trib 431269 Lanc 2304
# RB trib 431284 Lanc 2449
# init 431290  Lanc (2449 - 2297.26 + 2369.84 = 2521.58)

hoff_nodes <- nodes[nodes$NODE_ID %in% 431077:431290,]

hoff_trib_ids <- c(431077, 431104, 431116, 431146, 431155, 431174, 431179,
                   431194, 431219, 431228, 431240, 431252, 431269, 431284, 231290)
hoff_lanc_hip <- c(0, 342, 462, 856, 1013, 1203, 1302, 1468, 1736, 1857,
                   1971, 2103, 2304, 2449, 2522)

hoff_tribs <- nodes[nodes$NODE_ID %in% hoff_trib_ids, ]
hoff_tribs <- hoff_tribs[order(hoff_tribs$ToMouth_km),]
hoff_tribs$lanc_hip <- hoff_lanc_hip

hoff_tribs$out_hip <- 0

for (i in 2:nrow(hoff_tribs))	{
  hoff_tribs$out_hip[i] <- ((hoff_tribs$ToMouth_km[i] -
                               hoff_tribs$ToMouth_km[i-1]) * 1000) / (
                                 hoff_tribs$lanc_hip[i] - hoff_tribs$lanc_hip[i-1])
}


hoff_out <- hip_to_out(hf_vol$corr_hipchain_m, hoff_tribs$lanc_hip,
                       hoff_tribs$ToMouth_km * 1000, hoff_tribs$out_hip)

hoff_ids <- snap_to_node(hoff_out, hoff_nodes$ToMouth_km * 1000, hoff_nodes$NODE_ID)
hoff_ids[26] <- hoff_ids[26] - 1  # correct where two transects mapped to same node
hoff_ids[31] <- hoff_ids[31] - 1
hoff_ids[36] <- hoff_ids[36] + 1
hoff_ids[41] <- hoff_ids[41] + 1
hoff_ids <- c(hoff_ids, 431290)

hoff_trib_out <- nodes$ToMouth_km[nodes$NODE_ID %in% hoff_ids]

hoff_xsec_raw <- c(unlist(hf_vol[,13]), 0)
hoff_xsec_area <- interp_by_node(vals = hoff_xsec_raw,
                                trib_out = hoff_trib_out,
                                trib_id = hoff_ids,
                                chan_out = hoff_nodes$ToMouth_km,
                                 chan_id = hoff_nodes$NODE_ID)

hoff_contr_raw <- c(unlist(hf_vol[,5]), 10000)
hoff_contr_area <- interp_by_node(vals = hoff_contr_raw,
                                 trib_out = hoff_trib_out,
                                 trib_id = hoff_ids,
                                 chan_out = hoff_nodes$ToMouth_km,
                                 chan_id = hoff_nodes$NODE_ID)

hoff_valley_raw <- c(unlist(hf_vol[,7]), .5)
hoff_valley_width <- interp_by_node(vals = hoff_valley_raw,
                                    trib_out = hoff_trib_out,
                                    trib_id = hoff_ids,
                                    chan_out = hoff_nodes$ToMouth_km,
                                    chan_id = hoff_nodes$NODE_ID)

# local average bed slope
hoff_slp <- interp_by_node(vals = hf_vol$reach_slope,
                          trib_out = hoff_trib_out,
                          trib_id = hoff_ids,
                          chan_out = hoff_nodes$ToMouth_km,
                          chan_id = hoff_nodes$NODE_ID)


hoffman <- hoff_nodes
hoffman$xsec_area <- hoff_xsec_area
hoffman$contr_area <- hoff_contr_area
hoffman$valley_width <- hoff_valley_width
hoffman$slope <- hoff_slp

setwd('/home/crumplecup/work/muddier')
usethis::use_data(hoffman, overwrite = T)


# knowles

# mouth 383786
# RB trib 383799
# RB trib 383833
# RB trib 383854
# LB trib 383866 bear mouth
# LB trib 383892 lanc -724 K32
# RB trib 383907
# RB trib 383928
# LB trib 383955
# RB trib 383960
# RB trib 383965
# LB trib 383982
# RB trib 384005 lanc -426 K33
# RB trib 384009 lanc -340 K36
# RB trib 384024
# LB trib 384041 lanc 0 K38
# LB & RB trib 384111
# LB trib 384145 lanc 1239 'heron'
# RB trib 384162 lanc 1448
# LB trib 384190 lanc 1759 'boomer'
# RB trib 384214
# RB trib 384228
# RB trib 384327
# init 384864

knowles_nodes <- nodes[nodes$NODE_ID %in% 383786:384864,]
knowles_tribs <- nodes[nodes$NODE_ID %in% landmarks$nodeid,]

knowles_out <- hip_to_out(kn_vol$corrected_hipchain, landmarks$hip,
                       knowles_tribs$ToMouth_km * 1000, landmarks$outhip)

knowles_ids <- snap_to_node(knowles_out, knowles_nodes$ToMouth_km * 1000, knowles_nodes$NODE_ID)
knowles_ids[97] <- knowles_ids[97] + 1

kn_trib_out <- nodes$ToMouth_km[nodes$NODE_ID %in% knowles_ids]

kn_sub <- transects[,-c(11:13)]
kn_sub$slope <- interp_by_node(vals = kn_vol$local_slope,
                               trib_out = kn_trib_out,
                               trib_id = knowles_ids,
                               chan_out = kn_sub$ToMouth_km,
                               chan_id = kn_sub$NODE_ID)


# combine transect data for all creeks

bear$creek_name <- 'bear'
cedar$creek_name <- 'cedar'
hoffman$creek_name <- 'hoffman'
kn_sub$creek_name <- 'knowles'

names(kn_sub) <- names(bear)
kn_sub$contr_area <- kn_sub$contr_area * 1000000 # convert km2 to m2
creeks <- rbind(bear, rbind(cedar, rbind(hoffman, kn_sub)))
creeks$contr_area <- creeks$contr_area / 1000000 # convert m2 to km2

setwd('/home/crumplecup/work/muddier')
usethis::use_data(creeks, overwrite = T)


# plot spots
nodes@coords[nodes$NODE_ID == bear_id] %>% plot_pt(300)
nodes@coords[nodes$NODE_ID == ccc_id] %>% plot_pt(100)
nodes@coords[nodes$NODE_ID == ccf_id] %>% plot_pt(100)
nodes@coords[nodes$NODE_ID == grc_id] %>% plot_pt(50)
nodes@coords[nodes$NODE_ID == kn_id] %>% plot_pt(20)
nodes@coords[nodes$NODE_ID == kn2_id] %>% plot_pt(50)


# all creek transects
creek_transects <- nodes[nodes$NODE_ID %in% c(
  bear_ids, ced_ids_lng, hoff_ids, knowles_ids
), ]
creek_transects$creek_name <- 0
creek_transects$creek_name[creek_transects$NODE_ID %in% bear_ids] <- 'bear'
creek_transects$creek_name[creek_transects$NODE_ID %in% ced_ids_lng] <- 'cedar'
creek_transects$creek_name[creek_transects$NODE_ID %in% hoff_ids] <- 'hoffman'
creek_transects$creek_name[creek_transects$NODE_ID %in% knowles_ids] <- 'knowles'

setwd('/home/crumplecup/work/muddier')
usethis::use_data(creek_transects)
setwd('/home/crumplecup/work')
writeOGR(creek_transects, work_dir, 'creek_transects', driver = 'ESRI Shapefile')

# radiocarbon georefs

radio <- radio[,-5]

# bear creek subset
bear_radio <- radio[radio$location %in% c('BCL', 'BCU'), ]
bear_radio <- bear_radio[order(bear_radio$distance_m),]

bear_radio_out <- hip_to_out(bear_radio$distance_m, bear_tribs$lanc_hip,
                       bear_tribs$ToMouth_km * 1000, bear_tribs$out_hip)

bear_radio_ids <- snap_to_node(bear_radio_out, bear_nodes$ToMouth_km * 1000, bear_nodes$NODE_ID)

bear_radio$node_ids <- bear_radio_ids
bear_radio$creek_name <- 'bear'

# cedar creek subset

cedar_radio <- radio[radio$loc %in% c('CMS', 'CC1262', 'CC1282'), ]
cedar_radio <- cedar_radio[order(cedar_radio$distance_m),]

cedar_radio_out <- hip_to_out(cedar_radio$distance_m, ced_tribs$lanc_hip,
                             ced_tribs$ToMouth_km * 1000, ced_tribs$out_hip)

cedar_radio_ids <- snap_to_node(cedar_radio_out, ced_nodes$ToMouth_km * 1000,
                                ced_nodes$NODE_ID)

cedar_radio$node_ids <- cedar_radio_ids
cedar_radio$creek_name <- 'cedar'



# knowles creek subset
knowles_radio <- radio[radio$location %in% c('KCU', 'KCM', 'KCL'),]
knowles_radio <- knowles_radio[order(knowles_radio$distance_m),]

knowles_radio_out <- hip_to_out(knowles_radio$distance_m, landmarks$hip,
                                knowles_tribs$ToMouth_km * 1000, landmarks$outhip)

knowles_radio_ids <- snap_to_node(knowles_radio_out, knowles_nodes$ToMouth_km * 1000,
                                  knowles_nodes$NODE_ID)

knowles_radio$node_ids <- knowles_radio_ids
knowles_radio$creek_name <- 'knowles'


# combine subsets into creeks object

creeks_radio <- rbind(bear_radio, rbind(cedar_radio, knowles_radio))

setwd('/home/crumplecup/work/muddier')
usethis::use_data(creeks_radio, overwrite = T)
setwd('/home/crumplecup/work')

creeks_radiocarbon <- nodes[nodes$NODE_ID %in% creeks_radio$node_ids, ]
writeOGR(creeks_radiocarbon, work_dir, 'creeks_radiocarbon',
         driver = 'ESRI Shapefile')



