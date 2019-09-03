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

# Miller/Burnett model output
nodes <- readOGR('nodes_debrisflow.shp')

contr_delta <- function(vec) {
  len <- length(vec)
  delta <- vector(len, mode = 'numeric')
  for (i in 2:len) {
    delta[i] <- vec[i-1] - vec[i]
  }
  delta
}



# identify trib junctions by change in contr area
delta_contr_kn <- contr_delta(kn_vol$contr_area_km2)
delta_contr_br <- contr_delta(br_vol$corr_contr_area_m2)
delta_contr_cd <- contr_delta(cd_vol$corr_contr_area_m2)

bear_ref <- snap(bear_spot, nodes@coords)
ccc_ref <- snap(ccc_spot, nodes@coords)
ccf_ref <- snap(ccf_spot, nodes@coords)
grc_ref <- snap(grc_spot, nodes@coords)
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
kn_id <- nodes$NODE_ID[node_coords[,1] == kn_ref[1] &
                          node_coords[,2] == kn_ref[2] ]
kn2_id <- nodes$NODE_ID[node_coords[,1] == kn2_ref[1] &
                          node_coords[,2] == kn2_ref[2] ]



plot(nodes[nodes$NODE_ID >= 386220 & nodes$NODE_ID <= 386800,])
plot(delta_contr_kn)
plot(delta_contr_br[delta_contr_br > 40000])


plot_pt <- function(pt, mag = 1, so = nodes)  {
  ext <- extent(so)
  scale_x <- (ext[2] - ext[1]) / (2 * mag)
  scale_y <- (ext[4] - ext[3]) / (2 * mag)

  box <- 0
  box[1] <- pt[1] - scale_x
  box[2] <- pt[1] + scale_x
  box[3] <- pt[2] - scale_y
  box[4] <- pt[2] + scale_y

  pic <- raster::crop(so, box)
  plot(pic, pch = 20, col = 'slateblue')
  lines(circ(pt, rad = scale_x * 0.15))
  points(c(pt,pt) %>% matrix(ncol = 2) %>% t,
         pch = 19, col = 'coral3')
}

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

hip_to_out <- function(hipchain, trib_hip, trib_out, out_hip)  {

  out <- 0
  up <- 0
  down <- 0
  flag <- 0

  j <- 2
  for (i in 1:length(hipchain))	{

  while(hipchain[i] > trib_hip[j] & j < length(trib_hip)) j <- j + 1
  up <- trib_out[j] - ((trib_hip[j] - hipchain[i]) * out_hip[j])
  down <- trib_out[j-1] + ((hipchain[i] - trib_hip[j-1]) * out_hip[j])
  flag <- 0
  if((trib_out[j] - up) < (down - trib_out[j-1])) flag <- 1
  if(flag) out[i] <- up
  if(!flag) out[i] <- down

  }
  out
}

#snap transects to nearest node

snap_to_node <- function(dist, outlet, ids)  {
  id <- 0
  for (i in 1:length(dist))	{
    difs <- outlet - dist[i]
    id[i] <- ids[which.min(abs(difs))]
  }
  id
}

bear_out <- hip_to_out(br_vol$outlet_dist_m, bear_tribs$lanc_hip,
                       bear_tribs$ToMouth_km * 1000, bear_tribs$out_hip)

bear_ids <- snap_to_node(bear_out, bear_nodes$ToMouth_km * 1000, bear_nodes$NODE_ID)
bear_ids[20] <- bear_ids[20] + 1

bear_trib_out <- nodes$ToMouth_km[nodes$NODE_ID %in% bear_ids]



interp_by_node <-
  function(vals, trib_out, trib_id, chan_out, chan_id)  {
    res <- vector(length(chan_id), mode = 'numeric')
    j <- 1

    for (i in seq_along(res))  {
      if (chan_id[i] == trib_id[j])  {
        res[i] <- vals[j]
      }

      if (chan_id[i] > trib_id[j] &
          chan_id[i] < trib_id[j + 1])  {
        trib_len <- trib_out[j + 1] - trib_out[j]
        chan_len <- chan_out[i] - trib_out[j]
        pct_len <- chan_len / trib_len
        val_dif <- vals[j + 1] - vals[j]
        res[i] <- vals[j] + (val_dif * pct_len)
      }

      if (chan_id[i] == trib_id[j + 1])  {
        res[i] <- vals[j + 1]
        if ((j + 1) < length(trib_id))
          j <- j + 1
      }
    }
    res
  }



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

bear <- bear_nodes
bear$xsec_area <- bear_xsec_area
bear$contr_area <- bear_contr_area
bear$valley_width <- bear_valley_width

setwd('/home/crumplecup/work/muddier')
usethis::use_data(bear, overwrite = T)



# plot spots
nodes@coords[nodes$NODE_ID == bear_id] %>% plot_pt(300)
nodes@coords[nodes$NODE_ID == ccc_id] %>% plot_pt(100)
nodes@coords[nodes$NODE_ID == ccf_id] %>% plot_pt(100)
nodes@coords[nodes$NODE_ID == grc_id] %>% plot_pt(50)
nodes@coords[nodes$NODE_ID == kn_id] %>% plot_pt(20)
nodes@coords[nodes$NODE_ID == kn2_id] %>% plot_pt(50)


