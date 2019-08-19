#clean workspace
rm(list=ls(all=TRUE))


#assign workspace
#setwd('C:/Users/Erik Rose/Documents/')
setwd('/home/crumplecup/Documents/sediment')

source('crumplecup_3.5.0.R')
crs_ref <- '+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0'
pal <- brewer.pal(11,'Spectral')
pal <- pal[c(1,3,10,11)]


# cd C:\Program Files\R\R-3.4.3\bin\x64
# Rgui.ege --max-ppsize=500000

#IMPORT LIDAR DATA FOR MAPLETON#

#setwd('C:/Users/Erik Rose/Documents/sed lidar/')

uk <- raster('uk.tif')
ukh <- raster('ukh.tif')
#bear <- raster('bear.tif')
#bearh <- raster('bearh.tif')
#lk <- raster('lk.tif')
#lkh <- raster('lkh.tif')

#IMPORT VALLEY TRIBUTARY JUNCTION PTS

kn_vol <- read.csv('kn_vol.csv')
flowline <- read.csv('flowline.csv')
flowline <- flowline[,2:3] %>% as.matrix

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

gpsing[15,3]
gpsing[26,4] <- -724
gpsing[16,4] <- -426
gpsing[19,4] <- -340
gpsing[21,4] <- 0

#from lancaster 2006a
lanc <- c(4867939.89,441376.51,793)
lanc <- lanc %>% rbind(c(4866229.08,440854.59,2907))

#from lancaster 2006b
h631 <- c(440662,4867502)
lancs <- lanc[,c(2,1)] %>% rbind(h631)


# plot of lancaster 2006A 'gpsing knowles'  #
#png('gpsing_knowles.png')
frame <- extent(gpsing[,1:2]) + c(-1,1,-2,2)*200
croppy <- crop_uk %>% crop(frame)
plot(croppy)
gps_labs <- gpsing[,1:2] + c(-1.5,-1)*40
points(gpsing[,1:2],pch=19,col='orange')
text(gps_labs,labs)
#dev.off()




# plot lancaster gps refs for upper knowles and heron creek #
#png('lanc_uk_gps.png')
frame <- extent(lancs) + c(-2,2,-1,1)*200
crop_uk %>% crop(frame) %>% plot
points(lancs,pch=19,col='blue')
lanc_cords <- lancs + c(-2,1)*50
text(lanc_cords,c(lanc[,3],'Heron'))
#dev.off()

nodepath <- nodes
nodes <- nodepath
nodes <- nodes[nodes$NODE_ID>300000,]
nodes <- nodes[nodes$NODE_ID<384366,]
knodes <- nodes
nodes <- nodepath

density(knodes$DebrisFlow) %>% plot




#transect tools

#triangle length

tri_length <- function(pt,bear,dist,north=TRUE)	{
#given hypotenuse length of right triangle
#and coords of start point
#returns coords of end point along hypotenuse

	if(bear>=90 & bear<270)	{
		C <- abs(bear-180)
		north <- FALSE
		}

	if(bear<90) C <- bear
	if(bear>=270) C <- 360-bear
	a <- (dist*sin((90-C)*pi/180))/sin(90*pi/180)
	if(north) y <- pt[2] + a
	if(!north) y <- pt[2] - a

	c <- (dist*sin(C*pi/180))/sin(90*pi/180)
	if(bear<=180) x <- pt[1] + c
	if(bear>180) x <- pt[1] - c

	pt <- c(x,y)
	return(pt)
	}




#increment pathway

inc_path <- function(coords)	{
#given a matrix of coordinates representing a path or line
#take the bearing and length of each segment
#divide the segment into length-1 points
#for each point along the segment, calculate the x and y position
#return a new path consisting containing all old and new x and y values
#the purpose of this function is to increment along the channel line
#so I can take more transects and get a finer resolution of detail on the channel

	inc <- c(0,0) %>% matrix(ncol=2)
	for (i in 1:(nrow(coords)-1))	{
		inc <- inc %>% rbind(coords[i,])
		bear <- atan2(coords[i+1,1]-coords[i,1],coords[i+1,2]-coords[i,2])
		if(bear<0) bear <- 2*pi + bear
		bear <- bear * (180/pi)
		if(bear>360) bear <- bear - 360
		seg <- sqrt((coords[i+1,1]-coords[i,1])^2 +
					(coords[i+1,2]-coords[i,2])^2) %>% floor
		for (j in 1:(seg-1))	{
			pt <- coords[i,] %>% tri_length(bear,j)
			inc <- inc %>% rbind(pt)
			}
		}
	inc <- inc %>% rbind(coords[nrow(coords),])
	inc <- inc[-1,]
	return(inc)
	}






# take a transect

take_transect <- function(pt,bear,rad)	{
#given a point pt, a bearing bear and a radius rad
#take a transect centered on point pt, along bearing bear with radius rad

	start <- pt %>% tri_length(bear,rad)
	bear <- bear + 180
	if(bear>=360) bear <- bear - 360
	stop <- pt %>% tri_length(bear,rad)
	pt <- start %>% rbind(stop) %>% inc_path
	return(pt)
	}



# find the low path


low_path <- function(coords,rad=25,map=uk)	{
	path <- c(0) %>% matrix(nrow=1,ncol=2)
	for(i in 1:(nrow(coords)-1))	{
		bear <- atan2(coords[i+1,1]-coords[i,1],coords[i+1,2]-coords[i,2])
		if (bear<0) bear <- bear + 2*pi
		bear <- bear * (180/pi) + 90
		if(bear>360) bear <- bear-360
#		print(bear)
		seg <- coords[i:(i+1),] %>% inc_path
		for(j in 1:nrow(seg))	{
			low <- seg[j,] %>% take_transect(bear,rad)
			low <- low %>% cbind(extract(map,low))
#			print(low)
#			print(low[low[,3]==min(low[,3])])
			lowpt <- low[low[,3]==min(low[,3])][1:2]
			path <- path %>% rbind(lowpt)
			}
		}
	path <- path[-1,]
	return(path)
	}



# DRAW CIRCLE

circ <- function(pt,rad=15)	{
	mat <- c(0,0) %>% matrix(ncol=2)
	for (i in 1:360)	{
		mat <- mat %>% rbind(tri_length(pt,i,rad))
		}
	mat <- mat[-1,]
	return(mat)
	}


# HIPCHAIN #

hipchain <- function(pt,dist,inc=10,chan=flowline,us=T,filter=4000)	{
#travel along channel a set increment from a pt
#pt is a coord in crs_ref
#dist is total travel distance in meters along channel
#inc is in meters, the length in which to divide travel distance
#chan is a matrix of coords delineating channel from downstream to upstream
#if us=T distance is upstream, downstream when F
#filter removes coords off the channel path (errors from low_path I suspect)

	dif <- ((chan[,1] - pt[1])^2 + (chan[,2] - pt[2])^2) %>% sqrt
	pt <- chan[dif<=inc,]
#	print(pt%>%dim)
	dmax <- dif[dif<=inc]
#	print(dmax)

	if(us)	{
		pt <- pt[nrow(pt),]
		dmax <- dmax[length(dmax)]
		}
	if(!us)	{
		pt <- pt[1,]
		dmax <- dmax[1]
		}

#	print(dmax)
	packet <- list(pt,dist-dmax)

#	print(packet[[1]] %>% matrix(ncol=2))
	print(packet[[2]])
	while(packet[[2]]>inc)	{
		packet <- c(packet[[1]]) %>% hipchain(packet[[2]],inc,chan,us,filter)
		}
	if(packet[[2]]>1.5)	{
		packet <- c(packet[[1]]) %>% hipchain(packet[[2]],packet[[2]],chan,us,filter)
		}
	return(packet)
	}


#  BIRDSEYE  #


birdseye <- function(pt,mag=40,flow=flowline,nod=nodes,map=uk)	{
	frame <- map %>% zoom(pt,mag)
	crop <- map %>% crop(frame)
	plot(crop)
	points(flow,pch=20,col='blue')
	lines(pt %>% circ)
	points(nod,col='red')
	points(c(pt,pt)%>%matrix(ncol=2)%>%t,pch=19,col='red')
	}



est_hipchain <- function(pt,dist,chan=flowline,us=T)	{
	mat <- c(0,0) %>% matrix(ncol=2)
	for (i in 1:4)	{
		packet <- hipchain(pt,dist,inc=i*10,chan,us)
		mat <- mat %>% rbind(packet[[1]])
		}
	mat <- mat[-1,]
	return(mat)
	}


#  TRANSECTIFY  #

transectify <- function(pt,bear,rad=50,mat,map=uk)	{
	sbear <- bear + 90
	if(sbear>=360) sbear <- sbear - 360
	start <- pt %>% tri_length(sbear,10)
	ebear <- bear - 90
	if(ebear<0) ebear <- ebear + 360
	end <- pt %>% tri_length(ebear,10)
	path <- c(start,end) %>% matrix(ncol=2) %>% t %>% inc_path
	map %>% extract(path[11,] %>% take_transect(bear,rad)) %>% FtoM %>% plot(
		main='Transect',xlab='Valley Floor',ylab='Elevation(m)')
	for (i in 1:4)	{
		map %>% extract(path[i*6-5,] %>% take_transect(bear,rad)) %>% FtoM %>% lines(col=pal[i])
		}
	lines(mat,lwd=2)
	return(path)
	}

#  SNAP  #


snap <- function(pt,line=nodepath){
	dif <- ((line[,1] - pt[1])^2 + (line[,2] - pt[2])^2) %>% sqrt
	pt <- line[dif==min(dif),]
	return(pt)
	}



get_bear <- function(coords)	{
	bear <- atan2(coords[2,1]-coords[1,1],coords[2,2]-coords[1,2])
	if(bear<0) bear <- 2*pi + bear
	bear <- bear * (180/pi)
	if(bear>360) bear <- bear - 360
	return(bear)
	}



#  DIFFER  #

differ <- function(vec,delta=0)	{

#given a vector > length 2
#return a variable of differences between elements

	for (i in 2:length(vec))	{
		delta[i] <- (vec[i] - vec[i-1])
		}
	return(delta)
	}


# CUMULT #

cumult <- function(vec,dif,cum=0) {

# given a vector > length 2
# return a cumulative distribution curve

	for (i in 1:length(vec))	{
		cum[i] <- sum(vec[1:i]) / dif
		}
	return(cum)
	}






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


gps$hip[1] <- -724


gps[5,1:2] %>% as.numeric %>% birdseye(12)
knodes@coords[knodes$NODE_ID==384090,] %>% birdseye(12)

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


#calculate meters to outlet for transects

#based upon the 'outlet to hipchain' ratio
#estimate the distance from nearest landmark point for each transect

tout <- 0
up <- 0
down <- 0
flag <- 0

j <- 2
for (i in 1:length(kn_vol$corrected_hipchain))	{

	while(kn_vol$hipchain[i]>gps$hip[j]) j <- j+1
	up <- gps$kmtout[j] - ((gps$hip[j]-kn_vol$corrected_hipchain[i]) * gps$outhip[j]/1000)
	down <- gps$kmtout[j-1] + ((kn_vol$corrected_hipchain[i]-gps$hip[j-1]) * gps$outhip[j]/1000)
	flag <- 0
	if((gps$kmtout[j]-up)<(down-gps$kmtout[j-1])) flag <- 1
	if(flag) tout[i] <- up
	if(!flag) tout[i] <- down

	}

#snap transects to nearing node
ntout <- 0
for (i in 1:length(tout))	{
	dif <- knodes$ToMouth_km - tout[i]
	ntout[i] <- knodes$NODE_ID[which.min(abs(dif))]
	}

#subsect transect nodes
tranodes <- knodes[knodes$NODE_ID%in%ntout,]

knodes <- knodes[-nrow(knodes),]

frame <- extent(tranodes) + c(-2,2,-1,1)*300
pdel %>% crop(frame) %>% plot
lines(knodes@coords,lwd=2,col='purple')
points(tranodes@coords,pch=19,col='orange')
points(gps[,1:2],pch=20,col='red')

knodes$NODE_ID %>% summary

#subset nodes within transect length
trans <- knodes[knodes$NODE_ID>=min(tranodes$NODE_ID) & knodes$NODE_ID<=max(tranodes$NODE_ID),]

kn_vol %>% names

#record cross-sectional area at each transect
tranodes$area <- 0
for (i in 1:nrow(tranodes))  tranodes$area[i] <- kn_vol$total_area_m2[i]

#record contributing area at each transect
tranodes$contr <- 0
for (i in 1:nrow(tranodes))  tranodes$contr[i] <- kn_vol$contr_area_km2[i]

#interpolate cross-sectional area between transects
trans$area <- 0
trans$contr <- 0
seg <- 0

for (i in 1:(nrow(tranodes)-1))	{
	trans$area[trans$NODE_ID==tranodes$NODE_ID[i]]	<- tranodes$area[tranodes$NODE_ID==tranodes$NODE_ID[i]]
	trans$area[trans$NODE_ID==tranodes$NODE_ID[i+1]]	<- tranodes$area[tranodes$NODE_ID==tranodes$NODE_ID[i+1]]
	trans$contr[trans$NODE_ID==tranodes$NODE_ID[i]]	<- tranodes$contr[tranodes$NODE_ID==tranodes$NODE_ID[i]]
	trans$contr[trans$NODE_ID==tranodes$NODE_ID[i+1]]	<- tranodes$contr[tranodes$NODE_ID==tranodes$NODE_ID[i+1]]

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
		}
	}


plot(trans$area,trans$DebrisFlow)



# Locate tribs and add 'debris fan effect' to delivery prob

crop_uk %>% plot
points(nodes@coords, pch=19,col='red')
lines(trans@coords)

trans@coords[11,] %>% birdseye(200)
text(nodes@coords,labels = nodes$NODE_ID)

#variable for new debris flow probs
trans$dfp <- 0

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




trans$dfprob <- 0
trans$dfprob[trans$dfp>0] <- trans$dfp[trans$dfp>0]
trans$dfprob[trans$dfp==0] <- trans$DebrisFlow[trans$dfp==0]

#scatter plot debris flow prob vs. x-sec area
plot(trans$dfprob,trans$area,col='orange',pch=18,xlim=c(0,0.015),ylim=c(0,100),
	xlab='Delivery Probability',
	ylab='Cross-sectional Area (m2)')
points(trans$DebrisFlow,trans$area,col='blue')
points(trans$dfp[trans$dfp>0],trans$area[trans$dfp>0],col='red',pch=19)



#plot x-sec area of transects and nodes
plot(trans$area,type='l')
points(trans$area,col='red')
plot(tranodes$area,type='l')


#order transects by cross-sectional area
trans$delta_area <- 0
ao_trans <- trans[order(trans$area),]

#difference in max and min x-sec area
deltar <- ao_trans$area[nrow(ao_trans)] - ao_trans$area[1]

#incremental difference, use differ function
for (i in 2:nrow(ao_trans))	{
	ao_trans$delta_area[i] <- (ao_trans$area[i] - ao_trans$area[i-1])
	}


#order transects by debris flow probability
po_trans <- ao_trans[order(ao_trans$dfprob),]
po_trans$delta_prob <- po_trans$dfprob %>% differ

#difference btwn max and min prob
deltob <- po_trans$dfprob[nrow(po_trans)] - po_trans$dfprob[1]

#
cdf_ardf <- po_trans$delta_area %>% cumult(deltar)
cdf_dfar <- ao_trans$delta_prob %>% cumult(deltob)

plot(po_trans$dfprob,cdf_ardf,
	xlab='Debris Flow Probability',
	ylab='area change per debris flow prob change')


cdf_ar <- ao_trans$area %>% cdf
cdf_df <- ao_trans$dfprob %>% cdf

cdf_ar$y %>% length
cdf_df$y %>% length

plot(cdf_ar$x[1:321],cdf_df$x)


plot(trans$area ~ trans$dfprob)
abline(lm(trans$area ~ trans$dfprob))


probar <- trans$dfprob * trans$area
probar <- probar / sum(probar)

bit <- probar %>% cdf
bit[bit$y<.8] %>% plot
mat <- bit$x %>% cbind(bit$y)
mat <- mat[mat[,2]<1,]
mat %>% plot(type='l',col='orange',lwd=2.5,
	main='CDF of Area-Weighted Delivery Probabilities',
	xlab='Area-Weighted Delivery Probability',
	ylab='Cumulative Proportion of Sample At or Above Probability Value')
trans$dfprob %>% cdf %>% lines(col='purple',lwd=2.5)




pvec <- 0
adist <- 0

for (i in 1:nrow(trans))	{
	for (j in 1:(trans$area[i]%>%round))	{
		pvec <- c(pvec,trans$dfprob[i])
		adist <- c(adist,trans$ToMouth_km[i])
		}
	}

pvec <- pvec[-1]
adist <- adist[-1]

pvec %>% cdf %>% plot(type='l',lwd=2.5,col='purple',
	main='CDF of Delivery Probability',
	xlab='Delivery Probability',
	ylab='Proportion of Sample at or Below Probability')
trans$dfprob %>% cdf %>% lines(col='orange',lwd=2.5)
legend('bottomright',legend=c('Area-Weighted','Unweighted'),fill=c('purple','orange'))


# area-slope weighted sample space

asvec <- 0

trans %>% names
for (i in 1:nrow(trans))	{
	for (j in 1:((trans$contr[i]*trans$GRADIENT[i])%>%round))	{
		asvec <- c(asvec,trans$dfprob[i])
		}
	}

asvec <- asvec[-1]

# area-contr area weighted sample space

acvec <- 0

trans %>% names
for (i in 1:nrow(trans))	{
	for (j in 1:((trans$area[i]/trans$contr[i])%>%round))	{
		acvec <- c(acvec,trans$dfprob[i])
		}
	}

acvec <- acvec[-1]


# area-contr area weighted sample space prime

acvec1 <- 0

for (i in 1:nrow(trans))	{
	for (j in 1:((trans$area[i]*trans$contr[i])%>%round))	{
		acvec1 <- c(acvec1,trans$dfprob[i])
		}
	}

acvec1 <- acvec1[-1]

# area to contr area*slope ratio weighted sample space

aacvec <- 0

trans %>% names
for (i in 1:nrow(trans))	{
	for (j in 1:((trans$area[i]/(trans$contr[i]*trans$GRADIENT[i]))%>%round))	{
		aacvec <- c(aacvec,trans$dfprob[i])
		}
	}

aacvec <- aacvec[-1]





# relative probability of weighted spaces

prel <- 0	#area
asrel <- 0	#area-slope product
acrel <- 0	#area/contr-area ratio
acrel1 <- 0	#area*contr-area product
aacrel <- 0	#area/(contr-area*slope) ratio
uwrel <- 0	#unweighted

for (i in 1:nrow(trans))	{
	prel[i] <- trans$area[i]%>%floor * trans$dfprob[i] / sum(pvec)
	asrel[i] <- (trans$area[i]*trans$GRADIENT[i])%>%floor * trans$dfprob[i] / sum(asvec)
	acrel[i] <- (trans$area[i]/trans$contr[i])%>%round * trans$dfprob[i] / sum(acvec)
	acrel1[i] <- (trans$area[i]*trans$contr[i])%>%round * trans$dfprob[i] / sum(acvec1)
	aacrel[i] <- (trans$area[i]/(trans$contr[i]*trans$GRADIENT[i]))%>%round * trans$dfprob[i] / sum(aacvec)
	uwrel[i] <- trans$DebrisFlow[i] / sum(trans$DebrisFlow)
	}





# W

#png('weighted_cdfs.png')
pvec %>% cdf %>% plot(type='l',lwd=2.5,col='purple',
	main='Weighted CDFs of Delivery Probability',
	xlab='Delivery Probability',
	ylab='Proportion of Sample at or Below Probability')
trans$dfprob %>% cdf %>% lines(col='orange',lwd=2.5)
asvec %>% cdf %>% lines(col='blue',lwd=2.5)
acvec %>% cdf %>% lines(col='green',lwd=2.5)
legend('bottomright',legend=c('Area-Weighted','Area-Slope','Area-Contr Area','Unweighted'),
	fill=c('purple','blue','green','orange'))
#dev.off()

trans %>% names

#Scatter plot of X-Sec Area and Delivery Probability

#png('area_pdel.png')
plot(trans$area,trans$DebrisFlow,col= rgb(red=1, green=0, blue=.2, alpha=.5),pch=20,
	main='Cross-sectional Area vs. Delivery Probability',
	xlab='Cross-sectional Area (m2)',ylab='Delivery Probability' )
abline(lm(trans$DebrisFlow ~ trans$area),col='gray')
#dev.off()




po_trans$delta_prob %>% sum
ao_trans$delta_area %>% sum

kn_vol %>% names

ao_trans <- po_trans[order(po_trans$area),]
ao_delta_dfar <- ao_trans$delta_prob / ao_trans$delta_area
ao_cdf_dfar <- 0
ao_delta_dfar[ao_delta_dfar>1] <- 0
sum_dfar <- ao_delta_dfar %>% sum

for (i in 1:length(ao_delta_dfar))	{
	ao_cdf_dfar[i] <- ao_delta_dfar[1:i] %>% sum / sum_dfar
	}


ao_trans$area %>% cdf %>% plot


pdf_prob <- trans$dfprob %>% pdf


trans$normarea <- trans$area / max(trans$area)
trans$dfarea <- trans$dfprob * trans$normarea
trans$dfarea %>% cdf %>% plot(col='blue')

wt_prob <- trans$area / sum(trans$area)
wt_prob <- trans$dfprob * trans$area
wt_prob %>% cdf %>% plot
trans$dfprob %>% cdf %>% lines




plot(trans$area,trans$ToMouth_km)
points(tranodes$area,tranodes$ToMouth_km,col='red')



#plots relating gps locations to nodes

#png('gps_nodes_eg2.png')
frame <- crop_uk %>% zoom(gps[1,1:2] %>% as.numeric,25) + c(-2,2,-1,1)*30
crop_uk %>% crop(frame) %>% plot(main='Relating GPS locations to channel nodes')
lines(flowline,col='purple',lwd=2)		#flowline
points(gps[,1:2],col='blue',pch=19)	#gps locations
labmat <- gps[1,1:2]
labmat[,2] <- labmat[,2] + 12
text(labmat,labels=c('K-32'))	#labels
points(nodes@coords,pch=19)			#nodes
points(nodes@coords[nodes$NODE_ID%in%gps$nodeid[1]] %>% matrix(ncol=2),
	col='red',pch=19)
legend('topright',fill=c('blue','red','purple','black'),
	legend=c('GPS locations','Corresponding Node','LiDAR-derived Flowline','MB channel nodes'))
#dev.off()

#png('gps_nodes_eg1.png')
frame <- crop_uk %>% zoom(gps[3,1:2] %>% as.numeric + c(-20,20),
	25) + c(-2,2,-1,1)*30
crop_uk %>% crop(frame) %>% plot(main='Relating GPS locations to channel nodes')
lines(flowline,col='purple',lwd=2)		#flowline
points(gps[,1:2],col='blue',pch=19)	#gps locations
labmat <- gps[2:3,1:2]
labmat[,2] <- labmat[,2] + 12
text(labmat,labels=c('K-33','K-36'))	#labels
points(nodes@coords,pch=19)			#nodes
points(nodes@coords[nodes$NODE_ID%in%gps$nodeid[2:3]] %>% matrix(ncol=2),
	col='red',pch=19)
legend('bottomleft',fill=c('blue','red','purple','black'),
	legend=c('GPS locations','Corresponding Node','LiDAR-derived Flowline','MB channel nodes'))
#dev.off()



# estimate valley sediment volume #

#distance between transects
t_dist <- 0
seg_vol <- 0
for (i in 2:nrow(kn_vol))	{
	t_dist[i] <- kn_vol$corrected_hipchain[i] - kn_vol$corrected_hipchain[i-1]
	seg_vol[i] <- kn_vol$total_area_m2[(i-1):i] %>% sum / 2 * t_dist[i]
	}
plot(seg_vol)

tot_vol <- seg_vol %>% sum

# debris flow pct by bank height

bank_ht <- (knowles_samples$coarse_fluvial_m +
	knowles_samples$fine_alluvium_m +
	knowles_samples$df_deposit_m) %>% sum

df_ht <- knowles_samples$df_deposit_m %>% sum

# total debris flow volume (volume x pct debris flow)
df_vol <- tot_vol * df_ht / bank_ht

plot(kn_vol$total_volume_m3)
lines(seg_vol)


# mean debris flow residence time #

radio <- radio[-35,c(1,13:16)]
radio[,2] %>% as.character %>% as.numeric %>% mean


half_life <- radio[,2] %>% as.character %>% as.numeric %>% mean / 2
decay <- log(0.5) / half_life * -1

# annual total flux

tot_flux <- df_vol - df_vol * exp(-decay)

# expected flux per node
# per weighted space

pflux <- tot_flux * prel
asflux <- tot_flux * asrel
acflux <- tot_flux * acrel
acflux1 <- tot_flux * acrel1
aacflux <- tot_flux * aacrel
uwflux <- tot_flux * uwrel


#create weighted sample space of flux rates

af <- 0
asf <- 0
acf <- 0
acf1 <- 0
aacf <- 0

for (i in 1:nrow(trans))	{
	for (j in 1:(trans$area[i]%>%round))	{
		af <- c(af,pflux[i])
		}

	for (j in 1:((trans$contr[i]*trans$GRADIENT[i])%>%round))	{
		asf <- c(asf,asflux[i])
		}

	for (j in 1:((trans$area[i]/trans$contr[i])%>%round))	{
		acf <- c(acf,acflux[i])
		}
	for (j in 1:((trans$area[i]*trans$contr[i])%>%round))	{
		acf1 <- c(acf1,acflux1[i])
		}
	for (j in 1:((trans$area[i]/(trans$contr[i]*trans$GRADIENT[i]))%>%round))	{
		aacf <- c(aacf,aacflux[i])
		}
	}

af <- af[-1]
asf <- asf[-1]
acf <- acf[-1]
acf1 <- acf1[-1]
aacf <- aacf[-1]






# fit expected flux rates to delivery probabilities

dat <- data.frame(flux = acflux, dfprob = trans$dfprob)
mod <- lm(flux ~ dfprob, data=dat)
mod %>% summary

afdat <- data.frame(flux = af, dfprob = pvec)
asfdat <- data.frame(flux = asf, dfprob = asvec)
acfdat <- data.frame(flux = acf, dfprob = acvec)
acfdat1 <- data.frame(flux = acf1, dfprob = acvec1)
aacfdat <- data.frame(flux = aacf, dfprob = aacvec)
a2dat <- data.frame(flux=af, dfprob=pvec, df2=pvec^2)
a2ddat <- data.frame(flux=af, dfprob=pvec, df2=pvec^2, dist=adist)
afmod <- lm(flux ~ dfprob, data=afdat)
asfmod <- lm(flux ~ dfprob, data=asfdat)
acfmod <- lm(flux ~ dfprob, data=acfdat)
acfmod1 <- lm(flux ~ dfprob, data=acfdat1)
aacfmod <- lm(flux ~ dfprob, data=aacfdat)
a2mod <- lm(flux ~ dfprob + df2, data=a2dat)
a2dmod <- lm(flux ~ dfprob + df2 + dist, data=a2ddat)

safdat <- data.frame(flux = pflux, dfprob = trans$dfprob)
sasfdat <- data.frame(flux = asflux, dfprob = trans$dfprob)
sacfdat <- data.frame(flux = acflux, dfprob = trans$dfprob)
sacfdat1 <- data.frame(flux = acflux1, dfprob = trans$dfprob)
saacfdat <- data.frame(flux = aacflux, dfprob = trans$dfprob)
safmod <- lm(flux ~ dfprob, data=safdat)
sasfmod <- lm(flux ~ dfprob, data=sasfdat)
sacfmod <- lm(flux ~ dfprob, data=sacfdat)
sacfmod1 <- lm(flux ~ dfprob, data=sacfdat1)
saacfmod <- lm(flux ~ dfprob, data=saacfdat)

lafdat <- data.frame(flux = af, dfprob = (pvec+0.0001)%>%log)
lasfdat <- data.frame(flux = asf, dfprob = (asvec+0.0001)%>%log)
lacfdat <- data.frame(flux = acf, dfprob = (acvec+0.0001)%>%log)
lafmod <- lm(flux ~ dfprob, data=lafdat)
lasfmod <- lm(flux ~ dfprob, data=lasfdat)
lacfmod <- lm(flux ~ dfprob, data=lacfdat)

uwdat <- data.frame(flux = uwflux, dfprob = trans$DebrisFlow)
luwdat <- data.frame(flux = uwflux, dfprob = (trans$DebrisFlow+0.0001)%>%log)
uwmod <- lm(flux ~ dfprob, data=uwdat)
luwmod <- lm(flux ~ dfprob, data=luwdat)




# predict expected flux rates for nodes based upon model

nog <- readOGR('nodes_debrisflow.shp')
df <- data.frame(dfprob = nog$DebrisFlow)
df2 <- data.frame(dfprob=nog$DebrisFlow,df2=nog$DebrisFlow^2)
adf <- data.frame(dfprob=nog$DebrisFlow,df2=nog$DebrisFlow^2,dist=nog$ToMouth_km)
ldf <- data.frame(dfprob = (nog$DebrisFlow+0.0001)%>%log)
pred <- predict(mod, newdata=df, interval='confidence')

afpred <- predict(afmod, newdata=df, interval='confidence')	#area
asfpred <- predict(asfmod, newdata=df, interval='confidence')	#contr area * slope
acfpred <- predict(acfmod, newdata=df, interval='confidence')	#area / contr area
acfpred1 <- predict(acfmod1, newdata=df, interval='confidence')	#area * contr area
aacfpred <- predict(aacfmod, newdata=df, interval='confidence')	#area / (contr area * slope)
a2pred <- predict(a2mod,newdata=df2, interval='confidence')	#area + area^2
a2dpred <- predict(a2dmod, newdata=adf, interval='confidence')	#area + area^2 + dist

lafpred <- predict(lafmod, newdata=ldf, interval='confidence')	#log area
lasfpred <- predict(lasfmod, newdata=ldf, interval='confidence')	#log contr area * slope
lacfpred <- predict(lacfmod, newdata=ldf, interval='confidence')	#log area / contr area

safpred <- predict(safmod, newdata=df, interval='confidence')	#area
sasfpred <- predict(sasfmod, newdata=df, interval='confidence')	#contr area * slope
sacfpred <- predict(sacfmod, newdata=df, interval='confidence')	#area / contr area
sacfpred1 <- predict(sacfmod1, newdata=df, interval='confidence')	#area * contr area
saacfpred <- predict(saacfmod, newdata=df, interval='confidence')	#area / (contr area * slope)


#predict pdel values from relative pdel of weighted spaces

m_df <- data.frame(dfp=trans$DebrisFlow,pr=prel,dist=trans$ToMouth_km,slope=trans$GRADIENT)
m_pdel <- lm(pr ~ dfp + dist + slope, data=m_df)
#m_pdel <- lm(pr ~ dist + slope, data=m_df)

mrp_df <- data.frame(dfp=nog$DebrisFlow,dist=nog$ToMouth_km,slope=nog$GRADIENT)
mrp_pdel <- predict(m_pdel, newdata=mrp_df, interval='confidence')
mrp <- mrp_pdel[,1]
mrp[mrp<0] <- 0


lm_df <- data.frame(ldfp=(trans$DebrisFlow+0.0001)%>%log,lpr=(prel+0.0001)%>%log,
	dist=trans$ToMouth_km,slope=trans$GRADIENT)
lm_pdel <- lm(lpr ~ ldfp + dist + slope, data=lm_df)
lm_pdel <- lm(lpr ~ dist + slope, data=lm_df)

lmrp_df <- data.frame(ldfp=(nog$DebrisFlow+0.0001)%>%log,dist=nog$ToMouth_km,slope=nog$GRADIENT)
lmrp_pdel <- predict(lm_pdel, newdata=lmrp_df, interval='confidence')
lmrp <- exp(lmrp_pdel[,1])

#flux from lpdel + rel prob + gradient + dist
frp_df <- data.frame(flux=pflux, pdel=trans$dfprob, prel=prel, slope=trans$GRADIENT, dist=trans$ToMouth_km)
frp_mod <- lm(flux ~ pdel + prel, data=frp_df)
frpp_df <- data.frame(pdel=nog$DebrisFlow, prel=mrp, slope=nog$GRADIENT, dist=nog$ToMouth_km)
frpp <- predict(frp_mod, newdata=frpp_df, interval='confidence')
mfrp <- frpp[,1]
mfrp[mfrp<0] <- 0


frp_df <- data.frame(flux=pflux, lpdel=(trans$dfprob+0.0001)%>%log, prel=prel, slope=trans$GRADIENT, dist=trans$ToMouth_km)
frp_mod <- lm(flux ~ lpdel + prel, data=frp_df)
frpp_df <- data.frame(lpdel=(nog$DebrisFlow+0.0001)%>%log, prel=lmrp, slope=nog$GRADIENT, dist=nog$ToMouth_km)
frpp <- predict(frp_mod, newdata=frpp_df, interval='confidence')





maf <- afpred[,1]	#flux rates can not be negative, set to zero
maf[maf<0] <- 0
masf <- asfpred[,1]
masf[masf<0] <- 0
macf <- acfpred[,1]
macf[macf<0] <- 0
macf1 <- acfpred1[,1]
macf1[macf1<0] <- 0
maacf <- aacfpred[,1]
maacf[maacf<0] <- 0
ma2d <- a2dpred[,1]
ma2d[ma2d<0] <- 0


smaf <- safpred[,1]	#flux rates can not be negative, set to zero
smaf[smaf<0] <- 0
smasf <- sasfpred[,1]
smasf[smasf<0] <- 0
smacf <- sacfpred[,1]
smacf[smacf<0] <- 0
smacf1 <- sacfpred1[,1]
smacf1[smacf1<0] <- 0
smaacf <- saacfpred[,1]
smaacf[smaacf<0] <- 0


laf <- lafpred[,1]	#flux rates can not be negative, set to zero
laf[laf<0] <- 0
lasf <- lasfpred[,1]
lasf[lasf<0] <- 0
lacf <- lacfpred[,1]
lacf[lacf<0] <- 0

uwpred <- predict(uwmod, newdata=df, interval='confidence')	#unweighted
luwpred <- predict(luwmod, newdata=ldf, interval='confidence')	#log unweighted

uw <- uwpred[,1]
uw[uw<0] <- 0
luw <- luwpred[,1]
luw[luw<0] <- 0


#compare model differences

mods <- uw %>% matrix(ncol=1) %>% cbind(maf
	%>% cbind(masf
	%>% cbind(macf
	%>% cbind(luw
	%>% cbind(laf
	%>% cbind(lasf
	%>% cbind(lacf)))))))
sums <- mods %>% apply(2,sum)

difs <- (maf-uw) %>% matrix(ncol=1) %>% cbind((masf-uw)
	%>% cbind((macf-uw)
	%>% cbind((laf-luw)
	%>% cbind((lasf-luw)
	%>% cbind((lacf-luw))))))

pctmin <- 0
pctplus <- 0
pctzero <- 0
added <- 0
lost <- 0
for (i in 1:6)	{
	pctmin[i] <- length(difs[difs[,i]<0,i])/length(difs[,i])
	pctplus[i] <-  length(difs[difs[,i]>0,i])/length(difs[,i])
	pctzero[i] <-  length(difs[difs[,i]==0,i])/length(difs[,i])
	added[i] <-  difs[difs[,i]>0,i] %>% sum
	lost[i] <-  difs[difs[,i]<0,i] %>% sum

	}



#color ramp for intensity map
ramp <- colorRampPalette(c('orange','blue','green'))
nog$pred <- pred[,1]
nog$col <- ramp(30)[as.numeric(cut(nog$pred,breaks = 30))]

#flux intensity map
plot(nog,col=nog$col,pch=20,cex=0.25)
legend('topright',legend=c('low flux','mid flux','high flux'),fill=c('orange','blue','green'))

plot(nog,col=nog$col,pch=20,cex=0.25)

plot(acflux,col='purple',pch=20,cex=0.25,
	xlab='Channel Node',ylab='Flux in m3/yr')
points(asflux,col='orange',pch=20,cex=0.25)

points(pflux,col='yellow',pch=20,cex=0.5)

plot(nog$pred,col='forestgreen',pch=20,cex=0.25)



# siuslaw basin delivery probability heat map

nog$dfcol <- ramp(30)[as.numeric(cut(nog$DebrisFlow,breaks=30))]
plot(nog,col=nog$dfcol,pch=20,cex=0.25)
legend('topright',legend=c('low pdel','mid pdel','high pdel'),fill=c('orange','blue','green'))

# unweighted flux fit to delivery probabilities

uwdat <- data.frame(flux = uwflux, dfprob = trans$DebrisFlow)
uwmod <- lm(flux ~ dfprob, data=uwdat)
uwmod %>% summary

uwpred <- predict(uwmod, newdata=df, interval='confidence')
nog$uwpred <- uwpred[,1]
nog$uwcol <- ramp(30)[as.numeric(cut(nog$uwpred,breaks=30))]
plot(nog,col=nog$uwcol,pch=20,cex=0.25)
legend('topright',legend=c('low uw','mid uw','high uw'),fill=c('orange','blue','green'))

nog$pdif <- nog$uwpred - nog$pred
nog$pcol <- ramp(30)[as.numeric(cut(nog$pdif,breaks=30))]
plot(nog,col=nog$pcol,pch=20,cex=0.25)
legend('topright',legend=c('low dif','mid dif','high dif'),fill=c('orange','blue','green'))

# scatter plot with linear model fits

plot(acvec,acf,col= rgb(red=.2, green=.6, blue=.4, alpha=.5),pch=20,cex=.7)
points(pvec,af,col= rgb(red=.8, green=.5, blue=0, alpha=.5),pch=20,cex=.7)
points(asvec,asf,col= rgb(red=.2, green=.4, blue=.6, alpha=.5),pch=20,cex=.7)
abline(lm(acf~acvec),col= rgb(red=.2, green=.6, blue=.4, alpha=.5),lwd=2.5)
abline(lm(af~pvec),col= rgb(red=.8, green=.5, blue=0, alpha=.5),lwd=2.5)
abline(lm(asf~asvec),col= rgb(red=.2, green=.4, blue=.6, alpha=.5),lwd=2.5)

plot(acvec%>%log,acf,col= rgb(red=.2, green=.6, blue=.4, alpha=.5),pch=20,cex=.7)
points(pvec%>%log,af,col= rgb(red=.8, green=.5, blue=0, alpha=.5),pch=20,cex=.7)
points(asvec%>%log,asf,col= rgb(red=.2, green=.4, blue=.6, alpha=.5),pch=20,cex=.7)
abline(lm(acf~(acvec%>%log)),col= rgb(red=.2, green=.6, blue=.4, alpha=.5),lwd=2.5)
abline(lm(af~(pvec%>%log)),col= rgb(red=.8, green=.5, blue=0, alpha=.5),lwd=2.5)
abline(lm(asf~(asvec%>%log)),col= rgb(red=.2, green=.4, blue=.6, alpha=.5),lwd=2.5)

plot(acfpred[,1],col= rgb(red=.2, green=.4, blue=.6, alpha=.5),pch=20,cex=.25)
points(afpred[,1],col= rgb(red=.2, green=.6, blue=.4, alpha=.5),pch=20,cex=.25)
points(asfpred[,1],col= rgb(red=.8, green=.5, blue=0, alpha=.5),pch=20,cex=.25)
points(uwpred[,1],col= rgb(red=1, green=1, blue=1, alpha=.5),pch=20,cex=.25)

plot(lafpred[,1],col= rgb(red=.2, green=.6, blue=.4, alpha=.25),pch=20,cex=.05,ylim=c(0,4))
points(lasfpred[,1],col= rgb(red=.8, green=.5, blue=0, alpha=.25),pch=20,cex=.05)
points(lacfpred[,1],col= rgb(red=.2, green=.4, blue=.6, alpha=.25),pch=20,cex=.05)
points(uwpred[,1],col= rgb(red=0, green=0, blue=0, alpha=.05),pch=20,cex=.25)

plot(abs(lafpred[,1]-uwpred[,1]),col= rgb(red=.2, green=.6, blue=.4, alpha=.5),pch=20,cex=.25,ylim=c(0,4))
points(abs(lasfpred[,1]-uwpred[,1]),col= rgb(red=.8, green=.5, blue=0, alpha=.5),pch=20,cex=.25)
points(abs(lacfpred[,1]-uwpred[,1]),col= rgb(red=.2, green=.4, blue=.6, alpha=.5),pch=20,cex=.25)

plot((laf-uwpred[,1]),col= rgb(red=.2, green=.6, blue=.4, alpha=.05),pch=20,cex=.25,ylim=c(-1.4,.4))
points((lacf-uwpred[,1]),col= rgb(red=.2, green=.4, blue=.6, alpha=.05),pch=20,cex=.25)
points((lasf-uwpred[,1]),col= rgb(red=.8, green=.5, blue=0, alpha=.05),pch=20,cex=.25)

plot(aacvec,aacf,col= rgb(red=1, green=0, blue=0, alpha=.5),pch=20,cex=.7,
	xlab='Delivery Probability',ylab='Flux in m3/yr')
points(pvec,af,col= rgb(red=0, green=0, blue=1, alpha=.5),pch=20,cex=.7)
points(asvec,asf,col= rgb(red=.2, green=.9, blue=.2, alpha=.5),pch=20,cex=.7)
points(trans$DebrisFlow,uwflux,col= rgb(red=0, green=0, blue=0, alpha=.5),pch=20,cex=.7)
abline(lm(aacf~aacvec),col= rgb(red=1, green=0, blue=0, alpha=.5),lwd=2.5)
abline(lm(af~pvec),col= rgb(red=0, green=0, blue=1, alpha=.5),lwd=2.5)
abline(lm(asf~asvec),col= rgb(red=.2, green=.9, blue=.2, alpha=.5),lwd=2.5)
abline(lm(uwflux~trans$DebrisFlow),col= rgb(red=0, green=0, blue=0, alpha=.5),lwd=2.5)
legend('topleft',legend=c('Area','Area*Slope','Area:(Contr Area*Slope)','Unweighted'),fill=c('blue','green','red','black'))

plot(trans$dfprob,aacflux,col= rgb(red=1, green=0, blue=0, alpha=.5),pch=20,cex=.7,
	xlab='Delivery Probability',ylab='Flux in m3/yr')
points(trans$dfprob,pflux,col= rgb(red=0, green=0, blue=1, alpha=.5),pch=20,cex=.7)
points(trans$dfprob,asflux,col= rgb(red=.2, green=.9, blue=.2, alpha=.5),pch=20,cex=.7)
points(trans$dfprob,uwflux,col= rgb(red=0, green=0, blue=0, alpha=.5),pch=20,cex=.7)
abline(lm(aacflux~trans$dfprob),col= rgb(red=1, green=0, blue=0, alpha=.5),lwd=2.5)
abline(lm(pflux~trans$dfprob),col= rgb(red=0, green=0, blue=1, alpha=.5),lwd=2.5)
abline(lm(asflux~trans$dfprob),col= rgb(red=.2, green=.9, blue=.2, alpha=.5),lwd=2.5)
abline(lm(uwflux~trans$dfprob),col= rgb(red=0, green=0, blue=0, alpha=.5),lwd=2.5)
legend('topleft',legend=c('Area','Area*Slope','Area:(Contr Area*Slope)','Unweighted'),fill=c('blue','green','red','black'))



#difference plot between models, zoomed in on hot spot of change
plot((laf[99000:100000]-uwpred[99000:100000,1]),col= rgb(red=.2, green=.6, blue=.4, alpha=.5),pch=20,cex=.25,ylim=c(-1.4,.4))
points((lacf[99000:100000]-uwpred[99000:100000,1]),col= rgb(red=.2, green=.4, blue=.6, alpha=.5),pch=20,cex=.25)
points((lasf[99000:100000]-uwpred[99000:100000,1]),col= rgb(red=.8, green=.5, blue=0, alpha=.5),pch=20,cex=.25)

plot((laf[99000:100000]-luwpred[99000:100000,1]),col= rgb(red=.2, green=.6, blue=.4, alpha=.5),pch=20,cex=.25)
points((lacf[99000:100000]-luwpred[99000:100000,1]),col= rgb(red=.2, green=.4, blue=.6, alpha=.5),pch=20,cex=.25)
points((lasf[99000:100000]-luwpred[99000:100000,1]),col= rgb(red=.8, green=.5, blue=0, alpha=.5),pch=20,cex=.25)

#color ramp for intensity map
ramp <- colorRampPalette(c('orange','blue','green'))
nog$uw <- uw
nog$luw <- luw
nog$af <- maf
nog$asf <- masf
nog$acf <- macf
nog$acf1 <- macf1
nog$aacf <- maacf
nog$a2 <- a2pred[,1]
nog$a2d <- ma2d
nog$saf <- smaf
nog$sasf <- smasf
nog$sacf <- smacf
nog$sacf1 <- smacf1
nog$saacf <- smaacf
nog$frpp <- frpp
nog$laf <- laf
nog$lasf <- lasf
nog$lacf <- lacf
nog$afdif <- maf-uw
nog$asfdif <- masf-uw
nog$acfdif <- macf-uw
nog$lafdif <- laf-luw
nog$lasfdif <- lasf-luw
nog$lacfdif <- lacf-luw


nog$uwcol <- ramp(30)[as.numeric(cut(nog$uw,breaks = 30))]
nog$luwcol <- ramp(30)[as.numeric(cut(nog$luw,breaks = 30))]
nog$afcol <- ramp(30)[as.numeric(cut(nog$af,breaks = 30))]
nog$asfcol <- ramp(30)[as.numeric(cut(nog$asf,breaks = 30))]
nog$acfcol <- ramp(30)[as.numeric(cut(nog$acf,breaks = 30))]
nog$lafcol <- ramp(30)[as.numeric(cut(nog$laf,breaks = 30))]
nog$lasfcol <- ramp(30)[as.numeric(cut(nog$lasf,breaks = 30))]
nog$lacfcol <- ramp(30)[as.numeric(cut(nog$lacf,breaks = 30))]
nog$afdifcol <- ramp(30)[as.numeric(cut(nog$afdif,breaks = 30))]
nog$asfdifcol <- ramp(30)[as.numeric(cut(nog$asfdif,breaks = 30))]
nog$acfdifcol <- ramp(30)[as.numeric(cut(nog$acfdif,breaks = 30))]
nog$lafdifcol <- ramp(30)[as.numeric(cut(nog$lafdif,breaks = 30))]
nog$lasfdifcol <- ramp(30)[as.numeric(cut(nog$lasfdif,breaks = 30))]
nog$lacfdifcol <- ramp(30)[as.numeric(cut(nog$lacfdif,breaks = 30))]



#flux intensity map
dev.new()
plot(nog,col=nog$uwcol,pch=20,cex=0.25,main='unweighted')
dev.new()
plot(nog,col=nog$luwcol,pch=20,cex=0.25,main='log unweighted')
dev.new()
plot(nog,col=nog$afcol,pch=20,cex=0.25,main='area')
dev.new()
plot(nog,col=nog$asfcol,pch=20,cex=0.25,main='contr area*slope')
dev.new()
plot(nog,col=nog$acfcol,pch=20,cex=0.25,main='area/contr area')
dev.new()
plot(nog,col=nog$lafcol,pch=20,cex=0.25,main='log area')
dev.new()
plot(nog,col=nog$lasfcol,pch=20,cex=0.25,main='log contr area*slope')
dev.new()
plot(nog,col=nog$lacfcol,pch=20,cex=0.25,main='log area/contr area')
dev.new()
plot(nog,col=nog$afdifcol,pch=20,cex=0.25,main='area')
dev.new()
plot(nog,col=nog$asfdifcol,pch=20,cex=0.25,main='contr area*slope')
dev.new()
plot(nog,col=nog$acfdifcol,pch=20,cex=0.25,main='area/contr area')
dev.new()
plot(nog,col=nog$lafdifcol,pch=20,cex=0.25,main='log area')
dev.new()
plot(nog,col=nog$lasfdifcol,pch=20,cex=0.25,main='log contr area*slope')
dev.new()
plot(nog,col=nog$lacfdifcol,pch=20,cex=0.25,main='log area/contr area')
legend('topright',legend=c('low flux','mid flux','high flux'),fill=c('orange','blue','green'))
dev.new()



plot((lasf-luw),col= rgb(red=.8, green=.5, blue=0, alpha=.01),pch=20,cex=.25)
points((lacf-luw),col= rgb(red=.2, green=.4, blue=.6, alpha=.01),pch=20,cex=.25)
points(laf-luw,col= rgb(red=.2, green=.6, blue=.4, alpha=.01),pch=20,cex=.25)

plot(nog$lafdif[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.8, green=.5, blue=0, alpha=.8),pch=20,cex=.25)
points(nog$lasfdif[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.2, green=.4, blue=.6, alpha=.8),pch=20,cex=.25)
points(nog$lacfdif[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.2, green=.6, blue=.4, alpha=.8),pch=20,cex=.25)

plot(nog$laf[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.8, green=.5, blue=0, alpha=.8),pch=20,cex=.25)
points(nog$lasf[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.2, green=.4, blue=.6, alpha=.8),pch=20,cex=.25)
points(nog$lacf[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.2, green=.6, blue=.4, alpha=.8),pch=20,cex=.25)
points(acflux,col= rgb(red=.1, green=.1, blue=.1, alpha=.8),pch=20,cex=.5)

plot(nog$af[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.8, green=.5, blue=0, alpha=.8),pch=20,cex=.25)
points(nog$asf[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.2, green=.4, blue=.6, alpha=.8),pch=20,cex=.25)
points(nog$acf[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.2, green=.6, blue=.4, alpha=.8),pch=20,cex=.25)
points(acflux,col= rgb(red=.1, green=.1, blue=.1, alpha=.8),pch=20,cex=.5)
points(uwflux,col= rgb(red=.8, green=.1, blue=.1, alpha=.8),pch=20,cex=.5)


plot(nog$asf[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.0, green=.0, blue=.0, alpha=.8),pch=20,cex=.4)
points(asflux,col= rgb(red=.0, green=.0, blue=.8, alpha=.8),pch=20,cex=.4)
points(uwflux,col= rgb(red=.0, green=.8, blue=.0, alpha=.8),pch=4,cex=.5)

plot(nog$acf[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.0, green=.0, blue=.0, alpha=.8),pch=20,cex=.4)
points(acflux,col= rgb(red=.8, green=0, blue=0, alpha=.8),pch=20,cex=.4)
points(uwflux,col= rgb(red=.0, green=.8, blue=.0, alpha=.8),pch=4,cex=.5)

plot(uwflux,col= rgb(red=.0, green=.0, blue=.0, alpha=.8),pch=20,cex=.4,ylim=c(0,1.5))
points(pflux,col= rgb(red=1, green=.0, blue=.0, alpha=.8),pch=20,cex=.4)
points(asflux,col= rgb(red=0, green=1, blue=.0, alpha=.8),pch=20,cex=.4)
points(acflux,col= rgb(red=0, green=0, blue=1, alpha=.8),pch=20,cex=.4)
legend('topleft',legend=c('unweighted','area','contr area*slope','area/contr area'),
	fill=c('black','red','green','blue'))

dev.new()
plot(nog$af[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.0, green=.0, blue=1, alpha=.8),pch=20,cex=.4)
points(pflux,col= rgb(red=1, green=0, blue=0, alpha=.8),pch=20,cex=.4)
points(uwflux,col= rgb(red=.0, green=.0, blue=0, alpha=.8),pch=20,cex=.5)
legend('topleft',legend=c('unweighted','area est','area pred'),fill=c('black','red','blue'))

dev.new()
plot(nog$asf[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.0, green=.0, blue=1, alpha=.8),pch=20,cex=.4)
points(asflux,col= rgb(red=1, green=0, blue=0, alpha=.8),pch=20,cex=.4)
points(uwflux,col= rgb(red=.0, green=.0, blue=0, alpha=.8),pch=20,cex=.5)
legend('topleft',legend=c('unweighted','contr area*slope est','contr area*slope pred'),fill=c('black','red','blue'))

dev.new()
plot(nog$acf[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.0, green=.0, blue=1, alpha=.8),pch=20,cex=.4)
points(acflux,col= rgb(red=1, green=0, blue=0, alpha=.8),pch=20,cex=.4)
points(uwflux,col= rgb(red=.0, green=.0, blue=0, alpha=.8),pch=20,cex=.5)
legend('topleft',legend=c('unweighted','area/contr area est','area/contr area pred'),fill=c('black','red','blue'))

dev.new()
plot(nog$acf1[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.0, green=.0, blue=1, alpha=.8),pch=20,cex=.4)
points(acflux1,col= rgb(red=1, green=0, blue=0, alpha=.8),pch=20,cex=.4)
points(uwflux,col= rgb(red=.0, green=.0, blue=0, alpha=.8),pch=20,cex=.5)
legend('topleft',legend=c('unweighted','area*contr area est','area*contr area pred'),fill=c('black','red','blue'))

dev.new()
plot(nog$aacf[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.0, green=.0, blue=1, alpha=.8),pch=20,cex=.4)
points(aacflux,col= rgb(red=1, green=0, blue=0, alpha=.8),pch=20,cex=.4)
points(uwflux,col= rgb(red=.0, green=.0, blue=0, alpha=.8),pch=20,cex=.5)
legend('topleft',legend=c('unweighted','area/(contr area * slope) est','area/(contr area * slope) pred'),fill=c('black','red','blue'))

dev.new()
plot(aacflux,col= rgb(red=.0, green=.0, blue=1, alpha=.8),pch=20,cex=.4)
points(nog$saacf[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=1, green=0, blue=0, alpha=.8),pch=20,cex=.4)
points(uwflux,col= rgb(red=.0, green=.0, blue=0, alpha=.8),pch=20,cex=.5)
legend('topleft',legend=c('unweighted','area/(contr area * slope) est','area/(contr area * slope) pred'),fill=c('black','red','blue'))
points(nog$aacf[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.0, green=1, blue=.5, alpha=.8),pch=20,cex=.4)
points(nog$aacf[nog$NODE_ID%in%trans$NODE_ID]*(sum(aacflux)/sum(nog$aacf[nog$NODE_ID%in%trans$NODE_ID])),col= rgb(red=.0, green=1, blue=.5, alpha=.8),pch=20,cex=.4)

dev.new()
plot(nog$af[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.0, green=.0, blue=1, alpha=.8),pch=20,cex=.4)
points(pflux,col= rgb(red=1, green=0, blue=0, alpha=.8),pch=20,cex=.4)
points(uwflux,col= rgb(red=.0, green=.0, blue=0, alpha=.8),pch=20,cex=.5)
legend('topleft',legend=c('unweighted','area est','area pred'),fill=c('black','red','blue'))


dev.new()
plot(nog$laf[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.0, green=.0, blue=1, alpha=.8),pch=20,cex=.4)
points(pflux,col= rgb(red=1, green=0, blue=0, alpha=.8),pch=20,cex=.4)
points(nog$luw[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.0, green=.0, blue=0, alpha=.8),pch=20,cex=.5)
legend('topleft',legend=c('unweighted','area est','area pred'),fill=c('black','red','blue'))

plot(nog$luwf[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.0, green=.0, blue=0, alpha=.8),pch=20,cex=.5)


dev.new()
plot(nog$af[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.0, green=.0, blue=1, alpha=.8),pch=20,cex=.4)
points(nog$a2[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.2, green=.8, blue=.3, alpha=.8),pch=20,cex=.4)
points(nog$a2d[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=.8, green=.8, blue=0, alpha=.8),pch=20,cex=.4)

points(pflux,col= rgb(red=1, green=0, blue=0, alpha=.8),pch=20,cex=.4)
points(uwflux,col= rgb(red=.0, green=.0, blue=0, alpha=.8),pch=20,cex=.5)
legend('topleft',legend=c('unweighted','area pred', 'area + area^2 pred','area est'),fill=c('black','blue','green','red'))


dev.new()
#plot(nog$a2d[nog$NODE_ID%in%trans$NODE_ID],col= rgb(red=1, green=0, blue=.2, alpha=.8),pch=20,cex=.4)
plot(nog$a2d[nog$NODE_ID%in%trans$NODE_ID]* (sum(pflux)/sum(a2dpred[nog$NODE_ID%in%trans$NODE_ID,1])),
	col= rgb(red=.5, green=1, blue=0, alpha=.8),pch=20,cex=.4,xlab='channel node',ylab='flux in m3/yr')
points(mfrp[nog$NODE_ID%in%trans$NODE_ID]* (sum(pflux)/sum(mfrp[nog$NODE_ID%in%trans$NODE_ID])),
	col= rgb(red=1, green=0, blue=1, alpha=.8),pch=20,cex=.4)


points(pflux,col= rgb(red=0, green=0, blue=1, alpha=.8),pch=20,cex=.4)
points(uwflux,col= rgb(red=.0, green=.0, blue=0, alpha=.8),pch=20,cex=.5)
legend('topleft',legend=c('unweighted','area + area^2 + dist pred','prel pred','area est')
	,fill=c('black','green','purple','blue'))


plot(nog$frpp[nog$NODE_ID%in%trans$NODE_ID], col= rgb(red=1, green=0, blue=1, alpha=.8),pch=20,cex=.4)

nfrp <- mfrp* (sum(pflux)/sum(mfrp[nog$NODE_ID%in%trans$NODE_ID]))

max_contr <- trans$contr %>% max	#in km2
sius_contr <- (504000%>%ACtoHC)*.01	#in km2

contrat <- sius_contr / max_contr
sius_flux <- sum(pflux) * contrat



