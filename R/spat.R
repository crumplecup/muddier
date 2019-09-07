
#'  birdseye
#'
#'  convenience plotting function
#'
#'  @param pt is a coordinate pair (x,y)
#'  @param mag is an integer of magnification
#'  @param flow is the river network flowline
#'  @param nod is the river node output of the M/B model
#'  @param map is the background raster
#'  @return a plot centered on `pt` at zoom `mag`
#'  @export
#'  @importFrom magrittr %>%
#'  @seealso circ

birdseye <- function(pt,mag=40,flow=flowline,nod=nodes,map=uk)	{
  frame <- map %>% zoom(pt,mag)
  crop <- map %>% crop(frame)
  plot(crop)
  points(flow,pch=20,col='blue')
  lines(pt %>% circ)
  points(nod,col='red')
  points(c(pt,pt)%>%matrix(ncol=2)%>%t,pch=19,col='red')
}



#' contr delta
#'
#' utility function for identifying tributaries
#'
#' @param vec is a numeric vector of contributing areas
#' @return a vector of differences in contributing area per step
#' @export

contr_delta <- function(vec) {
  len <- length(vec)
  delta <- vector(len, mode = 'numeric')
  for (i in 2:len) {
    delta[i] <- vec[i-1] - vec[i]
  }
  delta
}



#' differ
#'
#' given a vector of length > 1, returns a vector differences between elements
#'
#' @param vec is a numeric vector of length > 1
#' @return a vector of differences between elements in `vec`
#' @export

differ <- function(vec)	{
  dif <- vector(length(vec)-1, mode = 'numeric')
  for (i in 1:(length(vec)-1))	{
    dif[i] <- (vec[i+1] - vec[i])
  }
  return(dif)
}



#' draw circle
#'
#'
#'
#' @param pt is a coordinate pair
#' @param rad is the numeric radius of the circle
#' @return matrix of coordinates for a circle of radius `rad` with center `pt`
#' @importFrom magrittr %>%
#' @seealso tri_length
#' @export

circ <- function(pt,rad=15)	{
  mat <- c(0,0) %>% matrix(ncol=2)
  for (i in 1:360)	{
    mat <- mat %>% rbind(tri_length(pt,i,rad))
  }
  mat <- mat[-1,]
  return(mat)
}


#' get bearing
#'
#' given a matrix of two coordinate pairs with rows (x,y), returns the
#' bearing of the angle of the line cross through `coords` in degrees
#'
#' @param coords is a matrix of two coordinate pairs with rows (x,y)
#' @return bearing of angle of line crossing through `coords` in degrees

get_bear <- function(coords)	{
  bear <- atan2(coords[2,1]-coords[1,1],coords[2,2]-coords[1,2])
  if(bear<0) bear <- 2*pi + bear
  bear <- bear * (180/pi)
  if(bear>360) bear <- bear - 360
  return(bear)
}



#' hipchain to outlet
#'
#' converts hipchain lengths to outlet distances given the outlet to hipchain ratio
#'
#' @param hipchain is a numeric vector of hipchain distances
#' @param trib_hip is a numeric vector of tributary hipchain distances
#' @param trib_out is a numeric vector of tributary outlet distances
#' @param out_hip is a numeric vector of outlet to hipchain ratios
#' @return outlet distances converted from hipchain distances
#' @export
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



#' hipchain
#'
#' Given a matrix representation of the channel `chan` with cols (x,y) and rows running
#' upstream to downstream, pulls all pts within `inc` of `pt` along `chan`.
#' If `us` is true, then the last row of the subset matrix is `inc` units upstream
#' If `us` is false, then the first row of the subset is `inc` units downstream
#' The function recursively calls itself using the new point position and
#' distance remaining in place of `pt` and `dist` until having fully traveled
#' the original `dist` by increments of `inc`.
#'
#' @param pt is a coordinate pair
#' @param dist is a numeric hipchain distance (meters along channel)
#' @param inc is a numeric step increment
#' @param chan is a matrix of coords delineating channel from downstream to upstream
#' @param us is a boolean indicating up or downstream travel
#' @param filter is a numeric buffer distance defining off-channel pts
#' @return coordinate position `dist` along `chan` from `pt`, and distance remaining in list
#' @importFrom magrittr %>%

hipchain <- function(pt,dist,inc=10,chan=flowline,us=T,filter=4000)	{

  dif <- ((chan[,1] - pt[1])^2 + (chan[,2] - pt[2])^2) %>% sqrt
  pt <- chan[dif <= inc, ]
  dmax <- dif[dif <= inc]

  if(us)	{
    pt <- pt[nrow(pt),]
    dmax <- dmax[length(dmax)]
  }
  if(!us)	{
    pt <- pt[1,]
    dmax <- dmax[1]
  }

  packet <- list(pt,dist-dmax)

  print(packet[[2]])
  while(packet[[2]]>inc)	{
    packet <- c(packet[[1]]) %>% hipchain(packet[[2]],inc,chan,us,filter)
  }
  if(packet[[2]]>1.5)	{
    packet <- c(packet[[1]]) %>% hipchain(packet[[2]],packet[[2]],chan,us,filter)
  }
  return(packet)
}





#' increment pathway
#'
#' given a matrix of coordinates represented a line, returns a line path
#' incremented into segments of approximately unit length.
#' not for use on branching line paths.
#' given a matrix of coordinates representing a path or line
#' take the bearing and length of each segment
#' divide the segment into length-1 points
#' for each point along the segment, calculate the x and y position
#' return a new path consisting containing all old and new x and y values
#' the purpose of this function is to increment along the channel line
#' so I can take more transects and get a finer resolution of detail on the channel
#'
#' @param coords is a matrix of (x,y) coordinate pairs
#' @return line path incremented into segments appr unit value in length
#' @importFrom magrittr %>%
#' @seealso tri_length
#' @export

inc_path <- function(coords)	{

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


#' interpolate by node
#'
#' interpolate values based on outlet distances of nodes
#'
#' @param vals is a vector of values at tributary nodes
#' @param trib_out is a vector of tributary outlet distances
#' @param trib_id is a vector of tributary node ids
#' @param chan_out is a vector of channel outlet distances
#' @param chan_id is a vector of channel node ids
#' @return interpolated `vals` based on outlet distances of tribs
#' @export

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




#' low path
#'
#' Given a matrix of coordinates representing a line, takes transects of
#' radius `rad` at a bearing perpendicular to the angle of the line,
#' returns the low point along each transect in a matrix representing
#' the low path.
#'
#' @param coords is a matrix of (x,y) coordinate pairs
#' @param rad is a radius in meters
#' @param map is elevation data (DEM or lidar)
#' @return the coordinates of the low path within radius of coords
#' @importFrom magrittr %>%
#' @seealso inc_path
#' @seealso take_transect
#' @export


low_path <- function(coords,rad=25,map=uk)	{
  path <- c(0) %>% matrix(nrow=1,ncol=2)
  for(i in 1:(nrow(coords)-1))	{
    bear <- atan2(coords[i+1,1]-coords[i,1],coords[i+1,2]-coords[i,2])
    if (bear<0) bear <- bear + 2*pi
    bear <- bear * (180/pi) + 90
    if(bear>360) bear <- bear-360
    seg <- coords[i:(i+1),] %>% inc_path
    for(j in 1:nrow(seg))	{
      low <- seg[j,] %>% take_transect(bear,rad)
      low <- low %>% cbind(extract(map,low))
      lowpt <- low[low[,3]==min(low[,3])][1:2]
      path <- path %>% rbind(lowpt)
    }
  }
  path <- path[-1,]
  return(path)
}


#' plot point
#'
#' produces a plot centered around point `pt` zoomed in at factor `mag`
#'
#' @param pt is a coordinate pair
#' @param mag is a magnification factor (numeric)
#' @param so is a spatial object to plot the point upon
#' @return a plot of `so` centered on `pt` at zoom factor `mag`
#' @export

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



#' snap point
#'
#' finds the nearest point on a line path to a given point.
#' does not interpolate betweens points on line to find a local minimum.
#'
#' @param pt is a coordinate pair
#' @param line is a matrix with rows (x,y) specifying a line path
#' @return point on `line` closest in position to `pt`
#' @importFrom magrittr %>%
#' @export

snap <- function(pt, line = nodepath){
  dif <- ((line[,1] - pt[1])^2 + (line[,2] - pt[2])^2) %>% sqrt
  pt <- line[dif == min(dif), ]
  return(pt)
}



#' snap to node
#'
#' return id of node nearest to given outlet distance
#'
#' @param dist is a vector of outlet distances to place on the channel
#' @param outlet is a vector of outlet distances from the channel node network
#' @param ids is a vector of ids from the channel node network
#' @return vector of ids from `ids` nearest to distances in `dist`
#' @export

snap_to_node <- function(dist, outlet, ids)  {
  id <- 0
  for (i in 1:length(dist))	{
    difs <- outlet - dist[i]
    id[i] <- ids[which.min(abs(difs))]
  }
  id
}





#' take transect
#'
#' given a point pt, a bearing bear and a radius rad
#' take a transect centered on point pt, along bearing bear with radius rad
#'
#' @param pt is a coordinate pair (centerpoint of transect)
#' @param bear is a bearing in degrees (angle of transect)
#' @param rad is a radius specifying (width of transect)
#' @return matrix of coordinate pairs along transect using inc_path()
#' @importFrom magrittr %>%
#' @seealso tri_length
#' @seealso inc_path
#' @export

take_transect <- function(pt,bear,rad)	{

  start <- pt %>% tri_length(bear,rad)
  bear <- bear + 180
  if(bear>=360) bear <- bear - 360
  stop <- pt %>% tri_length(bear,rad)
  pt <- start %>% rbind(stop) %>% inc_path
  return(pt)
}





#' triangle length
#'
#' Given a coordinate pair, and the length of the hypotenuse
#' of a right triangle, returns the coordinate pair at the
#' other end of the segment.
#'
#' @param pt is a coordinate pair
#' @param bear is a bearing in degrees
#' @param dist is a numeric distance (hypotenuse length)
#' @param north is a boolean indicating north or south
#' @return coords of end point along hypotenuse
#' @export

tri_length <- function(pt,bear,dist,north=TRUE)	{

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



#' zoom
#'
#' Zoom in on a raster like a microscope
#' draws a new extent m(x-,x+,y-,y+)
#' d(x,y) is centerpoint of raster
#' z is zoom factor scaled by ((x+ - x-)|(y+ - y-))/(2 * z)
#'
#' @param x is a spatial object
#' @param d is coordinate pair (centerpoint)
#' @param z is numeric zoom factor
#' @return spatial object scaled to new extent by zoom factor
#' @export

zoom <- function(x,d,z=2)	{
  m <- 0

  if (class(x) != 'Extent' &
      class(x) != 'numeric')	{
    x <- raster::extent(x)
  }

  m[1] <- d[1] - ( (x[2] - x[1]) / (2 * z) )
  m[2] <- d[1] + ( (x[2] - x[1]) / (2 * z) )
  m[3] <- d[2] - ( (x[4] - x[3]) / (2 * z) )
  m[4] <- d[2] + ( (x[4] - x[3]) / (2 * z) )

  return(m)
}


































