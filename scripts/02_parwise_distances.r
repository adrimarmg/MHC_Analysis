# distance.africa() calculates the geographic distance between all pairs of 
# populations, but forces the route to go through specific waypoint locations 
# when populations are on different continents. 

distance.africa <- function(geoinfo)
{
  # SUPPORT FUNCTIONS
  
  # Km between two pointscon the globe using the haversine. 
  geodist.hav <- function(loc1, loc2)
  { 
    ct <- pi/180 
    radius <- 6371
    
    #the latitudes and longitudes
    lat1 <- loc1[1]
    lon1 <- loc1[2]
    lat2 <- loc2[1]
    lon2 <- loc2[2]
    
    a <- sin( ((lat2 - lat1)*ct)/2)^2 +
      cos(lat1*ct) * cos(lat2*ct) * sin(((lon2 - lon1)*ct)/2)^2 
    c <- 2 * asin(min(1,sqrt(a)))
    d <- radius * c
  }
  
  
  # distances for interwaypoint travel, passing through other waypoints
  # where sensible. 
  interway.km <- function(waypoints, waypassage)
  {
    interway.dist <- rep(0,nrow(waypassage))
    # number of waypoints for each possible pair
    nlinks <-  rowSums(!is.na(waypassage))
    
    for(i in 1:nrow(waypassage))
    {
      #with one or fewer links, interway.dist = 0
      if(nlinks[i]==0|1) {interway.dist[i] <- 0}
      
      #with 2 or more links, interway.dist is calculated jump by jump
      if(nlinks[i]>1) {
        link.vec <- waypassage[i, which(!is.na(waypassage[i,]))]
        link.dists <- rep(0, length(link.vec)-1)
        
        for(j in 1:length(link.dists))
        {
          #calculate distance for each link
          way1 <- waypoints[which(rownames(waypoints)==names(link.vec[which(link.vec==j)])),]
          way2 <- waypoints[which(rownames(waypoints)==names(link.vec[which(link.vec==j+1)])),]
          link.dists[j] <- geodist.hav(way1, way2)
        }
        interway.dist[i] <- sum(link.dists)
        names(interway.dist) <- rownames(waypassage)
      }}
    
    interway <- data.frame(wayStartEnd=paste(region.clust, region.seq, sep=""),
                           wayPoints=nlinks, wayKm=interway.dist)
  }
  
  ############################################################
  
  # SUPPORT DATA (values do not depend on geoinfo)
  
  # waypoints: matrix of waypoints w/ latitude and longitude
  Cairo <- c(30, 31)
  Istanbul <- c(41, 29)
  PhnomPenh <- c(11.5, 105)
  Anadyr <- c(64.5, 177.5)
  Nome <- c(64.5, -165.4)
  MexicoCity <- c(19.5, -99)
  waypoints <- rbind(Cairo, Istanbul, PhnomPenh, Anadyr, Nome, MexicoCity)
  colnames(waypoints) <- c("lat", "long")
  
  # waypassage: for all possible pairs of interregional travel, the waypoints
  regions <- c("Africa", "MidEast", "Asia", "Europe", "Oceania", "NAmerica", "SAmerica")
  region.seq <- rep(regions, 7)
  region.clust <- rep(regions, each=7)
  waypassage.rownames  <- paste(region.clust, region.seq, sep=" > ")
  waypassage <- matrix(NA, nrow=length(waypassage.rownames), ncol=nrow(waypoints))
  rownames(waypassage) <- waypassage.rownames
  colnames(waypassage) <- rownames(waypoints)
  
  #Africa endpoints, in order
  waypassage[2:7,1] <- 1
  waypassage[4,2] <- 2
  waypassage[5,3] <- 2
  waypassage[6:7,4] <- 2
  waypassage[6:7, 5] <- 3
  waypassage[7,6] <- 4
  waypassage[1:7,]
  #Middle East endpoints (no waypoints between Asia)
  waypassage[8,1] <- waypassage[11,2]<- waypassage[12,3] <- 1 
  waypassage[13,] <- c(NA, NA, NA, 1, 2, NA)
  waypassage[14,] <- c(NA, NA, NA, 1, 2, 3)
  #Asia enpoints
  waypassage[15,1] <- waypassage[18,2] <- waypassage[19,3] <- 1
  waypassage[20,] <- c(NA, NA, NA, 1, 2, NA)
  waypassage[21,] <- c(NA, NA, NA, 1, 2, 3)
  #Europe enpoints
  waypassage[22,1:2] <- c(2,1)
  waypassage[23,2] <- waypassage[24,2] <- 1
  waypassage[26, 2:3] <- c(1,2)
  waypassage[27,] <- c(NA, 1, NA, 2, 3, NA)
  waypassage[28,] <- c(NA, 1, NA, 2, 3, 4)
  #Oceania endpoints
  waypassage[c(29:32, 34:35), 3] <- 1
  waypassage[29,1] <- waypassage[32,2] <- waypassage[34:35, 4] <- 2
  waypassage[34:35,5] <- 3
  waypassage[35, 6] <- 4
  #North America endpoints
  waypassage[36:40,5] <- 1
  waypassage[36:40,4] <- 2
  waypassage[36,1] <- waypassage[39,2] <- waypassage[40,3] <- 3
  waypassage[42,6] <- 1
  #South America endpoints
  waypassage[43,] <- waypassage[44,] <- waypassage[45,] <-
    waypassage[46,] <- waypassage[47,] <- c(NA, NA, NA, 3, 2, 1)
  waypassage[43,1] <- waypassage[46, 2] <- waypassage[47,3] <- 4
  waypassage[48,6] <- 1
  
  # ways.regions: for all possible regional jumps, the first and last waypoint
  single.way.rows <- which(rowSums(waypassage, na.rm=TRUE)==1)
  single.ways <- waypassage[single.way.rows,]
  start.ends <- paste(region.clust, region.seq, sep="")[single.way.rows]
  single.waypoints <- data.frame(waypoint=apply(single.ways, 1, function(x) names(which(x==1))),
                                 StartEnd=start.ends)
  rownames(single.waypoints) <- names(single.ways)
  
  with.way.rows <- which(apply(!is.na(waypassage), 1, any)==TRUE)
  with.way <- waypassage[with.way.rows,]
  with.way[is.na(with.way)]<- 0
  first.way <- apply(with.way, 1, function(x) names(which(x==1)))
  last.way <- apply(with.way, 1, function(x) names(which(x==max(x))))
  start.ends <- paste(region.clust, region.seq, sep="")[with.way.rows]
  ways.regions <- data.frame(firstWay=first.way, lastWay=last.way,
                             StartEnd=start.ends)
  rownames(ways.regions) <- names(single.ways)
  
  #for multiwaypoint travel, the distance between the first and last waypoint
  interway.distances <- interway.km(waypoints,waypassage)
  
  #####################################################################
  
  # FUNCTION OPERATIONS (values depend on geoinfo)
  
  # latlongs: data frame of latitudes, longitudes, and regions for sampled groups
  rownames(geoinfo) <- geoinfo[,3]
  latlongs <- data.frame(lat=geoinfo[,1], long=geoinfo[,2],
                         Region=geoinfo[,4])
  latlongnames <- as.character(geoinfo[,3])
  rownames(latlongs) <- latlongnames
  
  # geodists: total distance between pop pairs is the sum of the columns
  pairwise.region.mat <- matrix(NA, nrow(latlongs), nrow(latlongs))
  rownames(pairwise.region.mat)<-colnames(pairwise.region.mat) <- rownames(latlongs)
  pop.clust <- rep(rownames(latlongs), nrow(latlongs))
  pop.seq <- rep(rownames(latlongs), each=nrow(latlongs))
  Reg12 <- paste(rep(latlongs$Region, nrow(latlongs)), 
                 rep(latlongs$Region, each=nrow(latlongs)), sep="")
  
  
  wayPts <- interway.distances$wayPoints[match(Reg12, interway.distances$wayStartEnd)]
  pops.regs.ways <- data.frame(Pop1=pop.clust, Pop2=pop.seq, Reg12=Reg12, 
                               wayPts = wayPts)
  
  geodists <- matrix(NA, nrow(pops.regs.ways), ncol=4)
  colnames(geodists) <- c("wayKm", "toWayA", "toWayB", "Pop2Pop")
  rownames(geodists) <- 
    paste(pops.regs.ways$Pop1, pops.regs.ways$Pop2, sep=" > ")
  geodists[,1] <- 
    interway.distances$wayKm[match(Reg12, 
                                   interway.distances$wayStartEnd)]
  
  for(i in 1:nrow(geodists))
  {
    PopA <- pops.regs.ways$Pop1[i]
    latlongA <- as.numeric(latlongs[which(rownames(latlongs)==PopA),1:2])
    PopB <- pops.regs.ways$Pop2[i]
    latlongB <- as.numeric(latlongs[which(rownames(latlongs)==PopB),1:2])
    
    if(pops.regs.ways$wayPts[i]==0) {
      #toWayA
      geodists[i,2] <- 0 
      #toWayB
      geodists[i,3] <- 0
      #Pop2Pop
      geodists[i,4] <- geodist.hav(latlongA,latlongB)
    }
    
    if(pops.regs.ways$wayPts[i] > 0) {
      waypoints.row <- match(pops.regs.ways$Reg12[i], 
                             ways.regions$StartEnd)
      wayA <- ways.regions$firstWay[waypoints.row]
      wayB <- ways.regions$lastWay[waypoints.row]
      wayA.latlong <- waypoints[which(rownames(waypoints)==wayA),]
      wayB.latlong <- waypoints[which(rownames(waypoints)==wayB),]
      geodists[i,2] <- geodist.hav(latlongA, wayA.latlong) 
      geodists[i,3] <- geodist.hav(latlongB, wayB.latlong) 
      geodists[i,4] <- 0
    }
  }
  
  # convert to pairwise distance matrix
  pw.dists <- matrix(apply(geodists, 1, sum), 
                     nrow=nrow(latlongs), ncol=nrow(latlongs))
  rownames(pw.dists) <- colnames(pw.dists) <- rownames(latlongs)
  
  return(pw.dists)
}
setwd("C:/Users/adria/OneDrive/Desktop")
setwd("C:/Users/adria/OneDrive/Desktop/PostDoc/distance_africa")
datcoords<-read.csv("coords.csv",header=TRUE)
datcoords
r<-distance.overland(datcoords)
write.csv(r, file = "distance_africa.csv")

### graphic###
dim <- ncol(r)
image(1:dim, 1:dim, r, axes = FALSE, xlab="", ylab="")
axis(1, 1:dim, datcoords[1:54,3], cex.axis = 0.5, las=3)
axis(2, 1:dim, datcoords[1:54,3], cex.axis = 0.5, las=1)
text(expand.grid(1:dim, 1:dim), sprintf("%0.1f", r), cex=0.6)
