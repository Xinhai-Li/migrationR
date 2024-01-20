# Process and analyze satellite tracking data, analyze bird migration patterns.
# Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com); 
#============================================================================================================================


#############################################################################################################################
#' Format standard Movebank data and add basic movement parameters
#'
#' @description Import standard Movebank data, export data with basic movement parameters
#' such as distance (between two adjacent points), speed, moving direction, changes of moving direction,
#' as well as Year, Month, Julian day, Hour of every record.
#'
#' @author Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com)
#'
#' @param data A data.frame of standard Movebank data
#'
#' @return
#'
#' @examples
#'
#' attach(movebankdata)
#' trackdata = as_trackdata(data = movebankdata, min_time_interval = 6)
#'
#' @import lubridate
#' @importFrom chron month.day.year
#' @importFrom argosfilter distanceTrack
#'
#' @export

as_trackdata <- function(data, min_time_interval = 6) {
  ## data cleaning...
  data = data[, c("individual.local.identifier", "timestamp", "location.long", "location.lat")]
  names(data) = c("ID", "Time", "Lon", "Lat")
  data = data[!is.na(data$Lon), ]
  ind_rare <- names(table(data$ID))[table(data$ID)<10]
  data = data[!(data$ID %in% ind_rare), ] # remove individuals with < 10 points
  ## Origin date (the ealiest date)
  # origin <- as.Date(sprintf("%s-01-01", year(min(data$Time))))

  ## Generate new variables..
  data$Time      <- as.POSIXlt(data$Time)
  # data$Time_diff <- as.Date(data$Time) - origin
  data$Year      <- month.day.year(data$Time)[[3]]
  data$Month     <- month.day.year(data$Time)[[1]]
  data$Day       <- yday(data$Time)
  data$Hour      <- hour(data$Time)
  Minute         <- minute(data$Time)
  data$Day_fine  <- data$Day + (data$Hour + Minute/60)/24


  myFun <- function(id) {
    S <- data[data$ID == id, ]
    n <- nrow(S)
    S <- S[order(S$Time), ]

    ## Dist
    S$Dist <- c(NA, distanceTrack(S$Lat, S$Lon))

    ## Time interval
    S$Time_interval <- c(NA, as.numeric(diff(as.numeric(S$Time))/3600))

    ## Speed
    S$Speed <- c(NA, S$Dist[-1]/as.numeric(diff(as.numeric(S$Time))/3600))

    ## Direction
    Direction <- acos((S$Lat[-1] - S$Lat[-n]) * 39946.79/360/S$Dist[-1]) * 360/2/pi
    id_23 <- which((S$Lon[-1] - S$Lon[-n]) < 0)
    if (length(id_23) > 0) {Direction[id_23] <- 360 - Direction[id_23] }
    Direction[which(is.nan(Direction))] <- 0
    S$Direction <- c(NA, Direction)

    ## North index
    S$North = abs((S$Direction -180)/180) # North index

    ## Changes in directions
    S$Redirect <- c(NA, NA, diff(Direction))
    id1 <- which(S$Redirect >  180)
    id2 <- which(S$Redirect < -180)
    if (length(id1) > 0) {S$Redirect[id1] <- S$Redirect[id1] - 180}
    if (length(id2) > 0) {S$Redirect[id2] <- S$Redirect[id2] + 180}

    ## Reset these variables to NA if "Time_interval" great than the threshold value;
    id <- which(S$Time_interval > min_time_interval)
    if (length(id) > 0) {
      S$Dist[id]      <- NA
      S$Speed[id]     <- NA
      S$Direction[id] <- NA
      S$Redirect[id]  <- NA
    }
    return(S)
  }

  tmp <- do.call(rbind, (lapply(unique(data$ID), myFun)))
  return(tmp)
}





# Plot key areas
##############################################################################################################
#' Plot breeding and wintering areas
#'
#' @description Plot breeding and wintering areas based on kernel density maps. The starting and ending dates
#' of the breeding peroid and wintering period are needed. The areas can be adjusted by the percentage of kernel areas.
#'
#' @author Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com)
#'
#' @param trackdata A data.frame of processed satellite tracking data.
#' @param ext.par The buffer distance (degree) outside the species occurrences ranges
#' @param breed.start The starting day (Julian day) of breeding migration
#' @param breed.end The ending day (Julian day) of breeding migration
#' @param winter.start The starting day (Julian day) of wintering migration
#' @param winter.end The ending day (Julian day) of wintering migration
#' @param lat.min.b The minimum latitude in the breeding season
#' @param lat.max.w The maximum latitude in the wintering season
#' @param breed.percent A vector of three percentage values defining the kernel density areas in the breeding season
#' @param winter.percent A vector of three percentage values defining the kernel density areas in the wintering season
#'
#' @return
#'
#' @examples
#'
#'  plot_breeding_wintering(trackdata=trackdata, ext.par=4, breed.start=100,
#'    breed.end=210, winter.start=1, winter.end=60,
#'    lat.min.b = 45, lat.max.w = 25,
#'    breed.percent = c(90, 70, 50), winter.percent = c(99.9, 90, 85))
#'
#' @importFrom maps map
#'
#' @import adehabitatHR
#' @import sp
#'
#' @export

plot_breeding_wintering = function(trackdata=trackdata, ext.par=4, breed.start=90,
                                   breed.end=210, winter.start=1, winter.end=60,
                                   lat.min.b = -70, lat.max.w = 90,
                                   breed.percent = c(90, 70, 50), winter.percent = c(90, 70, 50)){

  xlim = range(trackdata$Lon, na.rm=T); ylim = range(trackdata$Lat, na.rm=T)
  xlim[1] = xlim[1]-ext.par; xlim[2] = xlim[2]+ext.par
  ylim[1] = ylim[1]-ext.par; ylim[2] = ylim[2]+ext.par
  library(maps)
  map('world', mar = c(4,4,0,0), xlim=xlim, ylim=ylim, fill=F)
  # plot(shape, mar = c(4,4,0,0), xlim=xlim, ylim=ylim, fill=T, col="grey90")
  LatLon = trackdata[trackdata$Day>breed.start & trackdata$Day<breed.end, c('Lon', 'Lat')]
  LatLon = LatLon[LatLon$Lat > lat.min.b, ]
  LatLon = round(LatLon, 1) # thinning the points
  LatLon = unique.data.frame(LatLon)
  library(adehabitatHR)
  xy <- SpatialPoints(LatLon)
  kud <- kernelUD(xy)  # h = href is the default - ad hoc method for determining h
  tryCatch({ver <- getverticeshr(kud, percent = breed.percent[1]); plot(ver, col=adjustcolor("orange", alpha.f=0.5), add=T, border="NA")}, error = function(e){cat("ERROR: ", conditionMessage(e),"\n")})
  tryCatch({ver <- getverticeshr(kud, percent = breed.percent[2]); plot(ver, col=adjustcolor("orange", alpha.f=0.5), add=T, border="NA")}, error = function(e){cat(".")})
  tryCatch({ver <- getverticeshr(kud, percent = breed.percent[3]); plot(ver, col=adjustcolor("orange", alpha.f=0.5), add=T, border="NA")}, error = function(e){cat(".")})

  LatLon = trackdata[trackdata$Day>winter.start & trackdata$Day<winter.end, c('Lon', 'Lat')]
  LatLon = LatLon[LatLon$Lat < lat.max.w, ]
  LatLon = round(LatLon, 2) # thinning the points
  LatLon = unique.data.frame(LatLon)
  xy <- SpatialPoints(LatLon)
  kud <- kernelUD(xy)  # h = href is the default - ad hoc method for determining h
  tryCatch({ver <- getverticeshr(kud, percent = winter.percent[1]); plot(ver, col=adjustcolor("blue", alpha.f=0.5), add=T, border="NA")}, error = function(e){cat(".")})
  tryCatch({ver <- getverticeshr(kud, percent = winter.percent[2]); plot(ver, col=adjustcolor("blue", alpha.f=0.5), add=T, border="NA")}, error = function(e){cat(".")})
  tryCatch({ver <- getverticeshr(kud, percent = winter.percent[3]); plot(ver, col=adjustcolor("blue", alpha.f=0.5), add=T, border="NA")}, error = function(e){cat(".")})

  for (yr in sort(unique(trackdata$Year))) {
    dat_yr <- trackdata[trackdata$Year == yr, ]
    for (id in unique(dat_yr$ID)) {
      dat_yr_id <- dat_yr[dat_yr$ID == id, ]
      color_id <- sample(1:20, 1)
      # with(dat_yr_id, points(Lon, Lat, col = color_id, cex = 0.5))
      with(dat_yr_id,  lines(Lon, Lat, col = color_id, lwd = 0.5))
    }
  }
}
##############################################################################################################











# Annual movement distance
##############################################################################################################
#' Calculate annual movement distance (km) for each individual
#'
#' @description Sum the annual total movement distance for each individual
#'
#' @author Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com)
#'
#' @param trackdata A data.frame of processed satellite tracking data.
#' @param n The minimum number of days in a year for calculating the total annual movement distance
#'
#' @return
#'
#' @examples
#'
#' Dist = dist_annual(trackdata=trackdata, n = 200); Dist
#'
#' @export
#'

dist_annual = function(trackdata, n = 10){
  tmp_Day  <- with(trackdata, aggregate(list(numDay = Day),    list(ID = ID, Year = Year), function(x) {length(unique(x))}))
  tmp_Dist <- with(trackdata, aggregate(list(Distance = Dist), list(ID = ID, Year = Year), sum, na.rm = TRUE))
  tmp <- merge(tmp_Dist, tmp_Day, all.x = TRUE)
  tmp <- tmp[which(tmp$numDay > n), ]
  return(tmp)
}
##############################################################################################################



# trajectory plot
##############################################################################################################
#' Plot trajectories of all individuals in all years
#'
#' @description Plot two types of trajectories of all individuals in all years. Type "ind" shows different individuals in different color. Type "time" shows locating points in different colors determined by Julian days.
#'
#' @author Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com)
#'
#' @param trackdata A data.frame of processed satellite tracking data.
#' @param ext.par The buffer distance (degree) outside the species occurrences ranges
#' @param transparent The tranparency of colors for points
#' @param p.size The size of points
#' @param l.size The width of lines
#' @param type There are two options for the type. "ind" shows different individuals in different colors. "time" shows the points from Julian day 1 to 365 by gradient colors red, yellow and green.
#'
#' @return
#'
#' @examples
#'
#' plot_traj(trackdata = data, type = "individual")
#' plot_traj(trackdata = data, type = "chronic")
#'
#' @importFrom maps map
#'
#' @export
#'

plot_traj = function(trackdata, type="individual", buffer = 1, transparent = 0.1, p.size = 0.2, l.size = 1) {
  xlim <- range(trackdata$Lon, na.rm=T)
  ylim <- range(trackdata$Lat, na.rm=T)
  xlim <- xlim + c(-buffer, buffer)
  ylim <- ylim + c(-buffer, buffer)

  map('world', mar = c(4,4,0,0), xlim=xlim, ylim=ylim, fill=FALSE)
  for (yr in sort(unique(trackdata$Year))) {
    dat_yr <- trackdata[trackdata$Year == yr, ]
    for (id in unique(dat_yr$ID)) {
      dat_yr_id <- dat_yr[dat_yr$ID == id, ]

      if (type == "individual") {
        color_pt <- sample(1:20, 1)
        color_ln <- color_pt
      } else if (type == "chronic") {
        color_pt <- colorRampPalette(c("red", "yellow", "green"))(max(dat_yr_id$Day))[dat_yr_id$Day]
        color_ln <- color_pt
      }

      with(dat_yr_id, points(Lon, Lat, col = color_pt, cex = p.size))
      with(dat_yr_id,  lines(Lon, Lat, col = color_ln, lwd = l.size))
    }
  }
}
##############################################################################################################



# Calculate tracking duration
##############################################################################################################
#' Calculate tracking duration (days) for all individuals
#'
#' @description Calculate tracking duration (from the first day to the last day (or the present day)) for all individuals
#'
#' @author Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com)
#'
#' @param trackdata A data.frame of processed satellite tracking data.
#'
#' @return
#'
#' @examples
#'
#' Dur = duration(trackdata); head(Dur)
#'
#' @export
#'

duration <- function(trackdata) {
  tmp_duration <- with(trackdata, aggregate(list(Duration = Time), list(Individual = ID),
                                       FUN=function(x){return(as.numeric(diff(as.numeric(range(as.Date(x), na.rm = TRUE)))))}))
  tmp_N_record <- with(trackdata, aggregate(list(Duration = Time), list(Individual = ID), length))
  tmp <- cbind(tmp_duration, No_record = tmp_N_record[,2])
  return(tmp)
}
##############################################################################################################



# Plot tracking duration
##############################################################################################################
#' Plot tracking duration (days) for all individuals
#'
#' @description Plot tracking duration (from the first day to the last day (or the present day)) for all individuals
#'
#' @author Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com)
#'
#' @param trackdata A data.frame of processed satellite tracking data.
#'
#' @return
#'
#' @examples
#'
#' plot_track_duration(trackdata)
#'
#' @export

plot_track_duration = function(trackdata=trackdata, cex.lab=1, cex.axis=1){

  Time <- strptime(trackdata$Time, "%Y-%m-%d %H:%M:%S") #capital Year for 4 digits  # :%S
  Date_D = as.numeric(difftime(Time, min(Time), units='days'))
  trackdata = cbind(trackdata, Date_D)
  trackdata = trackdata[order(trackdata$Date_D), ]
  plot(c(0, max(trackdata$Date_D)),c(1, length(unique(trackdata$ID))), col="white", xlab="Date", ylab="Individuals",
       xaxt='n', cex.lab=cex.lab, cex.axis=cex.axis)
  at = c(1, floor(nrow(trackdata)/2), nrow(trackdata))
  axis(side=1, at=trackdata$Date_D[ at ], labels = as.Date(trackdata$Date_D[ at ], origin = min(Time)), cex.axis=cex.axis)

  IDs <- unique(trackdata$ID)
  for (i in 1:length(IDs)){
    Ind <- trackdata[trackdata$ID == IDs[i], ]
    lines(c(min(Ind$Date_D, na.rm=T), max(Ind$Date_D, na.rm=T)), c(i,i), lwd=2)
  }
}
##############################################################################################################






# Unimode method for identifying nest sites
##############################################################################################################
#' Identify nest site for each individual
#'
#' @description Determie the nest site for each individual based on the most used site in the breeding season
#'
#' @author Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com)
#'
#' @param trackdata A data.frame of processed satellite tracking data.
#' @param breed.S The starting day of the breeding season
#' @param breed.E The ending day of the breeding season
#' @param minimum.rec The minimum number needed for estimating nest site.
#'
#' @return
#'
#' @examples
#'
#' nest_locating(trackdata, breed.S = 130, breed.E = 180, minimum.rec = 50)
#'
#' # Check a single individual
#' IDs = names(table(trackdata$ID))
#' i=1; YEAR=2018
#' Ind = trackdata[trackdata$ID == IDs[i] & trackdata$Year==YEAR, ]
#' breed.S = 140; breed.E = 180
#' Ind = Ind[Ind$Day > breed.S & Ind$Day < breed.E, ]
#' plot(Ind$Lon, Ind$Lat)
#' lines(Ind$Lon, Ind$Lat, col="grey", lwd=.5)
#' LAT = round(Ind$Lat, 4); LON = round(Ind$Lon, 4)
#' LATLON = as.character(LAT*LON*10^8) # numeric would cause no-match
#' frq = sort(table(LATLON), decreasing=T)
#' LAT = Ind$Lat[LATLON ==  names(frq [frq==max(frq)] )][1]
#' LON = Ind$Lon[LATLON ==  names(frq [frq==max(frq)] )][1]
#' points(LON, LAT, col=2, pch=16)
#'
#' @importFrom maps map
#'
#' @export
#'

nest_locating = function(trackdata=trackdata, breed.S = breed.S, breed.E = breed.E, minimum.rec = 100){
  IDs = unique(trackdata$ID) # all individuals of a species
  No = length(IDs)*length(unique(trackdata$Year))
  NEST = data.frame(ID=1:No, Ind = NA, Year=NA, Nest_Lat=NA, Nest_Lon=NA) # for holding all results
  n = 1
  for (yr in sort(unique(trackdata$Year))) {
    dat_yr <- trackdata[trackdata$Year == yr, ]
    for (id in unique(dat_yr$ID)) {
      dat_yr_id <- dat_yr[dat_yr$ID == id, ]
      dat_yr_id = dat_yr_id[dat_yr_id$Day > breed.S & dat_yr_id$Day < breed.E, ] # further restrict records
      if (nrow(dat_yr_id) >= minimum.rec) {
        LAT = round(dat_yr_id$Lat, 4); LON = round(dat_yr_id$Lon, 4)
        LATLON = as.character(LAT*LON*10^8) # numeric would cause no-match
        NEST$Ind[n] = dat_yr_id$ID[1]
        NEST$Year[n] = yr
        frq = sort(table(LATLON), decreasing=T)
        NEST$Nest_Lat[n] = dat_yr_id$Lat[LATLON ==  names(frq [frq==max(frq)] )][1]
        NEST$Nest_Lon[n] = dat_yr_id$Lon[LATLON ==  names(frq [frq==max(frq)] )][1]
      }
      n = n+1
      print(n)
    }
  }
  NEST = NEST[!is.na(NEST$Nest_Lat), ]
  return(NEST)
}
##############################################################################################################







# Daily distance
##############################################################################################################
#' Calculate daily movement distance (km) for each individual
#'
#' @description Sum the daily total movement distance for each individual
#'
#' @author Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com)
#'
#' @param trackdata A data.frame of processed satellite tracking data.
#'
#' @return
#'
#' @examples
#'
#' Daily.dist = dist_daily(trackdata); head(Daily.dist)
#' plot(Daily.dist$Day, Daily.dist$Dist2, col=as.numeric(as.factor(Daily.dist$Individual)), xlab="Julian day", ylab="Flying distance (km)")
#'
#' @export
#'

dist_daily <- function(trackdata) {
  with(trackdata, aggregate(list(Distance = Dist), list(Individual = ID, Year = Year, Day = Day, Month=Month), sum, na.rm = TRUE))
}

##############################################################################################################






# plot daily movement distance
##############################################################################################################
#' Plot daily movement distance (km) for each individual
#'
#' @description Plot the daily total movement distance for each individual in the range of Julian day (1-365)
#'
#' @author Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com)
#'
#' @param trackdata A data.frame returned from function daily_dist().
#'
#' @return
#'
#' @examples
#'
#' Daily.dist = dist_daily(trackdata)
#' plot_daily_dist(Daily.dist)
#'
#' @export
#'

plot_daily_dist = function(data=Daily.dist){
  IDs = names(table(data$Individual))
  plot(data$Day, data$Dist, xlab="Julian day", ylab="Flying distance (km)",type='n')
  n=0
  for (yr in sort(unique(data$Year))) {
    for (id in unique(data$Individual)) {
      n=n+1
      tryCatch({
        Ind = data[data$Individual == id & data$Year==yr, ]
        Ind = Ind[order(Ind$Day), ]
        if (nrow(Ind) >1) {
          points(Ind$Day, Ind$Dist, col=rainbow(length(unique(data$Individual)))[n])
          lines( Ind$Day, Ind$Dist, col=rainbow(length(unique(data$Individual)))[n])
        }
      }, error = function(e){cat(".")})
    }
  }
}
##############################################################################################################




# migration timing, using a minimum distance of one day's movement to determine a starting day
##############################################################################################################
#' Determine the starting day and ending day of migration
#'
#' @description If the daily movement distance of an individual is over a threshold (i.e. 100km), the day is defined as the starting day.
#' The last day of a series of continuous days that movement distance over the threshold is the ending day. Similarly, if the movement distance
#' in one hour on the migration day (starting or ending day) is over a threshold (i.e. 10km), the "hour" is the starting or ending hour of the movement.
#'
#' @author Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com)
#'
#' @param trackdata A data.frame of processed satellite tracking data.
#'
#' @param dist_min_day The threshold movement distance (e.g. 100km / day) for determining the start date of a migration.
#' @param dist_min_hour The threshold movement distance (e.g. 10km / day) for determining the start hour of a migration.
#' @param dist_outlier If the movement distance within one hour is over 150km, it is defined as an outlier (caused by missing records).
#'
#' @return
#'
#' @examples
#'
#' timing = mig_timing(trackdata=trackdata, dist_min_day = 100, dist_min_hour = 10, dist_outlier = 150)
#'
#' @export
#'

mig_timing <- function(trackdata, dist_min_day = 100, dist_min_hour = 10, dist_outlier = 150) {
  Dat <- unique(trackdata[, c("Year", "ID")])

  Dat$Breed_Start_Day  <- NA
  Dat$Breed_Start_Hour <- NA
  Dat$Breed_End_Day    <- NA
  Dat$Breed_End_Hour   <- NA
  Dat$Winter_Start_Day <- NA
  Dat$Winter_Start_Hour<- NA
  Dat$Winter_End_Day   <- NA
  Dat$Winter_End_Hour  <- NA

  for(i in 1:nrow(Dat)) {

    dat_org <- trackdata[which(trackdata$ID == Dat$ID[i] & trackdata$Year == Dat$Year[i]), ]
    dat_agg <- with(dat_org, aggregate(list(Dist = Dist), list(Day = Day), sum, na.rm = TRUE))
    dat <- dat_agg[dat_agg$Dist > dist_min_day, ]

    if (nrow(dat) > 1) {
      ## remove single day long distance movement
      if ((dat$Day[2] - dat$Day[1])>1) {dat <- dat[-1, ]}

      days_1 <- dat$Day[which(dat$Day <= 180)]
      days_2 <- dat$Day[which(dat$Day >  180)]

      ## Breed
      if (length(days_1) > 0) {
        ## Start
        Dat$Breed_Start_Day[i] <- days_1[1]
        ## Hour
        dd <- id <- NULL
        dd <- dat_org[dat_org$Day == Dat$Breed_Start_Day[i],]
        id <- which(dd$Dist > dist_min_hour & dd$Dist < dist_outlier)
        if (length(id) > 0) { Dat$Breed_Start_Hour[i] <- dd$Hour[min(id)] }
        ## End
        if (length(days_1) > 1) {
          Dat$Breed_End_Day[i] <- max(days_1)
          ## Hour
          dd <- id <- NULL
          dd <- dat_org[dat_org$Day == Dat$Breed_End_Day[i],]
          id <- which(dd$Dist > dist_min_hour & dd$Dist < dist_outlier)
          if (length(id) > 0) { Dat$Breed_End_Hour[i] <- dd$Hour[max(id)] }
        }
      }

      ## Winter
      if (length(days_2) > 0) {
        ## Start
        Dat$Winter_Start_Day[i] <- days_2[1]
        ## Hour
        dd <- id <- NULL
        dd <- dat_org[dat_org$Day == Dat$Winter_Start_Day[i],]
        id <- which(dd$Dist > dist_min_hour & dd$Dist < dist_outlier)
        if (length(id) > 0) { Dat$Winter_Start_Hour[i] <- dd$Hour[min(id)] }
        ## End
        if (length(days_2) > 1) {
          Dat$Winter_End_Day[i] <- max(days_2)
          ## Hour
          dd <- id <- NULL
          dd <- dat_org[dat_org$Day == Dat$Winter_End_Day[i],]
          id <- which(dd$Dist > dist_min_hour & dd$Dist < dist_outlier)
          if (length(id) > 0) { Dat$Winter_End_Hour[i] <- dd$Hour[max(id)] }
        }
      }

    }
  }
  return(Dat)
}
##############################################################################################################














# Plot moving directions
##############################################################################################################
#' Plot moving directions of an individual in a year
#'
#' @description Use calcualted moving directions by function format.data() to plot moving directions of an individual in a year
#'
#' @author Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com)
#'
#' @param ind Locating records of an individual in a year, from the processed satellite tracking data.
#'
#' @return
#'
#' @examples
#'
#' ind = trackdata[trackdata$ID=="蓑羽鹤06-683·BFU076" & trackdata$Year == 2018,]
#' plot_direction(ind)
#'
#' @importFrom plotrix polar.plot
#'
#' @export
#'

plot_direction = function(ind = ind){
  library(plotrix)
  ind = ind[!ind$Direction>360, ]
  oldpar <- polar.plot(log(ind$Speed+1), ind$Direction, main="Fly direction",
                       radial.lim = c(0, 5), start=90, clockwise=TRUE, lwd=3,
                       line.col=colorRampPalette(c("green","yellow","red"))(365)[ind$Day])
  par(oldpar) # reset everything
}
##############################################################################################################



# Plot movement trajectories of an individual at series periods
##############################################################################################################
#' Plot movement trajectories of an individual at any given number of periods
#'
#' @description Use calcualted moving directions by function format.data() to plot moving directions of an individual in a year
#'
#' @author Xinhai Li (xinhai_li_edu@126.com)
#'
#' @param ind Locating records of an individual in a year, from the processed satellite tracking data.
#' @param seg the number of trajectory segments
#' @param label a boolean value determining whether plot the labels of days
#'
#' @return
#'
#' @examples
#'
#' ind = trackdata[trackdata$ID==trackdata$ID[1] & trackdata$Year == trackdata$Year[1],]
#' plot_traj_segments(data=ind, seg=6, label=F)
#'
#' @importFrom calibrate textxy
#'
#' @export
#'

plot_traj_segments = function(ind=ind, seg=9, label=T){
  par(mfrow=c(round(sqrt(seg)+0.49), round(sqrt(seg))))
  duration <- as.numeric(diff(as.numeric(range(as.Date(ind$Time), na.rm = TRUE))))
  interval = round(duration/seg+0.5)
  segments = split(ind, factor(sort(rank(row.names(ind))%%seg))) # split dataframe

  for (i in 1:seg){
      segment = as.data.frame(segments[[i]])
      # names(segment) = c("ID","Time","Lon", "Lat", "Year","Month", "Day", "Hour","Day_fine", "Dist",
      #                  "Dist2", "Speed", "Direction", "Redirect")
      plot(segment$Lon, segment$Lat, pch=19, xlab='',ylab='',
           main=paste(as.Date(segment$Time[1]), as.Date(segment$Time[nrow(segment)]), sep=" to "),
           col = colorRampPalette(c("red","yellow","green"))(nrow(segment)))
      lines(segment$Lon, segment$Lat, col="grey",lwd=.5)
      if (label==T) textxy(segment$Lon, segment$Lat, labs=segment$Day, cex = 0.5, col = "black", m = c(0, 0)) #library(calibrate)
  }
}





# Compare movement speed between breeding season and wintering season
##############################################################################################################
#' Compare movement speed between breeding season and wintering season
#'
#' @description Apply a mixed effect model to compare movement speed at breeding season and wintering season, using hour (24h)
#'  as covariate (linear term, quadratic term and interaction term are included). The results show the t test table of fixed variables and terms,
#'  and Marginal R square (based on fixed effect) and Conditional R square (based on both fixed and random effects).
#'
#' @author Xinhai Li (xinhai_li_edu@126.com)
#'
#' @param trackdata A data.frame of processed satellite tracking data.
#' @param breed.start The starting day (Julian day) of breeding migration
#' @param breed.end The ending day (Julian day) of breeding migration
#' @param winter.start The starting day (Julian day) of wintering migration
#' @param winter.end The ending day (Julian day) of wintering migration
#'
#' @return
#'
#' @examples
#'
#' speed_mixedmodel(trackdata=trackdata, breed.start = 150, breed.end = 180, winter.start=1, winter.end=31)
#'
#' @importFrom MASS stepAIC
#' @importFrom nlme lme
#'
#' @export
#'
speed_mixedmodel = function(trackdata=trackdata, breed.start = 150, breed.end = 180, winter.start=1, winter.end=31){
  breed <- trackdata[trackdata$Day >= breed.start & trackdata$Day < breed.end, ]
  breed <- cbind(breed, Status="Breed")
  winter<- trackdata[trackdata$Day >= winter.start & trackdata$Day < winter.end, ]
  winter <- cbind(winter, Status="Winter")
  Dat <- rbind(breed, winter)
  Dat <- Dat[!is.na(Dat$Speed), ]
  Dat <- Dat[!is.infinite(Dat$Speed), ]
  fit <- lme(Speed ~ Status*Hour + I(Hour^2), random = ~1|ID, data = Dat, method="ML") # method="REML" does not support stepAIC()
  fit <- MASS:::stepAIC(fit)
  out <- list()
  out[[1]] <- summary(fit)$tTable
  out[[2]] <- r.squaredGLMM(fit)
  names(out) <- c("T test table for fixed variables and terms", "Marginal R square and Conditional R square")
  return(out)
}





# CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#' Provide pseudo-absent points and derive environmental variables
#'
#' @description This function provides pseudo-absent points using user defined sample size,
#'  derives values of environmental variables at occurrences and pseudo-absent points.
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param species A data.frame of species occurrences with columns "Lon", "Lat", "Time"
#'
#' @param buffer A value of distance (unit: degree) defining the width of buffer zone outside the occurrences
#'
#' @param absence The number of pseudo-absent points on the longitude side and the latitude side
#'
#' @param Envlayers A RasterBrick (a multi-layer raster object)
#'
#' @return Return a data.frame with values of environmental variables at occurrences and pseudo-absent points.
#'
#' @examples
#'  attach(kiang) # load occurrences data
#'  attach(BioClim)
#'  Data = getEnvDatakiang, buffer=0.5, absence=30, Envlayers=BioClim); head(Data)
#'
#' @importFrom raster extract
#' @importFrom raster nlayers
#' 
#' @export
#'

getEnvData = function(species, buffer, absence, Envlayers){
  species.loc = species[,  c('Lon','Lat')]
  lon = seq(min(species.loc[,1]) - buffer, max(species.loc[,1]) + buffer, length.out = absence)
  lat = seq(min(species.loc[,2]) - buffer, max(species.loc[,2]) + buffer, length.out = absence)
  back = expand.grid(lon, lat)
  bak <- raster::extract(Envlayers, back)
  bak = as.data.frame(bak)
  bak = cbind(Name='absent', Size=0, bak)
  bak = cbind(bak, back)
  no.layers = nlayers(Envlayers)
  names(bak)[(no.layers+3):(no.layers+4)] = c('Lon','Lat') 
  bak = bak[!is.na(bak$Elevation),]
  spe = raster::extract(Envlayers, species.loc) #Lon MUST be the first column
  spe = as.data.frame(spe)
  spe = cbind(Name="Crane", Size=1, spe)
  spe = cbind(spe, species.loc)
  ENV = rbind(spe, bak)
  return(ENV)
}




# Run hetero-occurrence species distribution models (HOSDMs)
##############################################################################################################
#' Run hetero-occurrence species distribution models (HOSDMs)
#'
#' @description HOSDM is a novel species distribution modeling approach that is specifically designed to distinguish 
#' between various types of occurrences. Since White-naped Cranes rely on two separate habitat types, croplands and water areas, 
#' for their daily activities, we categorized the occurrences into two groups and executed HOSDMs accordingly.
#'
#' @author Xinhai Li (xinhai_li_edu@126.com)
#' 
#' @param trackdata A data.frame of processed satellite tracking data with columns "Lon", "Lat", "Time"#'
#' @param species A subset of trackdata
#' @param prediction A boolean value (T/F) that signifies whether habitat suitability prediction is being performed
#' @param buffer define the buffer distance (in degree) for the extent of layers of environmental variables 
#' @param absence The number of pseudo-absence points are absence*absence
#' @param Envlayer A RasterBrick (a multi-layer raster object) for environmental variables
#'
#' @return
#'
#' @examples
#' # The number of points during the foraging stage (Points1) and the roosting stage (Points0), respectively.
#' # The R square values for species distribution models during the two stages (R0 and R1).
#' # The mean distances between adjacent points during the foraging stage (Dist1) and the roosting stage (Dist0).
#' # The distance between the centers of the forage area and the roost area (Dist).
#' results = HOSDM(trackdata = HOSDMdata, prediction=F, buffer=0.1, absence=30, Envlayer=BioClim)
#' results[[1]] 
#' 
#' # Predicted habitat suitability
#' plot(results[[2]][[1]], main=paste("Roosting habitat for individual", ind[1], sep=" " ))
#' plot(results[[3]][[1]], main=paste("Foraging habitat for individual", ind[1], sep=" " ))
#'
#' @importFrom raster extract
#' @importFrom raster nlayers
#' @import randomForest
#' @importFrom argosfilter distance
#'
#' @export
#'
HOSDM = function(trackdata = trackdata, prediction=F, buffer=0.1, absence=30, Envlayer=BioClim){
  
  getEnvData = function(species, buffer, absence, Envlayers){
    species.loc = species[,  c('Lon','Lat')]
    lon = seq(min(species.loc[,1]) - buffer, max(species.loc[,1]) + buffer, length.out = absence)
    lat = seq(min(species.loc[,2]) - buffer, max(species.loc[,2]) + buffer, length.out = absence)
    back = expand.grid(lon, lat)
    bak <- raster::extract(Envlayers, back)
    bak = as.data.frame(bak)
    bak = cbind(Name='absent', Size=0, bak)
    bak = cbind(bak, back)
    no.layers = nlayers(Envlayers)
    names(bak)[(no.layers+3):(no.layers+4)] = c('Lon','Lat') 
    bak = bak[!is.na(bak$Elevation),]
    spe = raster::extract(Envlayers, species.loc) #Lon MUST be the first column
    spe = as.data.frame(spe)
    spe = cbind(Name="Crane", Size=1, spe)
    spe = cbind(spe, species.loc)
    ENV = rbind(spe, bak)
    return(ENV)
  }
  
  ind = unique(trackdata$ID)
  N = 0
  out = data.frame(No=1:(length(ind)*length(unique(trackdata$Winter))), Ind=NA, Winter=NA, Points0=NA, Points1=NA, R0=NA, R1=NA, Dist0=NA, Dist1=NA, Dist=NA)
  pred0 =list(); pred1 =list()
  
  for (i in 1:length(ind)){
    for (j in min(trackdata$Winter):max(trackdata$Winter)){
      
      tryCatch({
        
        D0 = trackdata[trackdata$ID==ind[i] & trackdata$Winter==j & trackdata$Hetero==0, ]
        Data0 = getEnvData(D0, buffer=buffer, absence=absence, Envlayer=Envlayer)
        RF0 <- randomForest(Size ~ Landcover + Elevation + Footprint + Bio_1+ Bio_4 +Bio_12+ Bio_15 +Solar1 +Vapor1 +Wind1+Lat+Lon, 
                            data=Data0, importance=TRUE, ntree=500)
        
        D1 = trackdata[trackdata$ID==ind[i] & trackdata$Winter==j & trackdata$Hetero==1, ]
        Data1 = getEnvData(D1, buffer=0.1, absence=40, Envlayer=BioClim)
        RF1 <- randomForest(Size ~ Landcover + Elevation + Footprint + Bio_1+ Bio_4 +Bio_12+ Bio_15 +Solar1 +Vapor1 +Wind1+Lat+Lon, 
                            data=Data1, importance=TRUE, ntree=500)      
        
        N = N+1
        out$No[N] = N
        out$Ind[N] = ind[i]
        out$Winter[N] = j
        out$Points0[N] = nrow(D0)
        out$Points1[N] = nrow(D1)
        out$R0[N] = max(RF0$rsq)
        out$R1[N] = max(RF1$rsq)
        dist0 = aggregate(Dist ~ Day, sum, data=D0)
        out$Dist0[N] = mean(dist0$Dist)
        dist1 = aggregate(Dist ~ Day, sum, data=D1)
        out$Dist1[N] = mean(dist1$Dist)      
        out$Dist[N] <- argosfilter::distance(mean(D0$Lat), mean(D1$Lat), mean(D0$Lon), mean(D1$Lon)) # library(argosfilter)
        
        if (prediction){ # predict habitat suitability
          pred0[[N]] <- predict(BioClim, RF0, type="response")
          pred1[[N]] <- predict(BioClim, RF1, type="response")
        }
        
        print(paste(ind[i], j, sep=" "))
        
      }, error = function(e){cat("ERROR: ", conditionMessage(e),"\n")})  
    }
  }
  output = list(out, pred0, pred1)
  return(output)
}


