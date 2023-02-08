# Process and analyze satellite tracking data
# Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com)
#============================================================================================================================

# Import Movebank data
#############################################################################################################################
#' Format standard Movebank data and add basic movement parameters
#'
#' @description Import standard Movebank data, export data with basic movement parameters
#' such as distance (between two adjacent points), speed, moving direction, changes of moving direction,
#' as well as Year, Month, Julian day, Hour of every record.
#'
#' @author Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com)
#'
#' @param trackdata A data.frame of standard Movebank data
#'
#' @return
#'
#' @examples
#'
#' attach(trackdata)
#' data = as.trackdata(trackdata = trackdata)
#'
#' @import lubridate
#' @importFrom chron month.day.year
#'
#' @export

as.trackdata <- function(trackdata) {
  ## data cleaning...
  trackdata = trackdata[, c("individual.local.identifier", "timestamp", "location.long", "location.lat")]
  names(trackdata) = c("ID", "Time", "Lon", "Lat")
  trackdata = trackdata[!is.na(trackdata$Lon), ]
  ind_rare <- names(table(trackdata$ID))[table(trackdata$ID)<10]
  trackdata = trackdata[!(trackdata$ID %in% ind_rare), ] # remove rare individuals
  ## Origin date (the earliest date in the first year)
  # origin <- as.Date(sprintf("%s-01-01", year(min(trackdata$Time))))

  ## Generate new variables..
  trackdata$Time      <- as.POSIXlt(trackdata$Time)
  # trackdata$Time_diff <- as.Date(trackdata$Time) - origin
  trackdata$Year      <- month.day.year(trackdata$Time)[[3]]
  trackdata$Month     <- month.day.year(trackdata$Time)[[1]]
  trackdata$Day       <- yday(trackdata$Time)
  trackdata$Hour      <- hour(trackdata$Time)
  trackdata$Day_fine  <- trackdata$Day + trackdata$Hour/24 + minute(trackdata$Time)/60/24

  myFun <- function(id) {
    S <- trackdata[trackdata$ID == id, ] # track data of one individual
    n <- nrow(S)
    S <- S[order(S$Time), ]

    ## Distance between two points
    Dist <- acos(sin(S$Lat[-n] * pi/180) * sin(S$Lat[-1] * pi/180) +
                   cos(S$Lat[-n] * pi/180) * cos(S$Lat[-1] * pi/180) * cos(S$Lon[-1] * pi/180 - S$Lon[-n] * pi/180)) * 6371
    Dist[which(is.nan(Dist))] <- 0
    S$Dist <- c(NA, Dist)

    ## Dist2
    # S$Dist2 <- c(NA, distanceTrack(S$Lat, S$Lon))
    Dist2 <- (((S$Lat[-1] - S$Lat[-n]) *  39946.79/360)^2 + ((S$Lon[-1] - S$Lon[-n]) * pi * 12756.32/360 * cos(S$Lat[-n] * pi * 2/360))^2)^0.5
    S$Dist2 <- c(NA, Dist2)

    ## Speed
    S$Speed <- c(NA, S$Dist2[-1]/as.numeric(diff(as.numeric(S$Time))/3600))

    ## Direction, moving direction calculated from two adjacent points. The north direction is 0 or 360.
    Direction <- acos((S$Lat[-1] - S$Lat[-n]) * 39946.79/360/S$Dist2[-1]) * 360/2/pi
    id_23 <- which((S$Lon[-1] - S$Lon[-n]) < 0)
    if (length(id_23) > 0) {Direction[id_23] <- 360 - Direction[id_23] }
    Direction[which(is.nan(Direction))] <- 0
    S$Direction <- c(NA, Direction)

    ## Redirect, changes in directions calculated from two adjacent directions based on three adjacent points.
    S$Redirect <- c(NA, NA, diff(Direction))
    id1 <- which(S$Redirect >  180)
    id2 <- which(S$Redirect < -180)
    if (length(id1) > 0) {S$Redirect[id1] <- S$Redirect[id1] - 180}
    if (length(id2) > 0) {S$Redirect[id2] <- S$Redirect[id2] + 180}

    return(S)
  }

  ## Do parallel calculation
  # if(require(foreach) & require(doParallel2)) {
  #   registerDoParallel(cores = detectCores()-2)
  #   tmp <- bind_rows(foreach(id = unique(trackdata$ID)) %dopar% { myFun(id) })
  # } else {
  # tmp <- bind_rows(lapply(unique(trackdata$ID), myFun))
  # }

  tmp <- do.call(rbind, (lapply(unique(trackdata$ID), myFun)))

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
  tmp_Day  <- with(trackdata, aggregate(list(numDay = Day),     list(ID = ID, Year = Year), function(x) {length(unique(x))}))
  tmp_Dist <- with(trackdata, aggregate(list(Distance = Dist2), list(ID = ID, Year = Year), sum, na.rm = TRUE))
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

plot_traj = function(trackdata, type="individual", buffer = 4, transparent = 0.1, p.size = 0.2, l.size = 1) {
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

dist_daily = function(trackdata=trackdata){
  IDs = names(table(trackdata$ID)) # all individuals of a species
  DIST = list()
  n = 0
  # DAY = 1:365; DAY = as.data.frame(DAY)
  for (yr in sort(unique(trackdata$Year))) {
    for (id in unique(trackdata$ID)) {
      tryCatch({
      Ind = trackdata[trackdata$ID == id & trackdata$Year==yr, ]
        Dist = aggregate(Dist2 ~ Day, sum, data=Ind)
        Dist = data.frame(Dist, Individual=id, Year=yr)
        # DAYS = merge(DAY, Dist, by.x = "DAY", by.y = "Day", all.x=T)
        # DAYS = data.frame(DAYS, Year=YEAR)
      }, error = function(e){cat(".")})
      n = n+1
      DIST[[n]] = Dist
    }
  }
  DD = DIST[[1]]
  for (i in 2:length(DIST)){
    DD = rbind(DD, DIST[[i]])
  }
  DD = DD[!is.na(DD$Individual),]

  return(DD)
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
  plot(data$Day, data$Dist2, xlab="Julian day", ylab="Flying distance (km)",type='n')
  n=0
  for (yr in sort(unique(data$Year))) {
    for (id in unique(data$Individual)) {
      n=n+1
      tryCatch({
        Ind = data[data$Individual == id & data$Year==yr, ]
        Ind = Ind[order(Ind$Day), ]
        if (nrow(Ind) >1) {
          points(Ind$Day, Ind$Dist2, col=rainbow(length(unique(data$Individual)))[n])
          lines( Ind$Day, Ind$Dist2, col=rainbow(length(unique(data$Individual)))[n])
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
#' @description If the daily movement distance of an individual is over a threshold (i.e. 100km), the day is defined as the starting day. The last day of a series of continuous days that movement distance over the threshold is the ending day.
#'
#' @author Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com)
#'
#' @param trackdata A data.frame of processed satellite tracking data.
#'
#' @param  min.dist The threshold movement distance (e.g. 100km / day) for determining the start of a migration.
#'
#' @return
#'
#' @examples
#'
#' dates = mig_day(trackdata, min.dist=100) ## minimum distance 100 km
#'
#' @export
#'

mig_day = function(trackdata=trackdata, min.dist = min.dist){
  Mig.day = unique(trackdata[, c("ID", "Year")])
  Mig.day = data.frame(Mig.day, Breed_S_Day = NA, Breed_E_Day = NA, Winter_S_Day = NA, Winter_E_Day = NA)

  for (i in 1:nrow(Mig.day)){
    tryCatch({
      Ind = trackdata[trackdata$ID == Mig.day$ID[i] & trackdata$Year == Mig.day$Year[i], ]
      Dist = aggregate(Dist2 ~ Day, sum, data=Ind)
      move = Dist[Dist$Dist2 >100, ] # 100 km as a threshold
      if (nrow(move >1)) {
        if ((move$Day[2]-move$Day[1]) >1) move = move[-1, ] # remove single day long distance movement
        Mig.day$Breed_S_Day[i]  = min(move[move$Day <= 180, 1])
        Mig.day$Winter_S_Day[i] = min(move[move$Day >  180, 1])
        Mig.day$Breed_E_Day[i]  = max(move[move$Day <= 180, 1])
        Mig.day$Winter_E_Day[i] = max(move[move$Day >  180, 1])
      }
    }, error = function(e){cat(".")})
    print(i)
  }
  Mig.day$Breed_S_Day[is.infinite(Mig.day$Breed_S_Day)] <- NA
  Mig.day$Breed_E_Day[is.infinite(Mig.day$Breed_E_Day)] <- NA
  Mig.day$Winter_S_Day[is.infinite(Mig.day$Winter_S_Day)] <- NA
  Mig.day$Winter_E_Day[is.infinite(Mig.day$Winter_E_Day)] <- NA
  return(Mig.day)
}
##############################################################################################################





# migration timing - hour
##############################################################################################################
#' Determine the starting hour and ending hour on the starting and ending day of migration
#'
#' @description If the hourly movement distance of an individual is over a threshold (i.e. 5km), the hour is defined as the starting hour. The ending our is define as the last hour of a series of continuous movement hours in the last day.
#'
#' @author Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com)
#'
#' @param trackdata A data.frame of processed satellite tracking data.
#' @param Mig.day A data.frame returned from function mig_day().
#' @param outlier.dist The threshold distance (e.g. 150 km) for determining an outlier (caused by missing records). If the distance > 150 km, it is not likely the distance that the bird flied in an hour.
#' @param min.move The threshold distance (e.g. 10 km) for determining the start time of migration. Whenever a bird flies over 10 km in an hour, it starts migration.
#'
#' @return
#'
#' @examples
#'
#' DD = mig_hour(trackdata=trackdata, Mig.day = dates, outlier.dist = 150, min.move = 10)
#'
#' @export
#'

mig_hour = function(trackdata=trackdata, Mig.day=Mig.day, outlier.dist = outlier.dist, min.move = min.move){
  Mig.hour = data.frame(Mig.day$ID, Mig.day$Year, Breed_S_Hour=NA, Breed_E_Hour=NA,
                        Winter_S_Hour=NA, Winter_E_Hour=NA)
  Mig.s = list(nrow(Mig.day))
  hours = numeric(nrow(Mig.day))

  for (i in 1: nrow(Mig.day)){
    tryCatch({
      Mig.s[[i]] = trackdata[trackdata$ID==Mig.day$ID[i] & trackdata$Year==Mig.day$Year[i] & trackdata$Day==Mig.day$Breed_S_Day[i],]
      Mig.s[[i]] = Mig.s[[i]][Mig.s[[i]]$Dist < outlier.dist,] # points with distance over 150 km was not accumulated by one hour
      hours[i] = Mig.s[[i]][Mig.s[[i]]$Dist > min.move, ]$Hour[1] # Distance > 10 km means start to migrate
    }, error = function(e){cat(".")})
  }
  Mig.hour$Breed_S_Hour = hours

  Mig.s = list(nrow(Mig.day))
  hours = numeric(nrow(Mig.day))
  for (i in 1: nrow(Mig.day)){
    tryCatch({
      Mig.s[[i]] = trackdata[trackdata$ID==Mig.day$ID[i] & trackdata$Year==Mig.day$Year[i] & trackdata$Day==Mig.day$Breed_E_Day[i],]
      Mig.s[[i]] = Mig.s[[i]][Mig.s[[i]]$Dist < outlier.dist,] # points with distance over 150 km was not accumulated by one hour
      hours[i] = max(Mig.s[[i]][Mig.s[[i]]$Dist > min.move, ]$Hour) # Distance > 10 km means start to migrate
    }, error = function(e){cat(".")})
  }
  Mig.hour$Breed_E_Hour = hours

  Mig.s = list(nrow(Mig.day))
  hours = numeric(nrow(Mig.day))
  for (i in 1: nrow(Mig.day)){
    tryCatch({
      Mig.s[[i]] = trackdata[trackdata$ID==Mig.day$ID[i] & trackdata$Year==Mig.day$Year[i] & trackdata$Day==Mig.day$Winter_S_Day[i],]
      Mig.s[[i]] = Mig.s[[i]][Mig.s[[i]]$Dist < outlier.dist,] # points with distance over 150 km was not accumulated by one hour
      hours[i] = Mig.s[[i]][Mig.s[[i]]$Dist > min.move, ]$Hour[1] # Distance > 10 km means start to migrate
    }, error = function(e){cat(".")})
  }
  Mig.hour$Winter_S_Hour = hours

  Mig.s = list(nrow(Mig.day))
  hours = numeric(nrow(Mig.day))
  for (i in 1: nrow(Mig.day)){
    tryCatch({
      Mig.s[[i]] = trackdata[trackdata$ID==Mig.day$ID[i] & trackdata$Year==Mig.day$Year[i] & trackdata$Day==Mig.day$Winter_E_Day[i],]
      Mig.s[[i]] = Mig.s[[i]][Mig.s[[i]]$Dist < outlier.dist,] # points with distance over 150 km was not accumulated by one hour
      hours[i] = max(Mig.s[[i]][Mig.s[[i]]$Dist > min.move, ]$Hour) # Distance > 10 km means start to migrate
    }, error = function(e){cat(".")})
  }
  Mig.hour$Winter_E_Hour = hours

  return(Mig.hour)
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

plot_traj_segments = function(data=ind, seg=9, label=T){
  par(mfrow=c(round(sqrt(seg)+0.49), round(sqrt(seg))))
  duration <- as.numeric(diff(as.numeric(range(as.Date(ind$Time), na.rm = TRUE))))
  interval = round(duration/seg+0.5)
  segments = split(ind, factor(sort(rank(row.names(ind))%%seg))) # split dataframe

  for (i in 1:seg){
      segment = as.data.frame(segments[[i]])
      names(segment) = c("ID","Time","Lon", "Lat", "Year","Month", "Day", "Hour","Day_fine", "Dist",
                        "Dist2", "Speed", "Direction", "Redirect")
      plot(segment$Lon, segment$Lat, pch=19, xlab='',ylab='',
           main=paste(as.Date(segment$Time[1]), as.Date(segment$Time[nrow(segment)]), sep=" to "),
           col = colorRampPalette(c("red","yellow","green"))(nrow(segment)))
      lines(segment$Lon, segment$Lat, col="grey",lwd=.5)
      if (label==T) textxy(segment$Lon, segment$Lat, labs=segment$Day, cex = 0.5, col = "black", m = c(0, 0)) #library(calibrate)
  }
}








