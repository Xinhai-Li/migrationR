# Process and analyze satellite tracking data
# Xinhai Li (xinhai_li_edu@126.com)
#============================================================================================================================


# Import Movebank data
#############################################################################################################################
#' Format standard Movebank data and add basic movement parameters
#'
#' @description Import standard Movebank data, export data with basic movement parameters such as distance (between two adjacent points), speed, moving direction, changes of moving direction, as well as Year, Month, Julian day, Hour of every record
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param D A data.frame of standard Movebank data
#'
#' @return
#'
#' @examples
#'
#'  attach(D)
#'  data = format.data(D = D, Time.s = '2000-1-1 0:0:0')
#'
#' @importFrom lubridate yday
#' @importFrom chron month.day.year
#'
#' @export
#'
format.data = function(D=D, Time.s=Time.s){
  D = D[!is.na(D$location.long), ]
  D = D[, c("timestamp","location.long", "location.lat", "individual.local.identifier")]
  Time <- strptime(D$timestamp, "%Y-%m-%d %H:%M:%S") #capital Year for 4 digits  # :%S
  Time.s = strptime(Time.s, "%Y-%m-%d %H:%M:%S") # a starting time
  Time_diff = as.numeric(difftime(Time, Time.s, units='days'))
  library(lubridate)
  library(chron)
  year = month.day.year(Time)[3]
  month = month.day.year(Time)[1] # library(chron) #Julian  day to month #0:364
  day = yday(D$timestamp)
  Day_fine = (Time_diff/365 - floor(Time_diff/365)) *365 -1 #start from Jan. 1  # Julian day: 0:364
  Hour = hour(Time)
  names(D) = c("Time","Lon", "Lat", "ID")
  D = data.frame(D, Time_d = Time_diff, Year=year, Month=month, Day=day, Day_fine=Day_fine, Hour=Hour)
  D = D[order(D$Time_d), ]
  pop = list()
  list = as.factor(unique(D$ID)) # list of Individuals

  for (K in 1:length(list)) {
    S = D[D$ID==list[K], ]; nrow(S) #replace for each ind.
    S = cbind(S, Dist=NA, Dist2=NA, Speed_ave=NA, Theta=NA, Delta_th=NA)

    N = nrow(S)
    for (i in 1:(N-1)){
      tryCatch({
      S$Dist2[i+1] = (((S$Lat[i+1]-S$Lat[i])*39946.79/360)^2+
                        ((S$Lon[i+1]-S$Lon[i])*pi*12756.32/360*cos(S$Lat[i]*pi*2/360))^2)^0.5
      S$Dist[i+1]  = acos(sin(S$Lat[i]*pi/180)*sin(S$Lat[i+1]*pi/180) + cos(S$Lat[i]*pi/180)*cos(S$Lat[i+1]*pi/180) * # Calculates the geodesic distance between two points specified
                            cos(S$Lon[i+1]*pi/180-S$Lon[i]*pi/180)) * 6371 # by radian latitude/longitude using the Spherical Law of Cosines (slc)
      S$Speed_ave[i+1] = S$Dist2[i+1] / as.numeric(difftime(strptime(S$Time[i+1], "%Y-%m-%d %H:%M:%S"), strptime(S$Time[i], "%Y-%m-%d %H:%M:%S"), units='hours'))  #km/h #"%Y/%m/%d %H:%M:%S"
      S$Theta[i+1] = acos((S$Lat[i+1]-S$Lat[i])*39946.79/360 / S$Dist2[i+1]) *360/2/pi # top/bottom right quadrant #Dist would be wrong
      if(S$Lon[i+1] < S$Lon[i]) S$Theta[i+1] = 360 - S$Theta[i+1] # top/bottom left quadrant
      if(S$Speed_ave[i+1] == 0) S$Theta[i+1] = 0 # fixing the NA problem when animal did not move
      }, error = function(e){cat(".")})
    }

    # changing in moving direction
    for (j in 2:(N-1)){
      tryCatch({
      S$Delta_th[j+1] = S$Theta[j+1] - S$Theta[j]
      if(S$Delta_th[j+1] > 180)   S$Delta_th[j+1] = S$Delta_th[j+1] -180
      if(S$Delta_th[j+1] < -180)  S$Delta_th[j+1] = S$Delta_th[j+1] +180
      }, error = function(e){cat(".")})
    }

    pop[[K]] <- S
    print(K)
  }

  POP = rbind(pop[[1]], pop[[2]])
  for (i in 3:length(list)){
    POP = rbind(POP, pop[[i]])
  }
  names(POP)[6:8] = c("Year","Month","Day")
  return(POP)
}
###################################################################################################




# Plot key areas
##############################################################################################################
#' Plot breeding and wintering areas
#'
#' @description Plot breeding and wintering areas based on kernel density maps. The starting and ending dates of the breeding peroid and wintering period are needed. The areas can be adjusted by the percentage of kernel areas.
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param D A data.frame of processed satellite tracking data.
#'
#' @return
#'
#' @examples
#'
#'  plot.breeding.wintering(D=data, ext.par=4, breed.start=100,
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
#'
plot.breeding.wintering = function(D=D, ext.par=4, breed.start=breed.start, # shape=shape,
                                   breed.end=breed.end, winter.start=winter.start, winter.end=winter.end,
                                   lat.min.b = lat.min.b, lat.max.w = lat.max.w,
                                   breed.percent = c(90, 70, 50), winter.percent = c(90, 70, 50)){

  xlim = range(D$Lon, na.rm=T); ylim = range(D$Lat, na.rm=T)
  xlim[1] = xlim[1]-ext.par; xlim[2] = xlim[2]+ext.par
  ylim[1] = ylim[1]-ext.par; ylim[2] = ylim[2]+ext.par
  library(maps)
  map('world', mar = c(4,4,0,0), xlim=xlim, ylim=ylim, fill=F)
  # plot(shape, mar = c(4,4,0,0), xlim=xlim, ylim=ylim, fill=T, col="grey90")
  LatLon = D[D$Day>breed.start & D$Day<breed.end, c('Lon', 'Lat')]
  LatLon = LatLon[LatLon$Lat > lat.min.b, ]
  LatLon = round(LatLon, 1) # thinning the points
  LatLon = unique.data.frame(LatLon)
  library(adehabitatHR)
  xy <- SpatialPoints(LatLon)
  kud <- kernelUD(xy)  # h = href is the default - ad hoc method for determining h
  tryCatch({ver <- getverticeshr(kud, percent = breed.percent[1]); plot(ver, col=adjustcolor("orange", alpha.f=0.5), add=T, border="NA")}, error = function(e){cat("ERROR: ", conditionMessage(e),"\n")})
  tryCatch({ver <- getverticeshr(kud, percent = breed.percent[2]); plot(ver, col=adjustcolor("orange", alpha.f=0.5), add=T, border="NA")}, error = function(e){cat(".")})
  tryCatch({ver <- getverticeshr(kud, percent = breed.percent[3]); plot(ver, col=adjustcolor("orange", alpha.f=0.5), add=T, border="NA")}, error = function(e){cat(".")})
  # ver <- getverticeshr(kud, percent = 50); plot(ver, col=adjustcolor("orange", alpha.f=0.5), add=T, border="NA")

  LatLon = D[D$Day>winter.start & D$Day<winter.end, c('Lon', 'Lat')]
  LatLon = LatLon[LatLon$Lat < lat.max.w, ]
  LatLon = round(LatLon, 2) # thinning the points
  LatLon = unique.data.frame(LatLon)
  xy <- SpatialPoints(LatLon)
  kud <- kernelUD(xy)  # h = href is the default - ad hoc method for determining h
  tryCatch({ver <- getverticeshr(kud, percent = winter.percent[1]); plot(ver, col=adjustcolor("blue", alpha.f=0.5), add=T, border="NA")}, error = function(e){cat(".")})
  tryCatch({ver <- getverticeshr(kud, percent = winter.percent[2]); plot(ver, col=adjustcolor("blue", alpha.f=0.5), add=T, border="NA")}, error = function(e){cat(".")})
  tryCatch({ver <- getverticeshr(kud, percent = winter.percent[3]); plot(ver, col=adjustcolor("blue", alpha.f=0.5), add=T, border="NA")}, error = function(e){cat(".")})
  # ver <- getverticeshr(kud, percent = 50); plot(ver, col=adjustcolor("blue", alpha.f=0.5), add=T, border="NA")   # }, error = function(e){cat(".")})

  list = names(table(D$ID))
  for (YEAR in min(D$Year,na.rm=T):max(D$Year,na.rm=T)){
    for (i in 1:length(list)){
      tryCatch({
        Ind = D[D$ID == list[i] & D$Year==YEAR, ]
        if (nrow(Ind) >300) {
          lines(Ind$Lon, Ind$Lat, col=rainbow(length(list))[i], lwd=0.5)
        }
      }, error = function(e){cat(".")})
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
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param D A data.frame of processed satellite tracking data.
#'
#' @return
#'
#' @examples
#'
#' Dist = dist.annual(D=D, min.day=10); Dist
#'
#' @export
#'
dist.annual = function(D=D, min.day=min.day){
  list = names(table(D$ID))
  Yr = names(table(D$Year))
  DIST = data.frame(ID=1:(length(list)*length(Yr)), Year=NA, Ind=NA, Distance=NA)
  N = 0
  for (YEAR in min(D$Year,na.rm=T):max(D$Year,na.rm=T)){
    for (i in 1:length(list)){
      N = N+1

      tryCatch({

        Ind = D[D$ID == list[i] & D$Year==YEAR, ]
        # Ind = Ind[(Ind$Day > 50 & Ind$Day <180) | (Ind$Day > 250 & Ind$Day <340), ] # migration period
        if (length(unique(Ind$Day)) > min.day) { # valid records for at least 300 days
          DIST$Year[N] = YEAR
          DIST$Ind[N] = list[i]
          DIST$Distance[N] = sum(Ind$Dist2, na.rm=T)
        }

      }, error = function(e){cat(".")})

      print(N)
    }
  }
  DIST = DIST[!is.na(DIST$Ind),]
  return(DIST)
}
##############################################################################################################



# trajectory plot
##############################################################################################################
#' Plot trajectories of all individuals in all years
#'
#' @description Plot two types of trajectories of all individuals in all years. Type "ind" shows different individuals in different color. Type "time" shows locating points in different colors determined by Julian days.
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param D A data.frame of processed satellite tracking data.
#'
#' @return
#'
#' @examples
#'
#' plot.traj(D=D, ext.par=4, transparent = 0.2, p.size = 0.5, l.size = 1, type="ind")
#' plot.traj(D=D, ext.par=4, transparent = 0.2, p.size = 0.5, l.size = 1, type="time")
#'
#' @importFrom maps map
#'
#' @export
#'
plot.traj = function(D=D, ext.par=ext.par, transparent = 0.1, p.size = 0.2, l.size = 1, type="ind"){
  xlim = range(D$Lon, na.rm=T); ylim = range(D$Lat, na.rm=T)
  xlim[1] = xlim[1]-ext.par; xlim[2] = xlim[2]+ext.par
  ylim[1] = ylim[1]-ext.par; ylim[2] = ylim[2]+ext.par
  library(maps)
  map('world', mar = c(4,4,0,0), xlim=xlim, ylim=ylim, fill=F)
  list = names(table(D$ID))

  if (type=="ind") {
    for (YEAR in min(D$Year,na.rm=T):max(D$Year,na.rm=T)){
      for (i in 1:length(list)){
        tryCatch({
          Ind = D[D$ID == list[i] & D$Year==YEAR, ]
          col=sample(1:20, 1)
          points(Ind$Lon, Ind$Lat, col=adjustcolor(col, alpha.f=transparent), cex=p.size)
          lines(Ind$Lon, Ind$Lat, col=col, lwd = l.size)
        }, error = function(e){cat(".")})
      }
    }
  }
  else if (type=="time"){
    for (YEAR in min(D$Year,na.rm=T):max(D$Year,na.rm=T)){
      for (i in 1:length(list)){
        tryCatch({
          Ind = D[D$ID == list[i] & D$Year==YEAR, ]
          col = colorRampPalette(c("red", "yellow", "green"))(max(Ind$Day))[Ind$Day]
          points(Ind$Lon, Ind$Lat, col=adjustcolor(col, alpha.f=transparent), cex=p.size)
          lines(Ind$Lon, Ind$Lat, col="grey", lwd = l.size)
        }, error = function(e){cat(".")})
      }
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
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param D A data.frame of processed satellite tracking data.
#'
#' @return
#'
#' @examples
#'
#' Dur = duration(D); head(Dur)
#'
#' @export
#'
duration = function(D=D){
  list = names(table(D$ID))
  DUR = data.frame(Individual=list, Duration=NA, No=NA)
  for (i in 1:length(list)){
    Ind = D[D$ID == list[i], ]
    DUR$Duration[i] = (max(Ind$Year) - min(Ind$Year))*365 + Ind$Day[nrow(Ind)] - Ind$Day[1]
    DUR$No[i] = nrow(Ind)
  }
  names(DUR) = c("Individual", "Duration", "No_record")
  return(DUR)
}
##############################################################################################################



# Plot tracking duration
##############################################################################################################
#' Plot tracking duration (days) for all individuals
#'
#' @description Plot tracking duration (from the first day to the last day (or the present day)) for all individuals
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param D A data.frame of processed satellite tracking data.
#'
#' @return
#'
#' @examples
#'
#' track.duration.plot(D)
#'
#' @export
#'
plot.track.duration = function(D=D){
  list = names(table(D$ID))
  Time <- strptime(D$Time, "%Y-%m-%d %H:%M:%S") #capital Year for 4 digits  # :%S
  Date_D = as.numeric(difftime(Time, min(Time), units='days'))
  D = cbind(D, Date_D)
  D = D[order(D$Date_D), ]
  plot(c(0, max(D$Date_D)),c(1, length(list)), col="white", xlab="Date", ylab="Individuals", xaxt='n')
  at = c(1, floor(nrow(D)/2), nrow(D))
  axis(side=1, at=D$Date_D[ at ], labels = as.Date(D$Date_D[ at ], origin = min(Time)))

  for (i in 1:length(list)){
    Ind = D[D$ID == list[i], ]
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
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param D A data.frame of processed satellite tracking data.
#'
#' @return
#'
#' @examples
#'
#' nest.locating(D, breed.S = 130, breed.E = 180, minimum.rec = 100)
#'
#' # Check a single individual
#' list = names(table(D$ID))
#' i=1; YEAR=2018
#' Ind = D[D$ID == list[i] & D$Year==YEAR, ]
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
nest.locating = function(D=D, breed.S = breed.S, breed.E = breed.E, minimum.rec = 100){
  list = names(table(D$ID)) # all individuals of a species
  N = length(list)*length(unique(D$Year))
  NEST = data.frame(ID=1:N, Ind = NA, Year=NA, Nest_Lat=NA, Nest_Lon=NA) # for holding all results
  N = 1
  for (YEAR in min(D$Year):max(D$Year)){
    for (i in 1:length(list)){
      Ind = D[D$ID == list[i] & D$Year==YEAR, ]
      Ind = Ind[Ind$Day > breed.S & Ind$Day < breed.E, ]
      if (nrow(Ind) >= minimum.rec) {
        LAT = round(Ind$Lat, 4); LON = round(Ind$Lon, 4)
        LATLON = as.character(LAT*LON*10^8) # numeric would cause no-match
        NEST$Ind[N] = list[i]
        NEST$Year[N] = YEAR
        frq = sort(table(LATLON), decreasing=T)
        NEST$Nest_Lat[N] = Ind$Lat[LATLON ==  names(frq [frq==max(frq)] )][1]
        NEST$Nest_Lon[N] = Ind$Lon[LATLON ==  names(frq [frq==max(frq)] )][1]
      }
      N = N+1
      print(N)
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
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param D A data.frame of processed satellite tracking data.
#'
#' @return
#'
#' @examples
#'
#' Daily.dist = dist.daily(D); head(Daily.dist)
#' plot(Daily.dist$Day, Daily.dist$Dist2, col=as.numeric(as.factor(Daily.dist$Individual)), xlab="Julian day", ylab="Flying distance (km)")
#'
#' @export
#'
dist.daily = function(D=D){
  list = names(table(D$ID)) # all individuals of a species
  DIST = list()
  N = 0
  DAY = 1:365; DAY = as.data.frame(DAY)
  for (YEAR in min(D$Year):max(D$Year)){
    for (i in 1:length(list)){
      Ind = D[D$ID == list[i] & D$Year==YEAR, ]
      if (nrow(Ind>1)) {
        Dist = aggregate(Dist2 ~ Day, sum, data=Ind)
        Dist = data.frame(Dist, Individual=list[i])
        DAYS = merge(DAY, Dist, by.x = "DAY", by.y = "Day", all.x=T)
        DAYS = data.frame(DAYS, Year=YEAR)
      }
      N = N+1
      DIST[[N]] = DAYS
    }
  }
  DD = DIST[[1]]
  for (i in 2:length(DIST)){
    DD = rbind(DD, DIST[[i]])
  }
  DD = DD[!is.na(DD$Individual),]
  names(DD)[1] = "Day"
  return(DD)
}
##############################################################################################################






# plot daily movement distance
##############################################################################################################
#' Plot daily movement distance (km) for each individual
#'
#' @description Plot the daily total movement distance for each individual in the range of Julian day (1-365)
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param D A data.frame returned from function daily.dist().
#'
#' @return
#'
#' @examples
#'
#' Daily.dist = daily.dist(D)
#' plot.daily.dist(Daily.dist)
#'
#' @export
#'
plot.daily.dist = function(D=D){
  list = names(table(D$Individual))
  plot(D$Day, D$Dist2, xlab="Julian day", ylab="Flying distance (km)",type='n')
  for (YEAR in min(D$Year):max(D$Year)){
    for (i in 1:length(list)){
      tryCatch({
        Ind = D[D$Individual == list[i] & D$Year==YEAR, ]
        Ind = Ind[order(Ind$Day), ]
        if (nrow(Ind) >1) {
          points(Ind$Day, Ind$Dist2, col=rainbow(24)[i])
          lines( Ind$Day, Ind$Dist2, col=rainbow(24)[i])
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
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param D A data.frame of processed satellite tracking data.
#'
#' @return
#'
#' @examples
#'
#' dates = mig.day(D, min.dist=100) ## minimum distance 100 km
#'
#' @export
#'
mig.day = function(D=D, min.dist = min.dist){
  Mig.day = unique(D[, c("ID", "Year")])
  Mig.day = data.frame(Mig.day, Breed_S_Day = NA, Breed_E_Day = NA, Winter_S_Day = NA, Winter_E_Day = NA)

  for (i in 1:nrow(Mig.day)){
    tryCatch({
      Ind = D[D$ID == Mig.day$ID[i] & D$Year == Mig.day$Year[i], ]
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
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param D A data.frame of processed satellite tracking data.
#' @param Mig.day A data.frame returned from function mig.day().
#'
#' @return
#'
#' @examples
#'
#' DD = mig.hour(D=D, Mig.day = dates, outlier.dist = 150, min.move = 5)
#'
#' @export
#'
mig.hour = function(D=D, Mig.day=Mig.day, outlier.dist = outlier.dist, min.move = min.move){
  Mig.hour = data.frame(Mig.day$ID, Mig.day$Year, Breed_S_Hour=NA, Breed_E_Hour=NA,
                        Winter_S_Hour=NA, Winter_E_Hour=NA)
  Mig.s = list(nrow(Mig.day))
  hours = numeric(nrow(Mig.day))

  for (i in 1: nrow(Mig.day)){
    tryCatch({
      Mig.s[[i]] = D[D$ID==Mig.day$ID[i] & D$Year==Mig.day$Year[i] & D$Day==Mig.day$Breed_S_Day[i],]
      Mig.s[[i]] = Mig.s[[i]][Mig.s[[i]]$Dist < outlier.dist,] # points with distance over 150 km was not accumulated by one hour
      hours[i] = Mig.s[[i]][Mig.s[[i]]$Dist > min.move, ]$Hour[1] # Distance > 10 km means start to migrate
    }, error = function(e){cat(".")})
  }
  Mig.hour$Breed_S_Hour = hours

  Mig.s = list(nrow(Mig.day))
  hours = numeric(nrow(Mig.day))
  for (i in 1: nrow(Mig.day)){
    tryCatch({
      Mig.s[[i]] = D[D$ID==Mig.day$ID[i] & D$Year==Mig.day$Year[i] & D$Day==Mig.day$Breed_E_Day[i],]
      Mig.s[[i]] = Mig.s[[i]][Mig.s[[i]]$Dist < outlier.dist,] # points with distance over 150 km was not accumulated by one hour
      hours[i] = max(Mig.s[[i]][Mig.s[[i]]$Dist > min.move, ]$Hour) # Distance > 10 km means start to migrate
    }, error = function(e){cat(".")})
  }
  Mig.hour$Breed_E_Hour = hours

  Mig.s = list(nrow(Mig.day))
  hours = numeric(nrow(Mig.day))
  for (i in 1: nrow(Mig.day)){
    tryCatch({
      Mig.s[[i]] = D[D$ID==Mig.day$ID[i] & D$Year==Mig.day$Year[i] & D$Day==Mig.day$Winter_S_Day[i],]
      Mig.s[[i]] = Mig.s[[i]][Mig.s[[i]]$Dist < outlier.dist,] # points with distance over 150 km was not accumulated by one hour
      hours[i] = Mig.s[[i]][Mig.s[[i]]$Dist > min.move, ]$Hour[1] # Distance > 10 km means start to migrate
    }, error = function(e){cat(".")})
  }
  Mig.hour$Winter_S_Hour = hours

  Mig.s = list(nrow(Mig.day))
  hours = numeric(nrow(Mig.day))
  for (i in 1: nrow(Mig.day)){
    tryCatch({
      Mig.s[[i]] = D[D$ID==Mig.day$ID[i] & D$Year==Mig.day$Year[i] & D$Day==Mig.day$Winter_E_Day[i],]
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
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param ind Locating records of an individual in a year, from the processed satellite tracking data.
#'
#' @return
#'
#' @examples
#'
#' ind = D[D$ID=="蓑羽鹤06-683·BFU076" & D$Year == 2018,]
#' plot.direction(ind)
#'
#' @importFrom plotrix polar.plot
#'
#' @export
#'
plot.direction = function(ind = ind){
  library(plotrix)
  ind = ind[!ind$Theta>360, ]
  oldpar <- polar.plot(log(ind$Speed_ave+1), ind$Theta, main="Fly direction",
                       radial.lim = c(0, 5), start=90, clockwise=TRUE, lwd=3,
                       line.col=colorRampPalette(c("green","yellow","red"))(365)[ind$Day])
  par(oldpar) # reset everything
}
##############################################################################################################



