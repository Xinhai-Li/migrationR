# migrationR: a package for processing satellite telemetry data and analyzing migratory patterns

## Introduction

The R package migrationR is designed to efficiently process satellite telemetry data, thereby enabling in-depth analysis of animal migration patterns. It facilitates the extraction, organization, and interpretation of movement data collected from satellite tags attached to migratory species. With migrationR, researchers can uncover valuable insights into migratory routes, timing, and behavior, thereby enhancing our understanding of ecological and environmental influences on animal migrations.

```{r}
# Install migrationR
library(devtools)
install_github("Xinhai-Li/migrationR", force = TRUE)
library(migrationR)
data(movebankdata) # 96976 occurrences of 27 Demoiselle Cranes from 2018/1/1 to 2019/9/9 at 1h interval.
head(movebankdata)
```
 
## Functions

The current version has 14 functions. 

### as_trackdata()

The inaugural function, as_trackdata(), serves to import Movebank data and enriches it by appending a set of pertinent variables: Year, Month, Day, Hour, Day_fine, Dist, Speed, Direction, North, and Redirect. Of these, Day_fine is a continuous variable that represents Julian Day with a precision down to the second. The Dist variable denotes Euclidean distance between sequential data points, while Speed calculates the minimum velocity between two adjacent locations. North is a continuous variable assigned 1 for movements in a northerly direction and 0 for those heading southward. Direction signifies the moving direction, and Redirect tracks changes in the animal's directional movement.

```{r}
# Change variable names and add parameters
library(chron)
library(lubridate)
library(argosfilter)
trackdata = as_trackdata(data = movebankdata, min_time_interval = 6)
trackdata[1:3,]
table(trackdata$ID)
```

### plot_breeding_wintering()

Plot breeding area, wintering area and migration routes of a population.

```{r}
plot_breeding_wintering(trackdata = trackdata, ext.par=3, breed.start=100,
                        breed.end=210, winter.start=1, winter.end=60,
                        lat.min.b = -90, lat.max.w = 90,
                        breed.percent = c(80, 60, 40), winter.percent = c(99, 95, 90))
```

![Figure 1](https://github.com/Xinhai-Li/migrationR_data/blob/main/Rplot01.png)

**Figure 1. Breeding area (yellow), wintering area (blue) and migration routes of the Demoiselle Cranes**


```{r}
# Use latitude to constrain the range
plot_breeding_wintering(trackdata = trackdata, ext.par=3, breed.start=100,
                        breed.end=210, winter.start=1, winter.end=60,
                        lat.min.b = 45, lat.max.w = 25,
                        breed.percent = c(80, 60, 40), winter.percent = c(99, 95, 90))
```

![Figure 2](https://github.com/Xinhai-Li/migrationR_data/blob/main/Rplot02.png)

**Figure 2. Breeding area (yellow), wintering area (blue) and migration routes of the Demoiselle Cranes with minimum latitude 45 for breeding area and maximum latitude 25 for wintering area.**

### dist_annual()

Annual movement distance (km) per individual per year.

```{r}
# Some individuals were not tracked all the year. Argument n is the minimum number of days in a year
# for calculating the distance
Dist = dist_annual(trackdata=trackdata, n = 200) 
Dist
```

### plot_traj()

Plot movement trajectories of all individuals.

```{r}
# Use color to distinguish individuals
plot_traj(trackdata = trackdata, type = "individual")
```

![Figure 3](https://github.com/Xinhai-Li/migrationR_data/blob/main/Rplot03.png)

**Figure 3. Flying trajectories of the Demoiselle Cranes. Different colors indicate different individuals.**

```{r}
# Use color to represent time in a year
plot_traj(trackdata = trackdata, type = "chronic")
```

![Figure 4](https://github.com/Xinhai-Li/migrationR_data/blob/main/Rplot04.png)

**Figure 4. Flying trajectories of the Demoiselle Cranes. Different colors indicate different time in a year.**

### plot_track_duration()

Plot tracking duration of all individuals.

```{r}
plot_track_duration(trackdata, cex.lab=0.9, cex.axis=0.8)
```

![Figure 5](https://github.com/Xinhai-Li/migrationR_data/blob/main/Rplot05.png)

**Figure 5. Tracking duration for all the individuals of the Demoiselle Cranes.**


### nest_locating()

Estimate nest sites based on the most used location in the breeding season.

```{r}
# The starting time (breed.S) and ending time (breed.E) of the breeding season should be adjusted for other species
nest_locating(trackdata, breed.S = 130, breed.E = 180, minimum.rec = 100)

# Check the nest location of a single individual
ind = trackdata[trackdata$ID==trackdata$ID[1] & trackdata$Year == trackdata$Year[1],]
breed.S = 140; breed.E = 180
ind = ind[ind$Day > breed.S & ind$Day < breed.E, ]
plot(ind$Lon, ind$Lat)
lines(ind$Lon, ind$Lat, col="grey", lwd=.5)
LAT = round(ind$Lat, 4); LON = round(ind$Lon, 4)
LATLON = as.character(LAT*LON*10^8) # numeric would cause no-match
frq = sort(table(LATLON), decreasing=T)
LAT = ind$Lat[LATLON ==  names(frq [frq==max(frq)] )][1]
LON = ind$Lon[LATLON ==  names(frq [frq==max(frq)] )][1]
points(LON, LAT, col=2, pch=16)
```

![Figure 6](https://github.com/Xinhai-Li/migrationR_data/blob/main/Rplot06.png)

**Figure 6. Movements trajectories and estimated nest site of one Demoiselle Crane.**

### dist_daily()

Calculate daily movement distance (km).

```{r}
Daily.dist = dist_daily(trackdata);
head(Daily.dist)
# plot(Daily.dist$Day, Daily.dist$Dist2, col=as.numeric(as.factor(Daily.dist$Individual)), 
#     xlab="Julian day", ylab="Flying distance (km)")
```

### plot_daily_dist()

Plot daily movement distance across a year.

```{r}
Daily.dist = dist_daily(trackdata)
plot_daily_dist(Daily.dist)
```

![Figure 7](https://github.com/Xinhai-Li/migrationR_data/blob/main/Rplot07.png)

**Figure 7. Daily movement distance of Demoiselle Cranes across a year. Different colors indicate different individuals**

### plot_direction()

Plot movement directions of an individual in a year.

```{r}
ind = trackdata[trackdata$ID==trackdata$ID[1] & trackdata$Year == trackdata$Year[1],]
plot_direction(ind)
```

![Figure 8](https://github.com/Xinhai-Li/migrationR_data/blob/main/Rplot08.png)

**Figure 8. The flying directions of one Demoiselle Crane in a year.**

### plot_traj_segments()

Plot time series segments of movement trajectories of an individual.

```{r}
par(mar=c(4,4,4,2))
# The colors of points from red to green indicate the locating time of the points is from old to new.
plot_traj_segments(ind=ind, seg=4, label=F)
```

![Figure 9](https://github.com/Xinhai-Li/migrationR_data/blob/main/Rplot09.png)

**Figure 9. The trajectories of a single Demoiselle Crane across four distinct periods.**

```{r}
par(mar=c(4,4,4,2))
plot_traj_segments(ind=ind, seg=6, label=T)
```

![Figure 10](https://github.com/Xinhai-Li/migrationR_data/blob/main/Rplot10.png)

**Figure 10. The trajectories of a single Demoiselle Crane across six distinct periods.**

### mig_timing()

Estimate the date (day) and time (hour) of starting and ending of migrations.

```{r}
timing = mig_timing(trackdata=trackdata, dist_min_day = 100, dist_min_hour = 10, dist_outlier = 150)
timing[1:7,]
Winter_End_Day = timing$Winter_End_Day
Winter_End_Day = Winter_End_Day[!is.na(Winter_End_Day)]
Winter_End_Day = Winter_End_Day[Winter_End_Day>260]
hist(Winter_End_Day, main="", nclass=30, xlab="Julian Day")
```

![Figure 11](https://github.com/Xinhai-Li/migrationR_data/blob/main/Rplot11.png)

**Figure 11. The starting and ending dates of wintering migration of the Demoiselle Cranes.**

To implement Hetero-occurrence Species Distribution Models (HOSDMs), acquiring environmental variable data is essential, and such information can be sourced from the repository designated as "migrationR_data" at https://github.com/Xinhai-Li/migrationR_data.
