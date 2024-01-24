# migrationR: a package for processing satellite telemetry data and analyzing migratory patterns

## Introduction

The R package migrationR is designed to efficiently process satellite telemetry data, thereby enabling in-depth analysis of animal migration patterns. It facilitates the extraction, organization, and interpretation of movement data collected from satellite tags attached to migratory species. With migrationR, researchers can uncover valuable insights into migratory routes, timing, and behavior, thereby enhancing our understanding of ecological and environmental influences on animal migrations.

```{r}
# Install migrationR
library(devtools)
install_github("Xinhai-Li/migrationR", force = TRUE)
library(migrationR)
data(movebankdata) # 96976 occurrences of 27 Demoiselle Cranes from 2018/1/1 to 2019/9/9 at 1h interval.
library(knitr)
ncol(movebankdata)
knitr::kable(head(movebankdata[, 1:5]), "pipe") # Table 1 upper panel
knitr::kable(head(movebankdata[, 6:9]), "pipe") # Table 1 lower panel
```

**Table 1. The first six rows of the satellite tracking data for Demoiselle Cranes using the Movebank data format.**

|event.id  |timestamp     | location.long| location.lat|sensor.type     |
|:---------|:-------------|-------------:|------------:|:---------------|
|38DCBA5A2 |2018/1/1 0:00 |      69.40862|     21.85748|GPS-transmitter |
|38DCBA5A2 |2018/1/1 1:00 |      69.40854|     21.85745|GPS-transmitter |
|38DCBA5A2 |2018/1/1 2:00 |      69.40860|     21.85746|GPS-transmitter |
|38DCBA5A2 |2018/1/1 3:00 |      69.40868|     21.85743|GPS-transmitter |
|38DCBA5A2 |2018/1/1 4:00 |      69.40872|     21.85745|GPS-transmitter |
|38DCBA5A2 |2018/1/1 5:00 |      69.40864|     21.85742|GPS-transmitter |

|individual.taxon.canonical.name |taxon.detail       |individual.local.identifier |study.name                 |
|:-------------------------------|:------------------|:---------------------------|:--------------------------|
|Gruidae                         |Anthropoides virgo |hooded06_683_BFU076         |Guo Yuming's field surveys |
|Gruidae                         |Anthropoides virgo |hooded06_683_BFU076         |Guo Yuming's field surveys |
|Gruidae                         |Anthropoides virgo |hooded06_683_BFU076         |Guo Yuming's field surveys |
|Gruidae                         |Anthropoides virgo |hooded06_683_BFU076         |Guo Yuming's field surveys |
|Gruidae                         |Anthropoides virgo |hooded06_683_BFU076         |Guo Yuming's field surveys |
|Gruidae                         |Anthropoides virgo |hooded06_683_BFU076         |Guo Yuming's field surveys |


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
knitr::kable(trackdata[1:3,], "pipe") # Table 2
knitr::kable(table(trackdata$ID), "pipe") # Table 3
```

**Table 2. The first three rows of satellite tracking data for Demoiselle Cranes formatted according to the migrationR data structure and augmented with additional variables.**

|ID                  |Time                |      Lon|      Lat| Year| Month| Day| Hour| Day_fine|      Dist| Time_interval|     Speed| Direction|     North|  Redirect|
|:-------------------|:-------------------|--------:|--------:|----:|-----:|---:|----:|--------:|---------:|-------------:|---------:|---------:|---------:|---------:|
|hooded06_683_BFU076 |2018-01-01 00:00:00 | 69.40862| 21.85748| 2018|     1|   1|    0| 1.000000|        NA|            NA|        NA|        NA|        NA|        NA|
|hooded06_683_BFU076 |2018-01-01 01:00:00 | 69.40854| 21.85745| 2018|     1|   1|    1| 1.041667| 0.0088982|             1| 0.0088982| 248.03073| 0.3779485|        NA|
|hooded06_683_BFU076 |2018-01-01 02:00:00 | 69.40860| 21.85746| 2018|     1|   1|    2| 1.083333| 0.0062866|             1| 0.0062866|  79.83361| 0.5564799| -168.1971|

**Table 3. Individual IDs of Demoiselle Cranes and the corresponding number of occurrences.**

|Var1                         |  Freq|
|:----------------------------|-----:|
|hooded06_683_BFU076          | 16029|
|hooded07_687_BFU077          |  3502|
|hooded21_adultHV4BFU069_1    | 12988|
|hooded29_140_BFU260_20180809 |  1343|
|hooded30_138_BFU262_20180811 |  1246|
|hooded31_133_BFU263_20180811 |  1316|
|hooded32_134_BFU264_20180811 |  3879|
|hooded33_135_BFU265_20180811 |  5297|
|hooded34_136_BFU266_20180812 | 12509|
|hooded35_137_BFU267_20180816 |  1548|
|hooded38_NJGF006_20190717    |  2031|
|hooded39_NJGF007_20190717    |  2073|
|hooded40_NJGF008_20190719    |  1658|
|hooded45_NJGF015_20190723    |  1896|
|hooded46_NJGF017_20190723    |  1691|
|hooded47_NJGF071_20190724    |  2715|
|hooded48_NJGF072_20190726    |  2466|
|hooded49_NJGF073_20190727    |  2039|
|hooded50_NJGF074_20190727    |  2038|
|hooded51_NJGF075_20190727    |  2960|
|hooded52_NJGF076_20190801    |  3142|
|hooded53_NJGF077_20190802    |  1684|
|hooded54_NJGF078_20190802    |  2016|
|hooded55_BFU282_20190802     |  1738|
|hooded56_BFU285_20190803     |  3516|
|hooded57_BFU286_20190803     |  3437|

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
knitr::kable(Dist, "pipe") # Table 4
```

**Table 4. Annual flying distance (km) of Demoiselle Crane individuals and number of days with valid records in the year.**

|   |ID                           | Year| Distance| numDay|
|:--|:----------------------------|----:|--------:|------:|
|1  |hooded06_683_BFU076          | 2018| 20574.94|    365|
|2  |hooded06_683_BFU076          | 2019| 17721.67|    321|
|3  |hooded07_687_BFU077          | 2018| 15871.27|    280|
|4  |hooded21_adultHV4BFU069_1    | 2018| 17209.10|    360|
|5  |hooded21_adultHV4BFU069_1    | 2019| 19856.54|    339|
|14 |hooded34_136_BFU266_20180812 | 2019| 20563.93|    365|

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

**Figure 8. The flying directions of one Demoiselle Crane in a year. The color gradient, transitioning from green to red, represents the progression of Julian day, ranging from 1 to 365.**

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
knitr::kable(head(timing), "pipe") # Table 5
Winter_End_Day = timing$Winter_End_Day
Winter_End_Day = Winter_End_Day[!is.na(Winter_End_Day)]
Winter_End_Day = Winter_End_Day[Winter_End_Day>260]
hist(Winter_End_Day, main="", nclass=30, xlab="Julian Day") # Figure 11
```

**Table 5. Estimated date (day) and time (hour) of starting and ending of migrations for Demoiselle Crane individuals.**

|      | Year|ID                           | Breed_Start_Day| Breed_Start_Hour| Breed_End_Day| Breed_End_Hour| Winter_Start_Day| Winter_Start_Hour| Winter_End_Day| Winter_End_Hour|
|:-----|----:|:----------------------------|---------------:|----------------:|-------------:|--------------:|----------------:|-----------------:|--------------:|---------------:|
|1     | 2018|hooded06_683_BFU076          |              80|               13|           126|             20|              240|                11|            310|              18|
|33467 | 2019|hooded06_683_BFU076          |              84|               15|           121|             15|              242|                11|            304|              21|
|8705  | 2018|hooded07_687_BFU077          |              76|               14|           128|             22|              231|                12|            281|              18|
|12426 | 2018|hooded21_adultHV4BFU069_1    |              86|               14|           125|             22|              234|                 6|            309|              18|
|40792 | 2019|hooded21_adultHV4BFU069_1    |              88|               14|           134|             11|              230|                12|            323|              19|
|17778 | 2018|hooded29_140_BFU260_20180809 |              NA|               NA|            NA|             NA|              268|                13|            278|               0|

![Figure 11](https://github.com/Xinhai-Li/migrationR_data/blob/main/Rplot11.png)

**Figure 11. The starting and ending dates of wintering migration of the Demoiselle Cranes.**

### HOSDM()

HOSDM() facilitates the implementation of Hetero-occurrence Species Distribution Models (HOSDMs). These are innovative individual-based species distribution models that specialize in discerning between different types of species occurrences with the ultimate goal of boosting model accuracy. Notably, HOSDMs have been specifically designed for optimal utilization of time-series telemetry data derived from satellites.

To implement Hetero-occurrence Species Distribution Models (HOSDMs), acquiring environmental variable data is essential, and such information can be sourced from the repository designated as "migrationR_data" at https://github.com/Xinhai-Li/migrationR_data.

```{r}
library(raster)
BioClim <- brick('Env_Poyang.grd')
data(HOSDMdata)
results = HOSDM(trackdata = HOSDMdata, prediction=F, buffer=0.1, absence=30, Envlayer=BioClim)
results[[1]]

# Predicted habitat suitability
plot(results[[2]][[1]], main=paste("Roosting habitat for individual", ind[1], sep=" " ))
plot(results[[3]][[1]], main=paste("Foraging habitat for individual", ind[1], sep=" " ))
```
