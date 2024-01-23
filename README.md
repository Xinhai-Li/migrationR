# migrationR: a package for processing satellite telemetry data and analyzing migratory patterns

## Introduction

The R package migrationR is designed to efficiently process satellite telemetry data, thereby enabling in-depth analysis of animal migration patterns. It facilitates the extraction, organization, and interpretation of movement data collected from satellite tags attached to migratory species. With migrationR, researchers can uncover valuable insights into migratory routes, timing, and behavior, thereby enhancing our understanding of ecological and environmental influences on animal migrations.

```{r}
# Install migrationR
library(devtools)
install_github("Xinhai-Li/migrationR", force = TRUE)
library(migrationR)
data(movebankdata)
head(movebankdata)
```
 
  event.id  timestamp location.long location.lat sensor.type individual.taxon.canonical.name taxon.detail individual.local.identifier study.name
1	38DCBA5A2	2018/1/1 0:00	69.40862	21.85748	GPS-transmitter	Gruidae	Anthropoides virgo	hooded06_683_BFU076	Guo Yuming's field surveys
2	38DCBA5A2	2018/1/1 1:00	69.40854	21.85745	GPS-transmitter	Gruidae	Anthropoides virgo	hooded06_683_BFU076	Guo Yuming's field surveys
3	38DCBA5A2	2018/1/1 2:00	69.40860	21.85746	GPS-transmitter	Gruidae	Anthropoides virgo	hooded06_683_BFU076	Guo Yuming's field surveys
4	38DCBA5A2	2018/1/1 3:00	69.40868	21.85743	GPS-transmitter	Gruidae	Anthropoides virgo	hooded06_683_BFU076	Guo Yuming's field surveys
5	38DCBA5A2	2018/1/1 4:00	69.40872	21.85745	GPS-transmitter	Gruidae	Anthropoides virgo	hooded06_683_BFU076	Guo Yuming's field surveys
6	38DCBA5A2	2018/1/1 5:00	69.40864	21.85742	GPS-transmitter	Gruidae	Anthropoides virgo	hooded06_683_BFU076	Guo Yuming's field surveys

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

To implement Hetero-occurrence Species Distribution Models (HOSDMs), acquiring environmental variable data is essential, and such information can be sourced from the repository designated as "migrationR_data" at https://github.com/Xinhai-Li/migrationR_data.
