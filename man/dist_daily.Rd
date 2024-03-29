% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/migrationR functions.R
\name{dist_daily}
\alias{dist_daily}
\title{Calculate daily movement distance (km) for each individual}
\usage{
dist_daily(trackdata)
}
\arguments{
\item{trackdata}{A data.frame of processed satellite tracking data.}
}
\description{
Sum the daily total movement distance for each individual
}
\examples{

Daily.dist = dist_daily(trackdata); head(Daily.dist)
plot(Daily.dist$Day, Daily.dist$Dist2, col=as.numeric(as.factor(Daily.dist$Individual)), xlab="Julian day", ylab="Flying distance (km)")

}
\author{
Huidong Tian (tienhuitung@gmail.com); Xinhai Li (xinhai_li_edu@126.com)
}
