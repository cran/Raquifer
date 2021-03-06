% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Raquifer_package.R
\name{aquifer_predict.ykh_rad_bottom}
\alias{aquifer_predict.ykh_rad_bottom}
\title{S3 method for class 'aquifer_predict'}
\usage{
\method{aquifer_predict}{ykh_rad_bottom}(aquifer_lst, time_lst)
}
\arguments{
\item{aquifer_lst}{a list object of class 'aquifer'}

\item{time_lst}{a list object of class 'time'}
}
\value{
a data frame of cumulative water influx estimates using the Yildiz-Khosravi un-steady state radial flow model, bottom-water-drive
}
\description{
Return a data frame of estimated cumulative water influx for the Yildiz-Khosravi un-steady state radial flow model, bottom-water-drive
}
