% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_comp_data.R
\name{plot_comp_data}
\alias{plot_comp_data}
\title{Plot simulated compositional data}
\usage{
plot_comp_data(sim)
}
\arguments{
\item{sim}{A list containing the simulated data. The \code{sim$data} element must be a \code{SpatRaster}
or \code{PackedSpatRaster} object. The \code{sim$d} element specifies the number of compositional
components to plot.}
}
\value{
A plot of the compositional data layers. The function does not return a value;
it generates a plot as a side effect.
}
\description{
This function generates a plot of simulated compositional data stored in a \code{SpatRaster} object.
It extracts the relevant layers from the raster object and visualizes them as a series of plots.
}
\details{
If \code{sim$data} is a \code{PackedSpatRaster} object, it is first unpacked using
\code{terra::unwrap()}. The function ensures that \code{sim$data} is a \code{SpatRaster} object
before proceeding. The layers to be plotted are identified by their names, which
are expected to follow the pattern \code{"p_sim_1"}, \code{"p_sim_2"}, ..., \code{"p_sim_d"},
where \code{d} is the number of components.
}
