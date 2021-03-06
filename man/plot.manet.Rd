% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.manet.R
\name{plot.manet}
\alias{plot.manet}
\title{Plotting the output from the multiple allocation clustering.}
\usage{
\method{plot}{manet}(x, seed = 1, layout = layout_nicely, ...)
}
\arguments{
\item{x}{A manet object.}

\item{seed}{Random seed. Default is 1.}

\item{layout}{Layout of the network from the igraph package. Default is layout_nicely.}

\item{...}{Additional inputs to the igraph function.}
}
\value{
An actor-event network with events as round circles and actors as squared circles with the different colours corresponding to the identified communities.
}
\description{
This function plots the output of the manet function.
}
\examples{
data(deepsouth)
ds<-manet(deepsouth,K=2,maxT=100)
plot(ds)
}
