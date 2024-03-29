% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eigenmap.R
\name{eigenmap}
\alias{eigenmap}
\alias{eigenmap.score}
\title{Spatial Eigenvector Maps}
\usage{
eigenmap(
  x,
  alt.coord = NA,
  weighting = wf.sqrd,
  boundaries,
  wpar,
  tol = .Machine$double.eps^0.5
)

eigenmap.score(emap, target)
}
\arguments{
\item{x}{A set of coordinates defined in one (numeric vector) or many (a
coordinate x dimension matrix) dimensions or, alternatively, a distance
matrix provided by \code{\link{dist}}.}

\item{alt.coord}{Coordinates to be used when a distance matrix is
provided as x. Used for plotting purposes.}

\item{weighting}{The function to obtain the edge weighting matrix (see
details).}

\item{boundaries}{When required by argument \code{weighting}, a two-element
numeric vector containing the lower and upper threshold values used to obtain
the connectivity matrix (see \link{weighting-functions}).}

\item{wpar}{Shape parameter for argument \code{weignting} (optional).}

\item{tol}{The smallest absolute eigenvalue for a spatial eigenfunctions to
be considered as a suitable predictor. Default:
\code{.Machine$double.eps^0.5} (a machine-dependent value).}

\item{emap}{An \link{eigenmap-class} object.}

\item{target}{A (generally rectangular) distance matrix between a set of
target locations for which spatially-explicit predictions are being made
(rows), and the reference locations given to function \code{eigenmap}
(columns). See example 2.}
}
\description{
Function to calculate spatial eigenvector maps of a set of locations
in a space with an arbitrary number of dimension.
}
\details{
When function \code{eigenmap} is given coordinates as its argument
\code{x}, they are treated as Cartesian coordinates and the distances between
them are assumed to be Euclidean. Otherwise (e.g., when geodesic distances
are used), distances have to be provided as the argument \code{x} and
plotting coordinates have to be supplied as argument \code{alt.coord}.

The weighting function (see \link{weighting-functions}) must have the
distances as its first argument, optionally an argument named
\code{boundaries} giving the boundaries within which locations are regarded
as neighbours and/or an argument \code{wpar} containing any other weighting
function parameters.

Default values for argument \code{boundaries} are 0 for the minimum value and
\code{NA} for the maximum. For weighting functions with an argument
\code{bounraries}, The upper value \code{NA} indicates the function to take
the minimum value that allow every locations to form a single cluster
following single linkage clustering as a maximum value (obtained internally
from a call to \code{\link{hclust}}.
}
\section{Functions}{
\itemize{
\item \code{eigenmap()}: Main function for generating an eigenmap-class object from Cartesian
coordinates or pairwise distances.

\item \code{eigenmap.score()}: Generate scores for arbitrary locations within the scope of an existing
eigenvector map.

}}
\examples{
### Example 1: A linear transect.
data(salmon)

## A warning is issued when no boundaries are provided for a function that
## requires them.
## Example:
map <- eigenmap(x = salmon[,"Position"], weighting = wf.binary)
map
## plot(map)

## In the following examples, boundaries are provided; they are needed by the
## functions.
map <- eigenmap(x = salmon[,"Position"], weighting = wf.binary,
                boundaries = c(0,20))
map
## plot(map)

map <- eigenmap(x = salmon[,"Position"], weighting = wf.Drayf1,
                boundaries = c(0,20))
map
## plot(map)

map <- eigenmap(x = salmon[,"Position"], weighting = wf.Drayf2,
                boundaries = c(0,20))
map
## plot(map)

map <- eigenmap(x = salmon[,"Position"], weighting = wf.Drayf3,
                boundaries = c(0,20), wpar = 2)
map
## plot(map)

map <- eigenmap(x = salmon[,"Position"], weighting = wf.PCNM,
                boundaries = c(0,20))
map
## plot(map)

map <- eigenmap(x = salmon[,"Position"], weighting = wf.sqrd)
map
## plot(map)

map <- eigenmap(x = salmon[,"Position"], weighting = wf.RBF, wpar = 0.001)
map
## plot(map)

### Example 2: Using predictor scores

smpl <- c(4,7,10,14,34,56,61,64)  # A sample to be discarded
map <- eigenmap(x = salmon[-smpl,"Position"], weighting = wf.sqrd)
scr <- eigenmap.score(
         map, target = as.matrix(dist(salmon[,"Position"]))[,-smpl]
       )
## Scores of sampling points are the eigenvectors
scr[smpl,]

wh <- 5L   # You can try with other vectors.
plot(map$U[,wh] ~ salmon[-smpl,"Position"], ylab = expression(U[5]),
     xlab = "Position along transect")
points(y = scr[smpl,wh], x = salmon[smpl,"Position"], pch = 21L,
       bg = "black")

map <- eigenmap(x = salmon[-smpl,"Position"], weighting = wf.binary,
                boundaries = c(0,20))
scr <- eigenmap.score(
         map, target = as.matrix(dist(salmon[,"Position"]))[smpl,-smpl])

## Plot the 8 prediction sites along particular eigenvectors, here
## eigenvector #1:
wh <- 1L   # One could try the other vectors.
plot(map$U[,wh] ~ salmon[-smpl,"Position"], ylab = expression(U[1L]),
     xlab = "Position along transect (m)")
points(y = scr[,wh], x = salmon[smpl,"Position"], pch=21L, bg = "black")

map <- eigenmap(x = salmon[-smpl,"Position"], weighting = wf.PCNM,
                boundaries = c(0,100))
scr <- eigenmap.score(
         map, target = as.matrix(dist(salmon[,"Position"]))[smpl,-smpl]
       )

wh <- 1L   # You can try with other vectors.
plot(map$U[,wh] ~ salmon[-smpl,"Position"], ylab = expression(U[1]),
     xlab = "Position along transect (m)")
points(y = scr[,wh], x = salmon[smpl,"Position"], pch = 21L, bg = "black")

### Example 3: A unevenly sampled surface.

data(mite)

## Example using the principal coordinates of the square root of the
## (Euclidean) distances:
map <- eigenmap(x = as.matrix(mite.geo), weighting = wf.sqrd)
map
## plot(map)

## Example using the radial basis functions (RBF):
map <- eigenmap(x = as.matrix(mite.geo), weighting = wf.RBF)
map
## plot(map)

}
\references{
Borcard, D. and Legendre, P. 2002. All-scale spatial analysis of ecological
data by means of principal coordinates of neighbour matrices. Ecol. Model.
153: 51-68

Dray, S.; Legendre, P. and Peres-Neto, P. 2006. Spatial modelling: a
comprehensive framework for principal coordinate analysis of neighbor
matrices (PCNM). Ecol. Modelling 196: 483-493

Legendre, P. and Legendre, L. 2012. Numerical Ecology, 3rd English edition.
Elsevier Science B.V., Amsterdam, The Netherlands.
}
\author{
\packageAuthor{codep}
Maintainer: \packageMaintainer{codep}
}
