## **************************************************************************
##
##    (c) 2018-2022 Guillaume Guénard
##        Department de sciences biologiques,
##        Université de Montréal
##        Montreal, QC, Canada
##
##    **Legendre and Gallagher transformation function**
##
##    This file is part of codep
##
##    codep is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    codep is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with codep. If not, see <https://www.gnu.org/licenses/>.
##
##    R source code file
##
## **************************************************************************
##
#' Transformation for Species Abundance Data
#' 
#' Calculates the transformed species abundances following Legendre and
#' Gallagher.
#' 
#' @param x A species abundance matrix (rows: sites, columns: species).
#' @param method The transformation method, one of chord (the default), chisq,
#' profile, Hellinger (see details).
#' 
#' @return A matrix of the transformed species abundances.
#' 
#' @details Multivariate least squares methods for ordination methods, such as
#' the principal component analysis, or modelling methods, such as the
#' multiscale codependence analysis (\code{\link{MCA}}), the canonical
#' redundancy analysis (RDA).
#' 
#' @author \packageAuthor{codep}
#' Maintainer: \packageMaintainer{codep}
#' 
#' @references
#' Legendre, P. & Gallagher E. D. 2001. Ecologically meaningful transformations
#' for ordination of species data. Oecologia 129: 271-280
#' doi: 10.1007/s004420100716
#' 
#' @importFrom stats dist
#' @importFrom graphics mtext
#' 
#' @examples
#' 
#' data(Doubs)
#' 
#' ## Removing any species that have not been not observed:
#' Doubs.fish -> x
#' x[rowSums(x)!=0,] -> x
#' 
#' ## Transforming the abundances
#' LGTransforms(x,"chord") -> chord
#' LGTransforms(x,"chisq") -> chisq
#' LGTransforms(x,"profile") -> sp_pr
#' LGTransforms(x,"Hellinger") -> Helli
#' 
#' dist(chord)
#' dist(chisq)
#' dist(sp_pr)
#' dist(Helli)
#' 
#' ## Legendre & Gallagher synthetic examples:
#' 
#' data(LGDat)
#' 
#' ## Diastemograms:
#' 
#' as.matrix(dist(LGDat[,1L])) -> geo
#' geo[upper.tri(geo)] -> geo
#' 
#' ## Raw Euclidean distances
#' par(mfrow=c(1,1), mar=c(5,5,4,2))
#' 
#' as.matrix(dist(LGDat[,-1L])) -> eco
#' eco[upper.tri(eco)] -> eco
#' 
#' plot(eco~geo, data=data.frame(geo=geo, eco=eco),
#'      xaxp=c(1,18,17), las=1, ylim=c(0,max(eco)),
#'      xlab="True geographic distance",
#'      ylab="Euclidean distance")
#' 
#' ## Euclidean distances on the transformed abundances:
#' par(mfrow=c(2,2), mar=c(5,5,4,2))
#' 
#' LGTransforms(LGDat[,-1L],"chord") -> chord
#' as.matrix(dist(chord)) -> eco
#' eco[upper.tri(eco)] -> eco
#' plot(eco~geo,data=data.frame(geo=geo,eco=eco),
#'      xaxp=c(1,18,17),las=1,xlab="",ylab="",
#'      main="Chord distance",ylim=c(0,max(eco)))
#' 
#' LGTransforms(LGDat[,-1L],"chisq") -> chisq
#' as.matrix(dist(chisq)) -> eco
#' eco[upper.tri(eco)] -> eco
#' plot(eco~geo,data=data.frame(geo=geo,eco=eco),
#'      xaxp=c(1,18,17),las=1,xlab="",ylab="",
#'      main="Chi-square distance",ylim=c(0,max(eco)))
#' 
#' LGTransforms(LGDat[,-1L],"profile") -> sp_pr
#' as.matrix(dist(sp_pr)) -> eco
#' eco[upper.tri(eco)] -> eco
#' plot(eco~geo,data=data.frame(geo=geo,eco=eco),
#'      xaxp=c(1,18,17),las=1,xlab="",ylab="",
#'      main="Dist. between profiles",ylim=c(0,max(eco)))
#' 
#' LGTransforms(LGDat[,-1L],"Hellinger") -> Helli
#' as.matrix(dist(Helli)) -> eco
#' eco[upper.tri(eco)] -> eco
#' plot(eco~geo,data=data.frame(geo=geo,eco=eco),
#'      xaxp=c(1,18,17),las=1,xlab="",ylab="",
#'      main="Hellinger distance",ylim=c(0,max(eco)))
#' 
#' mtext(text="True geographic distance", side=1, line=-1.5, outer=TRUE)
#' mtext(text="Ecological distance", side=2, line=-1.5, outer=TRUE)
#' 
#' @export
LGTransforms <- function(x, method = c("chord","chisq","profile","Hellinger")) {
  method <- match.arg(method)
  method <- match(method, c("chord","chisq","profile","Hellinger"))
  if(!is.matrix(x))
    x <- as.matrix(x)
  storage.mode(x) <- "double"
  return(
    .C(
      "LGTr_C",
      x,
      nrow(x),
      ncol(x),
      method
    )[[1L]]
  )
}
#' 
