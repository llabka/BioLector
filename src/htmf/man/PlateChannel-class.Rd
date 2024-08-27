\name{PlateChannel-class}
\Rdversion{1.1}
\docType{class}
\alias{PlateChannel-class}
\alias{initialize,PlateChannel-method}
\alias{measure<-,PlateChannel,matrix-method}
\alias{measure,PlateChannel-method}
\alias{time<-,PlateChannel,matrix-method}
\alias{time,PlateChannel-method}

\title{Class "PlateChannel"}
\description{Channel measured on a plate}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PlateChannel",
  name, measure, time)}, or use the constructor \code{PlateChannel}.

}
\section{Slots}{
	 \describe{
    \item{\code{measure}:}{Object of class \code{"matrix"} with the measures }
    \item{\code{time}:}{Object of class \code{"matrix"} with the times
      at which each measure was performed }
    \item{\code{name}:}{Object of class \code{"character"} with information }
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "PlateChannel")}: ... }
    \item{measure<-}{\code{signature(.Object = "PlateChannel", value =
	"matrix")}: set the measure matrix }
    \item{measure}{\code{signature(.Object = "PlateChannel")}: get the
      measure matrix }
    \item{time<-}{\code{signature(.Object = "PlateChannel", value =
	"matrix")}: set the time matrix }
    \item{time}{\code{signature(x = "PlateChannel")}: get the time matrix }
	 }
}
%\author{ <laurent@cbs.dtu.dk> }
\note{
  The dimensions for the matrixes \code{"measure"} and \code{"time"}
  are such as there is one row per cell/well, and one column per cycle
  (it is assumed that the measures are taken in the same order, and a
  cycle occurs each time all well have been measured the same number of
  times since the begining of the experiment).
}
%
% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	A \code{PlateChannel} will often be found in a \code{\linkS4class{PlateRun}}
}
\examples{
showClass("PlateChannel")
}
\keyword{classes}
