\name{getYOffsets}
\Rdversion{1.1}
\alias{getYOffsets}
\docType{data}
\title{
Estimate the signal intensity at a fixed time point
}
\description{
Using splines, interpolate the y-offset for each well in a PlateChannel given an common time and amplitude.
}
\usage{getYOffsets(
  x,                ## a PlateChannel object
  target = 0,       ## x coordinate (time) of the alignment point
  spar   = 0.8
)}
\format{
  The format is:
 chr "getYOffsets"
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(getYOffsets)
## maybe str(getYOffsets) ; plot(getYOffsets) ...
}
\keyword{datasets}
