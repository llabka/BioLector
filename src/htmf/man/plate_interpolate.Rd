\name{plate_interpolate}
\Rdversion{1.1}
\alias{plate_interpolate}
\docType{data}
\title{
Spline interpolate selected channels of a PlateRun.
}
\description{
Applies channel_interpolate to selected channels (filtersets) of a PlateRun.
}
\usage{plate_interpolate(x, filtersets = NULL, spar = 0.8, deriv = 0, times = "actual")}
\format{
  The format is:
 chr "plate_interpolate"
}
\details{
See channel_interpolate for more details.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
   TBA
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
 %%datSm <- plate_interpolate(dat, filtersets="Biomass")
## maybe str(plate_interpolate) ; plot(plate_interpolate) ...
}
\keyword{datasets}
