\name{channel_interpolate}
\Rdversion{1.1}
\alias{channel_interpolate}
\docType{data}
\title{
Spline interpolate all wells of a PlateChannel for a single time or time series.
}
\description{
Estimates well amplitudes for a specified set of times or possibly a 
single time point. If times is NULL, a single time series will be 
estimated from the channel times (currently min time over well samples).
This function can be used to interpolate all well amplitudes for a single time series
or used to estimate a single amplitude, or offset, for well correction.
Additionally, if deriv = 1, then first derivative estimates are returned. 
If first derivatives are to be used to estimate growth rates, then amplitudes
should be log transformed prior to using this function.
}
\usage{channel_interpolate(x, times = "actual", spar = 0.8, deriv = 0)}
\format{
  The format is:
 chr "channel_interpolate"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
This function can be used to estimate the amplitude of all wells for a particular time point or time series.
Currently does not have a time intervals parameter. Should be extended to handle a time matrix.
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
%%  ~~ possibly secondary sources and usages ~~
Christopher Workman
}
%\examples{
%%channel_interpolate(dat$Biomass, times=3)
%}
\keyword{datasets}
