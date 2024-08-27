\name{channel_regressions}
\Rdversion{1.1}
\alias{channel_regressions}
\docType{data}
\title{
Apply linear regression to each well of a PlateChannel
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Returns an array with 3 dimensions: Intervals, Wells, Stats. 
Function is intended to operate on any number of specified time intervals (list of time pairs). 
The rerurned regression statistics (Stats) are: Intercept, slope (x), 
slope at 2.5% confidence lower-bound (lb), 
and slope at 97.5% confidence upper-bound (ub). 
}
\usage{channel_regressions(x, intervals=NULL)}
\format{
  The format is:
  PlateChannel x
  list intervals
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
 interval <- list( c(6,8) )
 %%lrDat <- channel_regressions(dat$Biomass, interval)
}
\keyword{datasets}
