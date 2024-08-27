\name{plate_regressions}
\Rdversion{1.1}
\alias{plate_regressions}
\docType{data}
\title{
Apply linear regression to a PlateRun
}
\description{
Using channel_regressions, estimate regression statistics for a selected filterset in a PlateRun 
possibly over a set of specified time intervals.
}
\usage{plate_regressions(x, intervals = NULL, filterset = "Biomass")}
\format{
  The format is:
  PlateRun x
  list intervals
  chr filterset
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
Simple wrapper for the channel_regressions function. Currently only works to extract 
regression statistics from one filterset for PlateRun (x).
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
   TBA
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Christopher Workman
}
\examples{
  interval <- list( c(6,8) )
  %%lrDat <- plate_regressions(dat, interval)
}
\keyword{datasets}
