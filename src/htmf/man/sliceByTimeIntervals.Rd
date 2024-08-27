\name{sliceByTimeIntervals}
\Rdversion{1.1}
\alias{sliceByTimeIntervals}
\docType{data}
\title{
Trim cycles of a PlateRun
}
\description{
Based on a list of time pairs (list-of-lists), remove cycles not covered by the time intervals.
}
\usage{sliceByTimeIntervals(
    x,                ## a PlateRun
    intervals = NULL,
    cycleTrim = 0     ## number of cycles between intervals to trim
)
}
\format{
  The format is:
  chr "sliceByTimeIntervals"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
   TBA
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
   TBA
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(sliceByTimeIntervals)
## maybe str(sliceByTimeIntervals) ; plot(sliceByTimeIntervals) ...
}
\keyword{datasets}
