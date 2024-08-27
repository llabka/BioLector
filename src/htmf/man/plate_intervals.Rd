\name{define_intervals}
\Rdversion{1.1}
\alias{define_intervals}
\docType{data}
\title{
Determine the time intervals for a PLateRun
}
\description{
Creates a list of time pairs that exclude the run pauses.
}
\usage{plate_intervals(
  x,           ## a PlateRun
  thresh = 0.2 ## percent cycle time exceeded
)
}
\format{
  The format is:
 chr "define_intervals"
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
data(define_intervals)
## maybe str(define_intervals) ; plot(define_intervals) ...
}
\keyword{datasets}
