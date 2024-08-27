\name{plate_compensation}
\Rdversion{1.1}
\alias{plate_compensation}
\docType{data}
\title{
Signal compensation for a PlateRun
}
\description{
Apply channel compensation to selected channels of a PlateRun. Intended for use with the temperature reference.
}
\usage{plate_compensation(
  x,       ## a PlateRun
  ref,     ## a PlateRun
  intervals   = NULL, ## has no effect if only one interval
  filtersets  = NULL, ## c("Biomass"),
  spar        = 0.8,
  y.min       = NULL,
  plotting    = TRUE
)
}
\format{
  The format is:
 chr "plate_compensation"
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
data(plate_compensation)
## maybe str(plate_compensation) ; plot(plate_compensation) ...
}
\keyword{datasets}
