\name{plate_valign}
\Rdversion{1.1}
\alias{plate_valign}
\docType{data}
\title{
Vertically align wells of a PlateRun
}
\description{
Vertically align wells of a PlateRun for selected channles given a target time and intensity
}
\usage{ plate_valign(x, Talign=0, Yalign=NULL, filtersets=NULL, spar=0.8, interpolate = TRUE) }
\format{
  The format is:
 chr "plate_valign"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
   Align all wells to a fixed xy coordinate (Talign, Yalign)
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
 TBA
}
\references{
%%  ~~ possibly secondary sources and usages ~~
 TBA
}
\examples{
data(plate_valign)
## maybe str(plate_valign) ; plot(plate_valign) ...
}
\keyword{datasets}
