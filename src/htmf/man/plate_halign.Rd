\name{plate_halign}
\Rdversion{1.1}
\alias{plate_halign}
\docType{data}
\title{
Horizontally align a PlateRun
}
\description{
Horizontally align selected channels over time intervals for a PlateRun using align_channel.
}
\usage{ plate_halign(x, intervals, filtersets=c("Biomass"), spar=0.8, extrapolating=TRUE)}
\format{
  The format is:
 chr "plate_halign"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
  Align each well across the time intervals. 
  Motivation: System pauses cause random Y-offsets that need to be removed. 
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
data(halign_plate)
## maybe str(plate_halign) ; plot(plate_halign) ...
}
\keyword{datasets}
