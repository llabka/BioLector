\name{plate_calibrate}
\Rdversion{1.1}
\alias{plate_calibrate}
\docType{data}
\title{
Calibrate a PlateRun
}
\description{
Calibrate the measures of PlateRun. 
}
\usage{plate_calibrate(x,  ## a plateRun
  filterset = "Biomass", 
  replicates = NULL,  ## which wells will be used for each fit
  times = NULL,       ## times that correspond to the columns of calibDat
  calibDat = NULL,    ## matrix/df of ODs or other calibration data
                      ## with well-labels as rownames
  fctName = "replicates",
  plotting = FALSE,
  ...
)
}
\format{
  The format is:
 chr "plate_calibrate"
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
data(plate_calibrate)
## maybe str(plate_calibrate) ; plot(plate_calibrate) ...
}
\keyword{datasets}
