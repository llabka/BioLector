\name{plate_normalize}
\Rdversion{1.1}
\alias{plate_normalize}
\docType{data}
\title{
Normalize a PLateRun
}
\description{
Normalize the selected PlateChannels of a PlateRun. 
}
\usage{plate_normalize(x,  ## a PlateRun
    xRef = NULL,
    filtersets = NULL,
    intervals = NULL,
    coefNames = NULL,  ## c("slope", "offset")
    outliers = FALSE,  ## correct outliers
    baseline = FALSE,  ## fix the baseline instability  
    valign  = FALSE,   ## align wells vertically to target time
    halign  = FALSE,   ## align across intervals, if any
    targetTime = NULL,       ## needed if valign true
    targetAmp  = NULL,
    y.min      = NULL,
    ctrim   = 0,
    spar    = 0.8,     ## for plate_halign
    outliersThr = 10,  ## threshold z-score for outlier detection and correction
    baselineThr = 5,   ## threshold z-score for baseline detection and correction    
    filtersets.log = NULL,      ## list of channel names c("Biomass", "GFP")
    plotting = FALSE
)
}
\format{
  The format is:
 chr "plate_normalize"
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
data(plate_normalize)
## maybe str(plate_normalize) ; plot(plate_normalize) ...
}
\keyword{datasets}
