\name{plotBioLector}
\Rdversion{1.1}
\alias{plotBioLector}
\docType{data}
\title{
Plot a PlateRun
}
\description{
A wrapper for plot_helper and, in turn, plot_biolector.
}
\usage{plotBioLector(
  x,   ## a plateRun
  oFile = NULL,  ## pdf output
  sFile = NULL,  ## integration output
  muFile = NULL, ## muMax output table
  yFile  = NULL,  ## interpolated amp
  lmFile = NULL,
  signalT = "AMP",
  filterS = "Biomass",
  plotBy  = 2,
  colorBy = 1,
  colorFcn = rainbow,
  xRef = NULL,
  plotRef = FALSE, ## should we plot the refernce samples?
  outliers = FALSE,
  outliersThr = 10,
  baseline = FALSE,
  baselineThr = 10,
  targetTime = NULL,
  targetAmp  = NULL,
  valign = FALSE,
  halign = FALSE,
  ctrim  = 0,
  spar   = 0.8,
  filterS.log = NULL,
  xlim   = NULL,
  ylim   = NULL,
  prows  = 1,
  pcols  = 1,
  lm.show = FALSE,
  lm.intervals = NULL,
  time.samples = NULL,
  coefNames = NULL,      ## to perform linear transformation
  ...   ## will be passed down to the plot() in plot_biolector
) 
}
\format{
  The format is:
 chr "plotBioLector"
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
data(plotBioLector)
## maybe str(plotBioLector) ; plot(plotBioLector) ...
}
\keyword{datasets}
