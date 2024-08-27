\name{plot_biolector}
\Rdversion{1.1}
\alias{plot_biolector}
\docType{data}
\title{
XY-plots for a PlateRun
}
\description{
Plot selected signals for selected channels (filter sets) of a PlateRun.
}
\usage{plot_biolector(
  x,                     ## a PlateRun
  filterset     = "Biomass",
  replicates    = NULL,  ## list of column indexes into X and Y
  intervals     = NULL,  ## list of start/stop time pairs (lists)
  colorList     = NULL,  ## same structure as replicates
  spar          = 0.8,
  show.spline   = TRUE,
  show.points   = TRUE,
  pch           = 19,
  pcex          = 0.3,
  ptype         = 'p',
  show.linearfit = FALSE,
  lm.intervals = NULL,
  xlim = NULL, 
  ylim = NULL, 
  xlab = "time", 
  ylab = "signal",
  legend.loc = NULL, ## "bottomright", "topleft", etc.
  legend.title = NULL,
  main.title = NULL,
  ...)
}
\format{
  The format is:
 chr "plot_biolector"
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
data(plot_biolector)
## maybe str(plot_biolector) ; plot(plot_biolector) ...
}
\keyword{datasets}
