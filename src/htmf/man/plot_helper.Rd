\name{plot_helper}
\Rdversion{1.1}
\alias{plot_helper}
\docType{data}
\title{
PlateRun plot halpre function
}
\description{
Calls plot_biolector after performing a number of normalization steps.
}
\usage{plot_helper(
   x,                   ## a PlateRun
   intervals   = NULL,
   plotBy      = 0,     ## cData column name or index OR -1 or "Wells"
   colorBy     = 1,     ## cData column name or index
   colorFcn    = rainbow,
   filtersets  = NULL,  ## if NULL, all filtersets/channels will be plotted
   spar        = 0.8,
   plotName    = "",
   ...)
}
\format{
  The format is:
 chr "plot_helper"
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
data(plot_helper)
## maybe str(plot_helper) ; plot(plot_helper) ...
}
\keyword{datasets}
