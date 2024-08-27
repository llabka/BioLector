\name{buildColors}
\Rdversion{1.1}
\alias{buildColors}
\docType{data}
\title{
buildColors
}
\description{
A helper function to convert a study factor to a color scheme
}
\usage{buildColors(
	 wellTypes,
    wellColors = NULL, # rainbow(length(levels(as.factor(wellTypes))))
    colorFcn = rainbow
)}
\format{
  The format is:
 chr "buildColors"
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(buildColors)
## maybe str(buildColors) ; plot(buildColors) ...
}
\keyword{datasets}
