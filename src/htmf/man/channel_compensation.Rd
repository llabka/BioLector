\name{channel_compensation}
\Rdversion{1.1}
\alias{channel_compensation}
\docType{data}
\title{
Correct channel signal with a reference
}
\description{
Intended to correct a channel's signal with the temperature reference
}
\usage{channel_compensation(
  x,       ## a PlateChannel
  ref,     ## a PlateChannel
  intervals   = NULL, ## has no effect if only one interval
  spar        = 0.8,
  y.min       = NULL,
  plotting    = TRUE
)}
\format{
  The format is:
 chr "channel_compensation"
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
data(channel_compensation)
## maybe str(channel_compensation) ; plot(channel_compensation) ...
}
\keyword{datasets}
