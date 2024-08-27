\name{read_biolector2}
\Rdversion{1.1}
\alias{read_biolector2}
\docType{data}
\title{
Read BioLector and plate design table
}
\description{
Extends read_biolector() by reading and adding design informatio to the PlateRun
}
\usage{read_biolector2(
   file,
   designfile  = NULL,
   measureType = c("AMP", "REF", "TEMP", "HUM")
)
}
\format{
  The format is:
  chr "read_biolector2"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
  Read a BioLector csv file and a design file (metadata table) that describes 
  the contents of the wells and describes the experiment.
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
  TBA
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(read_biolector2)
## maybe str(read_biolector2) ; plot(read_biolector2) ...
}
\keyword{datasets}
