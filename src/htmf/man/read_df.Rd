\name{read_df}
\Rdversion{1.1}
\alias{read_df}
\docType{data}
\title{
Read microtiter plate data and plate design from tab-delimited text files
}
\description{
Reads plate data with multiple channels(filters) from a simple text table.
  file     tab-delimited-text table of plate information

Calls df2plateRun which requires character lists:
      chanLbls
      designLbls 

If time label is not "Time" then timeLbl must also be passed to df2plateRun
}
\usage{read_df(
  file,
  wellLbl = "Well",
  ...
)
}
\format{
  The format is:
  chr "read_df"
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
data(read_df)
## maybe str(read_biolector2) ; plot(read_biolector2) ...
}
\keyword{datasets}
