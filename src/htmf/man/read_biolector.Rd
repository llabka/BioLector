\name{read_biolector}
\Rdversion{1.1}
\alias{read_biolector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read data from a DASGIP biolector
}
\description{
Read data from a DASGIP biolector
}
\usage{
read_biolector(
    file,
    measureType = "AMP" # c("AMP", "REF", "TEMP", "HUM")
)
}
\arguments{
  \item{file}{A connection or a filename}
  \item{measureType}{Abbreviation for the measurement type to be read.}
}
\details{
 Read data in a .csv file produce by a biolector machine/software into a \code{PlateRun}
\code{measureType}
"AMP" = Signal for the specified filter ("AMPLITUDE") for each well and each filterset
"REF" = Signal of the reference measurement, one per cycle (used for temperature compensation)
"TEMP" = Temperature ("ACT TEMP [degC]") as measured for each well and each filterset
"HUM" = Relative humidity ("ACT HUMIDITY [rH]") as measured for each well and each filterset

}
\value{
An object of class \code{PlateRun}
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
%%  ~~who you are~~
   Christopher Workman & Laurent Gautier
}
\note{
%%  ~~further notes~~
   TBA
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\linkS4class{PlateRun}}
}
\examples{

fn <- system.file("exampleData", "biolector_data.csv.gz", package="htmf")
# Create a gzip connection
gzf <- gzfile(fn)

pr <- read_biolector(gzf)

# Close the connection
close(gzf)

# list channel names
names(pr)

# extract date from one of the filters
biomass <- pr$Biomass

library(reshape)
dataf <- melt(measure(biomass))

ggplot(dataf) + 
  aes(x=cycle, y=value, group=well, col=well) + 
  geom_line()

dataf <- merge(
      melt(measure(biomass)),
      melt(time(biomass)),
      by = c("well", "cycle"),
      suffixes = c(".biomass", ".time"))

ggplot(dataf) + 
  aes(x=value.time, y=value.biomass, group=well, col=well) + 
  geom_line()

ggplot(melt(measure(pr$GFP))) + 
  aes(x=cycle, y=value, group=well, col=well) + 
  geom_line()


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

