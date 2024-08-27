\name{htmf-package}
\Rdversion{1.1}
\alias{htmf-package}
\alias{htmf}
\docType{package}
\title{
Analysis of Highly-Parallel Micro-Fermentation Data at the Well Level
}
\description{
Collection of functions for reading, normalizing, plotting  and extracting 
features from growth profiles.
}
\details{
\tabular{ll}{
Package: \tab htmf\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-08-11\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Who wrote it

Maintainer: Who to complain to <yourfault@somewhere.net>
~~ The author and/or maintainer of the package ~~
}
%\references{
%~~ Literature or other references for background information ~~
%}
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
#library(ggplot2)
fn <- system.file("exampleData", "biolector_data.csv.gz", package="htmf")

pr <- read_biolector(gzfile(fn))

dataf <- melt(measure(pr$Biomass))
dataf <- merge(dataf, cData(pr), 
               by.x = "well", by.y = 0)
ggplot(dataf) + 
  aes(x = cycle, y = value, col = plate_col) + 
  geom_line() + facet_wrap( ~ plate_row)

dataf_ph <- melt(measure(pr$pH))
names(dataf_ph)[ncol(dataf_ph)] <- "pH"

dataf <- merge(dataf, dataf_ph, by=c("well", "cycle"))

% this code breaks the build...
%ggplot(dataf) +
%  aes(x = cycle, y = value / pH) +
%  geom_line(aes(group = well), alpha = 0.5) +
%  facet_wrap( ~ condition)


}
