\name{PlateRun-class}
\Rdversion{1.1}
\docType{class}
\alias{PlateRun-class}
\alias{initialize,PlateRun-method}
\alias{channels,PlateRun-method}
\alias{channels<-,PlateRun,container-method}
\alias{cData,PlateRun-method}
\alias{cData<-,PlateRun,data.frame-method}
\alias{$,PlateRun-method}
\alias{$<-,PlateRun-method}
\alias{info,PlateRun-method}
\alias{names,PlateRun-method}
\title{Experimental run an plate}
\description{ Store all channels measured during a run }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PlateRun",
  celldata, channels)}, or the constructor \code{PlateRun}.
}
\section{Slots}{
	 \describe{
    \item{\code{channels}:}{Object of class \code{"container"} of class \code{"PlateChannel"}}
    \item{\code{celldata}:}{Object of class \code{"AnnotatedDataFrame"}
      with information for each cell}
    \item{\code{info}:}{Object of class \code{"list"}
      with information on the plate and the run}
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "PlateRun")}: ... }
    \item{cData}{\code{signature(.Object = "PlateRun")}: Return the
         \code{pData} for the slot \code{celldata}.
	 }
    \item{cData<-}{\code{signature(.Object = "PlateRun", value = "data.frame")}: Replace the
         \code{pData} for the slot \code{celldata}.
	 }
    \item{channels}{\code{signature(.Object = "PlateRun")}: Extract the
         container of \code{PlateChannel}s }
    \item{channels<-}{\code{signature(.Object = "PlateRun", value = "container")}: Replace
         the container of \code{PlateChannel}s }
    \item{names}{\code{signature(.Object = "PlateRun")}: Return the
	 names for the channels measured
	 }
    \item{info}{\code{signature(.Object = "PlateRun")}: Return 
         information about the the plate / run
	 }
    \item{$}{\code{signature(.Object = "PlateRun", name)}: 
         Extract the channel of the given \code{"name"}
	 }
    \item{$<-}{\code{signature(.Object = "PlateRun", name)}: 
         Replace the channel of the given \code{"name"}
	 }
   }
}
%\author{}
\note{
  The cell-associated data are stored in an \code{AnnotatedDataFrame}.
  The method \code{cData} extracts directly the \code{data.frame}. 
}
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%
\seealso{
	\code{\linkS4class{PlateChannel}}
}
%\examples{
%  
%}
\keyword{classes}
