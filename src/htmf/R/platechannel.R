# Considers each measurement gotten from an experimental
# "run" on the biolecor as a "channel"
# (e.g., number of cells, CO2, green fluorescence protein, etc...)


# AssayData seems extremely convoluted (and close to impossible
# to use/extend without loosing sanity).
setClass("PlateChannel",
         representation(
           measure = "matrix",
           time = "matrix",
           ##measureRef = "vector",
           ##timeRef = "vector",
           name = "character"
         ),
         prototype = prototype(
           measure= matrix(),
           time = matrix(),
           ##measureRef = vector(),
           ##timeRef = vector(),
           name = character()
           )
         )

setMethod("initialize",
          signature(.Object = "PlateChannel"),
          function(.Object, name, measure, time) {
            if (! identical(dim(measure), dim(time)) ) {
              stop("Mismatching dim() when initializing a PlateChannel")
            }
            if (! identical(dimnames(measure), dimnames(time)) ) {
              stop("Mismatching dimnames() when initializing a PlateChannel")
            }
            if (! is.null(colnames(measure)) ) {
              stop("Column names correspond to cycles (and should be set to NULL)")
            }
            callNextMethod(.Object,
                           measure = measure,
                           time = time,
                           name = name)
          }
          )

PlateChannel <-
  function(name, measure_matrix, time_matrix)
{
  new("PlateChannel", name, measure_matrix, time_matrix)
}

setGeneric("time",
           function(x, ...) {
             standardGeneric("time")
           })
setMethod("time", signature(x="PlateChannel"),
          function(x, ...) {
            x@time
          })
setGeneric("time<-",
           function(.Object, ..., value) {
             standardGeneric("time<-")
           })
           
setReplaceMethod("time", signature("PlateChannel", "matrix"),
          function(.Object, value) {
            if (!identical(dim(.Object@time), dim(value))) {
              stop("The new time matrix must be of the same size as the old time matrix")
            }
            if (!identical(dimnames(.Object@time), dimnames(value))) {
              stop("The new time matrix must have the same dimnames() as the old time matrix")
            }
            .Object@time <- value
            .Object
          })
setGeneric("measure",
           function(.Object, ...) {
             standardGeneric("measure")
           })
setMethod("measure", signature("PlateChannel"),
          function(.Object) {
            .Object@measure
          })
setGeneric("measure<-",
           function(.Object, ..., value) {
             standardGeneric("measure<-")
           })
setReplaceMethod("measure", signature("PlateChannel", "matrix"),
          function(.Object, value) {
            if (!identical(dim(.Object@measure), dim(value))) {
              stop("The new measure matrix must be of the same size as the old measure matrix")
            }
            if (!identical(dimnames(.Object@time), dimnames(value))) {
              stop("The new measure matrix must have the same dimnames() as the old measure matrix")
            }

            .Object@measure <- value
            .Object
          })



setClass("PlateRun",
         representation(
           channels = "container",
           celldata = "AnnotatedDataFrame",
           info = "list"
                        ),
         prototype = prototype(
           channels = new("container", content = "PlateRun")
           )
         )


setMethod("initialize",
          signature(.Object = "PlateRun"),
          function(.Object,
                   channels,
                   celldata,
                   info = list()
                   ) {
            content = "PlateChannel"
            channels <- new("container",
                            x = channels,
                            content = content)
            # check dimensions
            mydim <- NULL
            mydimnames <- NULL
            for (elt_i in seq(along = channels)) {
              if (is.null(mydim)) {
                mydim <- dim(channels[[elt_i]]@time)
              } else if (! identical(mydim, dim(channels[[elt_i]]@time))) {
              	print(mydim)
              	print(dim(channels[[elt_i]]@time))
                stop("Mismatching dim() in the channels")
              }
              if (is.null(mydimnames)) {
                mydimnames <- dimnames(channels[[elt_i]]@time)
              } else if (! identical(mydimnames, dimnames(channels[[elt_i]]@time))) {
                stop("Mismatching dimnames() in the channels")
              }              
            }
            if (nrow(celldata) != nrow(time(channels[[1]]))) {
              stop("celldata should have as many rows as rows in the matrixes in channels")
            }
            .Object@celldata <- celldata
            .Object@channels <- channels
            .Object@info <- info
            .Object
          }
          )

PlateRun <-
  function(list_channels, adf, info)
{
  if (missing(adf)) {
     adf <- new("AnnotatedDataFrame",
                data = data.frame(wellnum = rownames(time(list_channels[[1]]))))
  }
  if (missing(info)) {
    info <- list()
  }
  new("PlateRun", list_channels, adf, info)
}

setMethod("show", "PlateRun",
          function(object) {
              cat("PlateRun object with info:\n")
              for (i in seq(along = info(object))) {
                cat(names(info(object))[i], ":\n")
                tmp <- info(object)[[i]]
                for (j in seq(along = tmp)) { 
                  cat("  ", names(tmp)[j], ": " )
                  cat(paste(tmp[[j]], collapse = ","), "\n")
                }
              }
          })

setMethod("names", "PlateRun",
          function(x) {
            n <- length(x@channels)
            res <- rep("", n)
            for (i in 1:n) {
              res[i] <- x@channels[[i]]@name
            }
            res
          })


setGeneric("channels",
           function(.Object, ...) {
             standardGeneric("channels")
           })

setMethod("channels", "PlateRun",
          function(.Object) {
            .Object@channels
          })

setGeneric("channels<-",
           function(.Object, value) {
             standardGeneric("channels<-")
           })

setReplaceMethod("channels", c("PlateRun", "container"),
          function(.Object, value) {
            .Object@channels <- value
	    .Object
          })


setMethod("$", "PlateRun",
	  function(x, name) {
            # container does not recognise "names" :/
            i <- which(names(x) == name)
            x@channels[[i]]
          })

setReplaceMethod("$", "PlateRun",
	 function(x, name, value) {
            x@channels[[name]] <- value
            x
          })

setGeneric("info",
           function(.Object, ...) {
             standardGeneric("info")
           })

setMethod("info", signature("PlateRun"),
          function(.Object) {
            .Object@info
          })


setGeneric("cData",
           function(.Object, ...) {
             standardGeneric("cData")
           })
setMethod("cData", signature("PlateRun"),
          function(.Object) {
            pData(.Object@celldata) 
          })
setGeneric("cData<-",
           function(.Object, ..., value) {
              standardGeneric("cData<-")
           }) 
setReplaceMethod("cData", signature("PlateRun", "data.frame"),
          function(.Object, value) {
            pData(.Object@celldata) <- value
            .Object
          })

