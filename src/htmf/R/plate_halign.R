##
## Horizontally align a PlateRun for all or selected filtersets
##
## Returns: PlateRun
## Requires: channel_halign
##
plate_halign = function(
  x,       ## a PlateRun
  intervals     = NULL, ## has no effect if only one interval
  filtersets    = c("Biomass"),
  spar          = 0.8,
  extrapolating = TRUE
)
{
  ## if(!is.null(x)) { # panic }
  if( is.null(intervals) )
    return(x)
  if(is.null(filtersets))
    filtersets <- names(x)

  k <- length(names(x))  
  pcs <- vector("list", length=k)
  for( i in 1:k ) {
    chName <- names(x)[i]
    if(sum(chName == filtersets) > 0) { ## a valid filterset
  
      print(c("horizontally aligning", chName))
      chIntervals <- intervals
      if(is.list(intervals[[chName]])) {
      	print("Using channel and well specific intervals")
        chIntervals <- intervals[[chName]]
      }
      pcs[[i]] <- channel_halign(channels(x)[[i]], chIntervals, spar, extrapolating)
    }
    else {
      print(c("skipping", chName))
      pcs[[i]] <- channels(x)[[i]]
    }
  }
  adf <- new("AnnotatedDataFrame", data = cData(x) )
  pr <- PlateRun(pcs, adf, info = x@info)
  return(pr)
}

