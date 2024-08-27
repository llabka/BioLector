##
## Temperature compensation of signals for a PlateRun 
## for all or selected filtersets
##
## Returns: PlateRun
## Requires: channel_compensation
##
plate_compensation = function(
  x,       ## a PlateRun
  ref,     ## a PlateRun
  intervals   = NULL, ## has no effect if only one interval
  filtersets  = NULL, ## c("Biomass"),
  spar        = 0.8,
  y.min       = NULL,
  plotting    = TRUE
)
{
  ## if(!is.null(x)) { # panic }

  if(is.null(filtersets))
    filtersets <- names(x)

  k <- length(names(x))  
  pcs <- vector("list", length=k)
  for( i in 1:k ) {
    if(sum(names(x)[i] == filtersets) > 0) {
  	  print(paste("temp compensation for", names(x)[i]))
      pcs[[i]] <- channel_compensation(x=channels(x)[[i]],
                                       ref=channels(ref)[[i]], 
                                       intervals=intervals, spar=spar, 
                                       y.min = y.min, plotting=plotting)
  	}
  	else
      pcs[[i]] <- channels(x)[[i]]
  }
  adf <- new("AnnotatedDataFrame", data = cData(x) )
  pr <- PlateRun(pcs, adf, info = x@info)
  return(pr)
}

