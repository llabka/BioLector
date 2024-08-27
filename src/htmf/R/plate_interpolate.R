##
## Spline interpolate a PlateRun for all or selected filtersets
##
## Returns: PlateRun
## Requires: channel_interpolate
##
## WARNING: intervals parameter not implemented
##
plate_interpolate = function(
  x,       ## a PlateRun
  filtersets    = NULL, ## c("Biomass"),
  spar          = 0.8,
  deriv         = 0,
  times         = "actual" ## string or numeric 1D-array of times
                           ## "min", "max", "mean", "actual" are valid
)
{
  ## if(!is.null(x)) { # panic }

  if(is.null(filtersets))
    filtersets <- names(x)

  k <- length(names(x))
  pcs <- vector("list", length=k)
  for( i in 1:k ) {
    if(sum(names(x)[i] == filtersets) > 0) {
  	  print(c("spline interpolating", names(x)[i]))
  	  pc <- channels(x)[[i]]
      pcs[[i]] <- channel_interpolate(pc, times=times, spar=spar, deriv=deriv)
  	}
  	#else
    #  pcs[[i]] <- channels(x)[[i]]
  }
  adf <- new("AnnotatedDataFrame", data = cData(x) )
  pr <- PlateRun(pcs, adf, info = x@info)
  return(pr)
}

