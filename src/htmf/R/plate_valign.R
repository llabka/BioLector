## Vertically alignment the wells of a PlateRun
##
## Value returned: PlateRun
## Requires: getYOffsets, channel_interpolate
##
plate_valign = function(
   x,                    ## PlateRun,
   Talign      = 0,      ## time points to align
   Yalign      = NULL,   ## amplitude to align to (at Talign)
   filtersets  = NULL,
   spar        = 0.4,
   interpolate = TRUE,   ## use values taken from spline at time Talign
   LOS         = FALSE,  ## use log-of-slope method to find background levels
   xMin        = 3,      ## use log-of-slope method to find background levels
   xMax        = 10,     ## use log-of-slope method to find background levels
   yMin        = 1e-2,   ## use log-of-slope method to find background levels
   fdThr       = 1e-2,   ## use log-of-slope method to find background levels
   sdThr       = 0,      ## use log-of-slope method to find background levels
   lspar       = 0.04,   ## use log-of-slope method to find background levels
   plotting    = FALSE
)
{
  if(is.null(filtersets))
    filtersets <- names(x)

  pcs <- vector("list", length = length(channels(x)))

  for( i in 1:length(channels(x)) ) {
  	pc <- channels(x)[[i]]
    if(sum(names(x)[i] == filtersets) > 0) {
      print(paste("plate_valign", names(x)[i]))
      pcs[[i]] = channel_valign(x=pc, Talign=Talign, Yalign=Yalign, spar=spar, interpolate=interpolate,
                                LOS=LOS, xMin=xMin, xMax=xMax, yMin=yMin, fdThr=fdThr, sdThr=sdThr, lspar=lspar,
                                plotting=plotting)
    }
    else
       pcs[[i]] <- pc
  }
  adf <- new("AnnotatedDataFrame", data = cData(x) )
  pr <- PlateRun(pcs, adf, info = x@info)
  return(pr)
}
