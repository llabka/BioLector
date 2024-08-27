plate_normalize <- function(
    x,  ## a PlateRun
    xRef = NULL,     # DEPRECATED
    filtersets = NULL,
    intervals = NULL,
    coefNames = NULL,  ## c("slope", "offset")
    outliers = FALSE,  ## correct outliers
    baseline = FALSE,  ## fix the baseline instability
    valign  = FALSE,   ## align wells vertically to target time
    halign  = FALSE,   ## align across intervals, if any
    targetTime = NULL, ## used/infered if valign true
    targetAmp  = NULL, ## used/infered if valign true
    LOSflag    = FALSE,
    x.min = 3,
    x.max = 10,
    fdThr = 1e-2,
    sdThr = 0,
    lspar = 0.4,
    y.min   = NULL,
    ctrim   = 0,
    spar    = 0.4,     ## for plate_halign
    outliersThr = 10,  ## threshold z-score for outlier detection and correction
    baselineThr = 5,   ## threshold z-score for baseline detection and correction
    filtersets.log = NULL,      ## list of channel names c("Biomass", "GFP")
    plotting = FALSE
)
{
  if(is.null(intervals))
    intervals <- plate_intervals(x, 1)

  ## if there is a temerature reference (DEPRECATED)
  if( !is.null(xRef) )
    x <- plate_compensation(x=x, ref=xRef, intervals=intervals,
            filtersets=filtersets, spar=spar, y.min=y.min, plotting=plotting)

  ## linear transformation to OD600 or other biomass units
  if(!is.null(coefNames)) {
    for(fs in filtersets) {
      for(j in which(fs == names(x))) {
        #j <- grep(fs, names(x))
        mi <- measure(channels(x)[[j]])
        ## NEED ERROR HANDLING IF coefNames NOT FOUND
        a <- cData(x)[, coefNames] ## the (slope, offset) pairs

        tmin <- apply(time(channels(x)[[j]]), 2, min, na.rm=TRUE)
        tmax <- apply(time(channels(x)[[j]]), 2, max, na.rm=TRUE)
        for(jj in 1:length(intervals)) {
      	  jjj <- 0
          if(length(coefNames)/2 == length(intervals))
            jjj <- (jj-1)*2

          inds <- which(intervals[[jj]][1] <= tmax & tmin <= intervals[[jj]][2])
          mi[,inds] <- mi[,inds] * a[,jjj+1] + a[,jjj+2]
        }
        measure(channels(x)[[j]]) <- mi
      }
    }
  }

  x <- sliceByTimeIntervals(x=x, intervals=intervals, cycleTrim=ctrim)

  ## correct outliers and/or baseline shift errors
  if(outliers)
    x <- plate_correct_outliers(x, intervals=intervals, thr=outliersThr)

  if(baseline)
    x <- plate_correct_baseline(x, intervals=intervals, thr=baselineThr)

  if(halign)
    x <- plate_halign(x=x, intervals=intervals, filtersets=filtersets,
                      spar=spar, extrapolating=TRUE)
  if(valign) {
    if(!LOSflag) {
  	  if(is.null(targetTime))
  	     targetTime <- 0.1 + intervals[[1]][1]
  	  print(paste("Target time for valign = ", targetTime, targetAmp))
      x <- plate_valign(x=x, Talign=targetTime, Yalign=targetAmp, filtersets=filtersets, spar=spar)
    }
    else {
      x = plate_valign(x=x, filtersets=filtersets, spar=spar, LOS=LOSflag,
                       xMin=x.min, xMax=x.max, yMin=y.min,
                       fdThr=fdThr, sdThr=sdThr, lspar=lspar, plotting=plotting)
    }
  }

  ## enforce a minimum y-value
  if(!is.null(y.min)) {
  	for(fs in filtersets) {
  	  for(j in which(fs == names(x))) {
        mi <- measure(channels(x)[[j]])
        mi[mi < y.min] <- y.min
        measure(channels(x)[[j]]) <- mi
      }
    }
  }

  if(!is.null(filtersets.log)) {
    for(cn in filtersets.log) {
      for(j in which(cn == names(x))) {
        mi <- measure(channels(x)[[j]])
        measure(channels(x)[[j]]) <- log2(mi)
      }
    }
  }
  return(x)
}



