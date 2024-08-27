sliceByTimeIntervals <- function(
    x,                ## a PlateRun
    intervals = NULL,
    cycleTrim = 0     ## number of cycles between intervals to trim
)
{
  if(is.null(intervals)) { ## do nothing
    return(x)
  }
  else {
  	k <- length( channels(x) )
  	## if(k < 1) { panic }
    pcs <- vector("list", length=k)
    tm <- time(channels(x)[[1]])
    if(k > 1) {
      for( i in 2:k ) ## There must be a more elegant way to do this...
        tm <- rbind(tm, time( channels(x)[[i]] ))
  	}
    Tmin <- apply(tm, 2, min, na.rm=TRUE)
    Tmax <- apply(tm, 2, max, na.rm=TRUE)

    for( i in 1:k ) {
      pc <- channels(x)[[i]]
      pcName <- pc@name
      inds <- list()
      for(j in 1:length(intervals)) {
        ivl <- intervals[[j]]
        ivlInds <- which( ivl[1] <= Tmin & Tmax <= ivl[2] )
        kk <- length(ivlInds)
        if(cycleTrim > 0 & cycleTrim < kk) {
          ivlInds <- ivlInds[(1+cycleTrim):(kk-cycleTrim)]
        }
        inds[[j]] <- ivlInds
      }
      ii <- unlist(inds)
      ni <- dim(measure(pc))[1]
      mi <- length(ii)
      newMeasure <- matrix(measure(pc)[,ii], ni, mi)
      newTime <- matrix(time(pc)[,ii], ni, mi)
      dimnames(newMeasure) <- dimnames(newTime) <- dimnames(measure(pc))
      pcs[[i]] <- PlateChannel(name=pcName,
                               measure_matrix = newMeasure,
                               time_matrix = newTime)
    }
    adf <- new("AnnotatedDataFrame", data = cData(x) )
    pr <- PlateRun(pcs, adf, info = x@info)
    return(pr)
  }
}
