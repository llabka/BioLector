##
## Removes cycles that have been interrupted by system pauses
##
plate_intervals <- function(
  x,           ## a PlateRun
  thresh = 0.2 ## percent cycle time exceeded
)
{
  k <- length( channels(x) )
  ##if( k < 1 ){ panic }
  
  T <- time(channels(x)[[1]])
  if(k > 1) {
    for( i in 2:k ) ## There must be a more elegant way to do this...
      T <- rbind(T, time(channels(x)[[i]]))
  }
  Tmin <- apply(T, 2, min, na.rm=TRUE)
  Tmax <- apply(T, 2, max, na.rm=TRUE)
  ##print(c(min(Tmin), max(Tmax)))
  m <- length(Tmin)
  dC <- Tmin[2:m] - Tmin[1:(m-1)] ## difference between cycle starts
  gaps <- which(dC > median(dC)*(1+thresh)) ## 50% more than the median

  if(length(gaps)>0){
  	eTimes <- c(range(T)[1], 
  	            array(rbind(Tmax[gaps-1], Tmin[gaps+1])),
  	            range(T)[2])
    kk <- length(eTimes)
    inds <- t(matrix(1:kk, 2, kk/2))
    intervals <- list()
    for(i in 1:dim(inds)[1]){
   	  intervals[[i]] <- c(eTimes[inds[i,]])
    }
  }
  else {
    intervals <- list(range(T))
  }
  return(intervals)
}
