##
## DEPRECATED
## use plate_baseline_correction()
##
plate_baseline_intervals <- function(
  x,  ## a PlateRun 
  k = 4,  ## window size for test
  thr = 1e-5  ## p-value threshold for text
)
{
  chNames <- names(x)
  #print(chNames)
  kk <- length(chNames)
  pIntervals <- list()
  for( i in 1:kk ) {
  	cName <- chNames[i]
  	#print(cName)
  	pIntervals[[i]] <- channel_baseline_intervals(x=channels(x)[[i]],
  	                                              windowK=k, thr=thr)
  }
  names(pIntervals) <- chNames
  return(pIntervals)
}