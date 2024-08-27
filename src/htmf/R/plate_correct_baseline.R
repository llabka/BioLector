plate_correct_baseline <- function(
  x,  ## a PlateRun 
  intervals = NULL,
  thr = 5  ## p-value threshold for text
)
{
  chNames <- names(x)
  kk <- length(chNames)
  for( i in 1:length(chNames) ) {
  	channels(x)[[i]] <- channel_correct_baseline(x=channels(x)[[i]],
  	                                  intervals=intervals, thr=thr)
  }
  return(x)
}