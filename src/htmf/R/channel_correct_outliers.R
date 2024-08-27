##
## corrects outliers
##
channel_correct_outliers <- function(
  x,            ## a PlateChannel
  intervals = NULL,
  thr = 5
)
{
  A <- measure(x)
  T <- time(x)
  
  if(is.null(intervals))
    intervals <- list(range(T))

  wells <- rownames(A)
  for(j in 1:length(wells)) {
    for(intrvl in intervals) {
  	  inds <- which(intrvl[1]<=T[j,] & T[j,]<=intrvl[2])
      xi <- T[j,inds]
      yi <- A[j,inds]
      n <- length(xi)
      dy <- array(0, n)
      vlgc <- is.finite(yi)
      dy[vlgc] <- yi[vlgc] - smooth(yi[vlgc], endrule="copy", do.ends=TRUE)
      if(sum(abs(dy)) > 0) {
        zy <- scale(dy, center=TRUE, scale=TRUE) #(dy-median(dy))/mad(dy)
        outliers <- abs(zy) > thr
        if(sum(outliers) > 0) {
          print(paste("Outliers detected: well=", wells[j], ", samples=",
                      paste(inds[which(outliers)],collapse=","), sep=""))
          yi[which(outliers)] <- smooth(yi[which(outliers)])
        }
      }
    }
  }
  return(x)
}
