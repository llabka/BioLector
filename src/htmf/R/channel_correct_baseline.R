##
## corrects baseline shifts
##
channel_correct_baseline <- function(
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
    for(jj in 1:length(intervals)) {
      intrvl <- intervals[[jj]]
  	  inds <- which(intrvl[1]<=T[j,] & T[j,]<=intrvl[2])
      xi <- T[j,inds]
      yi <- A[j,inds]
      n <- length(xi)

      yi1 <- yi[2:n] - yi[1:(n-1)]
      zi1 <- scale(yi1, center=TRUE, scale=TRUE)
      baselineShifts <- abs(zi1) > thr
      bInds <- which(baselineShifts)
   
      if(sum(baselineShifts, na.rm=TRUE) > 0)
        print(paste("Baseline shifts detected: well=", wells[j], 
                    ", samples=", paste(inds[bInds], collapse=","),
                    ", z=", paste(round(zi1[bInds],2), collapse=","), sep=""))

      #print(bInds)
      k <- length(bInds)
      dy <- array(0, n)
      for(ii in bInds) { ## single offset
      	yOffset <- yi[ii]
        ii <- ii + 1
        ##print(c(j,jj,i,ii))
        if(ii > 2) {
          preii <- max(1, ii-4):(ii-1)
          lmRes <- lm(y~x, data.frame(x=xi[preii], y=yi[preii]))
          yOffset <- (xi[ii]*lmRes$coeff[2]+lmRes$coeff[1])
        }
        dy[ii:n] <- dy[ii:n] + (yOffset - yi[ii])
      }
      measure(x)[j,inds] <- yi + dy
    }
  }
  return(x)
}

