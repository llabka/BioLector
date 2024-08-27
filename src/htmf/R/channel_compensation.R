channel_compensation <- function(
  x,    ## a PlateChannel
  ref,  ## a PlateChannel
  intervals = NULL,
  spar = 0.8,
  y.min = NULL,
  plotting = TRUE
)
{
  #if(!is.null(x)) { panic }

  X = time(x)  ## hours
  Y = measure(x) 
  Tref <- time(ref)
  Aref <- measure(ref)
  ##print(str(Tref))
  if(is.null(intervals))
    intervals <- list(range(X))

  ## adjust Aref
  Aref <- Aref - min(Aref)  

  n <- dim(X)[1]
  m <- dim(X)[2]
  t.min <- apply(X, 2, min, na.rm=TRUE)
  t.max <- apply(X, 2, max, na.rm=TRUE)
  ArefFit <- Aref  ## IS THIS NEEDED
  sspl <- list()
  k <- length(intervals)
  
  for(ii in 1:k) {
    tPoints <- intervals[[ii]]
    tInds <- which(tPoints[1] <= t.max & t.min <= tPoints[2])
    #print(range(tInds))
    sspl[[ii]] <- smooth.spline(Tref[1,tInds], Aref[1,tInds], spar=spar)
    ArefFit[1,tInds] <- sspl[[ii]]$y
  }
  
  if(plotting) {
  	#print(range(Aref))
    plot(Tref, Aref, pch=19, cex=0.5, main=x@name)
    for(ii in 1:k) {
      tPoints <- intervals[[ii]]
      tInds <- tPoints[1] <= t.max & t.min <= tPoints[2]
      lines(Tref[1,which(tInds)], ArefFit[1,which(tInds)], col="red")
      abline(v=tPoints[1])
      abline(v=tPoints[2])
    }
  }
  
  REF <- Y
  for(j in 1:n) {
  	for(ii in 1:k) {
  	  tPoints <- intervals[[ii]]
      tInds <- which(tPoints[1] <= t.max & t.min <= tPoints[2])
  	  REF[j,tInds] <- predict(sspl[[ii]], x=X[j, tInds])$y
    }
  }
  pc <- PlateChannel(x@name, Y-REF, X)
  if(!is.null(y.min))
    measure(pc)[measure(pc) < y.min] <- y.min
  
  return( pc )
}
