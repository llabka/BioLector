##
## Input: PlateChannel
## Value returned: array of y-offsets
##
## WARNING: intervals not implemented
##
channel_interpolate = function(
  x,                ## a PlateChannel object
  times     = "actual", ## c("actual","min","mean","max"),
                        ## may also be a 1D array of times
  spar      = 0.8,
  deriv     = 0
)
{
  X <- time(x)
  A <- measure(x)
  wellNames <- rownames(X)
  colNames <- colnames(X)
  m <- dim(X)[1]

  if(times[1] == "min")
    times <- apply(X, 2, min)
  if(times[1] == "mean")
    times <- apply(X, 2, mean)
  if(times[1] == "max")
    times <- apply(X, 2, max)

  if(times[1] == "actual") {
  	n = dim(X)[2]
  	new_time = X
  }
  else {
    n = length(times) ## assumes array context for times
    new_time = t( matrix(data=times, nrow=n, ncol=m,
                     dimnames=list(cycle = NULL, well = wellNames))
               )
  }
  new_measure = matrix(data=A[,1:n], nrow=m, ncol=n,
                        dimnames=list(well = wellNames, cycle = NULL))

  ## Foreach row (plate well)...
  for(j in 1:m) {
    xj <- X[j,]
    yj <- A[j,]
    vlgc <- is.finite(xj) & is.finite(yj)
    sspl = smooth.spline(xj[vlgc], yj[vlgc], spar=spar)
    if(times[1] == "actual")
      new_measure[j, vlgc] <- predict(sspl, x=xj[vlgc], deriv=deriv)$y
    else
      new_measure[j, ] <- predict(sspl, x=times, deriv=deriv)$y

  }

  return( PlateChannel(name = x@name,
                      measure_matrix = new_measure,
                      time_matrix = new_time)
        )
}