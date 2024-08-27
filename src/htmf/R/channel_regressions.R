channel_regressions <- function(
  x,  ## PlateChannel
  intervals = NULL
)
{
  X = time(x)
  Y = measure(x)

  ## LOG TRANSFORMATION SHOULD BE DONE PRIOR 
  ## Y <- log2( Y + 1-min(Y) )
  
  if(is.null(intervals))
    intervals = list(range(X))

  n <- dim(X)[2]
  m <- dim(X)[1]
  kk <- length(intervals)
  
  ## y-intercept, slope, lower-bound(2.5%), upper-bound(97.5%)
  lrDat <- array(NA, dim=c(kk, m, 4), 
                dimnames=list(Intervals=1:kk, Wells=colnames(X),
                              Stats=c("Intercept","x","lb","ub")))

  ## Foreach x,y signal pair 
  for(j in 1:m) {
    x.well <- X[j,]
    y.well <- Y[j,]

    for(ii in 1:kk) { ## Foreach interval defined
      t.points <- intervals[[ii]]
      inds <- which(t.points[1] <= x.well & x.well <= t.points[2])
      reg <- lm(B ~ A, data=data.frame(A=x.well[inds], B=y.well[inds]))

      lrDat[ii,j,1:2] <- reg$coefficients
      lrDat[ii,j,3] <- confint(reg)[2]
      lrDat[ii,j,4] <- confint(reg)[4]
    }
  }
  return(lrDat)
}