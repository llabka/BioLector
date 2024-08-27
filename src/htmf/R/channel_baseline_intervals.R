##
## DEPRECETED
## use channel_baseline_correction()
##
## makes intervals for each well
##
channel_baseline_intervals <- function(
  x,            ## a PlateChannel 
  windowK = 4,  ## size of the sliding window for testing
  thr = 1e-5
)
{
  A <- measure(x)
  T <- time(x)

  intrvls <- list()
  wells <- rownames(A)
  for(j in 1:length(wells)) {
    xi <- T[j,]
    yi <- A[j,]
    wellIntrvls <- list()
    n <- length(xi)
    yps <- array(1, n) ## Watch the default value here...
    for(i in windowK:(n-windowK)) {
      w1inds <- (i-windowK):i
      w2inds <- (i+1):(i+windowK)
      if(var(yi[(i-windowK):(i+windowK)]) > 0.01) {
        yps[i] <- t.test(yi[w1inds], yi[w2inds])$p.value
      }
      #yps[i] <- abs( median(yi[w1inds]) - median(yi[w2inds]) )
    }
    inds <- which(yps < thr)
    #inds <- which(yps > thr)
    
    #print(c(max(yps), which(yps==max(yps))))
    #print(yps[inds])
    #print(xi[inds])

    tms <- c(range(xi)[1], xi[t(cbind(inds, inds))], range(xi)[2])
    k <- length(tms)
    for(i in 1:(k/2)) {
      ii <- 1 + (i-1)*2
      wellIntrvls[[i]] <- tms[ii:(ii+1)]
    }
    intrvls[[wells[j]]] <- wellIntrvls
  }
  return(intrvls)
}

