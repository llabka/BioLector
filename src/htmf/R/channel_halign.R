##
## Dependencies: smooth.spline, predict
##
channel_halign <- function(
  x,                    ## a Plate Channel Object
  intervals     = NULL, ## has no effect if only one interval
  spar          = 0.8,
  extrapolating = TRUE
)
{
    X <- time(x)
    Y <- Y.adj <- measure(x)
    m <- dim(X)[1]
    wellNames <- rownames(Y)
        
    ## Foreach column (growth well)...
    for(j in 1:m) {
      well <- wellNames[j]
      x.j <- as.matrix(X[j,])
      y.j <- as.matrix(Y[j,])
      wIntervals <- intervals
      if(!is.null(names(intervals)) & is.list(intervals[[well]]))
        wIntervals <- intervals[[well]]
      running.offset <- this.min <- last.max <- 0
      ## Foreach time interval in j
      for(ii in 1:length(wIntervals)) {
      	##print(c("channel_halign spar=", spar))
        t.points <- wIntervals[[ii]]
        #t.inds <- which(t.points[1] < x.j & x.j < t.points[2]) ## IS THIS THE DEFINITION WE WANT?
        t.inds <- which(t.points[1] <= x.j & x.j <= t.points[2] & is.finite(x.j) & is.finite(y.j)) ## 
        sspl <- NULL
        this.min <- y.j[t.inds[1]]
        if(length(unique(x.j[t.inds])) > 3) {
          sspl <- smooth.spline(x.j[t.inds], y.j[t.inds], spar=spar)
          this.min <- sspl$y[1]
        }
        ## subtract the difference from the last max and the current min
        if(ii > 1) {
          this.offset <- (y.j[t.inds[1]-1] - this.min)
          if(!is.null(sspl)) {
            last.max <- last.sspl$y[length(last.sspl$y)]
            if(extrapolating)
              last.max <- predict(last.sspl, sspl$x[1])$y
            this.offset <- (last.max - this.min)
          }
          Y.adj[j,t.inds] <- y.j[t.inds] + this.offset + running.offset
          running.offset <- running.offset + this.offset
          ## print(c(j, ii, round(c(last.max, this.min, last.max-this.min),2)))
        }
        last.sspl <- sspl
        ###last.sspl <- smooth.spline(x[t.inds], Y.adj[t.inds,j], spar=spar)
      }
    }
    ## measure(x) <- t(Y.adj)
    return( PlateChannel(name = x@name,
                         measure_matrix = Y.adj,
                         time_matrix = time(x)) )
}

