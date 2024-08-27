plot_biolector <- function(
  x,                     ## a PlateRun
  filterset     = "Biomass",
  replicates    = NULL,  ## list of column indexes into X and Y
  intervals     = NULL,  ## list of start/stop time pairs (lists)
  colorList     = NULL,  ## same structure as replicates
  spar          = 0.8,
  show.spline   = TRUE,
  show.points   = TRUE,
  pch           = 19,
  pcex          = 0.3,
  ptype         = 'p',
  show.linearfit = FALSE,
  lm.intervals = NULL,
  xlim = NULL, 
  ylim = NULL, 
  xlab = "time", 
  ylab = "signal",
  legend.loc = NULL, ## "bottomright", "topleft", etc.
  legend.title = NULL,
  main.title = NULL,
  ...
)
{
  #if(!is.null(x)) { panic }
  k <- length(channels(x)) ## number of channels

  ## for(fSet in filtersets) {
  chNum <- which(filterset == names(x))[1] ## the first matching channel name
  #if(length(chNum) == 0 | k == 0) { panic }
  
  pc <- channels(x)[[chNum]]
  X = time(pc)     ## hours
  Y = measure(pc)

  ## need to assert that X *has* row names
  if( is.null(replicates) )
    replicates <- split(1:dim(X)[1], rownames(X))
  if( is.null(intervals) )
    intervals <- list(range(X))
  if( is.null(colorList) ) {
    colorList <- array("black", dim(X)[1])
    names(colorList) <- rownames(X)
    legend.title <- ""
  }

  if(is.null(xlim)) #!is.finite(xlim)) {
    xlim <- range(intervals) ##xlim = range(X)
  
  if(is.null(ylim)) ##!is.finite(ylim))
    ylim = range(Y, na.rm=TRUE)

  if(is.null(lm.intervals))
    lm.intervals <- intervals
 
  
  wellsVSintervals <- match(rownames(X), names(intervals))
  wellIntervals <- sum(!is.na(wellsVSintervals)) == dim(X)[1]
    
  for(i in 1:length(replicates)) {
    s.m <- s <- names(replicates)[i]
    inds <- replicates[[s]]
    rep.colors <- colorList[inds]
    k <- length(inds)
    
    if(!is.null(main.title))
      s.m <- paste(main.title, ": ", s, sep="")

    plot(NA, type='n', main=s.m, sub=filterset, 
         xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
    
    for(j in 1:k) { ## start, foreach replicate
      x.well <- as.matrix(X[inds[j],])
      y.well <- as.matrix(Y[inds[j],])
      well <- rownames(X)[inds[j]]
      wIntervals <- intervals
      lm.wIntervals <- lm.intervals
      if(wellIntervals) {
         wIntervals <- intervals[[well]]
         lm.wIntervals <- lm.intervals[[well]]
      }
       kk <- length(wIntervals)
      for(ii in 1:kk) { ## start, foreach interval
        t.points <- wIntervals[[ii]]
        #t.inds <- which(t.points[1] < x.well & x.well < t.points[2])
        t.inds <- which(t.points[1] <= x.well & x.well <= t.points[2])
        ## fit the smoothing spline
        if(show.spline) {
          xx <- x.well[t.inds]
          yy <- y.well[t.inds]
          vlgc <- is.finite(xx) & is.finite(yy)
          sspl <- smooth.spline(xx[vlgc], yy[vlgc], spar=spar)
          lines(sspl$x, sspl$y, col=rep.colors[j], lty=1.2)
        }
        if(show.points)
          points(x.well[t.inds], y.well[t.inds], pch=pch, col=rep.colors[j], cex=pcex, type=ptype)
      }
      if(show.linearfit) { ## lm.intervals may be different than intervals
      	for(ii in 1:length(lm.wIntervals)) {
      	  t.points <- lm.wIntervals[[ii]]
          t.inds <- which(t.points[1] < x.well & x.well < t.points[2])
          reg <- lm(B ~ A, data=data.frame(A=x.well[t.inds], B=y.well[t.inds]))
          ## abline(a=reg$coefficients[1], b=reg$coefficients[2], col=rep.colors[j])
          lines(x.well[t.inds], reg$fitted.values, col=rep.colors[j])
        }
      }
      for(ii in 1:kk){ ## start 2.1
        t.points <- wIntervals[[ii]]
        if(ii > 1)
          abline(v=t.points[1])
        if(ii < kk)
          abline(v=t.points[2])
      }
    }
    if( !is.null(legend.loc) ) {
      rep.uniq <- sort( unique(names(rep.colors)) )
      ##rep.uniq <- as.character(sort( as.numeric(unique(names(rep.colors))) ))
      legend(legend.loc, legend=rep.uniq, fill=rep.colors[rep.uniq], 
             title=legend.title, bg=par("bg"))
    }
  }
}

