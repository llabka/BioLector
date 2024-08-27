plotBioLector <- function(
  x,     ## a plateRun
  oFile  = NULL,  ## pdf output
  sFile  = NULL,  ## integration output (IN ONE RESULTS TABLE)
  muFile = NULL,  ## muMax output table (IN ONE RESULTS TABLE)
  yFile  = NULL,  ## interpolated amp   (IN ONE RESULTS TABLE)
  lmFile = NULL,
  signalT = "AMP",
  filterS = "Biomass",
  plotBy  = 2,
  colorBy = 1,
  colorFcn = rainbow,
  xRef = NULL,      ## a plateRun used for temp compensation (DEPRECATED)
  plotRef = FALSE, ## should we plot the refernce samples?
  outliers = FALSE,
  outliersThr = 10,
  baseline = FALSE,
  baselineThr = 10,
  targetTime = NULL,
  targetAmp  = NULL,
  valign = FALSE,
  halign = FALSE,
  ctrim  = 0,
  spar   = 0.8,
  filterS.log = NULL,
  xlim   = NULL,
  ylim   = NULL,
  prows  = 1,
  pcols  = 1,
  lm.show = FALSE,
  lm.intervals = NULL,
  time.samples = NULL,
  coefNames = NULL,      ## to perform linear transformation
  ...   ## will be passed down to the plot() in plot_biolector
) 
{
  ii <- which(filterS == names(x))[1]
  colNames <- colnames(cData(x))

  pdf(oFile)
  par(mfrow=c(prows, pcols))
  
  if(!is.null(xlim)) {
    x <- sliceByTimeIntervals(x, intervals=list(xlim))
  }

  intervals <- plate_intervals(x, 0.3)
  print(c("targetAmp", targetAmp, "intervals", intervals))

  ## Normalization
  if(!is.null(targetAmp) & is.character(targetAmp)) {
    targetAmp <- cData(x)[,targetAmp]
    targetAmp[is.na(targetAmp)] <- median(targetAmp, na.rm=TRUE)
    print(c("targetAmp:", targetAmp))
  }

  dat <- plate_normalize(x=x, xRef=xRef, 
          filtersets=filterS, coefNames=coefNames, 
          targetTime=targetTime, targetAmp=targetAmp,
          baseline=baseline, baselineThr=baselineThr, 
          outliers=outliers, outliersThr=outliersThr,
          valign=valign, halign=halign, 
          ctrim=ctrim, spar=spar, plotting=plotRef)
  
  if(!is.null(filterS.log)) {
    pos <- measure(channels(dat)[[ii]]) > 0
    print(c("PosValues=", sum(pos), "NonPosValues=", sum(!pos)))
    measure(channels(dat)[[ii]])[!pos] <- min(measure(channels(dat)[[ii]])[pos])
    dat <- plate_normalize(x=dat, filtersets.log=filterS.log)
  }

  if(!is.null(sFile)) {
     print(c("Cumulative estimates", ii))
     pc <- channels(dat)[[ii]]
     pcx <- measure(pc)
     pcxs <- apply(pcx, 1, sum, na.rm=TRUE)
     des <- cData(dat)
     k <- dim(des)[2]
     tmpDat <- data.frame(Well= rownames(des), des[,3:k], cumulative=round(pcxs,3))
     write.table(file=sFile, x=tmpDat, row.names=FALSE, quote=FALSE, sep="\t")
  }
  
  if(!is.null(yFile)) {
     print(c("Interpolate y at fixed times", ii))
     interp <- plate_interpolate(x=dat, spar=spar, times=time.samples)
     y.interp <- measure(channels(interp)[[ii]])
     yDat <- round(y.interp, 4)
     colnames(yDat) <- as.character(time.samples)
     interpDat <- cbind(cData(dat), yDat)
     write.table(file=yFile, x=interpDat, col.names=TRUE, sep="\t", quote=FALSE)
  }
  
  if(!is.null(lm.intervals) & !is.null(lmFile)) {
     lmRes <- plate_regressions(dat, intervals=lm.intervals, filterset=filterS)
     #print(lmRes)
     #lmFile <- paste(getFileRoot(cFile), "lmRes.txt", sep=".")
     lmDat <-  t(apply(lmRes, 2, cbind))
     k <- dim(lmRes)[1]
     colnames(lmDat) <- rep(names(lmRes[1,1,]), each=k)
     lmDat <- cbind(cData(dat), round(lmDat,4))
     write.table(file=lmFile, x=lmDat, col.names=TRUE, sep="\t", quote=FALSE)
  }
  
  if(!is.null(muFile)){
    print(c("muMax estimates by dydt", ii))
    pc0 <- channels(dat)[[ii]]
    pcx0 <- measure(pc0)
    n <- dim(pcx0)[1]

    interpTime <- targetTime
    if(is.null(interpTime)) ## then take first time of first sample
      interpTime <- time(pc0)[1,1]
    datt <- plate_interpolate(dat, spar=spar, times=interpTime)
    initialAmps <- measure( channels(datt)[[ii]] )[,1]  ## CTW: Fixed bug here (2012-06-20), was [1,]
    
    dydx <- plate_interpolate(dat, spar=spar, times="actual", deriv=1)
    pc1 <- channels(dydx)[[ii]]
    pcx1 <- measure(pc1)
    pct1 <- time(pc1)
    muMax <- apply(pcx1, 1, max)
    muMaxInds <- apply(pcx1 - muMax, 1, function(tx){ order(tx, decreasing=TRUE)[1] })
    
    muMaxTimes <- array(NA, n)
    muMaxAmps  <- array(NA, n)
    for(i in 1:n) {
      muMaxTimes[i] <- pct1[i,muMaxInds[i]]
      muMaxAmps[i]  <- pcx0[i,muMaxInds[i]]
    }
    ## lag-time = muMaxTimes -Y(t=0)/muMax
    ## lag = (initY - muMaxY + muMax*muMaxT)/muMax
    #lagTime = muMaxTimes - initialAmps/muMax
    lagTime = (initialAmps - muMaxAmps + muMax*muMaxTimes)/muMax

    lmDat <- cbind(cData(dat), muMaxT=round(muMaxTimes,2), 
                   muMaxY=round(muMaxAmps,5), muMax=round(muMax,5),
                   lagT=round(lagTime,2))
    
    write.table(file=muFile, x=lmDat, col.names=TRUE, sep="\t", quote=FALSE)
  }
  
  plot_helper(x=dat, plotBy=plotBy, colorBy=colorBy, colorFcn=colorFcn,
              filtersets=filterS, spar=spar, 
              pcex=0.2, xlim=xlim, ylim=ylim,
              show.linearfit=lm.show, lm.intervals=lm.intervals, ...)

  if(!is.null(muFile)){
    plot_helper(x=dydx, plotBy=plotBy, colorBy=colorBy, colorFcn=colorFcn,
              filtersets=filterS, spar=spar, show.spline=FALSE, show.points=TRUE, 
              ptype='l', pcex=0.2, xlim=xlim, ylim=ylim, ylab="dy/dt")
  }
              
  ##plot_channel_error_by_well(dat)
  dev.off()
  
  print(warnings())
  return(dat)
}



















