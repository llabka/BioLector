plate_calibrate <- function(
  x,  ## a plateRun
  filterset = "Biomass", 
  replicates = NULL,  ## which wells will be used for each fit
  times = NULL,       ## times that correspond to the columns of calibDat
  calibDat = NULL,    ## matrix/df of ODs or other calibration data
                      ## with well-labels as rownames
  fctName = "replicates",
  plotting = FALSE,
  ...
)
{	
  pFcn <- function(cf){ 
    paste("Rsq= ", round(cf[3],5), " ", round(cf[2],3),
          "x + ", round(cf[1],3), sep="")
  }

  ##cdat <- cData(x)
  fsi <- which(filterset == names(x))[1] ## take the first channel with matching name
  ## LSdat <- plate_interpolate(x=x, filtersets=filterset, spar=0.4, times=times)
  ## LS <- measure(channels(LSdat)[[1]])
  LSdat <- plate_interpolate(x=x, spar=0.4, times=times)
  LS <- measure(channels(LSdat)[[fsi]])
  wells <- rownames(LS)

  if(is.null(replicates))
    replicates <- split(1:length(wells), wells)
 
  LSi <- match(wells, rownames(LS))
  Dsi <- match(wells, rownames(calibDat))

  Y <- as.matrix(calibDat[Dsi,])
  X <- as.matrix(LS[LSi,])

  m <- dim(LS)[1]
  wellCoef <- matrix(NA, m, 2)
  colnames(wellCoef) <- c("offset", "slope")
  rownames(wellCoef) <- rownames(LS)

  ####################
  k <- length(replicates)
  lmCoef <- matrix(NA, k, 3)

  if(plotting)
    plot(NA, type='n', xlim=range(X), ylim=range(Y), ...)

  for( j in 1:k ) {
  	s <- names(replicates)[j]
  	jj <- replicates[[j]]
    sx <- as.numeric(X[jj,])
    sy <- as.numeric(Y[jj,])

    lmRes <- lm(y~x, data.frame(x=sx, y=sy))
    if(plotting) {
      points(sx, sy, pch=19, cex=0.5, col=rainbow(k)[j])
      lines(sx, predict(lmRes), col=rainbow(k)[j])
    }
    lmCoef[j, 1:2] <- lmRes$coefficients
    lmCoef[j,   3] <- summary(lmRes)$r.squared
    wellCoef[jj,] <- t(matrix(lmRes$coefficients, 2, length(jj)))
  }
  if(plotting) {
    linFits <- apply(lmCoef, 1, pFcn)
    legend("bottomright", fill=rainbow(k), 
           legend=names(replicates), title=fctName)
    legend("topleft", fill=rainbow(k), legend=linFits, title="linear fit")
  }
  return(wellCoef)
}
