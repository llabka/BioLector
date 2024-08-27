channel_valign <- function(
  x,                 ## a Plate Channel Object
  Talign     = 0,    ## time points to align
  Yalign     = NULL, ## amplitude to align to (at Talign)
  spar       = 0.8,
  interpolate = TRUE,
  LOS = FALSE,       ## apply log-of-slope method to find background
  xMin = 3,
  xMax = 10,
  yMin = 1e-2,
  fdThr = 1e-2,
  sdThr = 0,
  lspar = 0.4,
  plotting=FALSE
)
{
  yPoint = Yalign
  Y = measure(x)
  pcDim = dim(Y)
  tMeans <- apply(time(x), 2, mean)
  tMinizer <- order(abs(tMeans - Talign))[1]
  Yoffsets <- measure(x)[, tMinizer]
  if(interpolate) {
    pc_offset <- channel_interpolate(x=x, times=Talign, spar=spar)
    Yoffsets <- measure(pc_offset)[,1]
  }
  if(is.null(Yalign))
    yPoint <- median(Yoffsets)
  dYs <- as.numeric(yPoint) - as.numeric(Yoffsets)

  if(LOS)  ## There are many thresholds/params in this call to estimate_bg()
    dYs = estimate_bg(pc=x, spar=spar, xMin=xMin, xMax=xMax, yMin=yMin,
                      fdThr=fdThr, sdThr=sdThr, lspar=lspar, plotting=plotting)

  Yadj <- Y + as.numeric(dYs)
  return( PlateChannel(name = x@name,
                           measure_matrix = Yadj,
                           time_matrix = time(x))
  )
}

estimate_bg = function(
  pc,         ## a plate channel
  spar,       ## smoothing parameter for fitting raw measurements
  xMin,       ## min time threshold for fitting range
  xMax,       ## max time threshold for fitting range
  fdThr,      ## first derivivative threshold for fitting range
  sdThr,      ## second derivivative threshold for fitting range
  lspar,      ## smoothing parameter for fitting log measurements (used bu SOL)
  yMin,       ## minimum y allowed when taking the log (used bu SOL)
  plotting=FALSE)
{
  SOL = function(x, y) {
    y[y < yMin] = yMin
    spl = smooth.spline(x=x, y=log(y), spar=lspar)
    mu = predict(spl, x, deriv=1)$y  # SOL
    return(mu)
  }

  CostFcn <- function(P) {
    mu = SOL(x=obs$time, y=obs$y + P[1])
    out = data.frame(time=obs$time, mu=mu)
    colnames(out) <- c("time", "mu")
    return(modCost(out, obs[,c("time", "mu")]))
  }

  X = time(pc)
  Y = measure(pc)
  wells = rownames(X)
  bgOpt = array(NA, length(wells))
  names(bgOpt) = wells
  if(plotting) par(mfrow=c(2,2))

  for(well in wells) {
    yi = Y[well,]
    xi = X[well,]
    sspl = smooth.spline(x=xi, y=yi, spar=spar)
    yt = predict(sspl, xi)$y
    fder = predict(sspl, xi, deriv=1)$y
    sder = predict(sspl, xi, deriv=2)$y
    l = fdThr < fder & sdThr < sder & xMin <= xi & xi <= xMax

    splF = splinefun(xi[l], yt[l], method="fmm")
    fderF = splF(xi[l], deriv=1)
    sderF = splF(xi[l], deriv=2)
    LOS = sderF/fderF  ## log-of-slope method
    obs = data.frame(time=xi[l], y=yi[l], mu=LOS)
    fitBg <- modFit(f=CostFcn, p = c(bg=0), lower = c(bg=-100), upper = c(bg=100))
    bgOpt[well] = fitBg$par

    if(plotting) {
      plot(xi, yi, type='l', log="y", ylim=c(yMin*0.5, max(yi)), main=well)
      lines(xi, yi+bgOpt[well], col="red")
      points(xi[l], yi[l], col="black")
      points(xi[l], yi[l]+bgOpt[well], col="red")
    }
  }
  return(bgOpt)
}
