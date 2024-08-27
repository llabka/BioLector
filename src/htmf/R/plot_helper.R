##
## Dependencies:
## buildColors
## plate_intervals
## plot_biolector
##
## legend.loc  = "topleft", ## option removed

plot_helper <- function(
   x,                   ## a PlateRun
   intervals   = NULL,
   plotBy      = 0,     ## cData column name or index OR -1 or "Wells"
   colorBy     = 1,     ## cData column name or index
   colorFcn    = rainbow,
   filtersets  = NULL,  ## if NULL, all filtersets/channels will be plotted
   spar        = 0.8,
   plotName    = "",
   ...)
{
  if(is.null(filtersets))
    filtersets <- names(x)
  #print(filtersets)

  colorFactor <- cData(x)[,colorBy]
  #print(colnames(cData(x)))
  colorName <- colnames(cData(x))[colorBy]

  colorList <- buildColors(colorFactor, colorFcn=colorFcn)
  plotFactor <- as.factor(rownames(cData(x))) ## the default, by well name

  if(plotBy>0 & plotBy!="Wells") {
    plotFactor <- as.factor(cData(x)[,plotBy])
    plotName  <- paste(plotName, colnames(cData(x))[plotBy], sep=" ")
  }

  if(plotBy == -1 | plotBy=="Wells") {
    plotFactor <- as.factor(array("plate", dim(cData(x))[1]))
    plotName  <- paste(plotName, "Wells", sep=" ")
  }
  replicates <- split(1:dim(cData(x))[1], plotFactor)

  if(is.null(intervals))
    intervals  <- plate_intervals(x, thresh=0.3)

  for(i in 1:length(names(x))) { ## over channels
    filterName <- names(x)[i]
    if(sum(filterName == filtersets) > 0) {
      plot_biolector(x, filterset=filterName, replicates=replicates,
                     intervals=intervals, colorList=colorList, spar=spar,
                     legend.title=colorName, main.title=plotName, ...)
    }
  }
}
