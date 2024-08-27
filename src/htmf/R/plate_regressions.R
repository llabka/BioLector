plate_regressions <- function(
  x,  ## PlateRun
  intervals = NULL,
  filterset = "Biomass"
)
{
  i <- which(filterset == names(x))[1] ## take the first channel that matches
  pc <- channels(x)[[i]]
  return( channel_regressions(pc, intervals) )
}