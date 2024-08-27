##
## Input: PlateChannel
## Value returned: array of y-offsets 
##
getYOffsets = function(
  x,                ## a PlateChannel object
  target = 0,       ## x coordinate (time) of the alignment point
  spar   = 0.8
)
{
  pci <- channel_interpolate(x, times=target[1], spar=spar)
  offsets <- measure(pci)[,1]
  
  return(offsets)
}