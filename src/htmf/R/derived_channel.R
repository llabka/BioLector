derived_channel <- function(
  x,  ## plateRun object
  c1 = NULL, ## first channel/filter name
  c2 = NULL, ## second channel/filter name
  cnew = "new", ## name of the new channel
  binOp="-"
)
{
  fnames <- names(x)
  pc1 <- channels(pr)[[which(fnames==c1)]] ## should have error checking here
  pc2 <- channels(pr)[[which(fnames==c2)]] ## should have error checking here
  
  tnew <- (time(pc1) + time(pc2))/2   ## just average the times between the 2 channels 
  mnew <- do.call( binOp, args=list(measure(pc1), measure(pc2)) )
  channels(x)[[cnew]] <- PlateChannel(cnew, mnew, tnew)
  return(x)
}