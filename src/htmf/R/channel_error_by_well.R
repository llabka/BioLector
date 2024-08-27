channel_error_by_well <- function(x, fcn=mad) ## IQR, var, sd, min, max, etc... 
{	
   smx <- channel_interpolate(x)
   dA <- measure(x) - measure(smx)
   ## diff(apply(dA, 1, summary)[c(2,5),]) ## inter-quartile range...
   return( apply(dA, 1, fcn) )
}

atoi <- function(x) {
  alphabet <- "ABCDEFGHI"
  integers <- "123456789"
  as.numeric(chartr(alphabet, integers, x))
}

plot_channel_error_by_well <- function(x, filterset = "Biomass", ...) 
{
  chInd <- which(filterset == names(x))[1] ## take the first matching channel name
  ch <- channels(x)[[chInd]]
  error <- channel_error_by_well(ch)
  rDate <- info(x)$header$`DATE START`[1]
  maxr <- info(x)$plate_info$`MTP ROWS`[1]
  maxc <- info(x)$plate_info$`MTP COLUMNS`[1]
  rowi <- atoi(as.character(cData(x)$plate_row))
  coli <- as.numeric(cData(x)$plate_col)

  ##par(bg="white", omd=c(0, 1, 0, 1), mar=c(5,6,6,6)) ##usr=c(0, max(rowi)+1, 0, max(coli)+1)
  plot(rowi, coli, cex=error, pch=19, axes=FALSE, ann=TRUE, 
       sub=rDate, xlab="", ylab="", cex.main=0.8, cex.sub=0.8,
       xlim=c(0.2, maxr+0.2), ylim=c(0.2, maxc+0.2), col="red", asp=maxc/maxr, ...)
  rlbls <- as.character(levels(cData(x)$plate_row))
  clbls <- as.character(levels(cData(x)$plate_col))
  axis(3, at=unique(rowi), labels=rlbls, tick=FALSE, cex.axis=0.5)
  axis(2, at=unique(coli), labels=clbls, tick=FALSE, cex.axis=0.5, las=2)
}
