##
## Convert a PlateRun to a flat data frame
##
## write.table(file="example_data.txt", x=dataf[inds,], 
##             row.names=FALSE, sep="\t", quote=FALSE)
##
plateRun2df <- function(
   x, 
   filtersets = NULL,
   intervals = NULL,
   cinds = NULL,
   cond  = NULL   ## condition column or name. Need intervals for this
)
{
  datafl <- list()
  chNames <- names(x)
  if(is.null(filtersets))
    filtersets <- chNames

  if(is.null(cinds))
    cinds <- 1:dim(cData(x))[2]

  for(i in 1:length(chNames)) {
  	chName <- chNames[i]
    if(sum(chName == filtersets) > 0) { ## a valid filterset
      ch <- channels(x)[[i]]
      tmp <- melt(time(ch))
      colnames(tmp)[ncol(tmp)] <- "time"
      dataf <- data.frame(tmp, measure=melt(measure(ch))[,3])
      dataf <- merge(dataf, cData(x)[,cinds], by.x = "well", by.y = 0)
      mType <- array(chName, dim(dataf)[1])
      datafl[[i]] <- data.frame(type=mType, dataf)
    }
  }
 
  dataf <- do.call("rbind", args=datafl)
  
  tInterval <- array(1, dim(dataf)[1])
  tCond <- array(1, dim(dataf)[1])
  if(!is.null(intervals)) {
  	k <- length(intervals)
  	tm <- as.numeric(dataf$'time')
    for(ii in 1:k) {
      tInterval[which(tm >= intervals[[ii]][1])] <- ii
    }
    if(!is.null(cond)) {
      design <- cData(x)[,cinds]
      kk <- dim(design)[2]
      Xconds <- t(matrix(unlist(sapply(design[,cond], strsplit, ",")), k, dim(design)[1]))
      rownames(Xconds) <- rownames(design)
      ## need to have the same number of conditions as intervals
      if(k == dim(Xconds)[2]) {
        Xcondsdf <- Xconds[dataf$well,]
        Xcond <- array(NA, dim(Xcondsdf)[1])
        for(ii in 1:k) {
          rinds <- which(tInterval == ii)
          Xcond[rinds] <- Xcondsdf[rinds, ii]
        }
      }
      rownames(Xcond) <- dataf$'well'
      dataf[,cond] <- Xcond 
      ##cnd <- merge(dataf$well, Xcond, by.x = "well", by.y = 0)[,(kk-k):kk]
    }
  }

  dataf <- data.frame(dataf, interval=tInterval)[order(as.numeric(dataf$'time')),]
  return(dataf)
}
