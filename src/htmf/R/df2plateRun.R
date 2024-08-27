df2plateRun <- function(
  df,
  wells,   ## which wells to extract
  chanLbls,   ## column names of df for the measured channels
  designLbls, ## column names of df that populate the AnnotatedDataFrame
  timeLbl = "Time", 
  wellLbl = "Well", 
  na.value=NA  ## should be removed
  ) 
{

  ## Should find the information to populate these lists...
  header <- list("FILENAME" = NA,
                 "PROTOCOL" = NA,
                 "FILE_VERSION" = NA,
                 "DATE START" = NA,
                 "DATE END" = NA)

  plate_info <- list("DEVICE" = NA,
                     "USER" = NA,
                     "PLATETYPE" = NA,
                     "MTP ROWS" = NA,
                     "MTP COLUMNS" = NA,
                     "FILTERSETS" = NA)

  ## PlateChannel list
  pcl <- vector("list", length = length(chanLbls))
  names(pcl) <- chanLbls

  initPheno <- data.frame(plate_row = sub("^([A-Z])[0-9]+$", "\\1", wells),
      plate_col = sub("^[A-Z]([0-9]+)$", "\\1", wells),
      row.names = wells)
      
  if(length(designLbls)>0) {
    dfFirstInst <- df[match(wells, df[, wellLbl]),]
    designDat <- subset(dfFirstInst, select=designLbls)
    for(ii in 1:dim(designDat)[2]) ## Seems that all cData must be text (not factors or numeric)?
      designDat[,ii] <- as.character(designDat[,ii])
    initPheno <- data.frame(initPheno, designDat)
  }
  adf <- new("AnnotatedDataFrame", data=initPheno)                

  times <- extractMatrix(df, chanLbl=timeLbl, wells=wells, wellLbl=wellLbl)
  for(chnl in chanLbls) {
  	measure <- extractMatrix(df, chanLbl=chnl, wells=wells, wellLbl=wellLbl)
  	if(!is.na(na.value)) ## hack to handle NA's
  	  measure[is.na(measure)] <- na.value
    pcl[[chnl]] <- PlateChannel(chnl, measure, times)
  }
  pr <- PlateRun(pcl, adf, info=list(header = header, plate_info = plate_info))
  return(pr)
}
################################################################################################
extractMatrix <- function(df, wells=NULL, chanLbl="Time", wellLbl="Well") {
  N <- dim(df)[1]
  allWells <- as.character(df[,wellLbl])
  wellTbl <- table(allWells)
  n <- length(wellTbl)
  wellIndexLsts <- split(1:N, allWells)
  
  if( sum(wellTbl[1:(n-1)] - wellTbl[2:n]) != 0) { 
	print("WARNING: extractMatrix finds differnet cycle totals in differnt wells")
  }		
  m <- min(as.numeric(wellTbl))
  x <- matrix(NA, n, m, dimnames=list(well=wells, cycle=NULL))
  for(well in wells) {
    inds <- wellIndexLsts[[well]]
    if(m < length(inds)) inds <- inds[1:m]
    x[well,] <- df[inds, chanLbl]  
  }
  return(x)
}
