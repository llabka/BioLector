##
## Calls read_df which requires character lists:
##    chanLbls   = NULL 
##    designLbls = NULL
##
## If time label is not "Time" then timeLbl must also be passed to read_df
##    timeLbl = "Time" by default
##
read_df2 <- function(
  file,     ## tab-delimited-text table of plate information
  designfile = NULL,
  wellLbl = "Well",
  ...        ## passed to read_df and then to df2plateRun
)
{
  pr <- read_df(file=file, ...)
  if( !is.null(designfile) ) {
    des <- read.delim(file= designfile, header=TRUE, row.names=1, as.is=TRUE)
    wells <- rownames(time(channels(pr)[[1]]))
    print(wells)
    cData(pr) <- cbind(cData(pr)[,1:2], des[wells,])
  }
  return(pr)
}