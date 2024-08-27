##
## Calls df2plateRun which requires character lists:
##      chanLbls and designLbls 
## If time label is not "Time" then timeLbl must also be passed to df2plateRun
##
read_df <- function(
  file,     ## tab-delimited-text table of plate information 
  wellLbl = "Well",
  ...        ## passed to df2plateRun
)
{
  df <- read.delim(file=file, header=TRUE, as.is=TRUE)
	
  ## extract all wells for now...
  allWells <- as.character(df[,wellLbl])
  wellTbl <- table(allWells)
  wellLbls <- sort(names(wellTbl))
    
  pr <- df2plateRun(df=df, wells=wellLbls, wellLbl=wellLbl, ...)
  return(pr)
}