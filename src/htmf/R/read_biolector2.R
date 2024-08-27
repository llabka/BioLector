read_biolector2 <- function(
    file,
    designfile  = NULL,
    measureType = c("AMP", "REF", "TEMP", "HUM")[1],
    blVersion = 2
)
{
  dat <- read_biolector(file, measureType=measureType, blVersion=blVersion)

  if(!is.null(designfile)) {
    design <- read.table(designfile, row.names=1, header=TRUE, sep="\t", as.is=TRUE)
    mInds <- match(rownames(cData(dat)), rownames(design))

    if(sum(!is.na(mInds)) != length(rownames(cData(dat))))
      break

    df <- cbind(cData(dat), design[mInds,])
    colnames(df) <- c(colnames(cData(dat)), colnames(design))
    dat@celldata <- new("AnnotatedDataFrame", data=df)
  }
  return(dat)
}