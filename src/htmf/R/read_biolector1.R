read_biolector1 <- function(
    file,
    measureType = "AMP" ##c("AMP", "REF", "TEMP", "HUM")
)
{

  header <- list("FILENAME" = NA,
                 "PROTOCOL" = NA,
                 "FILE_VERSION" = NA, #FIXME: watch this one
                 "DATE START" = NA,
                 "DATE END" = NA)

  plate_info <- list("DEVICE" = NA,
                     "USER" = NA,
                     "PLATETYPE" = NA,
                     "MTP ROWS" = NA,
                     "MTP COLUMNS" = NA,
                     "FILTERSETS" = NA)

  ## had to remove the \deg symbol from ACT TEMP [\degC] as it was causing complaints
  measureStrings <- c("REFERENCE", "AMPLITUDE", "ACT TEMP [C]",
                      "ACT HUMIDITY [rH]", "ACT O2 [%]", "ACT CO2 [%]")
  measure_str <- measureStrings[ grep(measureType, measureStrings) ]
  if( length(measure_str) != 1 ) {
  	print(c("ERROR: problem with measureType", measureType, "not in", measureStrings))
    break
  }

  ## print("Parsing 1")
  tmp <- readLines(file, encoding="latin1", warn=FALSE)

  ## print("Parsing 2  (SLOW!)")
  ## Problem here: lines ending in ';' don't include terminal "" when using strsplit :(
  ##lg <- sapply(tmp, function(x){ substr(x, nchar(x), nchar(x)) == ';' } )
  tmp <- sapply(tmp, function(x){ sub(pattern=";$", replacement="; ", x) }) ## Lame.
  rows <- strsplit(tmp, ";")
  rm(tmp)
  ## print("Parsing done.")

  incr_one <- as.integer(1)
  row_i <- as.integer(1)
  date_start <- NA

  if(rows[[3]][1] == "DATE END") ## then old csv
    header <- list("FILENAME" = NA, "DATE START" = NA, "DATE END" = NA)

  for (i in seq(along = header)) {
    row <- rows[[row_i]]

    # bailout at the first blank line
    if (length(row) == 0)
      break

    if(i == 1)
      row[1] <- gsub(" ", "", row[1]) ## to handle version 3.3 change to 'FILE NAME'

    if (names(header)[i] != row[1]) {
      warning(paste("Expected", names(header)[i], "but got", row[1]))
    }
    header[[row[1]]] <- rows[[row_i]][2:length(row)]
    row_i <- row_i + incr_one
  }
  if (! all(is.na(header$`DATE START`))) {
  	dateString <- paste(header$`DATE START`, collapse=" ")
  	if(grepl("[.]", dateString)) ## to handle version 3.3 date strings?
  	  dateString <- strptime(dateString, '%d.%m.%Y %H:%M:%S')
    header$`DATE START` <- as.POSIXct(dateString)
  }
  if (length(rows[[row_i]]) == 0)
    row_i <- row_i + incr_one
  for (i in seq(along = plate_info)) {
    row <- rows[[row_i]]

    # bailout at the first blank line
    if (length(row) == 0)
      break

    if (names(plate_info)[i] != row[1]) {
      warning(paste("Expected", names(plate_info)[i], "but got", row[1]))
    }
    plate_info[[i]] <- rows[[row_i]][2:length(row)]
    row_i <- row_i + incr_one
  }
  {if (is.na(plate_info$`MTP ROWS`) | is.na(plate_info$`MTP COLUMNS`))
    stop("No 'MTP ROWS' or 'MTP COLUMNS'")
  else {
    plate_info$`MTP ROWS` <- as.integer(plate_info$`MTP ROWS`)
    plate_info$`MTP COLUMNS` <- as.integer(plate_info$`MTP COLUMNS`)
  }}
  if (is.na(plate_info$`FILTERSETS`[1]))
    stop("No 'FILTERSETS'")
  else
    numFilterSets <- plate_info$`FILTERSETS`[1] <- as.integer(plate_info$`FILTERSETS`[1]) ## 20120309

  if (length(rows[[row_i]]) == 0)
    row_i <- row_i + incr_one
  row <- rows[[row_i]]
  if (row[1] != "FILTERSET") { ## header line for filterset information
    stop("Expected FILTERSET and got ", row[1])
  }

  fsets_i <- (row_i+1):(row_i + numFilterSets)
  dataf_fsets <- do.call("rbind", args=rows[fsets_i])
  dataf_fsets <- as.data.frame(dataf_fsets, optional = FALSE,
                               stringsAsFactors = FALSE)
  colnames(dataf_fsets) <- rows[[row_i]]
  rownames(dataf_fsets) <- dataf_fsets$`FILTERSET`
  dataf_fsets$`FILTERSET` <- as.integer(dataf_fsets$`FILTERSET`)

  while (row_i <= length(rows)) {
    row <- rows[[row_i]]
    # bailout at the first blank line
    if (length(row) == 0)
      break
    row_i <- row_i + incr_one
  }

  ## Read the bulk of the data
  if (length(rows[[row_i]]) == 0)
    row_i <- row_i + incr_one
  data_cols <- rows[[row_i]]
  row_i <- row_i + incr_one

  bulk_i <- row_i:length(rows)
  fieldOne <- unlist(lapply(rows[bulk_i], "[[", 1))
  acq_C <- grepl("^C", fieldOne)  ## Plate acquisition
  acq_R <- grepl("^R$", fieldOne) ## References acquisition
  acq_K <- grepl("^K$", fieldOne) ## Events record (bad)

  if(sum(acq_C) + sum(acq_R) < 1)
    stop("No acquisition found")

  if(sum(!acq_C & !acq_R & !acq_K) > 0) {
  	unknownTypes <- fieldOne[!acq_C & !acq_R & !acq_K] ##,substr, start=1, stop=1)
    stop("Unknown acquisitions of type:", unique(unknownTypes))
  }
  measuresRet <- acq_C
  if(measure_str == "REFERENCE")
    measuresRet <- acq_R

  if( length(measuresRet) == 0 )
    stop(paste("ERROR no ", measure_str, "found", sep=" "))

  dataf <- do.call("rbind", args=rows[bulk_i][measuresRet])
  dataf <- as.data.frame(dataf, optional = FALSE, stringsAsFactors = FALSE)
  colnames(dataf) <- data_cols[-length(data_cols)]

  {if(measure_str == "REFERENCE") {
    n.fsets <- dim(dataf_fsets)[1]
    k <- dim(dataf)[1]
    dataf$WELLNUM <- factor(array("stage", dim(dataf)[1]))  ## HERE HERE!!
    dataf$READING <- rep(1:ceiling(dim(dataf)[1]/n.fsets), each=n.fsets)[1:k]
  }
  else {
    dataf$WELLNUM <- factor(dataf$WELLNUM)
    dataf$READING <- as.integer(sub("^C|R([0-9]*)$", "\\1", dataf$READING)) ## PROBLEM HERE!!
  }}
  ## maxCycle is last complete cycle over ALL filtersets
  ## don't decrement maxCycle if all cycles are complete
  maxCycle <- max(dataf$READING)
  if( sum(dataf$READING == maxCycle) != sum(dataf$READING == (maxCycle-1)) )
    maxCycle <- maxCycle - 1

  ## Since the first cycle may be >1...
  totalCycles <- maxCycle - min(dataf$READING) + 1
  if(totalCycles < 1)
    stop("No complete acquisition cycles found")

  dataf_l <- split(dataf, dataf$`FILTERSET`)

  ## PlateChannel list
  pc_l <- vector("list", length = length(dataf_l))
  names(pc_l) <- names(dataf_l)

  ## wellNames <- levels(as.factor(dataf$WELLNUM))
  wellNames <- levels(factor(dataf$WELLNUM, levels=unique(dataf$WELLNUM)))
  wellN <- length(wellNames)

  for (fset_i in seq(along=dataf_l)) {
    dataf <- dataf_l[[fset_i]]
    ## trim the last incomplete cycle if necessary
    last_j <- max(which(dataf$READING == maxCycle))
    dataf <- dataf[1:last_j,]

    hours <- as.numeric(dataf$`TIME [h]`)  ##
    #seconds <- as.numeric(hours * 3600)  ## Is there a benefit to seconds??
    ##time <- as.POSIXct(seconds, origin = header$`DATE START`)
    time <- hours

    #o <- 1:wellN ## the natural order IS ordered by time over the cycle...
    ##o <- order(as.integer(dataf$WELLNUM), as.integer(cycle)) ## not working for me
    ##o <- order(wellNames) ## order by well name
    colName <- measure_str
    if(measure_str == "REFERENCE")
       colName <- "AMPLITUDE"
    if(measure_str == "ACT TEMP [C]") ## need index here as degree symbol causes problems
       colName <- 9
    measure <- matrix(as.numeric(dataf[,colName]),
                      nrow = wellN, ncol = totalCycles,
                      dimnames = list(well = wellNames, cycle = NULL)) #[o,]
    #print(rownames(measure))
    time <- matrix(time, nrow = wellN, ncol = totalCycles,
                   dimnames = list(well = wellNames, cycle = NULL)) #[o,]

    fname <- subset(dataf_fsets, `FILTERSET` == fset_i)$`FILTERNAME`  ## Problem here??
    #if(measure_str == "REFERENCE") {
      #print("Not Reordered")
      pc_l[[fset_i]] <- PlateChannel(fname, measure, time)
    #}
    #else {
      ## print("Reordering alternate rows...")
    #  pc_l[[fset_i]] <- PlateChannel(fname, measure[o,], time[o,])
    #}

  }
  if(measure_str == "REFERENCE") {
    initPheno <- data.frame(plate_row="A", plate_col="00", row.names="reference")
  }
  else {
    initPheno <- data.frame(plate_row = sub("^([A-Z])[0-9]+$", "\\1", rownames(measure)),
                            plate_col = sub("^[A-Z]([0-9]+)$", "\\1", rownames(measure)),
                            row.names = rownames(measure))
  }
  adf <- new("AnnotatedDataFrame", data=initPheno)
  pr <- PlateRun(pc_l, adf, info=list(header = header, plate_info = plate_info))

  return(pr)
}
