read_biolector3 = function(
  file,
  measureType = "AmpRef_1" ##c("Amp_1", "AmpRef_1", REF", "Temp_up", "Humidity", "O2")
)
{
  #measureStrings = c("Amp_1", "Amp_2", "AmpRef_1", "AmpRef_2", "Phase", "Cal",
  #                   "Temp_up", "Temp_down", "Temp_water", "O2", "CO2",
  #                   "Humidity", "Shaker", "Service", "Temp_Ch4")

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

  process_fields = function(x) {
    x = sub(pattern="^[[]", replacement="", x=x)
    r = unlist(strsplit(x, split="]\\s?")[[1]])
    if(length(r)==1) r = c(r, "")
    return(r)
  }

  process_section_heading = function(x) {
    x = sub(pattern="^=+\\s?", replacement="", x=x, perl=TRUE)
    x = sub(pattern="\\s?=+$", replacement="", x=x, perl=TRUE)
    return(x)
  }

  dat = readLines(file, encoding="latin1", warn=FALSE)
  n = length(dat)

  # Extract the following sections
  # "process" (2), "main", "measurement channels" (1-6)
  # "data", "end_process"  (end time only)
  secLgc = grepl("^=", dat)

  sinds = c(which(secLgc), n+1)
  k = length(sinds)  # number of sections
  secLen = sinds[2:k] - sinds[1:(k-1)] - 1 # Section lengths
  slabels = as.character(sapply(dat[secLgc], process_section_heading))
  inds = which(secLgc)

  # ################################
  # Get process, main information
  procList = list()
  for(i in which("process" == slabels)){ # <------------------------
    sdat = dat[(sinds[i]+1):(sinds[i]+secLen[i])]
    skvp = do.call(rbind, lapply(sdat, process_fields))
    procList[[as.character(i)]] = skvp
  }
  procTable = do.call(rbind, procList)
  i = which("main" == slabels)           # <------------------------
  sdat = dat[(sinds[i]+1):(sinds[i]+secLen[i])]
  mainTable = do.call(rbind, lapply(sdat, process_fields))

  # ################################
  # Get channel/filter information
  chanList = list()
  for(i in which("measurement channels" == slabels)){ # <------------------------
    sdat = dat[(sinds[i]+1):(sinds[i]+secLen[i])]
    skvp = do.call(rbind, lapply(sdat, process_fields))
    chanName = skvp[grep("_name", skvp[,1]),2]
    chanList[[chanName]] = skvp
  }
  numFilterSets = length(chanList)

  header2 = as.list(procTable[,2])
  names(header2) = procTable[,1]
  header3 = as.list(mainTable[,2])
  names(header3) = mainTable[,1]

  header$PROTOCOL = header2$process
  header$FILENAME = file
  header$FILE_VERSION = header2$file_version_number
  header$`DATE START` = header2$start_date_time
  i = which("end_process" == slabels)                # <------------------------
  if(length(i) > 0)
    header$"DATE END" = process_fields(dat[(sinds[i]+1)])[2]
  #header$`DATE START` = as.POSIXct(header2$start_date_time)

  # ################################
  # Get channel/filter information
  chanList = list()
  for(i in which("measurement channels" == slabels)){ # <------------------------
    sdat = dat[(sinds[i]+1):(sinds[i]+secLen[i])]
    skvp = do.call(rbind, lapply(sdat, process_fields))
    chanName = skvp[grep("_name", skvp[,1]),2]
    chanList[[chanName]] = skvp
  }
  chanTable = do.call(rbind, chanList)
  lgc1 = grep("_name", chanTable[,1])
  chanLabels = chanTable[lgc1, 2]
  names(chanLabels) = sapply(chanTable[lgc1, 1], sub, pattern="_name", replace="")
  numFilterSets = length(chanList)

  plate_info$DEVICE     = header2$hardware_id
  plate_info$USER       = header3$user
  plate_info$PLATETYPE  = header3$mtp
  plate_info$`MTP ROWS` = header3$`mtp-rows`
  plate_info$`MTP COLUMNS` = header3$`mtp-cols`
  plate_info$FILTERSETS = header3$numFilterSets

  # ################################
  # Get well information
  wellList = list()
  for(i in which("fermentation" == slabels)){ # <------------------------
    sdat = dat[(sinds[i]+1):(sinds[i]+secLen[i])]
    skvp = do.call(rbind, lapply(sdat, process_fields))
    wellList[[as.character(i)]] = skvp
  }
  wellTable = do.call(rbind, wellList)
  wellLbls = wellTable[grep("_well", wellTable[,1]),2]
  wellKeys = wellTable[grep("_well", wellTable[,1]),1]
  names(wellLbls) = sapply(wellKeys, sub, pattern="_well", replacement="")

  # ################################
  # Get data entries
  i = which("data" == slabels)
  datLbls = unlist(strsplit(dat[(sinds[i]+1)], split=";")[[1]])
  m = length(datLbls)
  sdat = dat[(sinds[i]+2):(sinds[i]+secLen[i])]
  lgc = !grepl("^C", sdat) & !grepl("^P", sdat) & !grepl("^R", sdat) # No REFERNCE for now
  if(measureType=="REF")
    lgc = grepl("^R", sdat)
  ddat = sdat[lgc]
  n1 = sum(lgc)

  dtab = matrix("", nrow=n1, ncol=m, dimnames=list(1:n1, datLbls))
  for(i in 1:n1) {
    dtab[i,] = (strsplit(ddat[i], split=";")[[1]])
  }

  # Trim any incomplete cycles
  cycleIndex = as.numeric(dtab[,"Cycle"])
  maxCycle = max(cycleIndex)
  if( sum(cycleIndex == maxCycle) != sum(cycleIndex == (maxCycle-1)) ) {
    maxCycle = maxCycle - 1
    last_j = max(which(cycleIndex == maxCycle))
    dtab = dtab[1:last_j,]
  }
  n2 = dim(dtab)[1]
  filtIndList = split(1:n2, dtab[, "Filterset"])
  #names(filtIndList) = chanLabels[names(filtIndList)]

  colName = measureType
  if(measureType == "REF")
    colName = "Amp_1"

  cycles = as.numeric(unique(dtab[,"Cycle"]))
  totalCycles = diff(range(cycles)) + 1

  pc_l <- vector("list", length = numFilterSets)
  names(pc_l) = names(filtIndList)
  for(fset_i in seq(along=filtIndList)) {
    fname = as.character(chanLabels[names(filtIndList)[fset_i]])
    dtab_i = dtab[filtIndList[[fset_i]], ]

    wellN = length(unique(dtab_i[,"Well"]))
    wellNames = wellLbls[unique(dtab_i[,"Well"])]
    hours = as.numeric(dtab_i[,"Time"]) / 3600

    measure = matrix(as.numeric(dtab_i[,colName]),
                     nrow = wellN, ncol = totalCycles,
                     dimnames = list(well = wellNames, cycle = NULL))
    time = matrix(hours, nrow = wellN, ncol = totalCycles,
                  dimnames = list(well = wellNames, cycle = NULL))
    pc_l[[fset_i]] <- PlateChannel(fname, measure, time)
  }
  if(measureType == "REF") {
    initPheno <- data.frame(plate_row="A", plate_col="00", row.names="reference")
  }
  else {
    initPheno <- data.frame(plate_row = sub("^([A-Z])[0-9]+$", "\\1", wellLbls),
                            plate_col = sub("^[A-Z]([0-9]+)$", "\\1", wellLbls),
                            row.names = wellLbls)
  }
  adf <- new("AnnotatedDataFrame", data=initPheno)
  pr <- PlateRun(pc_l, adf, info=list(header = header, plate_info = plate_info))

  return(pr)
}
