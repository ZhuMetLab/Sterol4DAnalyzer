OutputValues <- function(args, values, checkname = TRUE, global = FALSE) {
  if (missing(args)) {
    args <- names(values)
  } else {
    args <- as.character(substitute(args)[-1])
  }
  if (checkname & !is.null(names(values))) {
    args0 <- args
    values <- values[args0]
    values <- values[!sapply(values, is.null)]
    args <- names(values)
    idx <- which(is.na(match(args0, args)))
    if (length(idx) > 0) {
      warning(paste("The following arguments are not present in values and will be ignored:\n",
                    paste(args0[idx], collapse = ",")))

    }
  }
  lenArgs <- length(args)
  lenValues <- length(values)
  if (lenArgs != lenValues) {
    warning("Arguments and values are of different length!")
    minlen <- min(c(lenArgs, lenValues))
    idxkeep <- seq(minlen)
    args <- args[idxkeep]
    values <- values[idxkeep]
  }
  mapply(assign, args, values,
         MoreArgs = list(envir = parent.frame()))
}

GetParamTable <- function(params) {
  params <- params[!sapply(params, is.null)]
  paramName <- names(params)
  paramTableList <- lapply(paramName, function(nm) {
    param <- params[[nm]]
    if (inherits(param, 'list')) {
      nmSlots <- slotNames(param[[1]])
      valSlots <- unname(t(sapply(nmSlots, function(nmSlot) {
        sapply(seq_along(param), function(idx) {
          val <- slot(param[[idx]], nmSlot)
          paste(val, collapse = " ")
        })
      })))
    } else {
      nmSlots <- slotNames(param)
      valSlots <- unname(sapply(nmSlots, function(nmSlot) {
        val <- slot(param, nmSlot)
        paste(val, collapse = " ")
      }))
    }
    data.frame("paramName" = nm, "parameter" = nmSlots, valSlots)
  })
  numcol <- max(sapply(paramTableList, ncol))
  res <- do.call(rbind, lapply(paramTableList, function(paramTable) {
    if (ncol(paramTable) < numcol) {
      tbl <- paramTable[, 1:2, drop = FALSE]
      value <- paramTable[, -c(1:2), drop = FALSE]
      valueAppend <- data.frame(matrix("", ncol = numcol - ncol(value) - 2, nrow = nrow(value)))
      value <- cbind(value, valueAppend)
    } else {
      value <- paramTable[, -c(1:2)]
    }
    colnames(value) <- paste("value", seq(ncol(value)), sep = "_")
    cbind(paramTable[, 1:2, drop = FALSE], value)
  }))
  return(res)
}


LoadParam <- function(wd = ".",
                      libTypes = c("standard", "predict"),
                      spectraType = c("msp")) {

  spectraType <- match.arg(spectraType)

  fileParam <- system.file("params", paste0("params_", spectraType), package = getPackageName())
  params <- readRDS(fileParam)
# browser()
  params$matchParam <- lapply(libTypes, function(libType) {
    fileParam <- system.file("params", paste0("matchParam_", libType), package = getPackageName())
    readRDS(fileParam)
  })

  params$searchParam <- lapply(libTypes, function(libType) {
    fileParam <- system.file("params", paste0("searchParam_", libType), package = getPackageName())
    readRDS(fileParam)
  })

  params$experimentParam@wd <- wd
  wd0 <- getwd()
  setwd(wd)
  .CheckDir <- function(dirpath) {
    if (!dir.exists(dirpath)) {
      dir.create(dirpath)
    }
  }

  .CheckDir(params$experimentParam@resDir)
  .CheckDir(params$experimentParam@tmpDir)
  save(params, file = file.path(params$experimentParam@resDir, "paramsOnload.Rda"),
       version = 2)
  write.csv(GetParamTable(params),
            file.path(params$experimentParam@resDir, "paramsOnload.csv"),
            row.names = FALSE)
  for (nm in names(params)) {
    assign(nm, params[[nm]], envir = parent.frame())
  }
  options(mc.cores = params$experimentParam@nSlaves)
  setwd(wd0)
}
