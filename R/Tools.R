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

.ShowSeperateLines <- function(info) {
  sepLine <- paste(rep('-', 40), collapse = '')
  cat(' ', sepLine, paste0('  ', info), sepLine, sep = '\n')
}

.CheckSkip <- function(fn, rerun) {
  # check if use the existed results. if rerun = FALSE, ignore this check
  a <- file.exists(fn) & {!rerun}
  if (a) {
    cat('Using previous results: ', fn, '\n')
  }
  a
}


.MakeUniqueNames <- function(name) {
  duplicated <- TRUE
  i <- 2
  while (duplicated) {
    idxDuplicated <- which(duplicated(name))
    if (length(idxDuplicated) > 0) {
      if (i == 2) {
        name[idxDuplicated] <- paste0(name[idxDuplicated], '_', i)
      } else {
        name[idxDuplicated] <- gsub(paste0('_', i - 1), paste0('_', i),
                                    name[idxDuplicated])
      }
      i <- i + 1
    } else {
      duplicated <- FALSE
    }
  }
  return(name)
}


.GetPpmRange <- function(mz, ppm, resDefineAt = 400) {
  mz + c(-1, 1) * max(prod(mz, ppm, 1e-06), prod(resDefineAt,
                                                 ppm, 1e-06))
}

