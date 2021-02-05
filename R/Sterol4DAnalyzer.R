SterolAnalyzer <- function(filePeaks,
                           spectraType = c("msp"),
                           libdata = NULL,
                           adjustAccuracy = TRUE,
                           fileRtCalibrate = NULL,
                           fileRtRef = NULL,
                           experimentParam = NULL,
                           parseParam = NULL,
                           assignMSMSParam = NULL,
                           searchParam = NULL,
                           matchParam = NULL,
                           combineParam = NULL) {

  wd0 <- getwd()
  setwd(experimentParam@wd)
  spectraType <- match.arg(spectraType)
  params <- list("experimentParam" = experimentParam,
                 "assignMSMSParam" = assignMSMSParam,
                 "parseParam" = parseParam,
                 "searchParam" = searchParam,
                 "matchParam" = matchParam,
                 "combineParam" = combineParam
  )
  save(params, file = file.path(experimentParam@resDir, "params.Rda"), version = 2)
  write.csv(GetParamTable(params),
            file.path(experimentParam@resDir, "params.csv"),
            row.names = FALSE)

  rtcalExp <- rtcalRef <- NULL
  if (!is.null(fileRtCalibrate)) {
    rtcalExp <- read.csv(fileRtCalibrate, stringsAsFactors = FALSE)
    if (assignMSMSParam@rtUnitMS1 == "m") {
      rtcalExp$rt <- round(rtcalExp$rt * 60, 0)
    }
    if (!is.null(fileRtRef)) {
      rtcalRef <- read.csv(fileRtRef, stringsAsFactors = FALSE)
    } else {
      stop('Reference RT calibration table must be provided!')
    }
  }

  metInfo <- read.csv(filePeaks, stringsAsFactors = FALSE)
  if (adjustAccuracy) {
    clCCS <- match("ccs", tolower(colnames(metInfo)))
    if (length(clCCS)) {
      metInfo[, clCCS] <- round(metInfo[, clCCS], 1)
    }
  }
  filesSpec = list.files(experimentParam@wd,
                         pattern = paste0("(?i)\\.", spectraType, "$"),
                         recursive = TRUE,
                         full.names = TRUE)
  expdata <- AssignMSMS(assignMSMSParam, metInfo, experimentParam,
                        parseParam, files = filesSpec)

  .ShowSeperateLines("Identifying lipids ...")

  cat("    Searching librarial candidates ...\n")
  expdata <- SpectraTools::FilterNULL(expdata, keepInfo = TRUE)

  if (class(libdata) != 'list') {
    libdataAll <- list(libdata)
  } else {
    libdataAll <- libdata
  }

  specSearched <- lapply(seq_along(libdataAll), function(idx) {
    SpectraTools::SearchSpectra(dataExp = expdata, dataRef = libdataAll[[idx]],
                                searchParam = searchParam[[idx]],
                                rtcalExp = rtcalExp, rtcalRef = rtcalRef)
  })

  cat("    Matching MSMS ...\n")
  scoreMatchAll <- lapply(seq_along(libdataAll), function(idx) {
    scoreMatch <- BiocParallel::bplapply(specSearched[[idx]], function(specData) {
      dataExp <- specData$dataExp
      dataRef <- specData$dataRef
      nameRef <- dataRef@info$name
      nameRep <- unique(nameRef)
      idxKeep <- match(nameRep, nameRef)
      dataRef <- SpectraTools::SpectraData(info = dataRef@info[idxKeep, , drop = FALSE],
                                           spectra = dataRef@spectra[idxKeep])
      dataRef@info$name <- nameRep
      SpectraTools::MatchSpectra(dataExp, dataRef, matchParam[[idx]])
    })
    scoreMatch <- scoreMatch[!sapply(scoreMatch, is.null)]
  })

  cat("    Plotting MSMS match figures ...\n")
  for (idx in seq_along(libdataAll)) {
    dirPlot <- file.path(experimentParam@resDir, "MSMSMatchPlot", libnames[idx])
    PlotMatchResult(scoreMatchAll[[idx]], expdata@info, dirPlot, addname = FALSE)
  }

  scTables <- lapply(seq_along(libdataAll), function(idx) {
    scTable <- do.call(rbind, lapply(scoreMatchAll[[idx]], SpectraTools::GenOutputScore,
                                     matchParam[[idx]]@cutoff, type = "sterol"))
    pkTable <- MergeResTable(expdata@info, scTable)
    write.csv(pkTable, file.path(experimentParam@resDir,
                                 paste0("result_MSMSmatch_", libnames[idx], ".csv")),
              row.names = FALSE)
    return(scTable)
  })

  cat("    Finalizing scores ...\n")
  scoreMatchAll <- lapply(seq_along(libdataAll), function(idx) {
    scoreMatch <- BiocParallel::bplapply(scoreMatchAll[[idx]], function(sc) {
      SpectraTools::CombineScore(sc, combineParam)
    })
    scoreMatch <- scoreMatch[!sapply(scoreMatch, is.null)]
  })

  scTables <- lapply(seq_along(libdataAll), function(idx) {
    scTable <- do.call(rbind, lapply(scoreMatchAll[[idx]], SpectraTools::GenOutputScore,
                                     type = "sterol"))
    pkTable <- MergeResTable(expdata@info, scTable)
    write.csv(pkTable, file.path(experimentParam@resDir,
                                 paste0("result_ScoreCombine_", libnames[idx], ".csv")),
              row.names = FALSE)
    return(scTable)
  })

  setwd(wd0)
  cat("    Congratulations! ALL WORK DONE!!\n")
}

PlotMatchResult <- function(scoreMatch, expinfo, dirPlot,
                            addname = FALSE, plotPNG = FALSE) {
  if (!dir.exists(dirPlot)) {
    dir.create(dirPlot, recursive = TRUE)
  }
  BiocParallel::bplapply(names(scoreMatch), function(nm) {
    require(ggplot2)
    pkname <- expinfo[nm, "name"]
    matchScore <- scoreMatch[[nm]]
    if (!is.null(matchScore)) {
      filePlot = file.path(dirPlot, pkname)
      SpectraTools::PlotMirror(matchScore, pkname, addname = addname,
                               plotPNG = plotPNG,
                               plotPDF = TRUE,
                               filePlot = filePlot,
                               direction = 'both')
    }
  })
}


OutputMatchResult <- function(scoreMatch, expinfo, dirSave) {
  if (!dir.exists(dirSave)) {
    dir.create(dirSave)
  }

  lapply(names(scoreMatch), function(nm) {
    sc <- scoreMatch[[nm]]
    info <- sc@info
    frag <- lapply(sc@matchedFragments, function(frag) {
      frag[!is.na(frag$mz), , drop = FALSE]
    })
    dtCheck <- do.call(rbind, frag)
    dtCheck$lipname <- unlist(mapply(rep,
                                     times = sapply(frag, nrow),
                                     info$name,
                                     SIMPLIFY = FALSE))
    fnSave <- file.path(dirSave, paste0(expinfo[nm, "name"], ".csv"))
    write.csv(dtCheck, fnSave, row.names = FALSE)
  })
}

MergeResTable <- function(info, res) {
  if (!is.data.frame(res)) {
    res <- SpectraTools:::.Col2Numeric(res)
  }
  tmp <- data.frame(matrix(ncol = ncol(res), nrow = nrow(info)))
  colnames(tmp) <- colnames(res)
  rownames(tmp) <- rownames(info)
  tmp[, sapply(res, is.character)] <- ""
  tmp[, sapply(res, is.numeric)] <- 0
  tmp[rownames(res), ] <- res
  info <- cbind(info, tmp)
  return(info)
}

