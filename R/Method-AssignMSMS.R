#' @param assignMSMSParam \code{AssignMSMSParam} Parameters for MSMS asignment
#' @param metInfo \code{data.frame} User provided peak table
#' @param experimentParam \code{ExperimentParam} Parameters for experimental setup
#' @param files \code{character} MSMS spectrum files
#' @rdname Method-AssignMSMS
#' @export
setMethod(
  "AssignMSMS",
  signature = c("AssignMSMSParam", "data.frame", "ExperimentParam", "ParseSpectraParam"),
  function(assignMSMSParam, metInfo, experimentParam, parseParam,
           files = NULL, ...) {
    .ShowSeperateLines("Assigning MSMS spectra to features ...")

    if (assignMSMSParam@rtUnitMS1 == "m") {
      metInfo$rt <- round(metInfo$rt * 60, 0)
    }
    ccsFound <- any(grepl("^ccs$", colnames(metInfo)))
    if (!any(grepl("^name$", colnames(metInfo)))) {
      if (ccsFound) {
        ftnames <- apply(round(metInfo[, c("mz", "rt", "ccs")], 0), 1, function(dr) {
          paste0(c("M", "T", "C"), dr, collapse = "")
        })
      } else {
        ftnames <- apply(round(metInfo[, c("mz", "rt")], 0), 1, function(dr) {
          paste0(c("M", "T"), dr, collapse = "")
        })
      }
      ftnames <- .MakeUniqueNames(ftnames)
      metInfo <- data.frame("name" = ftnames, metInfo, stringsAsFactors = FALSE)
    }

    fileType <- tolower(unique(tools::file_ext(files)))
    if (length(fileType) > 1) {
      stop("Spectra files are not of the same types!")
    }

    fnSave <- file.path(experimentParam@tmpDir,
                        paste0("ms2-assigned-",assignMSMSParam@methodAssign, ".rdata"))
    if (.CheckSkip(fnSave, assignMSMSParam@rerun)) {
      load(fnSave)
    } else {
      if (missing(parseParam)) {
        parseParam <- switch(
          fileType,
          "msp" = {
            if (ccsFound) {
              SpectraTools::ParseSpectraParam(type = "msp",
                                              labelKeep = c("NAME",
                                                            "PRECURSORMZ",
                                                            "RETENTIONTIME",
                                                            "CCS"),
                                              labelName = c("name",
                                                            "mz",
                                                            "rt",
                                                            "ccs"),
                                              denoise = TRUE,
                                              ppmPrecursorFilter = 20)
            } else {
              SpectraTools::ParseSpectraParam(type = "msp",
                                              labelKeep = c("NAME",
                                                            "PRECURSORMZ",
                                                            "RETENTIONTIME"),
                                              labelName = c("name",
                                                            "mz",
                                                            "rt"),
                                              denoise = TRUE,
                                              ppmPrecursorFilter = 20)
            }
          }
        )

        if (!is.null(experimentParam@ms2range)) {
          SpectraTools::ms2range(parseParam) <- experimentParam@ms2range
        }
      }

      ms2spec <- BiocParallel::bplapply(files, function(file, parseParam) {
        SpectraTools::ParseSpectra(parseParam, file)
      }, parseParam = parseParam)

      if (assignMSMSParam@rtUnitMS2 == "m") {
        ms2spec <- lapply(ms2spec, function(specdata) {
          specdata@info$rt <- specdata@info$rt * 60
          specdata
        })
      }

      # ms2info <- lapply(ms2spec, function(specdata) specdata@info)
      if (fileType == "msp" && length(files) == 1) {
        rownames(metInfo) <- metInfo$name
        ms2spec <- lapply(ms2spec, function(spec) {
          rownames(spec@info) <- names(spec@spectra) <- spec@info$name
          spec
        })
        ms2spec <- ms2spec[[1]]
        ms2assigned <- lapply(metInfo$name, function(nm) {
          SpectraTools::SpectraData(info = ms2spec@info[nm, , drop = FALSE],
                                    spectra = ms2spec@spectra[nm])
        })
      } else {
        ms2all <- lapply(seq(nrow(metInfo)), function(nr) {
          # cat(nr, "\t")
          pk <- metInfo[nr, ]
          mzTol <- assignMSMSParam@mzTol
          if (mzTol <= 1) {
            mzrange <- pk$mz + c(-1, 1) * mzTol
          } else {
            mzrange <- .GetPpmRange(pk$mz,
                                    ppm = assignMSMSParam@mzTol,
                                    resDefineAt = experimentParam@resDefineAt)
          }
          rtrange <- pk$rt + c(-1, 1) * assignMSMSParam@rtTol
          if (ccsFound) {
            ccsrange <- pk$ccs + c(-1, 1) * assignMSMSParam@ccsTol * pk$ccs/100
          }

          specmatched <- unlist(sapply(ms2spec, function(specdata) {
            mz <- specdata@info$mz
            rt <- specdata@info$rt

            mzmatched <- mz >= mzrange[1] & mz <= mzrange[2]
            rtmatched <- rt >= rtrange[1] & rt <= rtrange[2]
            ccsmatched <- rep(TRUE, length(mz))
            if (ccsFound && !is.na(match("ccs", colnames(specdata@info)))) {
              ccs <- specdata@info$ccs
              ccsmatched <- ccs >= ccsrange[1] & ccs <= ccsrange[2]
            }
            idx <- which(mzmatched & rtmatched & ccsmatched)

            if (length(idx) > 0) {
              ms1int <- sapply(specdata@spectra[idx], function(spec) {
                sum(spec[, 2])
              })
              info <- cbind(specdata@info[idx, , drop = FALSE], "ms1int" = ms1int)
              SpectraTools::SpectraData(info = info,
                                        spectra = specdata@spectra[idx])
            }
          }))
        })

        assignParamList <- as.list(assignMSMSParam)

        args <- lapply(ms2all, function(specdata) {
          c("specdata" = list(specdata), assignParamList)
        })
        ms2assigned <- BiocParallel::bplapply(args, function(arg) {
          do.call(paste(".AssignMSMS", assignMSMSParam@methodAssign, sep = "."),
                  arg)
        })
      }

      save(ms2assigned, file = fnSave)
    }

    spectra <- lapply(ms2assigned, function(specdata) {
      if (is.null(specdata)) {
        return(NULL)
      }
      specdata@spectra[[1]]
    })

    object <- SpectraTools::SpectraData(info = metInfo,
                                        spectra = spectra)

    if (assignMSMSParam@outputSpec) {
      fnSave <- file.path(experimentParam@resDir, "spectra.msp")
      sink(fnSave)
      for (idx in seq(nrow(object@info))) {
        spec <- object@spectra[[idx]]
        if (!is.null(spec)) {
          ftname <- object@info[idx, "name"]
          mz <- object@info[idx, "mz"]
          rt <- object@info[idx, "rt"]

          cat("NAME: ", ftname, "\n", sep = "")
          cat("PRECURSORMZ: ", mz, "\n", sep = "")
          cat("RETENTIONTIME: ", rt, "\n", sep = "")
          if (ccsFound) {
            ccs <- object@info[idx, "ccs"]
            cat("CCS: ", ccs, "\n", sep = "")
          }
          cat("Num Peaks: ", nrow(spec), "\n", sep = "")
          for (nr in seq(nrow(spec))) {
            cat(paste(spec[nr, ], collapse = " "), "\n", sep = "")
          }
          cat("\n")
        }
      }
      sink()
    }
    write.csv(object@info,
              file.path(experimentParam@resDir, "PeakTable-MS2Assigned.csv"),
              row.names = FALSE)
    cat("Assigning work done!")
    return(object)
  })
