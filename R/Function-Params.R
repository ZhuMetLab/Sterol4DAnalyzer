#' Experiment Setup
#'
#' Parameters for experimental setup
#' @param wd \code{character} Working directory for experiment to be processed
#' @param resDir \code{character} Folder for storing results
#' @param tmpDir \code{character} Folder for storing files when processing. By default,
#'    it is the same with \code{wd}
#' @param polarity \code{character} Ion mode for MS data acquisition
#' @param lc \code{character} LC column used for LC seperation
#' @param instrument \code{character} Instrument used for MS data acquisition
#' @param ms1range \code{numeric(2)} Mass range setup for MS1 acquisition
#' @param ms2range \code{numeric(2)} Mass range setup for MS2 acquisition
#' @param resDefineAt \code{numeric} Mz value for defining resolution
#' @param nSlaves \code{numeric} Number of threads to be used
#' @return an \code{ExperimentParam} object
#' @rdname ExperimentParam
#' @export
ExperimentParam <- function(
  wd = ".",
  resDir = "results",
  tmpDir = file.path(resDir, "tmp"),
  polarity = c("positive", "negative"),
  ce = "30",
  lc = c("RP", "HILIC"),
  instrument = c("Sciex6600", "Agilent6560"),
  ms1range = NULL,
  ms2range = NULL,
  resDefineAt = 400,
  nSlaves = 4) {

  options(mc.cores = nSlaves)

  polarity <- match.arg(polarity)
  lc <- match.arg(lc)
  instrument <- match.arg(instrument)

  if (!dir.exists(wd)) {
    stop(paste0("Working directory setup wrong!\n The file path '",
                wd, "' is not found!"))
  }
  if (!dir.exists(resDir)) {
    dir.create(resDir, recursive = TRUE)
  }
  if (!dir.exists(tmpDir)) {
    dir.create(tmpDir, recursive = TRUE)
  }

  return(new("ExperimentParam",
             wd = wd,
             resDir = resDir,
             tmpDir = tmpDir,
             polarity = polarity,
             ce = ce,
             lc = lc,
             instrument = instrument,
             ms1range = ms1range,
             ms2range = ms2range,
             resDefineAt = resDefineAt,
             nSlaves = nSlaves))
}

#' AssignMSMSParam
#'
#' Parameters for MSMS assignment
#' @param methodAssign \code{character} Method for MSMS asignment
#' \itemize{
#'   \item highest - select the spectrum with highest precursor intensity when
#'    fragmentation in all samples
#'   \item consensus - consensus of all found spectra
#'   \item combined" - consensus of spectra selected from each sample with
#'    highest precursor intensity when fragmentation
#' }
#' @param methodRT \code{character} Raw or corrected RT for MSMS assignment
#' @param rtTol \code{numeric} RT tolerance when matching MSMS spectra and MS1 peaks
#' @param rtExtend \code{numeric} Extended RT  tolerance  when matching MSMS
#'  spectra and MS1 peaks
#' @param mzTol \code{numeric} MZ tolerance when matching MSMS spectra and MS1 peaks
#' \itemize{
#'   \item mzTol < 1 tolerance in Da
#'   \item mzTol >= 1 tolerance in ppm
#' }
#' @param ccsTol \code{numeric} CCS tolerance (in percentage) when matching MSMS
#'  spectra and MS1 peaks
#' @param rtUnitMS1 \code{character} RT unit for user provided MS1 peak table
#' @param rtUnitMS2 \code{character} RT unit for user provided MSMS spectrum files
#' @param isolationWindow \code{numeric} Isolation window when selecting precursor
#'  ions for fragmentation
#' @param thrIsolationIntensity \code{numeric} Intensity threshold of the
#'  contaminat precorsor ions in isolation window to be extracted
#' @param thrIntMS1Predicted \code{numeric} Intensity threshold of precursor ions
#'  for fragmentation to be kept
#' @param snthreshMSMS \code{numeric} S/N threshold for MSMS fragments
#' @param ppmBinMZ \code{numeric} Minimal mz in ppm for fragment binning
#' @param minBinMZ \code{numeric} Minimal mz in Da for fragment binning
#' @param thrRingMz \code{numeric} MZ tolerance in Da for finding satelite fragments
#' @param thrRingIntRel \code{numeric} Relative intensity tolerance for finding
#'  satelite fragments
#' @param minfracVote \code{numeric} Minimal fraction for presentation of fragment
#'  to be kept in consensus spectrum
#' @param outputSpec \code{logical} If output spectra to the spectral MSP file
#' @param plotMSMSEIC \code{logical} If plot EIC with found MSMS specta
#' @param typeMSMSEIC \code{character} Image format (either in PNG or PDF) for plotting
#' @param rerun \code{logical} If rerun all data processing progress

#' @rdname AssignMSMSParam
#' @export
AssignMSMSParam <- function(methodAssign = c("highest"),
                            methodRT = c("raw", "corrected"),
                            rtTol = 0,
                            rtExtend = 90,
                            mzTol = 0.01,
                            ccsTol = 0.5,
                            rtUnitMS1 = c("s", "m"),
                            rtUnitMS2 = c("s", "m"),
                            isolationWindow = 1,
                            thrIsolationIntensity = 150,
                            thrIntMS1Predicted = 500,
                            snthreshMSMS = 10,
                            minfracVote = 0.25,
                            outputSpec = FALSE,
                            plotMSMSEIC = TRUE,
                            typeMSMSEIC = c("pdf", "png", "both"),
                            rerun = FALSE) {
  methodAssign <- match.arg(methodAssign)
  methodRT <- match.arg(methodRT)
  rtUnitMS1 = match.arg(rtUnitMS1)
  rtUnitMS2 = match.arg(rtUnitMS2)
  typeMSMSEIC <- match.arg(typeMSMSEIC)
  new("AssignMSMSParam",
      methodAssign = methodAssign,
      methodRT = methodRT,
      rtTol = rtTol,
      rtExtend = rtExtend,
      mzTol = mzTol,
      ccsTol = ccsTol,
      rtUnitMS1 = rtUnitMS1,
      rtUnitMS2 = rtUnitMS2,
      isolationWindow = isolationWindow,
      thrIsolationIntensity = thrIsolationIntensity,
      thrIntMS1Predicted = thrIntMS1Predicted,
      snthreshMSMS = snthreshMSMS,
      minfracVote = minfracVote,
      outputSpec = outputSpec,
      plotMSMSEIC = plotMSMSEIC,
      typeMSMSEIC = typeMSMSEIC,
      rerun = rerun)
}

## Functions related to the Param class and sub-classes.##
#' @description Extract all slot values and put them into a list, names being
#'     the slot names. If a slot \code{addParams} exist its content will be
#'     appended to the returned list.
#'
#' @param x A Param class.
#'
#' @author Johannes Rainer, modified by YD Yin
#'
#' @noRd
.param2list <- function(x) {
  ## Get all slot names, skip those matching the provided pattern.
  sNames <- slotNames(x)
  skipSome <- grep(sNames, pattern = "^\\.")
  if (length(skipSome) > 0) {
    sNames <- sNames[-skipSome]
  }
  ## handle a slot called "addParams" differently: this is thougth to contain
  ## ... arguments thus we have to skip this one too.
  if (any(sNames == "addParams")) {
    sNames <- sNames[sNames != "addParams"]
    addP <- x@addParams
  } else {
    addP <- list()
  }
  if (length(sNames) > 0) {
    resL <- vector("list", length(sNames))

    for (i in 1:length(sNames)) {
      resL[i] <- list(slot(x, name = sNames[i]))
    }
    names(resL) <- sNames
    resL <- c(resL, addP)
    return(resL)
  }else{
    return(list())
  }
}
