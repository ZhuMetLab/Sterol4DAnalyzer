setClass("ParamSterol", representation = "VIRTUAL")

setClassUnion("nullOrCharacter", c("NULL", "character"))
setClassUnion("nullOrNumeric", c("NULL", "numeric"))

#' @export
setClass("ExperimentParam",
         slots = c(
           wd = "character",
           resDir = "character",
           tmpDir = "character",
           polarity = "character",
           ce = "character",
           lc = "character",
           instrument = "character",
           ms1range = "nullOrNumeric",
           ms2range = "nullOrNumeric",
           resDefineAt = "numeric",
           nSlaves = "numeric"
         ),
         contains = c("ParamSterol")
)

#' @export
setClass("AssignMSMSParam",
         slots = c(
           methodAssign = "character",
           methodRT = "character",
           rtTol = "numeric",
           rtExtend = "numeric",
           mzTol = "numeric",
           ccsTol = "numeric",
           rtUnitMS1 = "character",
           rtUnitMS2 = "character",
           isolationWindow = "numeric",
           thrIsolationIntensity = "numeric",
           thrIntMS1Predicted = "numeric",
           snthreshMSMS = "numeric",
           minfracVote = "numeric",
           outputSpec = "logical",
           plotMSMSEIC = "logical",
           typeMSMSEIC = "character",
           rerun = "logical"
         ),
         contains = c("ParamSterol")
)
