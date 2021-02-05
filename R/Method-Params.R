#' @export
setMethod("as.list", signature(x = "ParamSterol"), function(x, ...) {
  return(.param2list(x))
})

# The 'setAs' method.
setAs("ParamSterol" ,"list", function(from){
  return(.param2list(from))
})
