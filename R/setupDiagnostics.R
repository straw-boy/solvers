#' Setup a data.frame of diagnostics
#'
#' @param res the result from calling the C++ routine used to fit a model
#'   in SLOPE
#'
#' @return A data.frame
#'
#' @keywords internal
setupDiagnostics <- function(res) {
  iteration_timings <- res$iteration_timings
  primals <- res$primals
  duals <- res$duals
  loss <- res$loss

  nl <- length(iteration_timings)
  nn <- lengths(iteration_timings)
  iteration_timings <- unlist(iteration_timings)
  primal <- unlist(primals)
  dual <- unlist(duals)
  loss <- unlist(loss)
  
  data.frame(iteration = unlist(lapply(nn, seq_len)),
             time_taken = iteration_timings,
             primal = primal,
             dual = dual,
             loss = loss)
}
