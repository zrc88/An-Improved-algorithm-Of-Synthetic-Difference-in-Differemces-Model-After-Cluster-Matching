#' Summarize a synthdid object
#' @param object The object to summarize
#' @param weight.digits The number of digits to use when displaying weights (omega, lambda)
#' @param fast Be fast but less accurate, e.g. jackknife instead of bootstrap se estimate
#' @param ... Additional arguments (currently ignored).
#' @method summary synthdid_estimate
#' @export
summary.synthdid_estimate = function(object, weight.digits=3, fast=FALSE, ...) {
  N0 = attr(object, 'setup')$N0
  T0 = attr(object, 'setup')$T0
  list(estimate = c(object),
       se = sqrt(if(fast) { vcov(object, method = 'jackknife') } else { vcov(object) }),
       controls = round(synthdid_controls(object, weight.type='omega'),  digits=weight.digits),
       periods  = round(synthdid_controls(object, weight.type='lambda'), digits=weight.digits),
       dimensions = c( N1 = nrow(Y(object))-N0, N0 = N0, N0.effective = round(1 / sum(omega(object)^2),  weight.digits),
                       T1 = ncol(Y(object))-T0, T0 = T0, T0.effective = round(1 / sum(lambda(object)^2), weight.digits)))
}
