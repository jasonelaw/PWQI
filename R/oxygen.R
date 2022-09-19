
kCelsiusToKelvin <- 273.15

#' Compute salinity as a function of conductivity
#'
#' This function computes salinity as a function of conductivity
#' using equation 4 from Office of Water Quality Technical Memorandum 2011.03
#' \url{https://water.usgs.gov/admin/memo/QW/qw11.03.pdf}
#' @param x conductivity
#' @return a vector of salinity values
#' @export
convertCondToSalinity <- function(x){
  5.572e-4 * x + 2.02e-9 * x^2
}

#' Compute oxygen solubility as a function of temperature
#'
#' This function computes the solubility of oxygen in water in freshwater
#' at standard pressure (760 mm). See Office of Water Quality Technical Memorandum 2011.03
#' (\url{https://water.usgs.gov/admin/memo/QW/qw11.03.pdf}) and
#' \url{http://water.usgs.gov/software/DOTABLES/}. The formula is:
#' \deqn{-139.34411 + \frac{1.575701e5}{t} - \frac{6.642308e7}{t^2} + \frac{1.243800e10}{t^3} - \frac{8.621949e11}{t^4}}
#' @param t numeric, temperature in degrees Celsius
#' @return numeric vector of oxygen solubility in water values
#' @export
calculateOxygenSolubility <- function(t){
  stopifnot(t >= 0, t <= 40, is.vector(t))
  t <- t + kCelsiusToKelvin #convert to Kelvin
  exp(-1.3934411e2 + 1.575701e5 * 1/t + -6.642308e7 * 1/t^2 + 1.243800e10 * 1/t^3 + -8.621949e11 * 1/t^4)
}

#' Convert dissolved oxygen between \% saturation and concentration
#'
#' These functions convert between dissolved oxygen measurements expressed
#' as \% saturation and concentration (mg / L). They use \link{calculateOxygenSolubility}
#' to calculate the solubility of oxygen and use that to convert.
#'
#' @param do a numeric vector of dissolved oxygen measurements, \% saturation for
#' \code{convertDOSaturationToConcentration} and concentration in mg/L for
#' \code{convertDOConcentrationToSaturation}.
#' @param t a numeric vector of temperature measurements, in degrees Celsisus
#' @return a vector of dissolved oxygen measurements, in mg/L for
#' \code{convertDOSaturationToConcentration} and \% saturation for
#' \code{convertDOConcentrationToSaturation}
#' @name convertDO
NULL

#' @describeIn convertDO Convert a dissolved oxygen expressed as percent saturation to concentration.
#' @export
convertDOSaturationToConcentration <- function(do, t){
  stopifnot(do > 1, do < 120)
  # do in percent saturation and t in degrees Celsius
  sol <- calculateOxygenSolubility(t) # calculate saturated concentration at temperature
  sol * do / 100           # calculate DO concentration based on % saturation and saturated concentration
}

#' @describeIn convertDO Convert a dissolved oxygen concentration to percent saturation.
#' @export
convertDOConcentrationToSaturation <- function(do, t){
  stopifnot(do > 0, t > 0, t < 40)
  sol <- calculateOxygenSolubility(t)
  do / sol * 100
}



