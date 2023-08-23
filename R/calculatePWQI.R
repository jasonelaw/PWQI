#TODO: Deal with Temperature spawning versus nonspawning.

# Constants
kWatersheds <- c('Fanno Creek', 'Johnson Creek', 'Tualatin Tributaries', 'Tryon Creek',
                 'Willamette River Tributaries', 'Columbia River', 'Columbia Slough', 'Willamette River')

#' Polynomial coefficients
#'
#' These objects are polynomial coefficients for subindex functions. There
#' are separate index functions for each analyte. Some analytes are split
#' into multiple index functions by watershed, when some watersheds have
#' watershed specific issues (TMDLs, etc) that need to be incorporated
#' into the index. The polynomial coefficients are from lower to higher
#' order, left to right.
#'
#' These objects are not exported so cannot be accessed directly by users.
#' use \code{PWQI:::} syntax to access an object. See the examples section
#' for an example
#'
#' @examples
#' PWQI:::kDORiver
#' @name poly-coef
#' @import polynom
NULL

#' @rdname poly-coef
kDORiver           <- c( -41.343,   29.802,    -2.2154,    0.0596)
#' @rdname poly-coef
kDOStream          <- c( -58.232,   27.31,     -1.6464,    0.0392)
#' @rdname poly-coef
kCopper            <- c( 102.17,   -25.887,     2.6551,   -0.1171)
#' @rdname poly-coef
kTSS               <- c( 104.04,    -2.7238,    0.0278,   -0.0001)
#' @rdname poly-coef
kPhosphorusFanno   <- c( 103.89,  -399.76,    527.6,    -240.3)
#' @rdname poly-coef
kPhosphorus        <- c( 102.79,  -331.24,    395.32,   -167.67)
#' @rdname poly-coef
kTemperature       <- c(  29.668,   18.751,    -1.3944,    0.0248)
#' @rdname poly-coef
kTemperatureRiver  <- c(-101.87,    38.843,    -2.2667,    0.0365)
#' @rdname poly-coef
kMercury           <- c(  98.426,  -23.898,     2.5939,   -0.114)
#' @rdname poly-coef
kMercuryWillamette <- c( 100.41,   -54.278,    12.318,    -1.0076)
#' @rdname poly-coef
kEcoli             <- c( 101.29,    -0.1314,    8e-05,    -2.1e-08)
#' @rdname poly-coef
kAmmoniaRatio      <- c( 100,     -40)
#' @rdname poly-coef
kDOsat             <- c(0.55, 0.0095, -5e-05)

#' Internal functions
#'
#' Internal functions not meant to be called directly by users.
#'
#' @name pwqi-internal
NULL

#' @rdname pwqi-internal
Curry   <- functional::Curry

#' @rdname pwqi-internal
Compose <- functional::Compose

#' @rdname pwqi-internal
clamp <- function(x, a, b){ pmax(a, pmin(x, b)) }

#' @rdname pwqi-internal
iclamp <- function(x){ pmax(10, pmin(x, 100)) }

#' @rdname pwqi-internal
is.within <- function(x, r){ (x > r[1]) & (x < r[2]) }

#' @rdname pwqi-internal
filterRoots <- function(p, k){
  s <- solve(p, k)
  if(is.complex(s) && length(s) > 1){
    s <- Re(Filter(function(x) Im(x) == 0, s))
  } else if(length(s) > 1) {
    sol.deriv <- solve(deriv(p))
    f <- Curry(is.within, r = sol.deriv)
    s <- Filter(f, s)
  }
  return(s)
}

#' Creates a function to map a variable to a subindex score via a polynomial
#'
#' Function to deal with annoying discontinuities in index; goal is
#' to get index function as close to continous as possible. Solve
#' for roots when p(x) = 10 and p(x) = 100. When there are 3 real root
#' choose correct real root (using deriv(p(x)) = 0 to bracket).
#' Use roots (correct 10 and correct 100) to clamp input to poly, use (10,100) to clamp output of poly.
#' Poly either has one real root so the real Filter gets it or has three real
#' roots and need to bracket using the derivative = 0.
#' See http://en.wikipedia.org/wiki/Cubic_function#The_nature_of_the_roots
#' @param x a vector of polynomial coefficients
createIndexFun <- function(x){
  p <- polynom::polynomial(x)
  s10  <- filterRoots(p, 10)
  s100 <- filterRoots(p, 100)
  if (s10 < s100){
    f <- Curry(clamp, a = s10, b = s100)
  } else {
    f <- Curry(clamp, a = s100, b = s10)
  }
  Compose(f, polynom:::as.function.polynomial(p), iclamp)
}

#' Subindex scoring functions
#'
#' These functions map subindex raw data values to subindex scores
#' on a scale of 10 to 100. They are all created by the \link{createIndexFun} function.
#' @param ... a vector of raw subindex scores (e.g., copper concentrations for \code{calcCU}).
#' @examples
#' calcCu(c(0,4,8))
#' curve(calcCu, from = 0, to = 8)
#' @name subindex
NULL

#' @export
#' @rdname subindex
calcCu       <- createIndexFun(kCopper)
#' @export
#' @rdname subindex
calcTSS      <- createIndexFun(kTSS)
#' @export
#' @rdname subindex
calcDO       <- createIndexFun(kDOStream)
#' @export
#' @rdname subindex
calcDORiv    <- createIndexFun(kDORiver)
#' @export
#' @rdname subindex
calcP        <- createIndexFun(kPhosphorus)
#' @export
#' @rdname subindex
calcPF       <- createIndexFun(kPhosphorusFanno)
#' @export
#' @rdname subindex
calcTemp     <- createIndexFun(kTemperature)
#' @export
#' @rdname subindex
calcTempRiv  <- createIndexFun(kTemperatureRiver)
#' @export
#' @rdname subindex
calcHg       <- createIndexFun(kMercury)
#' @export
#' @rdname subindex
calcHgWill   <- createIndexFun(kMercuryWillamette)
#' @export
#' @rdname subindex
calcEcoli    <- createIndexFun(kEcoli)
#' @export
#' @rdname subindex
calcNH4      <- createIndexFun(kAmmoniaRatio)

#' @rdname pwqi-internal
adjustDO <- function(index, do, cel){
  do.sat <- convertDOConcentrationToSaturation(do, cel)
  fp <- polynom:::as.function.polynomial(polynom::polynomial(kDOsat))
  fc <- Curry(clamp, a = 0.45, b = 1)
  do.adjust <- Compose(fp, fc)
  ret <- ifelse(do.sat > 100, index * do.adjust(do.sat), index)
  ret
}

#' Calculate the Portland Water Quality Index

#' Calculate the Portland Water Quality Index as used in the watershed report cards.
#' The function is vectorized. Each argument should be the same length, excepting pH which
#' may be length 1.
#'
#' The PWQI uses cubic polynomials to map subindex raw values to subindex scores, which are supposed to range from 10 - 100.
#' However, the polynomials extend beyond the subindex range for allowable values of the raw data (e.g, the value of the
#' copper polynomial at 0 is 102.17). This means that a lot of senseless acrobatics is necessary to constrain the index
#' between 10 - 100. Unfortunately, we can't do this by just constraining the output; the functions are cubic polynomials
#' so can start going back DOWN for high values of the raw data (e.g., for copper 7 is index 10.89, but at 7.5 is index 7.97)
#' Thus, we have to constrain both the input and the output. I do this by solving for the roots at 10 and 100, and constraining
#' the input to be between the roots, then constrain the output to 10, 100 to remove numerical spill-over at 10 or 100.
#'
#' There are several watersheds that have special subindex scoring functions for certain parameters. The Columbia River, Willamette River,
#' and Columbia Slough have a separate dissolved oxygen function. The Fanno Watershed has a separate phosphorus function. The
#' Willamette and Columbia have their own temperature function. Finally, the Willamette has its own mercury function.
#'
#' @param watershed name of the watershed
#' @param Cu dissolved copper concentration in ug / L
#' @param DO dissolved oxygen in mg / L
#' @param Ecoli E. coli concentration in 1 / 100 mL
#' @param Hg total mercury concentration in ug / L
#' @param NH4 ammonia-nitrogen concentration in mg / L
#' @param P phosphorus concentration in mg / L
#' @param TSS total suspended solids concentration in mg / L
#' @param Cel temperature in Celsius
#' @param pH pH
#' @return A data.frame of water quality index subindices and the overall index.
#' @export
#' @examples
#'
#' kParameters <- c('Dissolved Copper' = 'Cu', 'Dissolved Oxygen' = 'DO', 'E. coli' = 'Ecoli',
#' 'Total Mercury' = 'Hg', 'Ammonia' = 'NH4',
#' 'Total Phosphorus' = 'P', 'Total Suspended Solids' = 'TSS', 'Temperature' = 'Cel')
#' kWatersheds <- c('Fanno Creek', 'Johnson Creek', 'Tualatin Tributaries', 'Tryon Creek',
#' 'Willamette River Tributaries', 'Columbia River', 'Columbia Slou  gh', 'Willamette River')
#'
#' kTestVals <- c(1.783, 6.1, 108, 5.142, 0.284, 0.259, 38, 18.5)
#' d <- matrix(kTestVals, dimnames = list(NULL, kParameters), ncol = length(kTestVals), nrow = 8, byrow=T)
#' d <- as.data.frame(d)
#' d$watershed <- kWatersheds
#'
#' with(d, PWQI(watershed, Cu, DO, Ecoli, Hg, NH4, P, TSS, Cel))
#' PWQI(watershed = 'Tryon Creek', Cu = 7.2, TSS = 98, Ecoli = 1596, NH4 = 2.25, DO = 3, P = 0.61, Cel = 27, Hg = 9.1)
#' PWQI(watershed = 'Tryon Creek', Cu = 0.05, TSS = 1, Ecoli = 9, NH4 = 0, DO = 12, P = 0.008, Cel = 5, Hg = 0)
PWQI <- function(watershed, Cu, DO, Ecoli, Hg, NH4, P, TSS, Cel, pH = 7.8){
  stopifnot(is.element(watershed, kWatersheds))
  i.cu    <- calcCu(Cu)
  i.tss   <- calcTSS(TSS)
  i.ecoli <- calcEcoli(Ecoli)
  NH4.crit <- ORDEQWaterQualityCriteria::nh4Criteria(pH, Cel, toxicity = 'chronic')
  i.nh4   <- calcNH4(NH4 / NH4.crit)
  i.do    <- ifelse(watershed %in% kWatersheds[c(6:8)], calcDORiv(DO), calcDO(DO)) # Will, Slough, Col
  i.p     <- ifelse(watershed %in% kWatersheds[c(1,3)], calcPF(P), calcP(P)) # Fanno
  i.temp  <- ifelse(watershed %in% kWatersheds[c(6,8)], calcTempRiv(Cel), calcTemp(Cel)) # Will, Col
  i.hg    <- ifelse(watershed %in% kWatersheds[8], calcHgWill(Hg), calcHg(Hg)) # Will
  i.do    <- adjustDO(i.do, DO, Cel)
  indices <-  cbind(Cu = i.cu, TSS = i.tss, Ecoli = i.ecoli, NH4 = i.nh4,
                    DO = i.do, P = i.p, Temp = i.temp, Hg = i.hg)
  pwqi <- apply(indices, 1, function(x) 1/mean(1/x))
  ret <- data.frame(watershed, Cu = i.cu, TSS = i.tss, Ecoli = i.ecoli, NH4 = i.nh4,
                    DO = i.do, P = i.p, Temp = i.temp, Hg = i.hg, PWQI = pwqi)
  ret[-1] <- lapply(ret[-1], round, digits = 0)
  ret
}
