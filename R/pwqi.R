#
# scl <- function(x, min , max){
#   (x - min) / (max - min)
# }
#
# iscl <- function(x, min, max){
#   x * (max - min) + min
# }
#
# flog   <- function(x, r = 1, m = 0.5){
#   1 / (1 + exp(r * (x - m)))
# }
#
# fexp <- function(x, r){
#   exp(r * x)
# }
#
# # Index Functions
# iTSS       <- Curry(fexp, r = -0.03)
# iPhos      <- Curry(fexp, r = -4.00)
# iPhosFanno <- Curry(fexp, r = -4.5)
#
# # Spline Interpolation
# s_tss <- splinefun(c(0,10,20,43,100), c(100,80,60,30,10), method = "hyman")
#
# #TSS
# plot(polynomial(PWQI:::kTSS), xlim = c(0, 120), ylim = c(10, 100))
# curve(iscl(fexp(x, -0.03), 10, 100), add = TRUE)
# #Phosphorus
# plot(polynomial(PWQI:::kPhosphorus), xlim = c(0, 1), ylim = c(10, 100))
# curve(iscl(fexp(x, -4), 10, 100), add = TRUE)
# plot(polynomial(PWQI:::kPhosphorusFanno), xlim = c(0, 1), ylim = c(10, 100))
# curve(iscl(fexp(x, -4.5), 10, 100), add = TRUE)
# #Temperature
# plot(polynomial(PWQI:::kTemperature), xlim = c(0, 30), ylim = c(10, 100))
# curve(iscl(flog(x, .4, m = 20), 10, 100), add = TRUE)
