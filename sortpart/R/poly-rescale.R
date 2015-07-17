#' poly.rescale
#'
#' Applies quadratic regression to input series to rescale peak values between
#' -1 and 1 
#'
#' @param xm input series (as matrix for multiple rescaling)
#' @param rescale polynomial order of rescaling (can be > default of 2)
#' @return vector of rescaled series

poly.rescale <- function (xm, rescale=2) {
    require (quantreg) 

    ybounds <- c (0.99, 0.01)
    n <- dim (xm) [1]
    bounds.upper <- apply (xm, 2, function (x) {
        dfr <- data.frame (tt=1:n, x=x)
        if (rescale == 2) {
            mod <- nlrq (x ~ a * tt ^ 2 + b * tt + cc, 
                         data=data.frame (tt=1:n, x=x),
                         tau=ybounds [1], start=list(a=0, b=0, cc=mean(x)))
        } else {
            mod <- nlrq (x ~ a * tt ^ 3 + b * tt ^ 2 + cc * tt + dd, 
                         data=data.frame (tt=1:n, x=x),
                         tau=ybounds [1], start=list(a=0, b=0, cc=0, dd=mean(x)))
        }
        predict (mod, newdata=dfr$tt)   })
    bounds.lower <- apply (xm, 2, function (x) {
        dfr <- data.frame (tt=1:n, x=x)
        if (rescale == 2) {
            mod <- nlrq (x ~ a * tt ^ 2 + b * tt + cc, 
                         data=data.frame (tt=1:n, x=x),
                         tau=ybounds [2], start=list(a=0, b=0, cc=mean(x)))
        } else {
            mod <- nlrq (x ~ a * tt ^ 3 + b * tt ^ 2 + cc * tt + dd, 
                         data=data.frame (tt=1:n, x=x),
                         tau=ybounds [2], start=list(a=0, b=0, cc=0, dd=mean(x)))
        }
        predict (mod, newdata=dfr$tt)   })
    xm <- 2 * (xm - bounds.lower) / (bounds.upper - bounds.lower) - 1
    return (xm)
} # end function rescale.xm
