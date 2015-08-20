#' calc.pnc
#'
#' Calculates probabilities of peak heights as returned from the num.clusts ()
#' function.  To generate useful statistics, this function should be run with a 
#' VERY large number of repeats (for example, 100,000), which will take a 
#' LONG time!
#'
#' @param city nyc, washington, chicago, boston, london (case insensitive)
#' @param method complete, ward, or skater, to be compared to k-means
#' @param nprts number of random trials used to calculate probabilities 
#' @param rescale order of polynomial rescaling of peak heights
#' @return data.frame and also dumps results to text file

calc.pnc <- function (city = "nyc", method="complete", nrpts=100, rescale=2)
{
    if (tolower (substring (city, 1, 1)) == "n")
        city <- "nyc"
    else if (tolower (substring (city, 1, 1)) == "l")
        city <- "london"
    else if (tolower (substring (city, 1, 1)) == "w" |
                tolower (substring (city, 1, 1)) == "d")
        city <- "washingtondc"
    else if (tolower (substring (city, 1, 1)) == "c")
        city <- "chicago"
    else if (tolower (substring (city, 1, 1)) == "b")
        city <- "boston"
    else 
        stop ("city not valid")

    method <- tolower (substring (method, 1, 1))
    if (!method %in% c ("w", "s", "k"))
        method <- "complete"
    else if (method == "w")
        method <- "ward"
    else if (method == "s")
        method <- "skater"
    
    nrpts <- nrpts * 2 
    # Presumes M \approx 4.5, so about half of all simulated series will have the
    # right M-values, therefore nrpts is doubled.
    
    dat <- num.clusts (city=city, method=method, plot=FALSE)
    cat ("\rcalculating peak probabilities for ", city, " ", method, "\n", sep="")
    dat <- dat [[2]] 
    # [[1]] is just the names of the four groups (to-h, to-k, from-h, from-k)
    # dat has 13 columns of $nc (number of peaks) plus 4 groups of G-values,
    # peak heights ($h), and lengths of series. The latter four are the position
    # of the final peak plus one. This is the length needed to feed in to the
    # random simulations. There are some differences between the four data sets
    # in $n values, but an average value is taken for all, as specified below.
    #
    # The whole thing is run here as a loop over the 14 values of $nc = (2:15),
    # because the actual time consuming bit is the poly.rescale calls internal
    # to this loop, and it makes the code clearer.
    pr.pk.hts <- num.samples <- array (NA, dim=c(nrow (dat), 4))
    
    pb <- txtProgressBar (0, 100, char="=", style=3)
    for (i in 1:nrow (dat)) {
        # First calculate the "average" length of each series:
        nvals <- as.numeric (dat [i, 10:13])
        if (length (which (nvals == min (nvals))) == 
            length (which (nvals == max (nvals)))) {
            n <- max (nvals)
        } else {
            n <- round (mean (nvals))
        }
        # Then extract the four peak heights
        pk.heights.i <- as.numeric (dat [i, 6:9])
        simvals <- array (runif (n * nrpts, min=-1, max=1), dim=c(n, nrpts))
        pks <- apply (simvals, 2, function (x) {
            which (diff (sign (diff (x))) == -2) + 1    })
        gvals <- lapply (pks, function (x) mean (diff (x)))
        gvals <- unlist (gvals)
        gvals0 <- as.numeric (dat [i, 2:5])
        # First trim x down to only gvals >= min (gvals0):
        # If observed G-values are the same, run all together, else run
        # individually
        if (sum (abs (diff (gvals0)))  == 0) { 
            g <- gvals0 [1]
            indx <- which (gvals >= g)
            num.samples [i, ] <- rep (length (indx), 4)
            prob.g <- length (indx) / nrpts
            if (length (indx > 1)) {
                sm <- simvals [, indx]
                sm <- poly.rescale (sm, rescale=rescale)
                sml <- as.list (as.data.frame (sm))
            
                pks <- apply (sm, 2, function (x) {
                    which (diff (sign (diff (x))) == -2) + 1    })
                non.pks <- lapply (pks, function (x) {
                    which (!1:n %in% x) })
                pk.height <- mapply (function (x, y, z) {
                    mean (z [x]) - mean (z [y])}, x=pks, y=non.pks, z=sml)
                pr.pk.height <- lapply (pk.heights.i, function (x) {
                    length (pk.height [pk.height >= x]) / length (pk.height)  })
                pr.pk.hts [i, ] <- unlist (pr.pk.height) * prob.g
            } else {
                pr.pk.hts [i, ] <- 0
            }
        } else { 
            # G-values differ, so each is calculated separately. Note that this
            # could be just repeated on unique values, but the time consuming
            # bit is poly.rescale, so this wouldn't make any difference anyway.
            gh0 <- as.list (as.data.frame (rbind (gvals0, pk.heights.i)))
            # gh0 is a list of four pairs of [g, h], so x [1] = g & x [2] = h.
            indx <- which (gvals >= min (gvals0))
            gvals <- gvals [indx]
            simvals <- simvals [, indx]
            simvals <- poly.rescale (simvals, rescale=rescale)
            
            temp <- lapply (gh0, function (x) {
                indx <- which (gvals >= x [1])
                if (length (indx) > 1) {
                    sm <- simvals [, indx]
                    sml <- as.list (as.data.frame (sm))
                    prob.m <- length (indx) / nrpts
        
                    pks <- apply (sm, 2, function (z) {
                        which (diff (sign (diff (z))) == -2) + 1    })
                    non.pks <- lapply (pks, function (z) {
                        which (!1:n %in% z) })
                    pk.height <- mapply (function (i1, i2, z) {
                        mean (z [i1]) - mean (z [i2])}, i1=pks, i2=non.pks, z=sml)
                    result <- length (pk.height [pk.height >= x [2]]) * prob.m / 
                        length (pk.height)
                } else { result <- list (0, 0, 0, 0)   }
                    # next line is just a dummy line to return the value
                    result [[1]] <- result [[1]] })
            pr.pk.hts [i, ] <- unlist (temp)
            num.samples [i, ] <- unlist (lapply (gvals0, function (x) 
                             length (which (gvals >= x))))
        }
        
        setTxtProgressBar (pb, 100 * i / length (dat$nc))
    } # end for i over the nc values
    close (pb)
    
    nrpts <- rep (nrpts, dim (dat) [1])
    dat <- cbind (dat, pr.pk.hts, nrpts, num.samples)
    names (dat) [14:22] <- c ("p1", "p2", "p3", "p4", "nrpts", 
                              "ns1", "ns2", "ns3", "ns4")
    fname <- paste (city, "-results-prob-m-", method, ".txt", sep="")
    write.table (dat, file=fname, row.names=FALSE, sep=",")
    return (dat)
} # end function calc.pnc
