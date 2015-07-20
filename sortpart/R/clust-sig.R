#' clust.sig
#'
#' Uses the output of calc.pnc, which calculates the probability of observing a
#' given peak height for a given number of clusters. This routine in turn uses
#' the output of num.clusts. 
#'
#' @param city nyc, washington, chicago, boston, london (case insensitive)
#' @param method = (ward, k-means, complete)
#' @param xmax maximum number of clusters to be analysed and plotted
#' @return statistical summary as screen dump

clust.sig <- function (city="nyc", method="complete", xmax=36)
{
    require (quantreg) 
    
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
    else if (method == "k")
        method <- "k-means"
    
    fname <- paste ("./results/", city, "-results-actual-to-", 
                    method, ".txt", sep="")
    if (!file.exists (fname))
        stop (cat (fname, " does not exist!", sep=""))
    else
        d.to <- read.csv (fname, header=TRUE)
    fname <- paste ("./results/", city, "-results-actual-from-", 
                    method, ".txt", sep="")
    if (!file.exists (fname))
        stop (cat (fname, " does not exist!", sep=""))
    else
        d.from <- read.csv (fname, header=TRUE)

    fname <- paste ("./results/", city, "-results-actual-to-k-means.txt", sep="")
    if (!file.exists (fname))
        stop (cat (fname, " does not exist!", sep=""))
    else
        d.to.kmeans <- read.csv (fname, header=TRUE)
    fname <- paste ("./results/", city, "-results-actual-from-k-means.txt", sep="")
    if (!file.exists (fname))
        stop (cat (fname, " does not exist!", sep=""))
    else
        d.from.kmeans <- read.csv (fname, header=TRUE)

    fname <- paste ("./results/", city, "-results-neutral.txt", sep="")
    if (!file.exists (fname))
        stop (cat (fname, " does not exist!", sep=""))
    else
        d.neutral <- read.csv (fname, header=TRUE)

    # The neutral files just contain $(nsd,dmn,dsd), where the distances are
    # total intra-cluster distances which can be converted to proportions
    # through dividing by the total overall distances of:
    total.dist <- get.total.dist (city=city)
    nc <- d.neutral$nc

    # First plot is ratio of observed to expected distance ridden within clusters
    d.in <- cbind (d.to$d.in / d.neutral$dmn, 
        d.from$d.in / d.neutral$dmn)
    d.in.kmeans <- cbind (d.to.kmeans$d.in / d.neutral$dmn, 
        d.from.kmeans$d.in / d.neutral$dmn)

    # Second is T-values calculated from observed versus expected distances
    d0 <- cbind (d.to$d.in, d.from$d.in)
    t0 <- (d0 - d.neutral$dmn) / (d.neutral$dsd / sqrt (nc))
    d0 <- cbind (d.to.kmeans$d.in, d.from.kmeans$d.in)
    tk <- (d0 - d.neutral$dmn) / (d.neutral$dsd / sqrt (nc))
    t05 <- -qt (0.05, nc)
    t01 <- -qt (0.01, nc)

    xdat <- nc
    ydat1 <- ydat2 <- list ()
    ydat1 [[1]] <- d.in
    ydat1 [[2]] <- t0
    ydat2 [[1]] <- d.in.kmeans
    ydat2 [[2]] <- tk
    ytxts <- list ()
    ytxts [[1]] <- "Proportional distance ridden within clusters"
    ytxts [[2]] <- "T-values"

    x11 (width=8, height=8)
    #layout (matrix (c(1:4,5,5), 3, 2, byrow=TRUE), height=c(0.46,0.46,0.08))
    par (mfrow=c(2,2)) # needs to be (1, 3) if mutual information plot is included
    par (mar=c(2.5,2.5,2,1), mgp=c(1.3,0.4,0), ps=10, tcl=-0.2)

    cols <- c ("red", "blue")
    xlabs <- "Number of Clusters"
    legxpos <- xdat [1] + 0.6 * diff (range (xdat))
    legxpos <- c (0.15, 0.6) * xmax

    xlims <- range (xdat)
    xlims [2] <- xmax
    indx <- which (xdat <= xmax)

    ybounds <- c (0.99, 0.01) # Upper and lower bounds for nlqr regressions
    for (i in 1:2) {
        ylims <- range (c (ydat1 [[i]] [indx, ], ydat2 [[i]] [indx, ]))
        plot (xdat, ydat1 [[i]] [,1], "l", col=cols [1], lwd=1, 
              xlim=xlims, ylim=ylims, xlab="", ylab=ytxts [[i]], xaxt="n", yaxt="n")
        axis (1, at=1:10 * 4)
        axis (2, at=pretty (ylims))
        for (j in 1:2) {
            lines (xdat, ydat1 [[i]] [,j], col=cols [j])
            lines (xdat, ydat2 [[i]] [,j], col=cols [j], lty=2)
        }
        for (j in 1:25) {
            lines (rep (j*4, 2), ylims, col="gray", lty=2)
        }
        title (xlab=xlabs)
        if (i == 1) {
            ypos <- ylims [1] + 0.2 * diff (ylims) # 0.88 with MI plots
            legend (legxpos [1], ypos, lwd=1, col=cols, bty="n",
                legend=c("to (h-clust)", "from (h-clust)"))
            legend (legxpos [2], ypos, lwd=1, col=cols, bty="n", lty=2,
                legend=c("to (k-means)", "from (k-means)"))
        }
        plotbounds <- FALSE; if (plotbounds) {
            # Then plot upper and lower bounds of h-clust T-statistics
            if (i == 2) {
                indx <- which (xdat <= xmax)
                ulbounds <- array (NA, dim=c(length (indx), 2))
                for (j in 1:2) { # Over (to, from) data
                    dfr <- data.frame (x=xdat [indx], y=ydat1 [[i]] [indx, j])
                    for (k in 1:2) {
                        mod <- nlrq (y ~ a * x ^ 2 + b * x + cc, data=dfr, 
                            tau=ybounds [k], start=list(a=0, b=0, cc=mean(dfr$y)))
                        ulbounds [,k] <- predict (mod, newdata=dfr$x)
                        pdfr <- data.frame (x=seq(min (dfr$x), max (dfr$x), len=100))
                        pdfr <- within (pdfr, y <- predict (mod, newdata=pdfr))
                        lines (pdfr$x, pdfr$y, col=cols [j], lty=2, lwd=2)
                    } # end for k
                } # end for j
            } # end if i
        } # end if plotbounds
    } # end for i over the 2 plots

    # Then plot re-scaled versions of hclust results
    indx <- which (xdat <= xmax)
    ulbounds <- array (NA, dim=c(length (indx), 2))
    y.resc <- array (NA, dim=c(length (indx), 4))
    ytemp <- cbind (t0 [indx, 1], tk [indx, 1], t0 [indx, 2], tk [indx, 2])
    for (j in 1:4) { # (to-complete, to-k-means, from-complete, from-k-means)
        #dfr <- data.frame (x=xdat [indx], y=t0 [indx, j])
        dfr <- data.frame (x=xdat [indx], y=ytemp [, j])
        for (k in 1:2) {
            mod <- nlrq (y ~ a * x ^ 2 + b * x + cc, data=dfr, 
                tau=ybounds [k], start=list(a=0, b=0, cc=mean(dfr$y)))
            ulbounds [,k] <- predict (mod, newdata=dfr$x)
        } # end for k
        y.resc [, j] <- 2 * (ytemp [, j] - ulbounds [,2]) /
            (ulbounds [,1] - ulbounds [,2]) - 1
    } # end for j
    plot (xdat [indx], y.resc [,1], "l", col=cols [1], lwd=1, ylim=range (y.resc),
          ylab="Rescaled T-values", main="Rescaled T-values")
    for (i in 1:25)
        lines (rep (i*4, 2), range (y.resc), col="gray", lty=2)
    for (i in 1:4)
        lines (xdat [indx], y.resc [,i], col=cols [ceiling (i/2)], lty=2-(i%%2))


    # *************************************************************************
    # ************************   STATISTICAL ANALYSES  ************************    
    # *************************************************************************
    ttxt <- c (" to (hc) ", "to (km)", "from (hc)", "from (km)")
    tvals <- cbind (t0 [,1], tk [,1], t0 [,2], tk [,2])
    gvals <- dvals <- rep (NA, 4)
    rescale <- 2 # If 3, then O(3) bounds are calculated
    npeaks <- get.num.clusts () # returns $num.clusts and $pmin (each of len=4)
    # The following are for subsequent analyses of mean values
    gout.hc <- gout.km <- nout.hc <- nout.km <- mout.hc <- mout.km <- NULL 
    meth <- rep (c ("complete", "k-means"), 2) # for get.partition.neighbours
    dirs <- c (rep ("to", 2), rep ("from", 2)) # ditto

    # Calculate how distinct the peaks are. This is done by first rescaling the
    # T-statistics according to non-linear regressions of lower and upper
    # bounds, so that they lie between -1 and 1. Peak distinction is then
    # measured as the average difference to the 2 adjacent points. This measure
    # thus has a theoretical maximum of 2, and a minimum of 0.
    #
    # Values for s(M=1) in the loop that follows need adjusting, because the
    # calculated value is 1, yet any real value of M<1.5 will still round to 1
    # and produce this value. These s values are therefore replaced with the
    # corresponding ones for M=1.5, which are:
    s1 <- 1 / (1 + log (1.5))
    e1 <- exp (1)

    tstats <- list (NULL, NULL, NULL, NULL)
    # tstats collects T-statistics for HC only for differences between means and
    # values of both e and 2. These are plotted in final panel.
    for (i in 1:4) { # Over (to, from) data
        # Find the position of the final peak (determined by npeaks), and fit the
        # nlrq lines only to that point.
        pks <- which (diff (sign (diff (tvals [,i]))) == -2) + 1
        if (length (pks) > npeaks$num.clusts [i]) { 
            pks <- pks [1:npeaks$num.clusts [i]]   }
        indx <- 1:max (pks)
        ulbounds <- array (NA, dim=c(length (indx), 2))

        dfr <- data.frame (x=xdat [indx], y=tvals [indx, i])
        for (j in 1:2) {
            if (rescale == 2) {
                mod <- nlrq (y ~ a * x ^ 2 + b * x + cc, data=dfr, 
                    tau=ybounds [j], start=list(a=0, b=0, cc=mean(dfr$y)))
            } else if (rescale == 3) {
                mod <- nlrq (y ~ a * x ^ 3 + b * x ^ 2 + cc * x + dd, data=dfr, 
                    tau=ybounds [j], start=list(a=0, b=0, cc=0, dd=mean(dfr$y)))
            }
            ulbounds [,j] <- predict (mod, newdata=dfr$x)
            # Fitted values, not used here
            #pdfr <- data.frame (x=seq(min (dfr$x), max (dfr$x), len=100))
            #pdfr <- within (pdfr, y <- predict (mod, newdata=pdfr))
        } # end for j
        dfr$y <- 2 * (dfr$y - ulbounds [,2]) / (ulbounds [,1] - ulbounds [,2]) - 1

        # See calc.pm for explanation and interpretation of the following lines:
        np <- length (pks)
        pks2 <- cbind (pks [1:(np - 1)], pks [2:np])
        mins <- apply (pks2, 1, function (z) {
                z [1] + which.min (dfr$y [z [1]:z [2]]) - 1 })
        pks2 <- cbind (pks2, mins)
        pk.depth <- apply (pks2, 1, function (z) {
                mean (dfr$y [z [1:2]]) - dfr$y [z [3]]  })
        gvals <- diff (xdat [pks])
        g <- mean (gvals)
        svals <- 1 / log (e1 * gvals / (e1 - 1))
        s <- mean (svals)
        nstars <- length (pks) * 3 + 46
        cat (rep ("*", nstars), "\n", sep="")
        cat ("********** T-VALUE PEAKS", ttxt [i], "= ", xdat [pks], "**********\n")
        cat (rep ("*", nstars), "\n", sep="")
        cat ("\tG = ", formatC (g, format="f", digits=3), " +/- ",
             formatC (sd (gvals), format="f", digits=3), "; s = ",
             formatC (s, format="f", digits=3), " +/- ",
             formatC (sd (svals), format="f", digits=3), "\n", sep="")

        # Analyse partition hierarchy:
        cat ("Sizes of new partitions formed from",
             "contiguous (c) and all (a) neighbours:\n")
        cat (" npks\tstage\t|\tN(c/a)\tM(c/a)\tN / M\t\t[T,p(N/M=e)]\t\t[T,p(N/M=2)]\t\t|\n")
        cat (rep ("-", 105), "\n", sep="")
        # Note that partition size and merge calculations exclude the initial
        # partition from 1 to xdat [pks] [1].
        part.sizesC <- part.sizes <- mvals <- NULL
        for (j in 2:length (xdat [pks])) 
        {
            nci <- c (xdat [pks] [j-1], xdat [pks] [j])
            nbs <- get.partition.neighbours (city=city, method=meth [i],
                                             nc=nci, dir=dirs [i])
            # -----------------------------------------------
            # To analyse ALL new partitions, instead of the just the contiguous
            # ones, insert this following definition of npeaks (returned from
            # get.num.clusts) above:
            #npeaks$num.clusts <- rep (50, 4)
            # And then use the following definitions of part.sizesC and mvals:
            #part.sizesC <- c (part.sizesC, 
            #                  unlist (lapply (nbs, function (x) sum (x))))
            #mvals <- c (mvals, unlist (lapply (nbs, function (x) length (x))))
            # The result is that ALL series converge to estimates of N/M=2, rather
            # than the theoretically predicted value of 2.7. 
            # ------------------------------------------------

            # Only analyse as long as the biggest contiguous partition represents
            # the majority of all new partitions:
            nfrac <- sum (nbs [[1]]) / sum (unlist (nbs))
            mfrac <- length (nbs [[1]]) / length (unlist (nbs))
            if (nfrac >= 0.5 & mfrac >= 0.5) 
            {
                part.sizes <- c (part.sizes, sum (unlist (nbs)))
                part.sizesC <- c (part.sizesC, sum (nbs [[1]]))
                mvals <- c (mvals, length (nbs [[1]]))
                cat ("| ", j, "\t", nci [1], "->", nci [2], "\t|\t", sum (nbs [[1]]), 
                     " / ", sum (unlist (nbs)), "\t", length (nbs [[1]]),
                     " / ", length (unlist (nbs)), sep="")
                if (j == 2) 
                    cat ("\t\t\t\t\t\t\t\t\t|\n")
                else 
                {
                    nm <- part.sizesC / mvals
                    # nm values for k-means are sometimes all identical, which
                    # causes the t-test to crash, as does the odd occasion when
                    # only a single value can be extracted before the fractions
                    # drop below 0.5:
                    if (sd (nm) == 0)
                        stop ("\nERROR: Estimated values of N/M",
                             " are all identical---Just run again!\n")
                    else if (length (nm) < 2)
                        stop ("\nERROR: Insufficient values of N/M generated",
                              "---Just run again!\n")
                    else 
                    {
                        tt <- t.test (nm - e1)
                        tt2 <- t.test (nm - 2)
                        cat ("\t", formatC (mean (nm, na.rm=TRUE), format="f", digits=2),
                            "+/-", formatC (sd (nm, na.rm=TRUE), format="f", digits=2),
                            "\t(", formatC (tt$statistic, format="f", digits=4),
                            ", ", formatC (tt$p.value, format="f", digits=4),
                            ")\t(", formatC (tt2$statistic, format="f", digits=4),
                            ", ", formatC (tt2$p.value, format="f", digits=4),
                            ")\t|\n", sep="")
                        if (i == 1) 
                        {
                            tstats [[1]] <- c (tstats [[1]], tt$statistic)
                            tstats [[2]] <- c (tstats [[2]], tt2$statistic)
                        } else if (i == 3) 
                        {
                            tstats [[3]] <- c (tstats [[3]], tt$statistic)
                            tstats [[4]] <- c (tstats [[4]], tt2$statistic)
                        }
                    } # end else not stop
                } # end else j > 2
            } # end else nfrac, mfrac > 0.5
        } # end for j
        cat (rep ("-", 105), "\n", sep="")
        # gvals then needs to be shorted to length of N & M:
        gvals <- gvals [1:length (mvals)]
        prop <- formatC (100 * sum (part.sizesC) / sum(part.sizes),
                         format="f", digits=2)
        cat ("contiguous / total = ", sum (part.sizesC), " / ",
             sum (part.sizes), " = ", prop, "%\n\n", sep="")

        # The following are used for the final average summaries
        if (meth [i] == "complete") 
        {
            gout.hc <- c (gout.hc, gvals)
            nout.hc <- c (nout.hc, part.sizesC)
            mout.hc <- c (mout.hc, mvals)
        } else 
        {
            gout.km <- c (gout.km, gvals)
            nout.km <- c (nout.km, part.sizesC)
            mout.km <- c (mout.km, mvals)
        }

        # Values of s are now calculated from observations of Delta G; from Delta G
        # estimated as N - M; from N; and from M.
        sG <- 1 / log (e1 * gvals / (e1 - 1))
        g.est <- part.sizesC - mvals
        sGC <- 1 / log (e1 * g.est / (e1 - 1))
        sN <- 1 / log (part.sizesC)
        sM <- 1 / (1 + log (mvals))
        sM [mvals == 1] <- s1

        sG.mean <- 1 / log (e1 * mean (gvals) / (e1 - 1))
        sGC.mean <- 1 / log (e1 * mean (g.est) / (e1 - 1))
        sN.mean <- 1 / log (mean (part.sizesC))
        sM.mean <- 1 / (1 + log (mean (mvals)))

        cat ("Mean Sizes of Partitions and Merges:\n")
        cat ("Data\t\t|\tMean\tSD\t|\ts.mn +/- SD\ts.mn (direct)\n")
        cat (rep ("-", 81), "\n", sep="")
        cat ("| All Data (G)\t|\t\t\t|\t",
             formatC (mean (sG), format="f", digits=3), " +/- ",
             formatC (sd (sG), format="f", digits=3), "\t\t",
             formatC (sG.mean, format="f", digits=3), "\t|\n", sep="")
        cat ("| Contig. (Gc)\t|\t\t\t|\t",
             formatC (mean (sGC), format="f", digits=3), " +/- ",
             formatC (sd (sGC), format="f", digits=3), "\t\t",
             formatC (sGC.mean, format="f", digits=3), "\t|\n", sep="")
        cat ("| Contig. (N)\t|\t", 
            formatC (mean (part.sizesC), format="f", digits=2), "\t",
            formatC (sd (part.sizesC), format="f", digits=2), "\t|\t", 
            formatC (mean (sN), format="f", digits=3), " +/- ",
            formatC (sd (sN), format="f", digits=3), "\t\t",
            formatC (sN.mean, format="f", digits=3), "\t|\n", sep="")
        cat ("| Merges (M)\t|\t",
            formatC (mean (mvals), format="f", digits=2), "\t",
            formatC (sd (mvals), format="f", digits=2),
            "\t|\t", formatC (mean (sM), format="f", digits=3), " +/- ",
            formatC (sd (sM), format="f", digits=3), "\t\t",
            formatC (sM.mean, format="f", digits=3), "\t|\n", sep="")
        cat (rep ("-", 81), "\n", sep="")

        # Then calculate T-tests between the four series of s-values
        cat ("T-statistics for differences in s-values (Gc = N - M):\n")
        cat ("|\t\t|\tStatistic\t\t|\tp-value\t\t\t|\n")
        cat ("|\t\t|\tGc\tN\tM\t|\tGc\tN\tM\t|\n")
        cat (rep ("-", 81), "\n", sep="")
        vlist <- c ("G", "Gc", "N")
        sarr <- cbind (sG, sGC, sN, sM)
        for (j in 1:3) 
        {
            cat ("|\t", vlist [j], "\t|\t")
            if (j > 1) 
                for (k in 1:(j-1)) 
                    cat ("\t")

            dj <- mean (diff (sarr [,j]))
            for (k in (j+1):4) 
            {
                dk <- mean (diff (sarr [,k]))
                if (dj == 0 & dk == 0)
                    cat ("NA\t")
                else
                {
                    tt <- t.test (sarr [,j], sarr [,k], var.equal=TRUE)
                    cat (formatC (tt$statistic, format="f", digits=3), "\t")
                }
            }
            # Then p-values
            cat ("|")
            for (k in 1:j) 
                cat ("\t")
            for (k in (j+1):4) 
            {
                dk <- mean (diff (sarr [,k]))
                if (dj == 0 & dk == 0)
                    cat ("NA")
                else
                {
                    tt <- t.test (sarr [,j], sarr [,k], var.equal=TRUE)
                    cat (formatC (tt$p.value, format="f", digits=3))
                }
                if (k < 4) 
                    cat ("\t")
                else 
                    cat ("\t|\n")
            }
        } # end for j
        cat (rep ("-", 81), "\n", sep="")

        # pk.depths are calculated w.r.t. background up to last peak
        dfr <- dfr [1:max (pks), ]
        non.pks <- which (!1:length (dfr$y) %in% pks)
        peak.difference <- mean (dfr$y [pks]) - mean (dfr$y [non.pks])
        cat ("   Mean peak depth = ", mean (pk.depth), 
             "; mean height above background = ", peak.difference, "\n", sep="")
        #if (i == 2)  cat (rep ("-", 80), "\n\n", sep="")
        cat (rep ("=", 81), "\n\n", sep="")
    } # end for i

    # ********************* PLOT T-STATISTICS ***********************
    #
    # These are T-statistics for respective distances to e and 2. Statistics are
    # only generated as along as the largest contiguous component represents the
    # majority of all new clusters. If this condition fails by the second peak,
    # then only one T-value will be generated for that series and there will
    # thus be nothing to plot.
    if (max (sapply (tstats, length)) > 1)
    {
        tstats <- lapply (tstats, function (x) abs (x))
        ylims <- lapply (tstats, function (x) range (abs (x)))
        ylims <- range (unlist (ylims))
        xmax <- lapply (tstats, function (x) length (x))
        xlims <- c (1, max (unlist (xmax)))
        cols <- c ("red", "red", "blue", "blue")
        ltys <- c (1, 2, 1, 2)
        plot (1:length (tstats [[1]]), tstats [[1]], "l", col=cols [1], lty=ltys [1],
              xlim=xlims, ylim=ylims, xlab="Number of peaks", ylab="T(M)")
        for (i in 1:4)
            lines (1:length (tstats [[i]]), tstats [[i]],
                   col=cols [i], lty=ltys [i], lwd=2)
        legend (xlims [1], ylims [2], lwd=2, col=cols, lty=ltys, bty="n",
                legend=c("TO(e)", "TO(2)", "FROM(e)", "FROM(2)"))

        for (i in 1:2) 
        {
            indx <- 2 * i - 1
            # indx is then "TO/FROM(e)", while (indx+1) is "TO/FROM(2)"
            i0 <- max (which (tstats [[indx]] < tstats [[indx + 1]]))
            m1 <- tstats [[indx]] [i0 + 1] - tstats [[indx]] [i0]
            m2 <- tstats [[indx + 1]] [i0 + 1] - tstats [[indx + 1]] [i0 + 1]
            c1 <- tstats [[indx]] [i0] - m1 * i0
            c2 <- tstats [[indx + 1]] [i0] - m2 * i0
            x0 <- (c2 - c1) / (m1 - m2)
            y0 <- m1 * x0 + c1
            points (x0, y0, pch=19, col=cols [indx])
            lines (rep (x0, 2), c (0, y0), col=cols [indx], lty=3)
        }
    }

    # ********************** SUMMARIES OF ALL FOUR ANALYSES ********************
    cat (rep ("*", 50), "\n", sep="")
    cat (rep ("*", 18), " AVERAGE VALUES ", rep ("*", 18), "\n", sep="")
    cat (rep ("*", 50), "\n\n", sep="")

    # Then finally calculate statistics for average values:
    txts <- c ("H-CLUST", "K-MEANS")
    gout <- list (hc=gout.hc, km=gout.km)
    nout <- list (hc=nout.hc, km=nout.km)
    mout <- list (hc=mout.hc, km=mout.km)
    snm <- NULL
    for (i in 1:2) 
    {
        nm <- nout [[i]] / mout [[i]]
        tt <- t.test (nm - e1)
        tt2 <- t.test (nm - 2)
        cat ("*** ", txts [i], " Average N / M = ", 
             formatC (mean (nm, na.rm=TRUE), format="f", digits=2),
             " +/- ", formatC (sd (nm, na.rm=TRUE), format="f", digits=2),
             ": p (N/M=e) = ", formatC (tt$p.value, format="f", digits=4),
             ", p (N/M=2) = ", formatC (tt2$p.value, format="f", digits=4),
             "\n", sep="")

        sg <- 1 / log (e1 * gout [[i]]/ (e1 - 1))
        sn <- 1 / log (nout [[i]])
        sm <- 1 / (1 + log (mout [[i]]))
        sm [mout [[i]] == 1] <- s1
        gvals <- cbind (gout [[i]], nout [[i]], mout [[i]])
        svals <- cbind (sg, sn, sm)
        snm <- c (snm, sn, sm)
        txt <- c ("G", "N", "M")
        cat ("Values averaged over TO and FROM:\n")
        for (j in 1:3) 
        {
            cat (txt [j], " = ", formatC (mean (gvals [,j]), format="f", digits=3), 
                 " +/- ", formatC (sd (gvals [,j]), format="f", digits=3), 
                 "; s = ", formatC (mean (svals [,j], na.rm=TRUE), format="f", digits=3), 
                 " +/- ", formatC (sd (svals [,j], na.rm=TRUE), format="f", digits=3), 
                 "\n", sep="")
        }

        cat ("T-statistics for differences in s-values:\n")
        cat ("\t\t|\tStatistic\t|\tp-value\n")
        cat ("\t\t|\tN\tM\t|\tN\tM\n")
        cat (rep ("-", 60), "\n", sep="")
        vlist <- c ("G", "N")
        for (j in 1:2) 
        {
            cat ("\t", vlist [j], "\t|\t")
            if (j > 1) 
                cat ("\t")
            for (k in (j+1):3) 
            {
                tt <- t.test (svals [,j], svals [,k], var.equal=TRUE)
                cat (formatC (tt$statistic, format="f", digits=3))
                cat ("\t")
            }
            # Then p-values
            cat ("|")
            for (k in 1:j) 
                cat ("\t")
            for (k in (j+1):3) 
            {
                tt <- t.test (svals [,j], svals [,k], var.equal=TRUE)
                cat (formatC (tt$p.value, format="f", digits=3))
                if (k < 3) 
                    cat ("\t")
                else 
                    cat ("\n")
            }
        }
        cat ("\n")
    } # end for i over (hc, km)
    cat ("\n***** FINAL SORTING EFFICIENCY = ",
         formatC (mean (snm, na.rm=TRUE), format="f", digits=4), " +/- ",
         formatC (sd (snm, na.rm=TRUE), format="f", digits=4), "\n", sep="")

    #fname <- "fig_clust_sig.eps"
    #junk <- dev.print (device = postscript, file=fname,
    #    onefile = FALSE, width = fwd, height = fht, 
    #    paper = "special", horizontal = FALSE)
}
