# Contains the following functions:
#   1. clust.sig - the main plotting & analysis routine
#   2. get.members
#   3. get.partition.neighbours
#   4. num.clusts
#   5. calc.pnc
#   6. rescale.xm
#
# clust.sig uses the output of calc.pnc, which calculates the probability of
# observing a given peak height for a given number of clusters. This routine in turn
# uses the output of num.clusts, which is therefore the "primary" routine to be run
# prior to any others.
#
# num.clusts calculates the G-values and peak heights for a range of numbers of
# peaks from 2 to 15 for all 4 data series.

# *****************************************************************
# *************************   CLUST.SIG   *************************
# *****************************************************************

clust.sig <- function (method="complete", xmax=36)
{
    require (quantreg) # For upper bound regressions
    require (data.table) # For timetaken function

    st <- Sys.time ()
    cat (rep ("\n", 100))
    # NOTE: This seed is only to allow it to always generate appropriate k-means to
    # get through running the whole thing. Comment out to run it properly!!
    set.seed (18)
    # These are results for different random seeds for k-means clusters.
    # groups is contiguous / all. These are calculated by analysing only those
    # hierarchical levels for which the largest contiguous partition represents more
    # than half of all clusters.
    #
    # Direction = "to":
    # -----------------
    # seed  groups      N           M           sN              sM
    # -----------------------------------------------------------------
    #   1   22/24   7.33-3.51   3.33-1.53   0.551-0.155 0.483-0.104
    #   2   17/17   8.50-2.12   4.00-1.41   0.474-0.056 0.430-0.066
    #   6   22/22   7.33-2.08   3.33-1.15   0.529-0.089 0.476-0.099
    #   7   15/15   7.50-0.71   3.50-0.71   0.497-0.023 0.448-0.041
    #   8   13/16   6.50-0.71   3.00-0.00   0.536-0.031 0.477-0.000
    #   13  21/23   7.00-3.00   3.33-1.53   0.557-0.148 0.483-0.104
    #   15  13/13   6.50-0.71   3.00-0.00   0.536-0.031 0.477-0.000
    #   16  20/20   6.67-0.58   3.00-0.00   0.529-0.026 0.477-0.000
    #   17  15/17   7.50-0.71   3.50-0.71   0.497-0.023 0.448-0.041
    #   18  20/20   6.67-0.58   3.00-0.00   0.529-0.026 0.477-0.000
    # -----------------------------------------------------------------
    #   17.8/18.7   7.15-1.47   3.30-0.70   0.524-0.061 0.468-0.046
    # -----------------------------------------------------------------
    # The corresponding 10 values for M(contig/all) for the (4,8,11,15) clusters were
    # seed  |   M(c,4)   M(c,8) M(c,11) |   M(a,4)  M(a,8)  M(a,11)
    # ---------------------------------------------------------------
    #   1   |   3       2       5       |   3       3       5
    #   2   |   3       5               |   3       5
    #   6   |   4       2       4       |   4       2       4
    #   7   |   3       4               |   3       4
    #   8   |   3       1       3       |   3       2       4
    #   13  |   3       2       5       |   3       3       5
    #   15  |   3       3       2       |   3       3       4
    #   16  |   3       3       3       |   3       3       3
    #   17  |   3       1       4       |   3       2       5
    #   18  |   3       3       3       |   3       3       3
    # ---------------------------------------------------------------
    #       |   3.1     2.6     3.625   |   3.1     3       4.125
    # ---------------------------------------------------------------
    
    #
    # Direction = "from":
    # -------------------
    # seed  groups      N           M           sN              sM
    # -----------------------------------------------------------------
    #   1   31/39   6.20-1.48   3.00-0.71   0.566-0.092 0.488-0.063
    #   2   26/30   8.67-3.79   4.00-1.73   0.487-0.087 0.437-0.068
    #   6   24/24   8.00-2.65   3.67-1.15   0.496-0.072 0.445-0.054
    #   7   17/20   5.67-1.53   2.67-0.58   0.598-0.109 0.515-0.066
    #   8   24/26   8.00-3.61   3.67-2.08   0.513-0.109 0.475-0.116
    #   13  27/33   6.75-2.50   3.25-1.26   0.557-0.121 0.482-0.085
    #   15  29/31   7.25-1.50   3.25-0.96   0.516-0.070 0.476-0.081
    #   16  35/39   7.00-4.12   3.40-2.07   0.579-0.146 0.495-0.104
    #   17  33/33   8.25-2.63   3.75-0.96   0.489-0.066 0.439-0.046
    #   18  34/36   8.50-1.73   4.00-0.82   0.473-0.041 0.424-0.039
    # -----------------------------------------------------------------
    #     28/31.1   7.43-2.55   3.44-1.23   0.527-0.091 0.468-0.072
    # -----------------------------------------------------------------
    # The corresponding 10 values for M(contig/all) for the (4,8,11,14,20) clusters were
    # seed  |   M(c)            |   M(a)
    # --------------------------------------------
    #   1   |   (3,2,3,4,3)     |   (3,3,3,6,4)
    #   2   |   (3,1,6,3)       |   (3,2,6,5)
    #   6   |   (3,3,2,5)       |   (3,3,4,5)
    #   7   |   (3,3,2,2)       |   (3,3,3,4)
    #   8   |   (3,1,2,6,2)     |   (3,2,3,6,4)
    #   13  |   (3,2,3,5)       |   (3,3,3,7)
    #   15  |   (4,2,4,3)       |   (4,2,4,4)
    #   16  |   (3,3,2,7,2)     |   (3,3,3,7,3)
    #   17  |   (3,3,4,5,2)     |   (3,3,4,5,4)
    #   18  |   (3,1,4,5,4)     |   (3,2,4,5,5)
    # --------------------------------------------
    #       |   (3.1,2.1,3.2,4.5,2.6)   (3.1,2.6,3.7,5.4,4)
    # --------------------------------------------
    #
    # Total number of (contiguous, all) clusters for all four series 
    # = 56.8 / 68 = 83.5%


    if (method == "ward") {
        d.to <- read.csv ("./results/results_actual_ward_to.txt", 
                            sep=",", header=TRUE)
        d.from <- read.csv ("./results/results_actual_ward_from.txt", 
                            sep=",", header=TRUE)
    } else if (method == "complete") {
        d.to <- read.csv ("./results/results_actual_complete_to.txt", 
                            sep=",", header=TRUE)
        d.from <- read.csv ("./results/results_actual_complete_from.txt", 
                            sep=",", header=TRUE)
    } else if (method == "skater") {
        d.to <- read.csv ("./results/results_actual_skater_to.txt", 
                            sep=",", header=TRUE)
        d.from <- read.csv ("./results/results_actual_skater_from.txt", 
                            sep=",", header=TRUE)
    }
    d.to.kmeans <- read.csv ("./results/results_actual_kmeans_to.txt", 
                            sep=",", header=TRUE)
    d.from.kmeans <- read.csv ("./results/results_actual_kmeans_from.txt", 
                            sep=",", header=TRUE)
    d.neutral <- read.csv ("./results/results_neutral.txt", sep=",", header=TRUE)

    # The k-means values include 10 repeats of each value, which are first averaged:
    indx <- sort (unique (d.to.kmeans$nc))
    d.to.kmeans.mn <- d.to.kmeans [indx,] # Dummy matrix of right size.
    d.from.kmeans.mn <- d.from.kmeans [indx,]
    for (i in 1:length (indx)) {
        indx2 <- which (d.to.kmeans == indx [i])
        d.to.kmeans.mn [i, ] <- colMeans (d.to.kmeans [indx2,])
        d.from.kmeans.mn [i, ] <- colMeans (d.from.kmeans [indx2,])
    }
    d.to.kmeans <- d.to.kmeans.mn
    d.from.kmeans <- d.from.kmeans.mn

    # The neutral files just contain $(nsd,dmn,dsd), where the distances are total
    # intra-cluster distances which can be converted to proportions through dividing
    # by the total overall distances of:
    total.dist <- 26858.03
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
    y.resc <- ulbounds <- array (NA, dim=c(length (indx), 2))
    for (j in 1:2) { # Over (to, from) data
        dfr <- data.frame (x=xdat [indx], y=t0 [indx, j])
        for (k in 1:2) {
            mod <- nlrq (y ~ a * x ^ 2 + b * x + cc, data=dfr, 
                tau=ybounds [k], start=list(a=0, b=0, cc=mean(dfr$y)))
            ulbounds [,k] <- predict (mod, newdata=dfr$x)
        } # end for k
        y.resc [, j] <- 2 * (t0 [indx, j] - ulbounds [,2]) /
            (ulbounds [,1] - ulbounds [,2]) - 1
    } # end for j
    plot (xdat [indx], y.resc [,1], "l", col=cols [1], lwd=2, ylim=range (y.resc),
          ylab="Rescaled T-values", main="Rescaled T-values")
    for (i in 1:25) {
        lines (rep (i*4, 2), range (y.resc), col="gray", lty=2)
    }
    for (i in 1:2) {
        lines (xdat [indx], y.resc [,i], col=cols [i], lwd=2)
    }


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
    # T-statistics according to non-linear regressions of lower and upper bounds, so
    # that they lie between -1 and 1. Peak distinction is then measured as the
    # average difference to the 2 adjacent points. This measure thus has a
    # theoretical maximum of 2, and a minimum of 0.
    #
    # Values for s(M=1) in the loop that follows need adjusting, because the
    # calculated value is 1, yet any real value of M<1.5 will still round to 1 and
    # produce this value. These s values are therefore replaced with the
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
        for (j in 2:length (xdat [pks])) {
            nci <- c (xdat [pks] [j-1], xdat [pks] [j])
            nbs <- get.partition.neighbours (nci, dir=dirs [i], method=meth [i])
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
            if (nfrac >= 0.5 & mfrac >= 0.5) {
                part.sizes <- c (part.sizes, sum (unlist (nbs)))
                part.sizesC <- c (part.sizesC, sum (nbs [[1]]))
                mvals <- c (mvals, length (nbs [[1]]))
                cat ("| ", j, "\t", nci [1], "->", nci [2], "\t|\t", sum (nbs [[1]]), 
                     " / ", sum (unlist (nbs)), "\t", length (nbs [[1]]),
                     " / ", length (unlist (nbs)), sep="")
                if (j == 2) { cat ("\t\t\t\t\t\t\t\t\t|\n")    }
                else {
                    nm <- part.sizesC / mvals
                    # nm values for k-means are sometimes all identical, which
                    # causes the t-test to crash, as does the odd occasion when only
                    # a single value can be extracted before the fractions drop
                    # below 0.5:
                    if (sd (nm) == 0) {
                        stop ("\nERROR: Estimated values of N/M",
                             " are all identical---Just run again!\n")
                    } else if (length (nm) < 2) {
                        stop ("\nERROR: Insufficient values of N/M generated",
                              "---Just run again!\n")
                    } else {
                        tt <- t.test (nm - e1)
                        tt2 <- t.test (nm - 2)
                        cat ("\t", formatC (mean (nm, na.rm=TRUE), format="f", digits=2),
                            "+/-", formatC (sd (nm, na.rm=TRUE), format="f", digits=2),
                            "\t(", formatC (tt$statistic, format="f", digits=4),
                            ", ", formatC (tt$p.value, format="f", digits=4),
                            ")\t(", formatC (tt2$statistic, format="f", digits=4),
                            ", ", formatC (tt2$p.value, format="f", digits=4),
                            ")\t|\n", sep="")
                        if (i == 1) {
                            tstats [[1]] <- c (tstats [[1]], tt$statistic)
                            tstats [[2]] <- c (tstats [[2]], tt2$statistic)
                        } else if (i == 3) {
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
        if (meth [i] == "complete") {
            gout.hc <- c (gout.hc, gvals)
            nout.hc <- c (nout.hc, part.sizesC)
            mout.hc <- c (mout.hc, mvals)
        } else {
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
        for (j in 1:3) {
            cat ("|\t", vlist [j], "\t|\t")
            if (j > 1) { for (k in 1:(j-1)) { cat ("\t") }  }
            for (k in (j+1):4) {
                tt <- t.test (sarr [,j], sarr [,k], var.equal=TRUE)
                cat (formatC (tt$statistic, format="f", digits=3))
                cat ("\t")
            }
            # Then p-values
            cat ("|")
            for (k in 1:j) { cat ("\t") }
            for (k in (j+1):4) {
                tt <- t.test (sarr [,j], sarr [,k], var.equal=TRUE)
                cat (formatC (tt$p.value, format="f", digits=3))
                if (k < 4) { cat ("\t") }
                else { cat ("\t|\n")   }
            }
        }
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
    tstats <- lapply (tstats, function (x) abs (x))
    ylims <- lapply (tstats, function (x) range (abs (x)))
    ylims <- range (unlist (ylims))
    xmax <- lapply (tstats, function (x) length (x))
    xlims <- c (1, max (unlist (xmax)))
    cols <- c ("red", "red", "blue", "blue")
    ltys <- c (1, 2, 1, 2)
    plot (1:length (tstats [[1]]), tstats [[1]], "l", col=cols [1], lty=ltys [1],
          xlim=xlims, ylim=ylims, xlab="Number of peaks", ylab="T(M)")
    for (i in 1:4) {        
        lines (1:length (tstats [[i]]), tstats [[i]],
               col=cols [i], lty=ltys [i], lwd=2)
    }
    legend (xlims [1], ylims [2], lwd=2, col=cols, lty=ltys, bty="n",
            legend=c("TO(e)", "TO(2)", "FROM(e)", "FROM(2)"))

    for (i in 1:2) {
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
    for (i in 1:2) {
        nm <- nout [[i]] / mout [[i]]
        tt <- t.test (nm - e1)
        cat ("*** ", txts [i], " Average N / M = ", 
             formatC (mean (nm, na.rm=TRUE), format="f", digits=2),
             " +/- ", formatC (sd (nm, na.rm=TRUE), format="f", digits=2),
             ": p (N/M=e) = ", formatC (tt$p.value, format="f", digits=4),
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
        for (j in 1:3) {
            cat (txt [j], " = ", formatC (mean (gvals [,j]), format="f", digits=3), 
                 " +/- ", formatC (sd (gvals [,j]), format="f", digits=3), 
                 "; s = ", formatC (mean (svals [,j], na.rm=TRUE), format="f", digits=3), 
                 " +/- ", formatC (sd (svals [,j], na.rm=TRUE), format="f", digits=3), 
                 "\n", sep="")
        }
        nm <- gvals [,2] / gvals [,3]
        tt <- t.test (nm - e1)
        cat ("  ***** N / M = ", formatC (mean (nm), format="f", digits=3), " +/- ",
             formatC (sd (nm), format="f", digits=3), ": p(N/M=e) = ",
             formatC (tt$p.value, format="f", digits=4), "\n", sep="")

        cat ("T-statistics for differences in s-values:\n")
        cat ("\t\t|\tStatistic\t|\tp-value\n")
        cat ("\t\t|\tN\tM\t|\tN\tM\n")
        cat (rep ("-", 60), "\n", sep="")
        vlist <- c ("G", "N")
        for (j in 1:2) {
            cat ("\t", vlist [j], "\t|\t")
            if (j > 1) { cat ("\t") }
            for (k in (j+1):3) {
                tt <- t.test (svals [,j], svals [,k], var.equal=TRUE)
                cat (formatC (tt$statistic, format="f", digits=3))
                cat ("\t")
            }
            # Then p-values
            cat ("|")
            for (k in 1:j) { cat ("\t") }
            for (k in (j+1):3) {
                tt <- t.test (svals [,j], svals [,k], var.equal=TRUE)
                cat (formatC (tt$p.value, format="f", digits=3))
                if (k < 3) { cat ("\t") }
                else { cat ("\n")   }
            }
        }
    } # end for i over (hc, km)
    cat ("\n***** FINAL SORTING EFFICIENCY = ",
         formatC (mean (snm, na.rm=TRUE), format="f", digits=4), " +/- ",
         formatC (sd (snm, na.rm=TRUE), format="f", digits=4), "\n", sep="")

    st <- timetaken (st)
    cat ("Time taken = ", st, "\n", sep="")

    #fname <- "fig_clust_sig.eps"
    #junk <- dev.print (device = postscript, file=fname,
    #    onefile = FALSE, width = fwd, height = fht, 
    #    paper = "special", horizontal = FALSE)
}

#clust.sig ()


# *****************************************************************
# ***********************   GET.MEMBERS   *************************
# *****************************************************************


get.members <- function (nc=8, method="ward", details=FALSE)
{
    # Creates clusters and then applies spatial constraint to reallocate stray
    # points to neighbouring clusters. "method" can be ward, complete, or k-means.
    # All methods other than k-means are passed to hclust.
    require (tripack) # For Delaunay triangulation and neighbour lists

    fname <- "./results/results_st_latlons.txt"
    dat <- read.csv (fname, header=TRUE) 
    # Construct a neighbour list to spatially constrain clusters, as defined by
    # Delaunay triangulations
    xy <- data.frame (cbind (dat$lon, dat$lat))
    names (xy) <- c ("x", "y")
    npts <- nrow (xy)
    tm <- tri.mesh (xy)
    nbs <- neighbours (tm)

    membs <- list ()
    fnames <- c ("./results/results_r2from.txt", "./results/results_r2to.txt")
    for (i in 1:2) {
        r2 <- read.csv (fnames [i], header=FALSE)
        r2 <- as.matrix (r2)
        r2 [r2 < -1] <- NA
        r2 <- sign (r2) * sqrt (abs (r2))
        dmat <- 1 - r2
        if (method != "k-means") {
            hc <- hclust (as.dist (dmat), method=method)
        }
        nc.add <- nc.temp <- 0 # nc.temp is observed value after reallocation
        while (nc.temp < nc) {
            if (method == "k-means") {
                km <- kmeans (as.dist (dmat), centers = nc + nc.add, 
                              iter.max = 20, nstart = 5)
                membs [[i]] <- as.vector (km$cluster)
            } else {
                membs [[i]] <- cutree (hc, k=nc + nc.add)
            }
            # Then constain membs to spatial neighbours only. First find points where no
            # others in cluster are neighbours
            non.nbs <- lapply (1:npts, function (x) {
                               cl.i <- membs [[i]] [nbs [[x]]]
                               if (!membs [[i]] [x] %in% cl.i) { nlist = x    }
                               else { nlist = NULL  }
                               return (nlist) })
            non.nbs <- unlist (non.nbs)
            if (details) {
                cat ("reallocated ", length (non.nbs), " points.\n", sep="")    }
            # Then reallocate those points to the neighbourhing cluster with minimal dmat
            d.min <- lapply (non.nbs, function (x) {
                             dnbs <- dmat [nbs [[x]], x]
                             di <- which.min (dnbs)
                             cbind (dnbs [di], membs [[i]] [nbs [[x]] [di]]) })
            # d.min is then a list of minimal distances and cluster IDS over all neighbours of
            # each point in non.nbs. The following extracts the cluster IDs only:
            new.clIDs <- unlist (lapply (d.min, function (x) x[,2]))
            # which can then be used to reassign the points:
            membs [[i]] [non.nbs] <- new.clIDs
            nc.temp <- length (unique (membs [[i]]))
            if (nc.temp < nc) { nc.add <- nc.add + 1    }
        } # end while nc.temp < nc
        if (nc.add > 0) { # Then renumber clusters to within nc0
            tb <- table (membs [[i]])
            membs.renum <- lapply (1:length (tb), function (x) {
                                   indx <- which (membs [[i]] == names (tb) [x])
                                   if (length (indx) > 0) cbind (indx, x)  })
            renum.indx <- unlist (lapply (membs.renum, function (x) x [,1]))
            renum.indx <- as.vector (renum.indx)
            renum.ID <- unlist (lapply (membs.renum, function (x) x [,2]))
            renum.ID <- as.vector (renum.ID)
            membs [[i]] [renum.indx] <- renum.ID
        }
    }
    membs <- do.call (cbind, membs) 
    if (details) {
        hf <- table (membs [,1])
        ht <- table (membs [,2])
        nf <- length (unique (membs [,1]))
        nt <- length (unique (membs [,2]))
        cat ("Minimal cluster sizes (from, to) = (", min (hf), ", ", min (ht), 
             ") for clusters#(", which.min (hf), ", ", which.min (ht),
             "); actual numbers of clusters = (", nf, ", ", nt, ")\n", sep="")
    }
    return (membs)
}


# *****************************************************************
# ***********************   GET.CLUSTERS   ************************
# *****************************************************************

get.clusters <- function (method="ward")
{
    require (spatstat) # for ppp
    require (geosphere) # for areaPolygon

    max_clust_size <- 100
    r2from <- read.csv ("./results/results_r2from.txt", header=FALSE)
    r2from <- as.matrix (r2from)
    n <- dim (r2from) [1]
    clust.from <- clust.to <- array (NA, dim=c(n, max_clust_size))
    clust.diam.from <- clust.diam.to <- array (NA, dim=c(max_clust_size, max_clust_size))

    fname <- "./results/results_st_latlons.txt"
    dat <- read.csv (fname, header=TRUE) # has lats & lons in it
    pb <- txtProgressBar (style=3)
    ptm <- proc.time ()
    for (nc in 1:max_clust_size) {
        membs <- get.members (nc, method=method)
        clust.from [, nc] <- as.vector (membs [,1])
        clust.to [, nc] <- as.vector (membs [,2])

        # Then calculate diameters of clusters:
        for (i in 1:nc) {
            for (j in 1:2) {
                indx <- which (membs [,j] == i)
                if (length (indx) > 2) {
                    x <- dat$lon [indx]
                    y <- dat$lat [indx]
                    xy <- ppp (x, y, xrange=range (x), yrange=range (y))
                    ch <- convexhull (xy)
                    chb <- cbind (ch$bdry [[1]]$x, ch$bdry [[1]]$y)
                    chb <- rbind (chb, chb [1,])
                    diam <- 2 * sqrt (areaPolygon (chb)) / (1000 * pi)
                    if (j == 1) { clust.diam.from [nc, i] <- diam  }
                    else { clust.diam.to [nc, i] <- diam   }
                }
            } # end for j over (from, to)
        } # end for i over nc
        setTxtProgressBar (pb, nc/max_clust_size)
    } # end for nc
    close (pb)
    cat ("Calculation time = ", format (.POSIXct ((proc.time () - ptm) [3],
                                        tz="GMT"), "%H:%M:%S"), "\n", sep="")
    indx <- 2:max_clust_size
    clust.from <- clust.from [,indx]
    clust.to <- clust.to [,indx]
    clust.diam.from <- clust.diam.from [indx, ]
    clust.diam.to <- clust.diam.to [indx, ]
    fname <- paste ("./results/clust_from_members_", method, ".txt", sep="")
    write.table (clust.from, file=fname, sep=",", row.names=FALSE, col.names=FALSE)
    fname <- paste ("./results/clust_to_members_", method, ".txt", sep="")
    write.table (clust.to, file=fname, sep=",", row.names=FALSE, col.names=FALSE)
    fname <- paste ("./results/clust_from_diameters_", method, ".txt", sep="")
    write.table (clust.diam.from, file=fname, sep=",", row.names=FALSE,
                 col.names=FALSE)
    fname <- paste ("./results/clust_to_diameters_", method, ".txt", sep="")
    write.table (clust.diam.to, file=fname, sep=",", row.names=FALSE,
                 col.names=FALSE)
} # end function get.clusters

# *****************************************************************
# *****************   GET.PARTITION.NEIGHBOURS   ******************
# *****************************************************************

get.partition.neighbours <- function (nc=c(11,15), dir="from", method="complete",
                                      plot=FALSE) {
    # Identifies the clusters in min (nc) that split to form separate clusters in 
    # max (nc). For each cluster, a list element is returned with the number of
    # groups into which each connected group is subsequently partitioned.
    #
    # Plot enables visual inspection of the spatial arrangements.
    require (spatstat)
    nlo <- min (nc)
    nhi <- max (nc)
    membs.lo <- get.members (nlo, method=method)
    membs.hi <- get.members (nhi, method=method)
    if (dir == "from") { diri <- 1  }
    else { diri <- 2    }

    if (plot) {
        fname <- "./results/results_st_latlons.txt"
        dat <- read.csv (fname, header=TRUE) 
        xy <- data.frame (cbind (dat$lon, dat$lat))
        names (xy) <- c ("x", "y")
        xlims <- range (xy$x)
        ylims <- range (xy$y)
        x11 (width=10, height=5)
        par (mfrow=c(1,2), mar=c(0,0,0,0))
        cols <- c ("red", "orange", "yellow", "lawngreen", "blue", "purple", "magenta",
                   "gray")
        while (length (cols) < max (nc)) { cols <- c (cols, cols)   }
        plot (NULL, NULL, xlim=xlims, ylim=ylims,
              xaxt="n", yaxt="n", xlab="", ylab="", frame=FALSE)
        for (i in 1:nc [1]) {
            indx <- which (membs.lo [,diri] == i)
            xyp <- ppp (xy$x [indx], xy$y [indx], 
                        xrange=range (xy$x [indx]), yrange=range (xy$y [indx]))
            ch <- convexhull (xyp)
            plot (ch, border=cols [i], add=TRUE)
            text (mean (xy$x [indx]), mean (xy$y [indx]), labels=i)
        }
        plot (NULL, NULL, xlim=xlims, ylim=ylims,
              xaxt="n", yaxt="n", xlab="", ylab="", frame=FALSE)
        for (i in 1:nc [1]) {
            indx <- which (membs.lo [,diri] == i)
            xyp <- ppp (xy$x [indx], xy$y [indx], 
                        xrange=range (xy$x [indx]), yrange=range (xy$y [indx]))
            ch <- convexhull (xyp)
            plot (ch, border=cols [i], add=TRUE)
        }
        for (i in 1:nc [2]) {
            indx <- which (membs.hi [,diri] == i)
            xyp <- ppp (xy$x [indx], xy$y [indx], 
                        xrange=range (xy$x [indx]), yrange=range (xy$y [indx]))
            ch <- convexhull (xyp)
            plot (ch, border=cols [i], add=TRUE, lty=2, lwd=2)
            text (mean (xy$x [indx]), mean (xy$y [indx]), labels=i)
        }
    } # end if plot


    indx.lo <- lapply (1:nlo, function (x) as.vector (which (membs.lo [,diri] == x)))
    indx.hi <- lapply (1:nhi, function (x) as.vector (which (membs.hi [,diri] == x)))
    # Find which groups in indx.lo contain the groups in indx.hi, and add up how
    # many of the latter belong in each of the former:
    indx <- lapply (indx.hi, function (i) {
                    indx2 <- lapply (indx.lo, function (j) {
                                     length (which (j %in% i))  })
                    indx2 <- unlist (indx2)
                    which.max (indx2)   })
    indx <- unlist (indx)
    # indx then maps each of the nc [2] clusters onto the corresponding nc [1].
    # It has a length of nc [2], and ranges from 1 to nc [1].
    tindx <- as.vector (table (indx)) # has a length of nc [1]
    tindx2 <- which (tindx > 1)
    tindx <- tindx [tindx2]
    # tindx2 now indexes ththe original nc1 clusters that subsequently became
    # divided in forming the nc [2], while tindx holds the actual numbers of
    # clusters into which these become partitioned. All that remains is to discern
    # whether these clusters from within nc [1] are neighbours.

    require (tripack) 
    fname <- paste ("./results/results_st_latlons.txt")
    dat <- read.csv (fname, header=TRUE) 
    # Construct a neighbour list to spatially constrain clusters, as defined by
    # Delaunay triangulations
    xy <- data.frame (cbind (dat$lon, dat$lat))
    names (xy) <- c ("x", "y")
    npts <- nrow (xy)
    tm <- tri.mesh (xy)
    nbs <- neighbours (tm)

    nbmat <- array (FALSE, dim=c(npts, npts))
    for (i in 1:length (nbs)) {
        nbmat [i, nbs [[i]]] <- nbmat [nbs [[i]], i] <- TRUE    }

    nbmat.groups <- array (FALSE, dim=c(nc [1], nc [1]))
    # ngroups will always be small, so a loop is used here because it's easier to
    # interpret than an lapply.
    for (i in 1:(nc [1] - 1)) {
        indx.i <- indx.lo [[i]]
        for (j in (i + 1):nc [1]) {
            indx.j <- indx.lo [[j]]
            nbmat.sub <- nbmat [indx.i, indx.j]
            if (length (which (nbmat.sub)) > 0) {
                nbmat.groups [i, j] <- nbmat.groups [j, i] <- TRUE
            }
        }
    }
    # Then reduce nbmat.groups to only those groups which become partitioned:
    nbmat.groups <- nbmat.groups [tindx2, tindx2]

    # Then finally extract clusters of neighbours
    require (igraph)
    nbmat.graph <- graph.adjacency (nbmat.groups, mode="undirected")
    nbmat.cl <- clusters (nbmat.graph)
    # Clusters has $membership IDs, $csize, and $no, where the latter is simply the
    # number of spatially separate clusters. The following results has one list
    # element for each distinct cluster. Each list element holds the number of
    # groups into which each neighbouring larger group (that is, each member of nlo)
    # is partitioned. The list is constructed with the largest group first, although
    # the size of a group just depends here on the number of contiguous clusters
    # within nc [1], and not the total number into which that contiguous group is
    # subsequently divided in nc [2].
    cindx <- order (nbmat.cl$csize, decreasing=TRUE)
    results <- lapply (1:length (cindx), function (i) {
                       membs <- which (nbmat.cl$membership == cindx [i])
                       tindx [membs]    })
    return (results)
} # end function get.partition.neighbours



# *****************************************************************
# ************************   NUM.CLUSTS   *************************
# *****************************************************************

num.clusts <- function (plot=FALSE, method="complete")
{
    # Calculates peak height and G-values as a function of numbers of clusters, to
    # determine the number of clusters that should be used for each of the four data
    # series. For plot=TRUE, it also reads the table of the same values (from
    # "results_prob_m.txt") to plot the probabilities of the observed peak heights,
    # as calcualted from calc.pnc (). (The latter routine itself requires, however,
    # the initial output of num.clusts.)
    #
    # These p-values are then used to determine the appropriate number of (most
    # highly significant) peaks to be used for analysing each series in the main
    # clust.sig() routine.
    require (quantreg) # For upper bound regressions
    
    if (method == "ward") {
        fname <- "./results/results_actual_ward_to.txt"
        d.to <- read.csv (fname, sep=",", header=TRUE)
        fname <- "./results/results_actual_ward_from.txt"
        d.from <- read.csv (fname, sep=",", header=TRUE)
    } else if (method == "complete") {
        fname <- "./results/results_actual_complete_to.txt"
        d.to <- read.csv (fname, sep=",", header=TRUE)
        fname <- "./results/results_actual_complete_from.txt"
        d.from <- read.csv (fname, sep=",", header=TRUE)
    } else if (method == "skater") {
        fname <- "./results/results_actual_skater_to.txt"
        d.to <- read.csv (fname, sep=",", header=TRUE)
        fname <- "./results/results_actual_skater_from.txt"
        d.from <- read.csv (fname, sep=",", header=TRUE)
    } 
    fname <- "./results/results_actual_kmeans_to.txt"
    d.to.kmeans <- read.csv (fname, sep=",", header=TRUE)
    fname <- "./results/results_actual_kmeans_from.txt"
    d.from.kmeans <- read.csv (fname, sep=",", header=TRUE)
    fname <- "./results/results_neutral.txt"
    d.neutral <- read.csv (fname, sep=",", header=TRUE)
    
    # The k-means values include 10 repeats of each value, which are first averaged:
    indx <- sort (unique (d.to.kmeans$nc))
    d.to.kmeans.mn <- d.to.kmeans [indx,] # Dummy matrix of right size.
    d.from.kmeans.mn <- d.from.kmeans [indx,]
    for (i in 1:length (indx)) {
        indx2 <- which (d.to.kmeans == indx [i])
        d.to.kmeans.mn [i, ] <- colMeans (d.to.kmeans [indx2,])
        d.from.kmeans.mn [i, ] <- colMeans (d.from.kmeans [indx2,])
    }
    d.to.kmeans <- d.to.kmeans.mn
    d.from.kmeans <- d.from.kmeans.mn
    nc <- d.neutral$nc
    
    d0 <- cbind (d.to$d.in, d.from$d.in)
    t0 <- (d0 - d.neutral$dmn) / (d.neutral$dsd / sqrt (nc))
    d0 <- cbind (d.to.kmeans$d.in, d.from.kmeans$d.in)
    tk <- (d0 - d.neutral$dmn) / (d.neutral$dsd / sqrt (nc))
    
    ttxt <- c (" to (hc) ", "to (km)", "from (hc)", "from (km)")
    tvals <- cbind (t0 [,1], tk [,1], t0 [,2], tk [,2])
    np.lim <- (2:15)
    gvals <- hvals <- nmax <- array (NA, dim=c(length (np.lim), 4))
    rescale <- 2 # If 3, then O(3) bounds are calculated
    ybounds <- c (0.99, 0.01) # Upper and lower bounds for nlqr regressions
    
    for (i in 1:4) { # Over (to, from) data
        pks <- which (diff (sign (diff (tvals [,i]))) == -2) + 1

        # The following could be done as lapply, but time is not as issue, and the
        # function would end up very complex and hard to read.
        for (j in 1:length (np.lim)) { # loop over number of peaks
            pks.j <- pks [1:min (np.lim [j], length (pks))] 
            # pks.j are then the positions of the first np.lim[j] pks, with the
            # following line making a continuous indx.
            indx <- 1:max (pks.j)
            ulbounds <- array (NA, dim=c(length (indx), 2))
            
            dfr <- data.frame (x=nc [indx], y=tvals [indx, i])
            for (k in 1:2) {
                if (rescale == 2) {
                    mod <- nlrq (y ~ a * x ^ 2 + b * x + cc, data=dfr, 
                                 tau=ybounds [k], start=list(a=0, b=0, cc=mean(dfr$y)))
                } else if (rescale == 3) {
                    mod <- nlrq (y ~ a * x ^ 3 + b * x ^ 2 + cc * x + dd, data=dfr, 
                                 tau=ybounds [k], start=list(a=0, b=0, cc=0, dd=mean(dfr$y)))
                }
                ulbounds [,k] <- predict (mod, newdata=dfr$x)
            } # end for j
            dfr$y <- 2 * (dfr$y - ulbounds [,2]) / (ulbounds [,1] - ulbounds [,2]) - 1
            
            nmax [j, i] <- max (nc [pks.j]) + 1
            # The first of the following lines calculates G-values including the
            # difference between the first peak and one. While this is strictly
            # correct, it is a stricter requirement than what is placed upon the
            # simulated series used to calculate the probabilities. The actual
            # calculation is therefore based only on the mean inter-peak spacing,
            # without presuming that partitioning starts at one.
            #gvals [j, i] <- mean (c (diff (nc [pks.j]), nc [pks.j] [1] - 1))
            gvals [j, i] <- mean (c (diff (nc [pks.j])))
            dfr <- dfr [1:max (pks.j), ]
            non.pks <- which (!1:length (dfr$y) %in% pks)
            hvals [j, i] <- mean (dfr$y [pks.j]) - mean (dfr$y [non.pks])
        } # end for j
    } # end for i
    
    if (plot) {
        dat <- read.csv ("./results/results_prob_m.txt", header=TRUE)
        ydat <- list ()
        ydat [[1]] <- gvals
        ydat [[2]] <- hvals
        ydat [[3]] <- dat [,14:17] # The simulated probabilities
        cols <- c ("red", "red", "blue", "blue")
        ltys <- c (1, 2, 1, 2)
        ylabs <- c ("G", "Peak height", "Prob (pk. height)")
        
        x11 (width=14, height=5)
        par (mfrow=c(1,3)) 
        par (mar=c(2.5,2.5,2,1), mgp=c(1.3,0.4,0), ps=10, tcl=-0.2)
        for (i in 1:3) {
            ylims <- range (ydat [[i]])
            #if (i == 3) { ylims [2] <- 0.1 }
            plot (np.lim, ydat [[i]] [,1], "l", col=cols [1], lty=ltys [1], 
                  ylim=ylims, xlab="Number of peaks", ylab=ylabs [i])
            for (j in 1:4) {
                lines (np.lim, ydat [[i]] [,j], col=cols [j], lty=ltys [j])
                points (np.lim, ydat [[i]] [,j], col=cols [j], pch=19)
            }
            if (i == 3) { # Highlight minimal-probability points
                mini <- apply (ydat [[i]], 2, which.min)
                # Adjust from-hclust to exclude first point:
                mini [3] <- 1 + which.min (ydat [[i]] [2:dim (ydat [[i]])[1], 3])
                mini.indx <- (0:3) * dim (ydat [[i]]) [1] + mini
                ydat.mini <- as.matrix (ydat [[i]]) [mini.indx]
                points (np.lim [mini], ydat.mini, col=cols, pch=19, cex=1.5)
                points (np.lim [mini], ydat.mini, col=cols, pch=1, cex=2.5)
            }
        } # end for i
        ypos <- ylims [1] + 1 * diff (ylims) # 0.88 with MI plots
        xpos <- min (np.lim) + c (0.1, 0.5) * diff (range (np.lim))
        cols <- cols [c(1,3)]
        legend (xpos [1], ypos, lwd=1, col=cols, bty="n",
                legend=c("to (h-clust)", "from (h-clust)"))
        legend (xpos [2], ypos, lwd=1, col=cols, bty="n", lty=2,
                legend=c("to (k-means)", "from (k-means)"))
    } else { # do not plot
        dat <- list ()
        dat [[1]] <- ttxt
        dat [[2]] <- data.frame (cbind (np.lim, gvals, hvals, nmax))
        names (dat [[2]]) <- c ("nc", "g1", "g2", "g3", "g4", "h1", "h2", "h3", "h4",
                                "n1", "n2", "n3", "n4")
        return (dat)
    } # end else not plot
} # end function num.clusts


# *****************************************************************
# **********************   GET.NUM.CLUSTS   ***********************
# *****************************************************************

get.num.clusts <- function ()
{
    # Just loads up the output of calc.pnc(), and returns the final numbers of
    # clusters corresponding to minimal joint probabilities of peak spacings and
    # depths.
    dat <- read.csv ("./results/results_prob_m.txt", header=TRUE)
    nc <- dat$nc
    dat <- dat [,14:17] # The simulated probabilities
    mini <- apply (dat, 2, which.min)
    # Adjust from-hclust to exclude first point:
    mini [3] <- 1 + which.min (dat [2:dim (dat)[1], 3])
    mini.indx <- (0:3) * dim (dat) [1] + mini
    dat.mini <- as.matrix (dat) [mini.indx]

    results <- data.frame (cbind (nc [mini], dat.mini))
    names (results) <- c ("num.clusts", "pmin")
    rownames (results) <- c ("to-h", "to-k", "from-h", "from-k")
    return (results)
} # end function get. num.clusts


# *****************************************************************
# *************************   CALC.PNC   **************************
# *****************************************************************

calc.pnc <- function (dat=NA, nrpts=100)
{
    # Calculates probabilities of peak heights as returned from the num.clusts ()
    # function.  If is.na (dat), then the table is called internally. To generate
    # useful statistics, this function should be run with a VERY large number of
    # repeats (for example, 100,000), which will take a LONG time!
    require (data.table) # For timetaken function
    require (quantreg) 
    st0 <- Sys.time ()
    st.loop <- st0
    nrpts <- nrpts * 2 
    # Presumes M \approx 4.5, so about half of all simulated series will have the
    # right M-values, therefore nrpts is doubled.
    rescale <- 2 # Quadratic upper and lower boundaries
    
    if (is.na(dat)) {
        dat <- num.clusts (plot=FALSE)
    }
    dat <- dat [[2]] 
    # [[1]] is just the names of the four groups (to-h, to-k, from-h, from-k)
    # dat has 13 columns of $nc (number of peaks) plus 4 groups of G-values, peak
    # heights ($h), and lengths of series. The latter four are the position of the
    # final peak plus one. This is the length needed to feed in to the random
    # simulations. There are some differences between the four data sets in $n
    # values, but an average value is taken for all, as specified below.
    #
    # The whole thing is run here as a loop over the 14 values of $nc = (2:15),
    # because the actual time consuming bit is the processing internal to this loop,
    # and it makes the code clearer.
    pr.pk.hts <- num.samples <- array (NA, dim=c(length (dat$nc), 4))
    
    for (i in 1:length (dat$nc)) {
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
        # If observed G-values are the same, run all together, else run individually
        if (sum (diff (gvals0))  == 0) { 
            g <- gvals0 [1]
            indx <- which (gvals >= g)
            num.samples [i, ] <- rep (length (indx), 4)
            if (length (indx > 1)) {
                sm <- simvals [, indx]
                sm <- rescale.xm (sm, rescale=rescale)
                sml <- as.list (as.data.frame (sm))
                prob.g <- length (indx) / nrpts
            
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
        } else { # G-values differ, so each is calculated separately
            hvals0 <- as.numeric (dat [i, 6:9])
            gh0 <- as.list (as.data.frame (rbind (gvals0, hvals0)))
            # gh0 is a list of four pairs of [g, h], so x [1] = g & x [2] = h.
            indx <- which (gvals >= min (gvals0))
            gvals <- gvals [indx]
            simvals <- simvals [, indx]
            simvals <- rescale.xm (simvals, rescale=rescale)
            
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
        
        cat ("[", i, "]: p = [", formatC (pr.pk.hts [i,], format="f", digits=4), 
             "] after (loop, total) = (", timetaken (st.loop), ", ", 
             timetaken (st0), ")\n")
        st.loop <- Sys.time ()
    } # end for i over the nc values
    st <- timetaken (st0)
    cat ("Total calculation time = ", st, "\n", sep="")
    
    nrpts <- rep (nrpts, dim (dat) [1])
    dat <- cbind (dat, pr.pk.hts, nrpts, num.samples)
    names (dat) [14:22] <- c ("p1", "p2", "p3", "p4", "nrpts", 
                              "ns1", "ns2", "ns3", "ns4")
    write.table (dat, file="results_prob_m.txt", row.names=FALSE, sep=",")
    return (dat)
} # end function calc.pnc


# *****************************************************************
# ************************   RESCALE.XM  **************************
# *****************************************************************

rescale.xm <- function (xm, rescale=2) {
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

