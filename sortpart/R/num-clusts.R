#' num.clusts
#'
#' Calculates peak height and G-values as a function of numbers of clusters, to
#' determine the number of clusters that should be used for each of the four
#' data series. For plot=TRUE, it also reads the table of the same values (from
#' "results_prob_m.txt") to plot the probabilities of the observed peak heights,
#' as calcualted from calc.pnc (). (The latter routine itself requires, however,
#' the initial output of num.clusts.)
#'
#' These p-values are then used to determine the appropriate number of (most
#' highly significant) peaks to be used for analysing each series in the main
#' clust.sig() routine.
#'
#' @param city nyc, washington, chicago, boston, london (case insensitive)
#' @param method complete, ward, or skater, to be compared to k-means
#' @return data.frame

num.clusts <- function (city="nyc", plot=FALSE, method="complete")
{
    require (quantreg) # For upper bound regressions

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
    
    wd0 <- getwd ()
    count <- 0
    wd <- ""
    while (!"bike-clusters" %in% list.files ("."))
    {
        wd <- paste (wd, "../", sep="")
        setwd ("../")
        count <- count + 1
    }
    setwd (wd0)
    wd <- paste (wd, "bike-clusters/results/", sep="")
    wdk <- paste (wd, "k-means-clusters/", sep="")

    fname <- paste (wd, city, "-results-actual-to-", method, ".txt", sep="")
    d.to <- read.csv (fname, sep="," , header=TRUE)
    fname <- paste (wd, city, "-results-actual-from-", method, ".txt", sep="")
    d.from <- read.csv (fname, sep="," , header=TRUE)

    fname <- paste (wd, city, "-results-actual-to-k-means.txt", sep="")
    d.to.kmeans <- read.csv (fname, sep=",", header=TRUE)
    fname <- paste (wd, city, "-results-actual-from-k-means.txt", sep="")
    d.from.kmeans <- read.csv (fname, sep=",", header=TRUE)
    fname <- paste (wd, city, "-results-neutral.txt", sep="")
    d.neutral <- read.csv (fname, sep=",", header=TRUE)
    nc <- d.neutral$nc
    
    d0 <- cbind (d.to$d.in, d.from$d.in)
    t0 <- (d0 - d.neutral$dmn) / (d.neutral$dsd / sqrt (nc))
    d0 <- cbind (d.to.kmeans$d.in, d.from.kmeans$d.in)
    tk <- (d0 - d.neutral$dmn) / (d.neutral$dsd / sqrt (nc))
    
    ttxt <- c (" to (hc) ", "to (km)", "from (hc)", "from (km)")
    tvals <- cbind (t0 [,1], tk [,1], t0 [,2], tk [,2])

    # ****** Maximum Number of Peaks is set here to 15 *****
    np.lim <- (2:15)

    gvals <- hvals <- nmax <- array (NA, dim=c(length (np.lim), 4))
    rescale <- 2 # If 3, then O(3) bounds are calculated
    ybounds <- c (0.99, 0.01) # Upper and lower bounds for nlqr regressions
    
    for (i in 1:4) { # to-hc, to-km, from-hc, from-km
        pks <- which (diff (sign (diff (tvals [,i]))) == -2) + 1

        # The following is quick, so is done as a loop over number of peaks,
        # rather than lapply
        for (j in 1:length (np.lim)) { 
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
        p0 <- 0.05 # plots a horizontal reference line
        fname <- paste ("./results/", city, "-results-prob-m-", method,
                        ".txt", sep="")
        if (file.exists (fname))
            dat <- read.csv (fname, header=TRUE)
        else
            dat <- calc.pnc (city=city, method=method)
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
                lines (range (np.lim), rep (p0, 2), col="gray", lwd=2, lty=2)
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

