#' compare.dists
#'
#' Calculates mean ratio of observed to expected intra-cluster distances
#' calculated out to given number of clusters.
#'
#' @param city nyc, washington, chicago, boston, london (case insensitive)
#' @param method = (ward, k-means, complete)
#' @param xmax maximum number of clusters to be analysed and plotted
#' @return data.frame

compare.dists <- function (city="nyc", method="complete", xmax=36,
                           ycut=c(1.2,1.3))
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

    nc <- d.neutral$nc

    d.in <- cbind (d.to$d.in / d.neutral$dmn, 
        d.from$d.in / d.neutral$dmn)
    d.in.kmeans <- cbind (d.to.kmeans$d.in / d.neutral$dmn, 
        d.from.kmeans$d.in / d.neutral$dmn)

    indx <- which (nc <= xmax)
    d.in <- d.in [indx,]
    d.in.kmeans <- d.in.kmeans [indx,]
    nc <- nc [indx]
    d.in <- apply (d.in, 2, cumsum) / cbind (nc - 1, nc - 1)
    d.in.kmeans <- apply (d.in.kmeans, 2, cumsum) / cbind (nc - 1, nc - 1)

    ydat <- cbind (d.in, d.in.kmeans)
    dirs <- rep (c ("to", "from"), 2)
    meths <- c (rep (method, 2), rep ("k-means", 2))
    cols <- rep (c ("red", "blue"), 2)
    ltys <- c (1, 1, 2, 2)

    x11 ()
    par (mar=c(2.5,2.5,2,1), mgp=c(1.3,0.4,0), ps=10, tcl=-0.2)

    ylims <- range (c (d.in, d.in.kmeans))
    legxpos <- nc [1] + c (0.05, 0.4) * diff (range (nc))
    legypos <- ylims [1] + 0.1 * diff (ylims)

    plot (nc, ydat [,1], "l", col=cols [1], lwd=1, 
          ylim=ylims, xlab="", xaxt="n", yaxt="n",
          ylab="Proportional distance ridden within clusters")
    axis (1, at=1:10 * 4)
    axis (2, at=pretty (ylims))
    for (i in 1:4)
        lines (nc, ydat [,i], col=cols [i], lty=ltys [i])
    title (xlab="Number of Clusters")

    legend (legxpos [1], legypos, lwd=1, col=cols, bty="n",
        legend=c("to (h-clust)", "from (h-clust)"))
    legend (legxpos [2], legypos, lwd=1, col=cols, bty="n", lty=2,
        legend=c("to (k-means)", "from (k-means)"))

    # Then lines that intercept at ycut
    for (i in 1:length (ycut))
    {
        xmax <- apply (ydat, 2, function (x) {
                       if (max (x) < ycut [i])
                           NA
                       else
                           max (which (x > ycut [i]))
        })
        for (j in 1:4)
            if (!is.na (xmax [j]) & xmax [j] < length (nc))
            {
                y0 <- ydat [xmax [j], j]
                lines (rep (nc [xmax [j]], 2), c (0, y0), col=cols [j], lty=2)
                lines (c (0, nc [xmax [j]]), rep (y0, 2), col=cols [j], lty=2)
                points (nc [xmax [j]], y0, pch=19, col=cols [j])
                y0 <- formatC (y0, format="f", digits=3)
                cat ("mean dists (", dirs [j], ", ", meths [j], ") d > ",
                     y0, " for nc <= ", nc [xmax [j]], "\n", sep="")
            }
    }
}
