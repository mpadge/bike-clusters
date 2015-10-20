#' get-members
#'
#' Creates clusters and then applies spatial constraint to reallocate stray
#' points to neighbouring clusters. "method" can be ward, complete, or k-means.
#' All methods other than k-means are passed to hclust.
#'
#' @param city nyc, washington, chicago, boston, london (case insensitive)
#' @param nc number of clusters
#' @param method = (ward, k-means, complete)
#' @param details prints statistics of cluster sizes
#' @param plot enables visual inspection of results
#' @return a vector of 2 columns for (from, to) data, each containing cluster
#' numbers for every bike station.

get.members <- function (city="nyc", nc=8, method="complete", details=FALSE,
                         plot=FALSE)
{
    require (tripack) # For Delaunay triangulation and neighbour lists

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
    
    cordir <- "/data/Dropbox/mark/analyses/bike-correlations/"
    wdd <- paste (cordir, "data/", sep="")
    wdr <- paste (cordir, "results/", sep="")

    if (city == "nyc" | city == "london" | city == "washingtondc")
    {
        fname <- paste (wdd, "station_latlons_", city, ".txt", sep="")
        lls <- c ("long", "lat")
    } else if (city == "chicago") {
        fname <- paste (wdd, "Divvy_Stations_2014-Q3Q4.csv", sep="")
        lls <- c ("longitude", "latitude")
    } else if (city == "boston") {
        fname <- paste (wdd, "hubway_stations.csv", sep="")
        lls <- c ("lng", "lat")
    }
    dat <- read.csv (fname, header=TRUE) 

    # Construct a neighbour list to spatially constrain clusters, as defined by
    # Delaunay triangulations
    xy <- data.frame (cbind (dat [lls [1]], dat [lls [2]]))
    names (xy) <- c ("x", "y")
    npts.full <- dim (xy) [1]

    if (method == "ward")
        method <- "ward.D" # or ward.D2? see ?hclust

    if (city == "washingtondc" | city == "london")
    {
        fnames <- paste (wdr, "R2_", city, "_from_all.csv", sep="")
        fnames <- c (fnames, paste (wdr, "R2_", city, "_to_all.csv", sep=""))
    } else {
        fnames <- paste (wdr, "R2_", city, "_from_all_00.csv", sep="")
        fnames <- c (fnames, paste (wdr, "R2_", city, "_to_all_00.csv", sep=""))
    }

    # Missing data for particular stations produce NAs in R2 matrices, so an
    # index of non-missing stations is made
    #nn <- as.numeric (row.names (unique (xy)))
    #ni <- which ((1:length (nn) - nn) != 0)[1] # index of duplicates
    r2 <- read.csv (fnames [1], header=FALSE) # doesn't matter which
    r2 <- as.matrix (r2)
    r2 [r2 < -1] <- NA
    r2 <- sign (r2) * sqrt (abs (r2))
    dmat.full <- 1 - r2
    indx <- which (rowSums (dmat.full, na.rm=TRUE) > 0)
    xy <- xy [indx,]
    npts <- length (indx)
    # Note that row.names (xy) are potentially discontinuous because of removal
    # of NA stations, but tri.mesh (xy) indexes directly by row numbers, so
    # values from this point on *DO NOT* index back into the original xy or r2
    # data.
    tm <- tri.mesh (xy)
    nbs <- neighbours (tm)

    membs <- list ()

    for (i in 1:2) {
        if (i == 2)
        {
            r2 <- read.csv (fnames [i], header=FALSE)
            r2 <- as.matrix (r2)
            r2 [r2 < -1] <- NA
        }
        r2 <- sign (r2) * sqrt (abs (r2))
        dmat.full <- 1 - r2
        indx <- which (rowSums (dmat.full, na.rm=TRUE) > 0)
        if (length (indx) != npts)
            stop ("ERROR: dim (", fnames [i], ") != ", npts, "\n", sep="")
        # Missing stations not in indx also have to be removed from nbs
        dmat <- dmat.full [indx, indx]
        membs [[i]] <- rep (NA, npts)

        if (method != "k-means") {
            hc <- hclust (as.dist (dmat), method=method)
        }

        nc.add <- nc.temp <- 0 # nc.temp is observed value after reallocation
        while (nc.temp < nc & (nc + nc.add) < length (nbs)) {
            if (method == "k-means") {
                km <- kmeans (as.dist (dmat), centers = nc + nc.add, 
                              iter.max = 20, nstart = 5)
                membs [[i]] <- as.vector (km$cluster)
            } else {
                membs [[i]] <- cutree (hc, k=nc + nc.add)
            }
            # Then constain membs to spatial neighbours only. First find points
            # where no others in cluster are neighbours
            non.nbs <- lapply (1:npts, function (x) {
                               cl.i <- membs [[i]] [nbs [[x]]]
                               if (!membs [[i]] [x] %in% cl.i) { nlist = x    }
                               else { nlist = NULL  }
                               return (nlist) })
            non.nbs <- unlist (non.nbs)

            # Then reallocate those points to the neighbourhing cluster with
            # minimal dmat. 
            d.min <- lapply (non.nbs, function (x) {
                             dnbs <- dmat [nbs [[x]], x]
                             di <- which.min (dnbs)
                             cbind (dnbs [di], membs [[i]] [nbs [[x]] [di]]) })

            # d.min is then a list of minimal distances and cluster IDS over all
            # neighbours of each point in non.nbs. The following extracts the
            # cluster IDs only: (7,66)
            new.clIDs <- unlist (lapply (d.min, function (x) x[,2]))
            # which can then be used to reassign the points:
            membs [[i]] [non.nbs] <- new.clIDs
            nc.temp <- length (unique (membs [[i]]))
            if (nc.temp < nc) { nc.add <- nc.add + 1    }
        } # end while nc.temp < nc

        if (nc.add > 0) { # Then renumber clusters to within nc0
            tb <- table (membs [[i]])
            membs.renum <- lapply (1:length (tb), function (x) {
                                   mindx <- which (membs [[i]] == names (tb) [x])
                                   if (length (mindx) > 0) cbind (mindx, x)  })
            renum.indx <- unlist (lapply (membs.renum, function (x) x [,1]))
            renum.indx <- as.vector (renum.indx)
            renum.ID <- unlist (lapply (membs.renum, function (x) x [,2]))
            renum.ID <- as.vector (renum.ID)
            membs [[i]] [renum.indx] <- renum.ID
        }

        # Next clause may be met in small cities (boston) for k-means, where
        # resultant numbers of clusters will not be precisely nc, and so nc.add
        # may be incremented extensively.
        if ((nc + nc.add) >= length (nbs))
            membs [[i]] <- rep (NA, length (nbs))
    } # end for i over (from, to)
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

    if (plot)
    {
        require (spatstat) # for convex hulls

        x11 (width=10)
        par (mfrow=c(1,2), mar=c(2,2,2,0.5), mgp=c(1.3,0.7,0), ps=10)
        cols <- rainbow (nc)
        mts <- c ("from", "to")

        for (i in 1:2)
        {
            plot (NULL, NULL, xlim=range (xy$x), ylim=range(xy$y), 
                  xlab="", ylab="", main=mts [i])
            for (j in 1:nc)
            {
                jindx <- which (membs [,i] == j)
                x <- xy$x [jindx]
                y <- xy$y [jindx]
                jindx <- which (!is.na (x) & !is.na (y))
                x <- x [jindx]
                y <- y [jindx]
                points (x, y, col=cols [j], pch=1)
                xyp <- ppp (x, y, xrange=range (x), yrange=range (y))
                ch <- convexhull (xyp)
                plot (ch, border=cols [j], add=TRUE)
            }
        }
    }

    # Then finally fill members out to original
    membs2 <- array (NA, dim=c(npts.full, 2))
    membs2 [indx,] <- membs
    membs <- membs2
    return (membs)
}
