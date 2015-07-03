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
#' @return a vector of 2 columns for (from, to) data, each containing cluster
#' numbers for every bike station.

get.members <- function (city="nyc", nc=8, method="ward", details=FALSE)
{
    require (tripack) # For Delaunay triangulation and neighbour lists

    if (tolower (substring (city, 1, 2)) == "ny")
        city <- "nyc"
    else if (tolower (substring (city, 1, 2)) == "lo")
        city <- "london"
    else if (tolower (substring (city, 1, 2)) == "wa" |
                tolower (substring (city, 1, 2)) == "dc")
        city <- "washingtondc"
    else if (tolower (substring (city, 1, 2)) == "ch")
        city <- "chicago"
    else if (tolower (substring (city, 1, 2)) == "bo")
        city <- "boston"
    else 
        stop ("city not valid")

    wd0 <- getwd ()
    count <- 0
    wdd <- wdr <- ""
    while (!"bike-correlations" %in% list.files ("."))
    {
        wdd <- paste (wdd, "../", sep="")
        wdr <- paste (wdr, "../", sep="")
        setwd ("../")
        count <- count + 1
    }
    setwd (wd0)
    wdd <- paste (wdd, "bike-correlations/data/", sep="")
    wdr <- paste (wdr, "bike-correlations/results/", sep="")

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
    # London has a duplicate xy for some reason
    nn <- as.numeric (row.names (unique (xy)))
    ni <- which ((1:length (nn) - nn) != 0)[1] # index of duplicates
    xy <- unique (xy) 
    names (xy) <- c ("x", "y")
    npts <- nrow (xy)
    tm <- tri.mesh (xy)
    nbs <- neighbours (tm)

    if (method == "ward")
        method <- "ward.D" # or ward.D2? see ?hclust

    membs <- list ()

    if (city == "washingtondc" | city == "london")
    {
        fnames <- paste (wdr, "R2_", city, "_from_all.csv", sep="")
        fnames <- c (fnames, paste (wdr, "R2_", city, "_to_all.csv", sep=""))
    } else {
        fnames <- paste (wdr, "R2_", city, "_from_all_00.csv", sep="")
        fnames <- c (fnames, paste (wdr, "R2_", city, "_to_all_00.csv", sep=""))
    }

    for (i in 1:2) {
        r2 <- read.csv (fnames [i], header=FALSE)
        r2 <- as.matrix (r2)
        r2 [r2 < -1] <- NA
        r2 <- sign (r2) * sqrt (abs (r2))
        dmat.full <- 1 - r2
        indx <- which (rowSums (dmat.full, na.rm=TRUE) > 0)
        # r2 matrices sometimes have entire missing columns, which are
        # ultimately returned with NA cluster memberships. These missing columns
        # also have to be removed from nbs. Also, the duplicate xy for London
        # has to be removed from the indx
        indx <- indx [which (!indx %in% ni)]
        dmat <- dmat.full [indx, indx]
        membs [[i]] <- rep (NA, npts)
        nbs2 <- lapply (nbs, function (i) i [which (i %in% indx)])

        if (method != "k-means") {
            hc <- hclust (as.dist (dmat), method=method)
        }

        nc.add <- nc.temp <- 0 # nc.temp is observed value after reallocation
        while (nc.temp < nc) {
            if (method == "k-means") {
                km <- kmeans (as.dist (dmat), centers = nc + nc.add, 
                              iter.max = 20, nstart = 5)
                membs [[i]] [indx] <- as.vector (km$cluster)
            } else {
                membs [[i]] [indx] <- cutree (hc, k=nc + nc.add)
            }
            # Then constain membs to spatial neighbours only. First find points
            # where no others in cluster are neighbours
            non.nbs <- lapply (1:npts, function (x) {
                               cl.i <- membs [[i]] [nbs2 [[x]]]
                               if (!membs [[i]] [x] %in% cl.i) { nlist = x    }
                               else { nlist = NULL  }
                               return (nlist) })
            non.nbs <- unlist (non.nbs)
            non.nbs <- non.nbs [which (non.nbs %in% indx)]

            # Then reallocate those points to the neighbourhing cluster with
            # minimal dmat. nbs2 maps on to the original npts structure, and so
            # must be related to dmat.full
            d.min <- lapply (non.nbs, function (x) {
                             dnbs <- dmat.full [nbs2 [[x]], x]
                             di <- which.min (dnbs)
                             cbind (dnbs [di], membs [[i]] [nbs2 [[x]] [di]]) })

            # d.min is then a list of minimal distances and cluster IDS over all
            # neighbours of each point in non.nbs. The following extracts the
            # cluster IDs only:
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
