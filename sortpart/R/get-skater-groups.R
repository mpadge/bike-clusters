#' get.skater.groups
#'
#' Produces .csv files with cluster memberships across a range of numbers of
#' clusters (1-max.groups=50). 
#'
#' @param city nyc, washington, chicago, boston, london (case insensitive)
#' @param dir = ("from", "to")
#' @param max.groups only seems to work with skater up to around 50
#' @return nothing (dumps files to dir)

get.skater.groups <- function (city="nyc", dirf="from", max.groups=50)
{
    require (sp)
    require (spdep)
    require (tripack)

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

    if (dirf != "from")
        dirf <- "to"

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

    # First read in correlations used to weight connections in the mstree. The
    # skater routine constructs *MINIMAL* variance sub-trees, so the data
    # submitted to skater thus have to be distances, rather than correlations.
    if (city == "washingtondc" | city == "london")
        fname <- paste (wdr, "R2_", city, "_", dirf, "_all.csv", sep="")
    else
        fname <- paste (wdr, "R2_", city, "_", dirf, "_all_00.csv", sep="")
    r2 <- read.csv (fname, header=FALSE)

    r2 <- as.matrix (r2)
    r2 [r2 < -1] <- NA
    r2 <- sign (r2) * sqrt (abs (r2))
    dfull <- 1 - r2

    # Then read in the coordinates
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
    n <- dim (dat) [2]

    # Construct a neighbour list to spatially constrain clusters, as defined by
    # Delaunay triangulations
    xy <- data.frame (cbind (dat [lls [1]], dat [lls [2]]))
    # London has a duplicate xy for some reason
    nn <- as.numeric (row.names (unique (xy)))
    ni <- which ((1:length (nn) - nn) != 0)[1] # index of duplicates
    names (xy) <- c ("x", "y")

    # r2 matrices sometimes have entire missing columns, which are ultimately
    # returned with NA cluster memberships. These missing columns also have to
    # be removed from nbs. Also, the duplicate xy for London has to be removed
    # from the indx
    indx <- which (rowSums (dfull, na.rm=TRUE) > 0)
    indx <- indx [which (!indx %in% ni)]
    d <- dfull [indx, indx]
    npts <- length (indx)
    xy <- xy [indx,]
    
    tm <- tri.mesh (xy)
    nbs <- neighbours (tm)
    nbs.max <- unlist (lapply (nbs, function (x) length (x)))
    nbs.max <- max (nbs.max) # = 10

    # The "mstree" function in spdep requires a neighbour object of class "nb".
    # This can be make by calculating k-nearest neighbours, and then converting
    # with the knn2nb.  In this case, however, rather than having k-nearest
    # neighbours, the neighbours must be defined by the triangular mesh. The
    # next lines thus replace the matrix of k-nearest neighbours in knn$nn with
    # the results of the neighbours call above.
    knn <- knearneigh (coordinates (xy), k=nbs.max)
    # Then replace the knn list with the triangulation neighbours
    nn <- lapply (nbs, function (x) {
                  c (x, rep (NA, nbs.max - length (x))) })
    knn$nn <- do.call (rbind, nn)

    dlist <- lapply (1:npts, function (i) {
                    indxd <- knn$nn [i,]
                    indxd <- indxd [!is.na (indxd)]
                    c (d [i, indxd], rep (NA, nbs.max - length (indxd)))  })
    dlist <- data.frame (do.call (rbind, dlist))

    nb <- knn2nb (knn)
    lcosts <- nbcosts (nb, dlist)
    nb.w <- nb2listw (nb, lcosts, style="B")
    mst.sp <- mstree (nb.w)

    # Then apply the actual skater grouping. The data are padded out here with
    # NAs for those stations not in indx.
    #
    # Note that the 3 in the skater call is the min group size constraint, which
    # can be replaced with c(3,max) for a min and max constraint, but the skater
    # routine appears to have problems, as can be seen for example with c(3,50).
    #
    # Note also that the "i" parameter is the number of times the tree is
    # divided, so dividing it N times creates (N+1) groups.
    if (dirf == "from")
        fname <- paste (city, "_skater_from_members.txt", sep="")
    else
        fname <- paste (city, "_skater_to_members.txt", sep="")

    sk.groups <- array (NA, dim=c(nrow (dfull), max.groups - 1))
    require (data.table)
    st0 <- Sys.time ()
    pb <- txtProgressBar (0, 100, char="=", style=3)
    for (i in 1:(max.groups - 1)) {
        res <- skater (mst.sp [,1:2], dlist, i, crit=3, vec.crit=rep (1, npts))
        sk.groups [indx,i] <- res$groups
        write.table (sk.groups, file=fname, sep=",",
                     row.names=FALSE, col.names=FALSE)

        setTxtProgressBar (pb, 100 * i / (max.groups - 1))
    }
    close (pb)
    cat ("finished in ", timetaken (st0), "\n", sep="")
} # end function
