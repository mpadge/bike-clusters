#' get.partition.neighbours
#'
#' Identifies the clusters in min (nc) that split to form separate clusters in
#' max (nc). For each cluster, a list element is returned with the number of
#' groups into which each connected group is subsequently partitioned.
#'
#' @param city nyc, washington, chicago, boston, london (case insensitive)
#' @param method complete, ward, or skater, to be compared to k-means
#' @param nc list of 2 elements with nc[1]>nc[2], so that routine can identify
#' clusters in nc[1] that split to form clusters in nc[2]
#' @param dir to or from
#' @param plot enables visual inspection of the spatial arrangements.
#' @return list of neighbouring clusters in nc[1] that subsequently split in
#' nc[2]

get.partition.neighbours <- function (city="nyc", method="complete",
                                      nc=c(11,15), dir="from", plot=FALSE)
{
    require (spatstat)
    require (tripack) 
    require (igraph)

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
    
    nlo <- min (nc)
    nhi <- max (nc)
    membs.lo <- get.members (city=city, nc=nlo, method=method)
    membs.hi <- get.members (city=city, nc=nhi, method=method)

    if (dir == "from") 
        diri <- 1
    else 
        diri <- 2

    wd0 <- getwd ()
    wd <- ""
    while (!"bike-correlations" %in% list.files ("."))
    {
        wd <- paste (wd, "../", sep="")
        setwd ("../")
    }
    setwd (wd0)
    wd <- paste (wd, "bike-correlations/data/", sep="")

    if (city == "nyc" | city == "london" | city == "washingtondc")
    {
        fname <- paste (wd, "station_latlons_", city, ".txt", sep="")
        lls <- c ("long", "lat")
    } else if (city == "chicago") {
        fname <- paste (wd, "Divvy_Stations_2014-Q3Q4.csv", sep="")
        lls <- c ("longitude", "latitude")
    } else if (city == "boston") {
        fname <- paste (wd, "hubway_stations.csv", sep="")
        lls <- c ("lng", "lat")
    }
    dat <- read.csv (fname, header=TRUE) 
    xy <- data.frame (cbind (dat [lls [1]], dat [lls [2]]))
    names (xy) <- c ("x", "y")

    if (plot) {
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
    # clusters into which these become partitioned. All that remains is to
    # discern whether these clusters from within nc [1] are neighbours.

    npts <- nrow (xy)
    # London has a duplicate xy for some reason
    nn <- as.numeric (row.names (unique (xy)))
    ni <- which ((1:length (nn) - nn) != 0)[1] # index of duplicates
    xy <- unique (xy) 
    tm <- tri.mesh (xy)
    nbs <- neighbours (tm)
    nb.indx <- which (!1:npts %in%ni)

    nbmat <- array (FALSE, dim=c(npts, npts))
    for (i in 1:length (nbs)) 
    {
        nbmat [nb.indx [i], nb.indx [nbs [[i]] ] ] <- TRUE
        nbmat [nb.indx [nbs [[i]] ], nb.indx [i] ] <- TRUE
    }

    nbmat.groups <- array (FALSE, dim=c(nc [1], nc [1]))
    # ngroups will always be small, so a loop is used here because it's easier
    # to interpret than an lapply.
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
    nbmat.graph <- graph.adjacency (nbmat.groups, mode="undirected")
    nbmat.cl <- clusters (nbmat.graph)
    # Clusters has $membership IDs, $csize, and $no, where the latter is simply
    # the number of spatially separate clusters. The following results has one
    # list element for each distinct cluster. Each list element holds the number
    # of groups into which each neighbouring larger group (that is, each member
    # of nlo) is partitioned. The list is constructed with the largest group
    # first, although the size of a group just depends here on the number of
    # contiguous clusters within nc [1], and not the total number into which
    # that contiguous group is subsequently divided in nc [2].
    cindx <- order (nbmat.cl$csize, decreasing=TRUE)
    results <- lapply (1:length (cindx), function (i) {
                       membs <- which (nbmat.cl$membership == cindx [i])
                       tindx2 [membs]    })
    return (results)
} # end function get.partition.neighbours
