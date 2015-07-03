get.skater.groups <- function (dir="FROM", max.groups=50)
{
    require (sp)
    require (spdep)
    require (tripack)

    setwd ("/data/analyses/bikes_and_cities/cluster_dynamics/results/")
    fname <- paste ("results_kvals.txt")
    dat <- read.csv (fname, header=TRUE) 
    xy <- data.frame (cbind (dat$lon, dat$lat))
    names (xy) <- c ("x", "y")
    npts <- nrow (xy)
    tm <- tri.mesh (xy)
    nbs <- neighbours (tm)
    nbs.max <- unlist (lapply (nbs, function (x) length (x)))
    nbs.max <- max (nbs.max) # = 10

    # The "mstree" function in spdep requires a neighbour object of class "nb". This can
    # be make by calculating k-nearest neighbours, and then converting with the knn2nb.
    # In this case, however, rather than having k-nearest neighbours, the neighbours
    # must be defined by the triangular mesh. The next lines thus replace the matrix of
    # k-nearest neighbours in knn$nn with the results of the neighbours call above.
    knn <- knearneigh (coordinates (xy), k=nbs.max)
    # Then replace the knn list with the triangulation neighbours
    nn <- lapply (nbs, function (x) {
                  c (x, rep (NA, nbs.max - length (x))) })
    knn$nn <- do.call (rbind, nn)

    # Then read in correlations and use these to weight connections in the mstree. The
    # skater routine constructs *MINIMAL* variance sub-trees, so the data submitted to
    # skater thus have to be distances, rather than correlations.
    if (dir == "FROM") {
        r2 <- read.csv ("results_r2from.txt", header=FALSE)
        fname <- "skater_from_members.txt"
    } else {
        r2 <- read.csv ("results_r2to.txt", header=FALSE)
        fname <- "skater_to_members.txt"
    }
    r2 <- as.matrix (r2)
    r2 [r2 < -1] <- NA
    r2 <- sign (r2) * sqrt (abs (r2))
    d <- 1 - r2
    dlist <- lapply (1:npts, function (i) {
                    indx <- knn$nn [i,]
                    indx <- indx [!is.na (indx)]
                    c (d [i, indx], rep (NA, nbs.max - length (indx)))  })
    dlist <- data.frame (do.call (rbind, dlist))

    nb <- knn2nb (knn)
    lcosts <- nbcosts (nb, dlist)
    nb.w <- nb2listw (nb, lcosts, style="B")
    mst.sp <- mstree (nb.w)

    # Note that the 3 in the skater call is the min group size constrain, which can
    # be replaced with c(3,max) for a min and max constraint, but the skater routine
    # appears to have problems, as can be seen for example with c(3,50).
    #
    # Note also that the "i" parameter is the number of times the tree is divided,
    # so dividing it N times creates (N+1) groups.
    setwd ("/data/analyses/bikes_and_cities/cluster_dynamics/")
    sk.groups <- array (NA, dim=c(npts, max.groups - 1))
    require (data.table)
    st0 <- Sys.time ()
    for (i in 1:(max.groups - 1)) {
        res <- skater (mst.sp [,1:2], dlist, i, crit=3, vec.crit=rep (1, npts))
        sk.groups [,i] <- res$groups
        st <- timetaken (st0)
        cat ("[", i + 1, " / ", max.groups, "] ", st, "\n", sep="")
        write.table (sk.groups, file=fname, sep=",",
                     row.names=FALSE, col.names=FALSE)
    }
    cat ("finished. Total time taken = ", timetaken (st0), "\n", sep="")
} # end function
