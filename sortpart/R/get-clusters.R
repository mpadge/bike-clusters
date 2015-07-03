#' get-clusters
#'
#' Produces .csv files with cluster memberships across a range of numbers of
#' clusters (1-max_clust_size=100), along with corresponding average cluster
#' diameters (in km).
#'
#' @param city nyc, washington, chicago, boston, london (case insensitive)
#' @param method = (ward, k-means, complete)
#' @return nothing (dumps files to dir)

get.clusters <- function (city="nyc", method="ward")
{
    require (spatstat) # for ppp
    require (geosphere) # for areaPolygon

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

    max_clust_size <- 100

    wd0 <- getwd ()
    count <- 0
    wd <- ""
    while (!"bike-correlations" %in% list.files ("."))
    {
        wd <- paste (wd, "../", sep="")
        setwd ("../")
        count <- count + 1
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
    n <- dim (dat) [1]
    clust.from <- clust.to <- array (NA, dim=c(n, max_clust_size))
    arr <- array (NA, dim=c(max_clust_size, max_clust_size))
    clust.diam.from <- clust.diam.to <- arr

    pb <- txtProgressBar (0, 100, char="=", style=3)
    for (nc in 1:max_clust_size) {
        membs <- get.members (city=city, nc, method=method)
        clust.from [, nc] <- as.vector (membs [,1])
        clust.to [, nc] <- as.vector (membs [,2])

        # Then calculate diameters of clusters:
        for (i in 1:nc) {
            for (j in 1:2) {
                indx <- which (membs [,j] == i)
                if (length (indx) > 2) {
                    x <- dat [[lls [1]]] [indx]
                    y <- dat [[lls [2]]] [indx]
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
        setTxtProgressBar (pb, 100 * i / max_clust_size)
    } # end for nc
    close (pb)
    setwd (wd0)
    indx <- 2:max_clust_size
    clust.from <- clust.from [,indx]
    clust.to <- clust.to [,indx]
    clust.diam.from <- clust.diam.from [indx, ]
    clust.diam.to <- clust.diam.to [indx, ]
    fname <- paste (city, "-clust-from-members-", method, ".txt", sep="")
    write.table (clust.from, file=fname, sep=",", row.names=FALSE, col.names=FALSE)
    fname <- paste (city, "-clust-to-members-", method, ".txt", sep="")
    write.table (clust.to, file=fname, sep=",", row.names=FALSE, col.names=FALSE)
    fname <- paste (city, "-clust-from-diameters-", method, ".txt", sep="")
    write.table (clust.diam.from, file=fname, sep=",", row.names=FALSE,
                 col.names=FALSE)
    fname <- paste (city, "-clust-to-diameters-", method, ".txt", sep="")
    write.table (clust.diam.to, file=fname, sep=",", row.names=FALSE,
                 col.names=FALSE)
} # end function get.clusters
