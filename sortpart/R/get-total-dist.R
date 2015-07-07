#' get.total.dist
#'
#' Calculates total distance ridden for given city
#'
#' @param city nyc, washington, chicago, boston, london (case insensitive)
#' @return total distance in km

get.total.dist <- function (city = "nyc")
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

    wd0 <- getwd ()
    wd <- ""
    while (!"bike-correlations" %in% list.files ("."))
    {
        wd <- paste (wd, "../", sep="")
        setwd ("../")
    }
    setwd (wd0)
    wd <- paste (wd, "bike-correlations/results/", sep="")

    dname <- paste (wd, "stationDistsMat-", city, ".csv", sep="")
    if (!file.exists (dname))
        stop (paste (dname, " does not exist!", sep=""))
    dists <- read.csv (dname, header=FALSE)
    if (city %in% c ("boston", "chicago", "nyc"))
        nname <- paste (wd, "NumTrips_", city, "_00.csv", sep="")
    else
        nname <- paste (wd, "NumTrips_", city, ".csv", sep="")
    if (!file.exists (nname))
        stop (paste (nname, " does not exist!", sep=""))
    ntrips <- read.csv (nname, header=FALSE)

    return (sum (dists * ntrips))
}
