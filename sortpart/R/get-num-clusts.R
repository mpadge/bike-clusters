#' get.num.clusts
#'
#' Just loads up the output of calc.pnc(), and returns the final numbers of
#' clusters corresponding to minimal joint probabilities of peak spacings and
#' depths.
#'
#' @param city nyc, washington, chicago, boston, london (case insensitive)
#' @param method = (ward, k-means, complete)
#' @return data frame

get.num.clusts <- function (city="nyc", method="complete")
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
    
    fname <- paste (wd, city, "-results-prob-m-", method, ".txt", sep="")
    dat <- read.csv (fname, header=TRUE)
    nc <- dat [,10:13]
    index <- dat$nc
    dat <- dat [,14:17] # The simulated probabilities

    if (city == "dummy")
    {
        # These lines equate numbers of clusters with minimal probabilities
        # (excluding the first peak, which is sometimes the minimum).
        mini <- 1 + apply (dat [2:dim(dat)[1],], 2, which.min)
    } else {
        # Whereas more peaks are included if the number is selected from the
        # largest peak with p < p0:
        p0 <- 0.05
        mini <- apply (dat, 2, function (x) {
                       if (length (which (x < p0)) == 0)
                           which.min (x)
                       else
                           which (x < p0)
                })
        mini <- sapply (mini, max)

        if (city == "nyc")
            mini [1] <- max (which (dat [,1] < 0.1))
    }
    dat.mini <- mapply (function (x, i) x[i], x=dat, i=mini)
    nc.mini <- mapply (function (x, i) x[i], x=nc, i=mini)

    # Then get cluster diameters
    diams <- rep (NA, 4)
    ttxt <- c ("to", "from")
    for (i in 1:2)
    {
        fname <- paste (wd, city, "-clust-", ttxt [i], 
                        "-diameters-complete.txt", sep="")
        dtemp <- read.csv (fname, header=FALSE)
        dtemp <- as.numeric (rowMeans (dtemp, na.rm=TRUE))
        diams [i * 2 - 1] <- dtemp [nc.mini [i * 2 - 1]]
    }
    

    results <- data.frame (cbind (index [mini], nc.mini, dat.mini, diams))
    names (results) <- c ("num.pks", "num.clusts", "pmin", "diam")
    rownames (results) <- c ("to-h", "to-k", "from-h", "from-k")
    return (results)
} # end function get. num.clusts
