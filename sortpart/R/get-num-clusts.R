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
    nc <- dat$nc
    dat <- dat [,14:17] # The simulated probabilities

    if (city == "chicago")
    {
        # These lines equate numbers of clusters with minimal probabilities
        # (excluding the first peak, which is sometimes the minimum).
        mini <- 1 + apply (dat [2:dim(dat)[1],], 2, which.min)
    } else {
        # Whereas more peaks are included if the number is selected from the
        # largest peak with p < p0 (=0.05).
        p0 <- 0.05
        # with a couple of manual tweaks
        # TODO: Check whether these are necessary with the new results!
        if (city == "nyc" | city == "boston")
            p0 <- 0.057
        else if (city == "london")
            p0 <- 0.065
        mini <- apply (dat, 2, function (x) {
                       if (length (which (x < p0)) == 0)
                           which.min (x)
                       else
                           which (x < p0)
                })
        mini <- sapply (mini, max)
    }
    mini.indx <- (0:3) * dim (dat) [1] + mini
    dat.mini <- as.matrix (dat) [mini.indx]

    results <- data.frame (cbind (nc [mini], dat.mini))
    names (results) <- c ("num.clusts", "pmin")
    rownames (results) <- c ("to-h", "to-k", "from-h", "from-k")
    return (results)
} # end function get. num.clusts
