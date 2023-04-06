library(matrixTests)

#--- functions -----------------------------------------------------------------

car_levene <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  ot <- og <- dft <- dfr <- fst <- p <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- mat[i,!bad]
    grp <- factor(groups[!bad])

    res <- car::leveneTest(vec ~ grp, center="mean")

    # if p-vallue is NA turn df to NA as well
    if(is.na(res[["Pr(>F)"]][1])) res[["Df"]] <- NA

    dft[i] <- res[["Df"]][1]
    dfr[i] <- res[["Df"]][2]
    fst[i] <- res[["F value"]][1]
    p[i]   <- res[["Pr(>F)"]][1]
    ot[i]  <- length(vec)
    og[i]  <- length(unique(grp))
  }

  data.frame(obs.tot=ot, obs.groups=og, df.between=dft, df.within=dfr,
             statistic=fst, pvalue=p
             )
}


#--- montecarlo ----------------------------------------------------------------

# two groups
x <- matrix(rnorm(10000), ncol=10)
g <- sample(letters[1:2], 6, replace=TRUE)
g <- sample(c("a", "a", "b", "b", g))  # ensure both groups have at least 2 obs
res1 <- car_levene(x, factor(g))
res2 <- row_levene(x, g)
stopifnot(all.equal(res1, res2))

# lots of groups
x <- matrix(rnorm(100000), ncol=100)
g <- sample(letters[1:15], 100, replace=TRUE)
res1 <- car_levene(x, factor(g))
res2 <- row_levene(x, g)
stopifnot(all.equal(res1, res2))


#--- extreme numbers -----------------------------------------------------------

# big numbers
x <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000,
       100000000000003, 100000000000002, 100000000000003, 100000000000000
       )
g <- c(rep("a", 4), rep("b", 4))
res1 <- car_levene(x, g)
res2 <- row_levene(x, g)
stopifnot(all.equal(res1, res2))

# small numbers
x <- c(1.00000000000004, 1.00000000000002, 1.00000000000003, 1.00000000000000,
       1.00000000000003, 1.00000000000002, 1.00000000000003, 1.00000000000000
       )
g <- c(rep("a", 4), rep("b", 4))
res1 <- car_levene(x, g)
res2 <- suppressWarnings(row_levene(x, g))
stopifnot(all.equal(res1, res2))

# large sample
x <- rnorm(3 * 10^6)
g <- rep(letters[1:3], each=10^6)
res1 <- car_levene(x, g)
res2 <- row_levene(x, g)
stopifnot(all.equal(res1, res2))

# TODO: add tests for Inf and -Inf values once decided how to handle them.


#--- minimal sample size -------------------------------------------------------

# two groups one with three values and another with one value
x <- c(1,2,3,4)
g <- c("a","a","a","b")
res1 <- car_levene(x, g)
res2 <- row_levene(x, g)
stopifnot(all.equal(res1, res2))

# two groups one with three values and another with two values
x <- c(1,2,3,4,5)
g <- c("a","a","a","b","b")
res1 <- car_levene(x, g)
res2 <- row_levene(x, g)
stopifnot(all.equal(res1, res2))

# five groups one with three values and others with one
x <- c(rnorm(3), 1, 2, 3, 4)
g <- c("a","a","a","b","c","d","e")
res1 <- car_levene(x, g)
res2 <- row_levene(x, g)
stopifnot(all.equal(res1, res2))


#--- constant values -----------------------------------------------------------

# all values are zero
x <- c(0,0,0,0,0,0)
g <- c("a","a","a","b","b","b")
res1 <- car_levene(x, g)
res2 <- suppressWarnings(row_levene(x, g))
stopifnot(all.equal(res1, res2))

# all values are constant
x <- c(1,1,1,1,1,1)
g <- c("a","a","a","b","b","b")
res1 <- car_levene(x, g)
res2 <- suppressWarnings(row_levene(x, g))
stopifnot(all.equal(res1, res2))

# within group values are constant
x <- c(1,1,1,2,2,2)
g <- c("a","a","a","b","b","b")
res1 <- car_levene(x, g)
res2 <- suppressWarnings(row_levene(x, g))
stopifnot(all.equal(res1, res2))

# one group's values are constant
x <- c(1,1,1,2,3,4)
g <- c("a","a","a","b","b","b")
res1 <- car_levene(x, g)
res2 <- row_levene(x, g)
stopifnot(all.equal(res1, res2))

# all groups have equal variances
x <- c(1:3,2:4,3:5,4:6,5:7)
g <- rep(1:5, each=3)
res1 <- car_levene(x, g)
res2 <- row_levene(x, g)
stopifnot(all.equal(res1, res2))


#--- constant residuals --------------------------------------------------------

# two groups - one has all constant values after residuals
x <- c(rnorm(4), 3, 3, 2, 2)
g <- rep(letters[1:2], each=4)
res1 <- car_levene(x, g)
res2 <- row_levene(x, g)
stopifnot(all.equal(res1, res2))

# many groups - all but one have constant values after residuals
x <- c(rnorm(3), 2, 0,3, 4,4,4, 6,6,5,5)
g <- rep(letters[1:5], c(3,1,2,3,4))
res1 <- car_levene(x, g)
res2 <- row_levene(x, g)
stopifnot(all.equal(res1, res2))

# many groups - all but one have constant values after residuals plus NAs
x <- c(rnorm(3), 2,NA, 0,3,NA, 4,4,4,NA, 6,6,5,5,NA)
g <- rep(letters[1:5], c(3,2,3,4,5))
res1 <- car_levene(x, g)
res2 <- row_levene(x, g)
stopifnot(all.equal(res1, res2))

