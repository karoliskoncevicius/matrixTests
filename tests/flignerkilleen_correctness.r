library(matrixTests)

#--- functions -----------------------------------------------------------------

base_fligner <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad    <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  ng <- nt <- ks <- p <- df <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- na.omit(mat[i,!bad])
    grp <- na.omit(groups[!bad])

    res <- fligner.test(vec, factor(grp))

    # if p-value is NA turn df to NA as well
    if(is.na(res$p.value)) res$parameter <- NA

    ng[i] <- length(unique(grp))
    nt[i] <- length(vec)
    ks[i] <- res$statistic
    p[i]  <- res$p.value
    df[i] <- res$parameter
  }

  data.frame(obs.tot=nt, obs.groups=ng, df=df, statistic=ks, pvalue=p,
             stringsAsFactors=FALSE
             )
}


#--- montecarlo ----------------------------------------------------------------

# two groups
x <- matrix(rnorm(10000), ncol=10)
g <- sample(letters[1:2], 10, replace=TRUE)
res1 <- base_fligner(x, factor(g))
res2 <- row_flignerkilleen(x, g)
stopifnot(all.equal(res1, res2))

# lots of groups
x <- matrix(rnorm(100000), ncol=100)
g <- sample(letters[1:15], 100, replace=TRUE)
res1 <- base_fligner(x, factor(g))
res2 <- row_flignerkilleen(x, g)
stopifnot(all.equal(res1, res2))

# TODO: test why the following returns NA: row_flignerkilleen(c(1,2,4,3), c(1,1,2,2))
#       it returns NAs for both base and matrixTests version, so probably not an error

#--- extreme numbers -----------------------------------------------------------

# big numbers
x <- c(100000000000004, 100000000000002, 100000000000003, 100000000000000,
       100000000000003, 100000000000002, 100000000000003, 100000000000000
       )
g <- c(rep("a", 4), rep("b", 4))
res1 <- base_fligner(x, g)
res2 <- row_flignerkilleen(x, g)
stopifnot(all.equal(res1, res2))

# small numbers
x <- c(1.00000000000004, 1.00000000000002, 1.00000000000003, 1.00000000000000,
       1.00000000000003, 1.00000000000002, 1.00000000000003, 1.00000000000000
       )
g <- c(rep("a", 4), rep("b", 4))
res1 <- base_fligner(x, g)
res2 <- suppressWarnings(row_flignerkilleen(x, g))
stopifnot(all.equal(res1, res2))

# TODO: add tests for Inf and -Inf values once decided how to handle them.


#--- minimal sample size -------------------------------------------------------

# two groups one with 2 samples and another with 1
x <- rnorm(3)
g <- c("a", "a", "b")
res1 <- base_fligner(x, g)
res2 <- row_flignerkilleen(x, g)
stopifnot(all.equal(res1, res2))

# two groups both with 2 values
x <- rnorm(4)
g <- c("a", "a", "b", "b")
res1 <- base_fligner(x, g)
res2 <- row_flignerkilleen(x, g)
stopifnot(all.equal(res1, res2))

# many groups one with 2 values and the rest with 1
x <- rnorm(15)
g <- c("a", letters[1:14])
res1 <- base_fligner(x, g)
res2 <- row_flignerkilleen(x, g)
stopifnot(all.equal(res1, res2))

# two groups one with 3 samples and another with 1
x <- rnorm(4)
g <- c("a", "a", "a", "b")
res1 <- base_fligner(x, g)
res2 <- row_flignerkilleen(x, g)
stopifnot(all.equal(res1, res2))

# two groups one with 2 samples and another with 1 plus NAs
x <- c(rnorm(2), NA, NA, NA, rnorm(1))
g <- c("a","a","a","b","b","b")
res1 <- base_fligner(x, g)
res2 <- row_flignerkilleen(x, g)
stopifnot(all.equal(res1, res2))


#--- constant values -----------------------------------------------------------

# one group's values are constant
x <- c(1,1,1,2,3,4)
g <- c("a","a","a","b","b","b")
res1 <- base_fligner(x, g)
res2 <- row_flignerkilleen(x, g)
stopifnot(all.equal(res1, res2))

# all groups have equal variances
x <- c(1:3,2:4,3:5,4:6,5:7)
g <- rep(1:5, each=3)
res1 <- base_fligner(x, g)
res2 <- row_flignerkilleen(x, g)
stopifnot(all.equal(res1, res2))

