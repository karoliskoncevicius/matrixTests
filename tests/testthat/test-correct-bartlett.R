context("correctness of bartlett")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_bartlett <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  badgr <- names(which(table(groups)==1))
  bad   <- groups %in% badgr
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  ng <- nt <- vt <- ks <- p <- df <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- na.omit(mat[i,!bad])
    grp <- na.omit(groups[!bad])
    res <- bartlett.test(vec, grp)

    ng[i] <- length(unique(grp))
    nt[i] <- length(vec)
    vt[i] <- sum((tapply(vec, grp, length)-1) * tapply(vec, grp, var, na.rm=TRUE)) / (nt[i]-ng[i])
    ks[i] <- res$statistic
    p[i]  <- res$p.value
    df[i]  <- res$parameter
  }

  data.frame(var.tot=vt, obs.tot=nt, obs.groups=ng, ksq.statistic=ks,
             p.value=p, df=df, stringsAsFactors=FALSE
             )
}

################################################################################
##################### TEST CONSISTENCY WITH BARTLETT.TEST ######################
################################################################################

test_that("monte-carlo random testing gives equal results", {
  set.seed(14)
  X <- matrix(rnorm(10000), ncol=100)
  X[sample(length(X), 100)] <- NA
  groups <- sample(c("a","b","c","d"), 100, replace=TRUE)

  t1 <- base_bartlett(X, groups)
  t2 <- bartlett(X, groups)

  expect_equal(t1, t2)
})

################################## EDGE CASES ##################################

test_that("weird numbers give equal results", {
  x <- rnorm(12, sd=0.000001); g <- rep(c("a","b"), each=6)
  expect_equal(base_bartlett(x, g), bartlett(x, g))
})

test_that("minumum allowed sample sizes give equal results", {
  x <- rnorm(4); g <- c("a", "b", "a", "b")
  expect_equal(base_bartlett(x, g), bartlett(x, g))
  expect_equal(base_bartlett(c(x, NA), c(g, "b")), bartlett(c(x, NA), c(g, "b")))
})

test_that("groups with one element remaining are dropped correctly", {
  x <- rnorm(12)
  g <- rep(letters[1:4], each=3)
  gg <- g; gg[1:2] <- NA
  expect_equal(base_bartlett(x, gg), suppressWarnings(bartlett(x, gg)))
  expect_equal(bartlett(x[-c(1:3)], gg[-c(1:3)]), suppressWarnings(bartlett(x, gg)))
  gg <- g; gg[c(1:2,4:5)] <- NA
  expect_equal(base_bartlett(x, gg), suppressWarnings(bartlett(x, gg)))
  expect_equal(bartlett(x[-c(1:6)], gg[-c(1:6)]), suppressWarnings(bartlett(x, gg)))
})

test_that("constant values give equal results", {
  x <- c(1,1,1,1); g <- c("a","b","a","b")
  expect_equal(base_bartlett(x, g), suppressWarnings(bartlett(x, g)))
  x <- c(1,1,2,2); g <- c("a","a","b","b")
  expect_equal(base_bartlett(x, g), suppressWarnings(bartlett(x, g)))
  x <- c(1,1,2,3); g <- c("a","a","b","b")
  expect_equal(base_bartlett(x, g), suppressWarnings(bartlett(x, g)))
})

