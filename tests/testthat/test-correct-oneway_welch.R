context("correctness of oneway_equalvar")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_oneway_welch <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  ot <- og <- dft <- dfr <- fst <- p <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- na.omit(mat[i,!bad])
    grp <- na.omit(groups[!bad])
    res <- oneway.test(vec ~ grp)

    dft[i] <- res$parameter[1]
    dfr[i] <- res$parameter[2]
    fst[i] <- res$statistic
    p[i]   <- res$p.value
    ot[i]  <- length(vec)
    og[i]  <- length(unique(grp))
  }

  data.frame(obs.tot=ot, obs.groups=og, df.treatment=dft, df.residuals=dfr,
             F.statistic=fst, p.value=p
             )
}

################################################################################
##################### TEST CONSISTENCY WITH aov() ##############################
################################################################################

test_that("monte-carlo random testing gives equal results", {
  set.seed(14)
  X <- matrix(rnorm(10000), ncol=100)
  X[sample(length(X), 100)] <- NA
  groups <- sample(c("a","b","c","d"), 100, replace=TRUE)
  groups[sample(length(groups), 10)] <- NA

  t1 <- base_oneway_welch(X, groups)
  t2 <- suppressWarnings(oneway_welch(X, groups))

  expect_equal(t1, t2)
})

################################## EDGE CASES ##################################

test_that("weird numbers give equal results", {
  x <- rnorm(12, sd=0.000001); g <- rep(c("a","b"),6)
  expect_equal(base_oneway_welch(x, g), oneway_welch(x, g))
})

test_that("minimum allowed sample size gives equal results", {
  x <- rnorm(12); g <- rep(letters[1:4], each=3)
  gg <- g; gg[seq(1,12,3)] <- NA
  expect_equal(base_oneway_welch(x, gg), suppressWarnings(oneway_welch(x, gg)))
  gg <- g; gg[seq(1,12,3)] <- NA; x[1:3] <- NA
  expect_equal(base_oneway_welch(x, gg), suppressWarnings(oneway_welch(x, gg)))
})

test_that("constant values give equal results", {
  # x <- c(1,1,1,1); g <- c("a","b","a","b")
  # expect_equal(base_oneway_welch(x, g), suppressWarnings(oneway_welch(x, g)))
  # x <- c(1,1,2,2); g <- c("a","a","b","b")
  # expect_equal(base_oneway_welch(x, g), suppressWarnings(oneway_welch(x, g)))
})

