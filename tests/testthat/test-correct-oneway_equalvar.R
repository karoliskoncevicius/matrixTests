context("correctness of oneway_equalvar")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_oneway <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  st <- sr <- mt <- mr <- ot <- og <- dft <- dfr <- fst <- p <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- na.omit(mat[i,!bad])
    grp <- na.omit(groups[!bad])
    res <- summary(aov(vec ~ grp))[[1]]

    st[i]  <- res[1,2]
    sr[i]  <- res[2,2]
    mt[i]  <- res[1,3]
    mr[i]  <- res[2,3]
    dft[i] <- res[1,1]
    dfr[i] <- res[2,1]
    fst[i] <- res[1,4]
    p[i]   <- res[1,5]
    ot[i]  <- length(vec)
    og[i]  <- length(unique(grp))
  }

  data.frame(sum.sq.treatment=st, sum.sq.residuals=sr, mean.sq.treatment=mt,
             mean.sq.residuals=mr, obs.tot=ot, obs.groups=og, df.treatment=dft,
             df.residuals=dfr, F.statistic=fst, p.value=p
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

  t1 <- base_oneway(X, groups)
  t2 <- suppressWarnings(oneway_equalvar(X, groups))

  expect_equal(t1, t2)
})

################################## EDGE CASES ##################################

test_that("weird numbers give equal results", {
  x <- rnorm(12, sd=0.000001); g <- rep(c("a","b"),6)
  expect_equal(base_oneway(x, g), oneway_equalvar(x, g))
})

test_that("groups with one remaining member give equal results", {
  x <- rnorm(12); g <- rep(letters[1:4], each=3)
  gg <- g; gg[1:2] <- NA
  expect_equal(base_oneway(x, gg), suppressWarnings(oneway_equalvar(x, gg)))
  gg <- g; gg[c(1:2,4:5,7,10)] <- NA
  expect_equal(base_oneway(x, gg), suppressWarnings(oneway_equalvar(x, gg)))
  gg <- g; gg[c(1:2,4:5,7:8,10)] <- NA
  expect_equal(base_oneway(x, gg), suppressWarnings(oneway_equalvar(x, gg)))
})

test_that("constant values give equal results", {
  x <- c(1,1,1,1); g <- c("a","b","a","b")
  expect_equal(base_oneway(x, g), suppressWarnings(oneway_equalvar(x, g)))
  # x <- c(1,1,2,2); g <- c("a","a","b","b")
  # expect_equal(base_oneway(x, g), suppressWarnings(oneway_equalvar(x, g)))
})

