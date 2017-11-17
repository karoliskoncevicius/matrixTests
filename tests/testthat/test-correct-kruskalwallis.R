context("correctness of kruskalwallis")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_kruskal <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)

  bad <- is.na(groups)
  mat    <- mat[,!bad,drop=FALSE]
  groups <- groups[!bad]

  ot <- og <- df <- chs <- p <- numeric(nrow(mat))
  for(i in 1:nrow(mat)) {
    bad <- is.na(mat[i,])
    vec <- na.omit(mat[i,!bad])
    grp <- na.omit(groups[!bad])
    res <- kruskal.test(vec, grp)

    ot[i]  <- length(vec)
    og[i]  <- length(unique(grp))
    df[i]  <- res$parameter
    chs[i] <- res$statistic
    p[i]   <- res$p.value
  }

  data.frame(obs.tot=ot, obs.groups=og, df, chsq.statistic=chs, p.value=p)
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

  t1 <- base_kruskal(X, factor(groups))
  t2 <- suppressWarnings(kruskalwallis(X, groups))

  expect_equal(t1, t2)
})

################################## EDGE CASES ##################################

test_that("weird numbers give equal results", {
  x <- rnorm(12, sd=0.000001); g <- rep(c("a","b"),6)
  expect_equal(base_kruskal(x, factor(g)), kruskalwallis(x, g))
})

test_that("groups with one remaining member give equal results", {
  x <- rnorm(12); g <- rep(letters[1:4], each=3)
  gg <- factor(g); gg[1:2] <- NA
  expect_equal(base_kruskal(x, gg), suppressWarnings(kruskalwallis(x, gg)))
  gg <- factor(g); gg[c(1:2,4:5,7,10)] <- NA
  expect_equal(base_kruskal(x, gg), suppressWarnings(kruskalwallis(x, gg)))
  gg <- factor(g); gg[c(1:2,4:5,7:8,10)] <- NA
  expect_equal(base_kruskal(x, gg), suppressWarnings(kruskalwallis(x, gg)))
  gg <- factor(LETTERS[1:12])
  expect_equal(base_kruskal(x, gg), suppressWarnings(kruskalwallis(x, gg)))
})

test_that("minimal allowed sample size gives equal results", {
  x <- c(1,2); g <- factor(c("a","b"))
  expect_equal(base_kruskal(x, g), suppressWarnings(kruskalwallis(x, g)))
  x <- c(1,2,NA,3); g <- factor(c("a","b","a",NA))
  expect_equal(base_kruskal(x, g), suppressWarnings(kruskalwallis(x, g)))
})

test_that("constant values give equal results", {
  x <- c(1,1,1,1); g <- factor(c("a","b","a","b"))
  expect_equal(base_kruskal(x, g), suppressWarnings(kruskalwallis(x, g)))
  x <- c(1,1,2,2); g <- factor(c("a","a","b","b"))
  expect_equal(base_kruskal(x, g), suppressWarnings(kruskalwallis(x, g)))
})

