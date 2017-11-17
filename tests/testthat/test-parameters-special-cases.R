context("Parameter special cases")

################################################################################
######################### ALLOWED DATA MATRIX FORMATS ##########################
################################################################################

test_that("x and y can be numeric vectors", {
  x <- 1:10; X <- matrix(x, nrow=1)
  y <- 10:1; Y <- matrix(y, nrow=1)
  grp <- sample(c(1,0), 10, replace=TRUE)
  expect_equal(ttest_onegroup(x=x), ttest_onegroup(x=X))
  expect_equal(ttest_equalvar(x=x, y=y), ttest_equalvar(x=X, y=Y))
  expect_equal(ttest_welch(x=x, y=y), ttest_welch(x=X, y=Y))
  expect_equal(ttest_paired(x=x, y=y), ttest_paired(x=X, y=Y))
  expect_equal(oneway_equalvar(x=x, groups=grp), oneway_equalvar(x=X, groups=grp))
  expect_equal(oneway_welch(x=x, groups=grp), oneway_welch(x=X, groups=grp))
  expect_equal(kruskalwallis(x=x, groups=grp), kruskalwallis(x=X, groups=grp))
  expect_equal(bartlett(x=x, groups=grp), bartlett(x=X, groups=grp))
  expect_equal(ievora(x=x, groups=grp), ievora(x=X, groups=grp))
})

test_that("x and y can be numeric data.frames", {
  x <- t(iris[1:75,-5]); X <- data.frame(x)
  y <- t(iris[76:150,-5]); Y <- data.frame(y)
  grp <- iris$Species[1:75]=="setosa"
  expect_equal(ttest_onegroup(x=x), ttest_onegroup(x=X))
  expect_equal(ttest_equalvar(x=x, y=y), ttest_equalvar(x=X, y=Y))
  expect_equal(ttest_welch(x=x, y=y), ttest_welch(x=X, y=Y))
  expect_equal(ttest_paired(x=x, y=y), ttest_paired(x=X, y=Y))
  expect_equal(oneway_equalvar(x=x, groups=grp), oneway_equalvar(x=X, groups=grp))
  expect_equal(oneway_welch(x=x, groups=grp), oneway_welch(x=X, groups=grp))
  expect_equal(kruskalwallis(x=x, groups=grp), kruskalwallis(x=X, groups=grp))
  expect_equal(bartlett(x=x, groups=grp), bartlett(x=X, groups=grp))
  expect_equal(ievora(x=x, groups=grp), ievora(x=X, groups=grp))
})

test_that("x and y can have 0 rows and 0 columns", {
  x <- matrix(0, nrow=0, ncol=0)
  y <- matrix(0, nrow=0, ncol=0)
  grp <- numeric()
  expect_equal(nrow(ttest_onegroup(x=x)), 0)
  expect_equal(nrow(ttest_equalvar(x=x, y=y)), 0)
  expect_equal(nrow(ttest_welch(x=x, y=y)), 0)
  expect_equal(nrow(ttest_paired(x=x, y=y)), 0)
  expect_equal(nrow(oneway_equalvar(x=x, groups=grp)), 0)
  expect_equal(nrow(oneway_welch(x=x, groups=grp)), 0)
  expect_equal(nrow(kruskalwallis(x=x, groups=grp)), 0)
  expect_equal(nrow(bartlett(x=x, groups=grp)), 0)
  expect_equal(nrow(ievora(x=x, groups=grp)), 0)
})

test_that("NA and NaN are treated the same", {
  Xna <- Xnan <- matrix(rnorm(100), nrow=10)
  Yna <- Ynan <- matrix(rnorm(100), nrow=10)
  Xna[sample(length(Xna), 20)] <- NA
  Xnan[is.na(Xna)] <- NaN
  Yna[sample(length(Yna), 20)] <- NA
  Ynan[is.na(Yna)] <- NaN
  grp <- c(rep(0,5), rep(1,5))
  expect_equal(ttest_onegroup(x=Xna), ttest_onegroup(x=Xnan))
  expect_equal(ttest_equalvar(x=Xna, y=Yna), ttest_equalvar(x=Xnan, y=Ynan))
  expect_equal(ttest_welch(x=Xna, y=Yna), ttest_welch(x=Xnan, y=Ynan))
  expect_equal(ttest_paired(x=Xna, y=Yna), ttest_paired(x=Xnan, y=Ynan))
  expect_equal(oneway_equalvar(x=Xna, groups=grp), oneway_equalvar(x=Xnan, groups=grp))
  expect_equal(oneway_welch(x=Xna, groups=grp), oneway_welch(x=Xnan, groups=grp))
  expect_equal(kruskalwallis(x=Xna, groups=grp), kruskalwallis(x=Xnan, groups=grp))
  expect_equal(bartlett(x=Xna, groups=grp), bartlett(x=Xnan, groups=grp))
  expect_equal(ievora(x=Xna, groups=grp), ievora(x=Xnan, groups=grp))
})

################################################################################
############################ ALLOWED GROUP FORMATS #############################
################################################################################

test_that("NA and NaN produce the same result", {
  set.seed(14)
  x <- t(iris[,-5])
  grp1 <- grp2 <- grp3 <- as.numeric(iris$Species=="setosa")
  inds <- sample(length(grp1), 10)
  grp1[inds] <- NA
  grp2[inds] <- NaN
  grp3[inds[1:5]] <- NA
  grp3[inds[6:10]] <- NaN
  expect_equal(suppressWarnings(oneway_equalvar(x=x, groups=grp1)), suppressWarnings(oneway_equalvar(x=x, groups=grp2)))
  expect_equal(suppressWarnings(oneway_equalvar(x=x, groups=grp1)), suppressWarnings(oneway_equalvar(x=x, groups=grp3)))
  expect_equal(suppressWarnings(oneway_welch(x=x, groups=grp1)), suppressWarnings(oneway_welch(x=x, groups=grp2)))
  expect_equal(suppressWarnings(oneway_welch(x=x, groups=grp1)), suppressWarnings(oneway_welch(x=x, groups=grp3)))
  expect_equal(suppressWarnings(kruskalwallis(x=x, groups=grp1)), suppressWarnings(kruskalwallis(x=x, groups=grp2)))
  expect_equal(suppressWarnings(kruskalwallis(x=x, groups=grp1)), suppressWarnings(kruskalwallis(x=x, groups=grp3)))
  expect_equal(suppressWarnings(bartlett(x=x, groups=grp1)), suppressWarnings(bartlett(x=x, groups=grp2)))
  expect_equal(suppressWarnings(bartlett(x=x, groups=grp1)), suppressWarnings(bartlett(x=x, groups=grp3)))
  expect_equal(suppressWarnings(ievora(x=x, groups=grp1)), suppressWarnings(ievora(x=x, groups=grp2)))
  expect_equal(suppressWarnings(ievora(x=x, groups=grp1)), suppressWarnings(ievora(x=x, groups=grp3)))
})

test_that("groups can all be NA", {
  x <- t(iris[,-5])
  grp <- rep(NA, ncol(x))
  expect_equal(suppressWarnings(oneway_equalvar(x=x, groups=grp)$obs.tot), rep(0, nrow(x)))
  expect_equal(suppressWarnings(oneway_welch(x=x, groups=grp)$obs.tot), rep(0, nrow(x)))
  expect_equal(suppressWarnings(kruskalwallis(x=x, groups=grp)$obs.tot), rep(0, nrow(x)))
  expect_equal(suppressWarnings(bartlett(x=x, groups=grp)$obs.tot), rep(0, nrow(x)))
  expect_equal(suppressWarnings(ievora(x=x, groups=grp)$obs.0), rep(0, nrow(x)))
})

test_that("groups can all be NaN", {
  x <- t(iris[,-5])
  grp <- rep(NaN, ncol(x))
  expect_equal(suppressWarnings(oneway_equalvar(x=x, groups=grp)$obs.tot), rep(0, nrow(x)))
  expect_equal(suppressWarnings(oneway_welch(x=x, groups=grp)$obs.tot), rep(0, nrow(x)))
  expect_equal(suppressWarnings(kruskalwallis(x=x, groups=grp)$obs.tot), rep(0, nrow(x)))
  expect_equal(suppressWarnings(bartlett(x=x, groups=grp)$obs.tot), rep(0, nrow(x)))
  expect_equal(suppressWarnings(ievora(x=x, groups=grp)$obs.0), rep(0, nrow(x)))
})

test_that("all observations can be in the same group", {
  x <- t(iris[,-5])
  grp <- rep(0, ncol(x))
  expect_equal(suppressWarnings(oneway_equalvar(x=x, groups=grp)$obs.groups), rep(1, nrow(x)))
  expect_equal(suppressWarnings(oneway_welch(x=x, groups=grp)$obs.groups), rep(1, nrow(x)))
  expect_equal(suppressWarnings(kruskalwallis(x=x, groups=grp)$obs.groups), rep(1, nrow(x)))
  expect_equal(suppressWarnings(bartlett(x=x, groups=grp)$obs.groups), rep(1, nrow(x)))
  expect_equal(suppressWarnings(ievora(x=x, groups=grp)$obs.0), rep(ncol(x), nrow(x)))
})


test_that("groups can be character", {
  x <- t(iris[,-5])
  grp1 <- sample(c("A", "B"), ncol(x), replace=TRUE)
  grp2 <- ifelse(grp1==unique(grp1)[1], 0, 1)
  expect_equal(oneway_equalvar(x=x, groups=grp1), oneway_equalvar(x=x, groups=grp2))
  expect_equal(oneway_welch(x=x, groups=grp1), oneway_welch(x=x, groups=grp2))
  expect_equal(kruskalwallis(x=x, groups=grp1), kruskalwallis(x=x, groups=grp2))
  expect_equal(bartlett(x=x, groups=grp1), bartlett(x=x, groups=grp2))
  expect_equal(ievora(x=x, groups=grp1), ievora(x=x, groups=grp2))
})

test_that("groups can be factor", {
  x <- t(iris[,-5])
  grp1 <- sample(c("A", "B"), ncol(x), replace=TRUE)
  grp2 <- factor(grp1)
  expect_equal(oneway_equalvar(x=x, groups=grp1), oneway_equalvar(x=x, groups=grp2))
  expect_equal(oneway_welch(x=x, groups=grp1), oneway_welch(x=x, groups=grp2))
  expect_equal(kruskalwallis(x=x, groups=grp1), kruskalwallis(x=x, groups=grp2))
  expect_equal(bartlett(x=x, groups=grp1), bartlett(x=x, groups=grp2))
  expect_equal(ievora(x=x, groups=grp1), ievora(x=x, groups=grp2))
})

test_that("groups can be logical", {
  x <- t(iris[,-5])
  grp1 <- sample(c("A", "B"), ncol(x), replace=TRUE)
  grp2 <- grp1==unique(grp1)[2]
  expect_equal(oneway_equalvar(x=x, groups=grp1), oneway_equalvar(x=x, groups=grp2))
  expect_equal(oneway_welch(x=x, groups=grp1), oneway_welch(x=x, groups=grp2))
  expect_equal(kruskalwallis(x=x, groups=grp1), kruskalwallis(x=x, groups=grp2))
  expect_equal(bartlett(x=x, groups=grp1), bartlett(x=x, groups=grp2))
  expect_equal(ievora(x=x, groups=grp1), ievora(x=x, groups=grp2))
})

test_that("groups can be complex", {
  x <- t(iris[,-5])
  grp1 <- sample(c("A", "B"), ncol(x), replace=TRUE)
  grp2 <- ifelse(grp1=="A", complex(1,2,3), complex(1,3,2))
  expect_equal(oneway_equalvar(x=x, groups=grp1), oneway_equalvar(x=x, groups=grp2))
  expect_equal(oneway_welch(x=x, groups=grp1), oneway_welch(x=x, groups=grp2))
  expect_equal(kruskalwallis(x=x, groups=grp1), kruskalwallis(x=x, groups=grp2))
  expect_equal(bartlett(x=x, groups=grp1), bartlett(x=x, groups=grp2))
  expect_equal(ievora(x=x, groups=grp1), ievora(x=x, groups=grp2))
})

test_that("groups can be a list", {
  x <- t(iris[,-5])
  grp1 <- sample(c("A", "B"), ncol(x), replace=TRUE)
  grp2 <- ifelse(grp1=="A", list("A"), list("B"))
  expect_equal(oneway_equalvar(x=x, groups=grp1), oneway_equalvar(x=x, groups=grp2))
  expect_equal(oneway_welch(x=x, groups=grp1), oneway_welch(x=x, groups=grp2))
  expect_equal(kruskalwallis(x=x, groups=grp1), kruskalwallis(x=x, groups=grp2))
  expect_equal(bartlett(x=x, groups=grp1), bartlett(x=x, groups=grp2))
  expect_equal(ievora(x=x, groups=grp1), ievora(x=x, groups=grp2))
})

test_that("groups can be a matrix", {
  x <- t(iris[,-5])
  grp1 <- sample(c("A", "B"), ncol(x), replace=TRUE)
  grp2 <- matrix(grp1, ncol=1)
  expect_equal(oneway_equalvar(x=x, groups=grp1), oneway_equalvar(x=x, groups=grp2))
  expect_equal(oneway_welch(x=x, groups=grp1), oneway_welch(x=x, groups=grp2))
  expect_equal(kruskalwallis(x=x, groups=grp1), kruskalwallis(x=x, groups=grp2))
  expect_equal(bartlett(x=x, groups=grp1), bartlett(x=x, groups=grp2))
  expect_equal(ievora(x=x, groups=grp1), ievora(x=x, groups=grp2))
})

test_that("groups can be infinite", {
  x <- t(iris[,-5])
  grp1 <- sample(c("A", "B"), ncol(x), replace=TRUE)
  grp2 <- ifelse(grp1=="A", Inf, -Inf)
  expect_equal(oneway_equalvar(x=x, groups=grp1), oneway_equalvar(x=x, groups=grp2))
  expect_equal(oneway_welch(x=x, groups=grp1), oneway_welch(x=x, groups=grp2))
  expect_equal(kruskalwallis(x=x, groups=grp1), kruskalwallis(x=x, groups=grp2))
  expect_equal(bartlett(x=x, groups=grp1), bartlett(x=x, groups=grp2))
  expect_equal(ievora(x=x, groups=grp1), ievora(x=x, groups=grp2))
})

################################################################################
########################## ALLOWED ALTERNATIVE VALUES ##########################
################################################################################

test_that("alternative can be partially completed", {
  x <- t(iris[1:75,-5])
  y <- t(iris[76:150,-5])
  alt1 <- c("greater", "less", "two.sided", "two.sided")
  alt2 <- c("g", "l", "t", "t")
  # NA
  expect_equal(ttest_onegroup(x=x, alternative=alt1), ttest_onegroup(x=x, alternative=alt2))
  expect_equal(ttest_equalvar(x=x, y=y, alternative=alt1), ttest_equalvar(x=x, y=y, alternative=alt2))
  expect_equal(ttest_welch(x=x, y=y, alternative=alt1), ttest_welch(x=x, y=y, alternative=alt2))
  expect_equal(ttest_paired(x=x, y=y, alternative=alt1), ttest_paired(x=x, y=y, alternative=alt2))
})

################################################################################
############################## ALLOWED MU VALUES ###############################
################################################################################

test_that("MU can be Infinite", {
  x <- t(iris[1:75,-5])
  y <- t(iris[76:150,-5])
  # Inf
  expect_equal(ttest_onegroup(x=x, mu=Inf)$t.statistic, rep(-Inf, nrow(x)))
  expect_equal(ttest_equalvar(x=x, y=y, mu=Inf)$t.statistic, rep(-Inf, nrow(x)))
  expect_equal(ttest_welch(x=x, y=y, mu=Inf)$t.statistic, rep(-Inf, nrow(x)))
  expect_equal(ttest_paired(x=x, y=y, mu=Inf)$t.statistic, rep(-Inf, nrow(x)))
  # -Inf
  expect_equal(ttest_onegroup(x=x, mu=-Inf)$t.statistic, rep(Inf, nrow(x)))
  expect_equal(ttest_equalvar(x=x, y=y, mu=-Inf)$t.statistic, rep(Inf, nrow(x)))
  expect_equal(ttest_welch(x=x, y=y, mu=-Inf)$t.statistic, rep(Inf, nrow(x)))
  expect_equal(ttest_paired(x=x, y=y, mu=-Inf)$t.statistic, rep(Inf, nrow(x)))
})

