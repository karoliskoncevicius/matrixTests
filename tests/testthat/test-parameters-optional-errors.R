context("Errors in optional parameters")

################################################################################
######################## VALUES MUST HAVE CORRECT TYPE #########################
################################################################################

test_that("alternative must be a character", {
  er <- '"alternative" must be a character vector with length 1 or nrow\\(x\\)'
  mat <- matrix(1:12, nrow=3)
  # NULL
  expect_error(ttest_onegroup(x=mat, alternative=NULL), er)
  expect_error(ttest_equalvar(x=mat, y=mat, alternative=NULL), er)
  expect_error(ttest_welch(x=mat, y=mat, alternative=NULL), er)
  expect_error(ttest_paired(x=mat, y=mat, alternative=NULL), er)
  # NA
  expect_error(ttest_onegroup(x=mat, alternative=NA), er)
  expect_error(ttest_equalvar(x=mat, y=mat, alternative=NA), er)
  expect_error(ttest_welch(x=mat, y=mat, alternative=NA), er)
  expect_error(ttest_paired(x=mat, y=mat, alternative=NA), er)
  # numeric
  expect_error(ttest_onegroup(x=mat, alternative=1), er)
  expect_error(ttest_equalvar(x=mat, y=mat, alternative=1), er)
  expect_error(ttest_welch(x=mat, y=mat, alternative=1), er)
  expect_error(ttest_paired(x=mat, y=mat, alternative=1), er)
  # complex
  expect_error(ttest_onegroup(x=mat, alternative=complex(1)), er)
  expect_error(ttest_equalvar(x=mat, y=mat, alternative=complex(1)), er)
  expect_error(ttest_welch(x=mat, y=mat, alternative=complex(1)), er)
  expect_error(ttest_paired(x=mat, y=mat, alternative=complex(1)), er)
  # in list
  expect_error(ttest_onegroup(x=mat, alternative=list("less")), er)
  expect_error(ttest_equalvar(x=mat, y=mat, alternative=list("less")), er)
  expect_error(ttest_welch(x=mat, y=mat, alternative=list("less")), er)
  expect_error(ttest_paired(x=mat, y=mat, alternative=list("less")), er)
  # data frame
  expect_error(ttest_onegroup(x=mat, alternative=data.frame("less")), er)
  expect_error(ttest_equalvar(x=mat, y=mat, alternative=data.frame("less")), er)
  expect_error(ttest_welch(x=mat, y=mat, alternative=data.frame("less")), er)
  expect_error(ttest_paired(x=mat, y=mat, alternative=data.frame("less")), er)
})

test_that("mu must be numeric", {
  er <- '"mu" must be a numeric vector with length 1 or nrow\\(x\\)'
  mat <- matrix(1:12, nrow=3)
  # NULL
  expect_error(ttest_onegroup(x=mat, mu=NULL), er)
  expect_error(ttest_equalvar(x=mat, y=mat, mu=NULL), er)
  expect_error(ttest_welch(x=mat, y=mat, mu=NULL), er)
  expect_error(ttest_paired(x=mat, y=mat, mu=NULL), er)
  # NA
  expect_error(ttest_onegroup(x=mat, mu=NA), er)
  expect_error(ttest_equalvar(x=mat, y=mat, mu=NA), er)
  expect_error(ttest_welch(x=mat, y=mat, mu=NA), er)
  expect_error(ttest_paired(x=mat, y=mat, mu=NA), er)
  # character
  expect_error(ttest_onegroup(x=mat, mu="1"), er)
  expect_error(ttest_equalvar(x=mat, y=mat, mu="1"), er)
  expect_error(ttest_welch(x=mat, y=mat, mu="1"), er)
  expect_error(ttest_paired(x=mat, y=mat, mu="1"), er)
  # complex
  expect_error(ttest_onegroup(x=mat, mu=complex(1)), er)
  expect_error(ttest_equalvar(x=mat, y=mat, mu=complex(1)), er)
  expect_error(ttest_welch(x=mat, y=mat, mu=complex(1)), er)
  expect_error(ttest_paired(x=mat, y=mat, mu=complex(1)), er)
  # in list
  expect_error(ttest_onegroup(x=mat, mu=list(1)), er)
  expect_error(ttest_equalvar(x=mat, y=mat, mu=list(1)), er)
  expect_error(ttest_welch(x=mat, y=mat, mu=list(1)), er)
  expect_error(ttest_paired(x=mat, y=mat, mu=list(1)), er)
  # data frame
  expect_error(ttest_onegroup(x=mat, mu=data.frame(1)), er)
  expect_error(ttest_equalvar(x=mat, y=mat, mu=data.frame(1)), er)
  expect_error(ttest_welch(x=mat, y=mat, mu=data.frame(1)), er)
  expect_error(ttest_paired(x=mat, y=mat, mu=data.frame(1)), er)
})


test_that("conf.level must be numeric", {
  er <- '"conf.level" must be a numeric vector with length 1 or nrow\\(x\\)'
  mat <- matrix(1:12, nrow=3)
  # NULL
  expect_error(ttest_onegroup(x=mat, conf.level=NULL), er)
  expect_error(ttest_equalvar(x=mat, y=mat, conf.level=NULL), er)
  expect_error(ttest_welch(x=mat, y=mat, conf.level=NULL), er)
  expect_error(ttest_paired(x=mat, y=mat, conf.level=NULL), er)
  # NA
  expect_error(ttest_onegroup(x=mat, conf.level=NA), er)
  expect_error(ttest_equalvar(x=mat, y=mat, conf.level=NA), er)
  expect_error(ttest_welch(x=mat, y=mat, conf.level=NA), er)
  expect_error(ttest_paired(x=mat, y=mat, conf.level=NA), er)
  # character
  expect_error(ttest_onegroup(x=mat, conf.level="1"), er)
  expect_error(ttest_equalvar(x=mat, y=mat, conf.level="1"), er)
  expect_error(ttest_welch(x=mat, y=mat, conf.level="1"), er)
  expect_error(ttest_paired(x=mat, y=mat, conf.level="1"), er)
  # complex
  expect_error(ttest_onegroup(x=mat, conf.level=complex(1)), er)
  expect_error(ttest_equalvar(x=mat, y=mat, conf.level=complex(1)), er)
  expect_error(ttest_welch(x=mat, y=mat, conf.level=complex(1)), er)
  expect_error(ttest_paired(x=mat, y=mat, conf.level=complex(1)), er)
  # in list
  expect_error(ttest_onegroup(x=mat, conf.level=list(1)), er)
  expect_error(ttest_equalvar(x=mat, y=mat, conf.level=list(1)), er)
  expect_error(ttest_welch(x=mat, y=mat, conf.level=list(1)), er)
  expect_error(ttest_paired(x=mat, y=mat, conf.level=list(1)), er)
  # data frame
  expect_error(ttest_onegroup(x=mat, conf.level=data.frame(1)), er)
  expect_error(ttest_equalvar(x=mat, y=mat, conf.level=data.frame(1)), er)
  expect_error(ttest_welch(x=mat, y=mat, conf.level=data.frame(1)), er)
  expect_error(ttest_paired(x=mat, y=mat, conf.level=data.frame(1)), er)
})


test_that("p-value cutoffs must be numeric", {
  er1 <- '"cutT" must be a numeric vector with length 1'
  er2 <- '"cutBfdr" must be a numeric vector with length 1'
  mat <- matrix(1:12, nrow=3)
  grp <- c(1,1,0,0)
  # NULL
  expect_error(ievora(x=mat, groups=grp, cutT=NULL), er1)
  expect_error(ievora(x=mat, groups=grp, cutBfdr=NULL), er2)
  # NA
  expect_error(ievora(x=mat, groups=grp, cutT=NA), er1)
  expect_error(ievora(x=mat, groups=grp, cutBfdr=NA), er2)
  # character
  expect_error(ievora(x=mat, groups=grp, cutT="0.05"), er1)
  expect_error(ievora(x=mat, groups=grp, cutBfdr="0.05"), er2)
  # complex
  expect_error(ievora(x=mat, groups=grp, cutT=complex(1)), er1)
  expect_error(ievora(x=mat, groups=grp, cutBfdr=complex(1)), er2)
  # in list
  expect_error(ievora(x=mat, groups=grp, cutT=list(1)), er1)
  expect_error(ievora(x=mat, groups=grp, cutBfdr=list(1)), er2)
  # data frame
  expect_error(ievora(x=mat, groups=grp, cutT=data.frame(1)), er1)
  expect_error(ievora(x=mat, groups=grp, cutBfdr=data.frame(1)), er2)
})

################################################################################
##################### VALUES MUST HAVE CORRECT DIMENSIONS ######################
################################################################################

test_that("alternative has correct dimensions", {
  er <- '"alternative" must be a character vector with length 1 or nrow\\(x\\)'
  mat <- matrix(1:12, ncol=3)
  # too short
  expect_error(ttest_onegroup(x=mat, alternative=c("less", "greater")), er)
  expect_error(ttest_equalvar(x=mat, y=mat, alternative=c("less", "greater")), er)
  expect_error(ttest_welch(x=mat, y=mat, alternative=c("less", "greater")), er)
  expect_error(ttest_paired(x=mat, y=mat, alternative=c("less", "greater")), er)
  # too long
  expect_error(ttest_onegroup(x=mat, alternative=rep("less", 5)), er)
  expect_error(ttest_equalvar(x=mat, y=mat, alternative=rep("less", 5)), er)
  expect_error(ttest_welch(x=mat, y=mat, alternative=rep("less", 5)), er)
  expect_error(ttest_paired(x=mat, y=mat, alternative=rep("less", 5)), er)
  # matrix format
  alt <- matrix(rep("less", 4), ncol=2)
  expect_error(ttest_onegroup(x=mat, alternative=alt), er)
  expect_error(ttest_equalvar(x=mat, y=mat, alternative=alt), er)
  expect_error(ttest_welch(x=mat, y=mat, alternative=alt), er)
  expect_error(ttest_paired(x=mat, y=mat, alternative=alt), er)
})

test_that("mu has correct dimensions", {
  er <- '"mu" must be a numeric vector with length 1 or nrow\\(x\\)'
  mat <- matrix(1:12, ncol=3)
  # too short
  expect_error(ttest_onegroup(x=mat, mu=c(1,2)), er)
  expect_error(ttest_equalvar(x=mat, y=mat, mu=c(1,2)), er)
  expect_error(ttest_welch(x=mat, y=mat, mu=c(1,2)), er)
  expect_error(ttest_paired(x=mat, y=mat, mu=c(1,2)), er)
  # too long
  expect_error(ttest_onegroup(x=mat, mu=rep(0, 5)), er)
  expect_error(ttest_equalvar(x=mat, y=mat, mu=rep(0, 5)), er)
  expect_error(ttest_welch(x=mat, y=mat, mu=rep(0, 5)), er)
  expect_error(ttest_paired(x=mat, y=mat, mu=rep(0, 5)), er)
  # matrix format
  mus <- matrix(rep(1, 4), ncol=2)
  expect_error(ttest_onegroup(x=mat, mu=mus), er)
  expect_error(ttest_equalvar(x=mat, y=mat, mu=mus), er)
  expect_error(ttest_welch(x=mat, y=mat, mu=mus), er)
  expect_error(ttest_paired(x=mat, y=mat, mu=mus), er)
})

test_that("conf.level has correct dimensions", {
  er <- '"conf.level" must be a numeric vector with length 1 or nrow\\(x\\)'
  mat <- matrix(1:12, ncol=3)
  # too short
  expect_error(ttest_onegroup(x=mat, conf.level=c(1,2)), er)
  expect_error(ttest_equalvar(x=mat, y=mat, conf.level=c(1,2)), er)
  expect_error(ttest_welch(x=mat, y=mat, conf.level=c(1,2)), er)
  expect_error(ttest_paired(x=mat, y=mat, conf.level=c(1,2)), er)
  # too long
  expect_error(ttest_onegroup(x=mat, conf.level=rep(0, 5)), er)
  expect_error(ttest_equalvar(x=mat, y=mat, conf.level=rep(0, 5)), er)
  expect_error(ttest_welch(x=mat, y=mat, conf.level=rep(0, 5)), er)
  expect_error(ttest_paired(x=mat, y=mat, conf.level=rep(0, 5)), er)
  # matrix format
  cfs <- matrix(rep(1, 4), ncol=2)
  expect_error(ttest_onegroup(x=mat, conf.level=cfs), er)
  expect_error(ttest_equalvar(x=mat, y=mat, conf.level=cfs), er)
  expect_error(ttest_welch(x=mat, y=mat, conf.level=cfs), er)
  expect_error(ttest_paired(x=mat, y=mat, conf.level=cfs), er)
})

test_that("p-value cutoffs have the right dimensions", {
  er1 <- '"cutT" must be a numeric vector with length 1'
  er2 <- '"cutBfdr" must be a numeric vector with length 1'
  mat <- matrix(1:12, ncol=3)
  grp <- c(0,0,1)
  # too long
  expect_error(ievora(x=mat, groups=grp, cutT=c(1,2)), er1)
  expect_error(ievora(x=mat, groups=grp, cutBfdr=c(1,2)), er2)
  # matrix format
  cts <- matrix(rep(1, 4), ncol=2)
  expect_error(ievora(x=mat, groups=grp, cutT=cts), er1)
  expect_error(ievora(x=mat, groups=grp, cutBfdr=cts), er2)
})

################################################################################
######################## VALUES MUST BE IN CORRECT SET #########################
################################################################################

test_that("alternative is in: less, greater, two-sided)", {
  er <- 'all "alternative" values must be in: two\\.sided, less, greater'
  mat <- matrix(1:12, nrow=3)
  # one value
  expect_error(ttest_onegroup(x=mat, alternative="ga"), er)
  expect_error(ttest_equalvar(x=mat, y=mat, alternative="ga"), er)
  expect_error(ttest_welch(x=mat, y=mat, alternative="ga"), er)
  expect_error(ttest_paired(x=mat, y=mat, alternative="ga"), er)
  # for each row and one incorrect
  expect_error(ttest_onegroup(x=mat, alternative=c("t","l","c")), er)
  expect_error(ttest_equalvar(x=mat, y=mat, alternative=c("t","l","c")), er)
  expect_error(ttest_welch(x=mat, y=mat, alternative=c("t","l","c")), er)
  expect_error(ttest_paired(x=mat, y=mat, alternative=c("t","l","c")), er)
})

test_that("conf.level is in: 0-1)", {
  er <- 'all "conf.level" values must be between: 0 and 1'
  mat <- matrix(1:12, nrow=3)
  # slightly below
  expect_error(ttest_onegroup(x=mat, conf.level=-0.001), er)
  expect_error(ttest_equalvar(x=mat, y=mat, conf.level=-0.001), er)
  expect_error(ttest_welch(x=mat, y=mat, conf.level=-0.001), er)
  expect_error(ttest_paired(x=mat, y=mat, conf.level=-0.001), er)
  # slightly above
  expect_error(ttest_onegroup(x=mat, conf.level=1.001), er)
  expect_error(ttest_equalvar(x=mat, y=mat, conf.level=1.001), er)
  expect_error(ttest_welch(x=mat, y=mat, conf.level=1.001), er)
  expect_error(ttest_paired(x=mat, y=mat, conf.level=1.001), er)
  # special values
  expect_error(ttest_onegroup(x=mat, conf.level=NA_integer_), er)
  expect_error(ttest_equalvar(x=mat, y=mat, conf.level=NaN), er)
  expect_error(ttest_welch(x=mat, y=mat, conf.level=Inf), er)
  expect_error(ttest_paired(x=mat, y=mat, conf.level=-Inf), er)
})

test_that("mu is in: -Inf:Inf)", {
  er <- 'all "mu" values must be between: -Inf and Inf'
  mat <- matrix(1:12, nrow=3)
  # NA
  expect_error(ttest_onegroup(x=mat, mu=NA_integer_), er)
  expect_error(ttest_equalvar(x=mat, y=mat, mu=NA_integer_), er)
  expect_error(ttest_welch(x=mat, y=mat, mu=NA_integer_), er)
  expect_error(ttest_paired(x=mat, y=mat, mu=NA_integer_), er)
  # NaN
  expect_error(ttest_onegroup(x=mat, mu=NaN), er)
  expect_error(ttest_equalvar(x=mat, y=mat, mu=NaN), er)
  expect_error(ttest_welch(x=mat, y=mat, mu=NaN), er)
  expect_error(ttest_paired(x=mat, y=mat, mu=NaN), er)
})

test_that("p-value cut-offs must be in 0-1", {
  er1 <- 'all "cutT" values must be between: 0 and 1'
  er2 <- 'all "cutBfdr" values must be between: 0 and 1'
  mat <- matrix(1:12, nrow=3)
  grp <- c(0,0,1,1)
  # slightly below
  expect_error(ievora(x=mat, groups=grp, cutT=-0.001), er1)
  expect_error(ievora(x=mat, groups=grp, cutBfdr=-0.001), er2)
  # slightly above
  expect_error(ievora(x=mat, groups=grp, cutT=1.001), er1)
  expect_error(ievora(x=mat, groups=grp, cutBfdr=1.001), er2)
  # special values
  expect_error(ievora(x=mat, groups=grp, cutT=NA_integer_), er1)
  expect_error(ievora(x=mat, groups=grp, cutBfdr=NaN), er2)
  expect_error(ievora(x=mat, groups=grp, cutT=Inf), er1)
  expect_error(ievora(x=mat, groups=grp, cutBfdr=-Inf), er2)
})


