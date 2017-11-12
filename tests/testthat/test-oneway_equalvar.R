context("oneway_equalvar")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_oneway <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)
  res <- Map(function(x, g) summary(aov(x~g))[[1]], split(mat, row(mat)), list(groups))

  gs <- matrix(groups, ncol=length(groups), nrow=nrow(mat), byrow=TRUE)
  gs[is.na(mat)] <- NA

  nGroups  <- apply(gs, 1, function(x) length(unique(na.omit(x))))
  nSamples <- apply(gs, 1, function(x) sum(table(x)))

  data.frame(sum.sq.treatment=sapply(res, `[`, 1, 2),
             sum.sq.residuals=sapply(res, `[`, 2, 2),
             mean.sq.treatment=sapply(res, `[`, 1, 3),
             mean.sq.residuals=sapply(res, `[`, 2, 3),
             obs.tot=nSamples,
             obs.groups=nGroups,
             df.treatment=sapply(res, `[`, 1, 1),
             df.residuals=sapply(res, `[`, 2, 1),
             F.statistic=sapply(res, `[`, 1, 4),
             p.value=sapply(res, `[`, 1, 5),
             row.names=rownames(mat)
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

test_that("edge cases give equal results", {
  x <- rnorm(4, sd=0.000001); g <- c("a","b","a","b")
  expect_equal(base_oneway(x, g), oneway_equalvar(x, g))
  x <- c(1,1,1,1); g <- c("a","b","a","b")
  expect_equal(base_oneway(x, g), suppressWarnings(oneway_equalvar(x, g)))
})

test_that("missing values in x are dealt with", {
  x <- t(iris[,1:4])
  g <- iris$Species
  x[1,c(1:10,51:60,101:110)] <- NA
  x[2,1:50] <- NA
  x[3,1:98] <- NA
  x[4,c(3:50,53:100,103:150)] <- NA
  expect_equal(base_oneway(x, g), oneway_equalvar(x, g))
})

test_that("missing values in groups are dealt with", {
  x <- rnorm(9)
  g <- rep(letters[1:3], each=3)
  expect_equal(base_oneway(x, g), suppressWarnings(oneway_equalvar(x, g)))
  gg <- ifelse(g=="a", NA, g)
  expect_equal(base_oneway(x, gg), suppressWarnings(oneway_equalvar(x, gg)))
  gg <- g; gg[c(1,4,7)] <- NA
  expect_equal(base_oneway(x, gg), suppressWarnings(oneway_equalvar(x, gg)))
})

test_that("missing values in both are dealth with", {
  x <- t(iris[,1:4]); x[sample(length(x), 100)] <- NA
  g <- iris$Species; g[sample(length(g), 100)] <- NA
  expect_equal(base_oneway(x, g), suppressWarnings(oneway_equalvar(x, g)))
})

test_that("groups with one element remaining are dealth with correctly", {
  x <- rnorm(12)
  g <- rep(letters[1:4], each=3)
  gg <- g; gg[1:2] <- NA
  expect_equal(base_oneway(x, gg), suppressWarnings(oneway_equalvar(x, gg)))
  gg <- g; gg[c(1:2,4:5,7,10)] <- NA
  expect_equal(base_oneway(x, gg), suppressWarnings(oneway_equalvar(x, gg)))
  gg <- g; gg[c(1:2,4:5,7:8,10)] <- NA
  expect_equal(base_oneway(x, gg), suppressWarnings(oneway_equalvar(x, gg)))
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("unsolvable situations produce warnings", {
  er <- "1 rows had less than 2 groups with enough observations"
  expect_warning(oneway_equalvar(NA_integer_, "a"), er)
  expect_warning(oneway_equalvar(rep(NA_integer_, 4), c("a", "a", "b", "b")), er)
  expect_warning(oneway_equalvar(1:2, c("a", "a")), er)
  er <- "1 rows had one observation per group"
  expect_warning(oneway_equalvar(1:3, c("a", "b", "c")), er)
  expect_warning(oneway_equalvar(1:3, c("a", NA, NA)), er)
  er <- "1 rows had essentially perfect fit"
  expect_warning(oneway_equalvar(c(1,1:3), c("a", "a", "b", "c")), er)
  er <- "1 columns dropped due to missing group information"
  expect_warning(oneway_equalvar(1:3, c("a", "b", NA)), er)
  er <- "2 columns dropped due to missing group information"
  expect_warning(oneway_equalvar(1:6, c("a", "a", "b", "b", NA, NA)), er)
})

################################################################################
############################### TEST PARAMETERS ################################
################################################################################

###################################### X #######################################

test_that("x cannot be missing", {
  er <- 'argument "x" is missing, with no default'
  expect_error(oneway_equalvar(), er)
})

test_that("invalid x values produce errors", {
  er <- "'x' must be a numeric matrix or a vector"
  expect_error(oneway_equalvar(iris), er)
  expect_error(oneway_equalvar(matrix(LETTERS[1:10], ncol=2)), er)
  expect_error(oneway_equalvar(list(1:10)), er)
  expect_error(oneway_equalvar(NA), er)
  expect_error(oneway_equalvar(NULL), er)
  expect_error(oneway_equalvar(matrix(nrow=0, ncol=5)), er)
  expect_error(oneway_equalvar(matrix(nrow=5, ncol=0)), er)
})

test_that("data.frame values pass when they are all numeric", {
  expect_error(oneway_equalvar(t(iris[,-5]), iris$Species), NA)
})

test_that("zero rows matrix pass", {
  expect_error(oneway_equalvar(matrix(NA_integer_, nrow=0, ncol=2), c("a", "a")), NA)
})

#################################### GROUPS ####################################

test_that("groups cannot be missing", {
  er <- 'argument "groups" is missing, with no default'
  expect_error(oneway_equalvar(1:5), er)
})

test_that("invalid groups values produce errors", {
  er <- "'groups' must be a character vector with length ncol\\(x\\)"
  expect_error(oneway_equalvar(1:5, NULL), er)
  expect_error(oneway_equalvar(1:5, rep("a",4)), er)
  expect_error(oneway_equalvar(1:5, list(rep("a", 5))), er)
  expect_error(oneway_equalvar(matrix(1:8, nrow=2), matrix(rep(c("a", "b"), 4), ncol=4)), er)
})

test_that("group values convertable to character pass", {
  expect_error(oneway_equalvar(1:4, factor(rep(c("a", "b"), 2))), NA)
  expect_error(oneway_equalvar(1:4, as.list(rep(c("a", "b"), 2))), NA)
  expect_error(oneway_equalvar(1:4, as.list(rep(c(1, 2), 2))), NA)
  expect_error(oneway_equalvar(1:4, matrix(rep(c(1, 2), 2), nrow=1)), NA)
  expect_error(oneway_equalvar(1:4, matrix(rep(c(1, 2), 2), ncol=1)), NA)
})

