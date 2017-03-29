context("bartlett")

################################################################################
############################### HELPER FUNCTIONS ###############################
################################################################################

base_bartlett <- function(mat, groups) {
  if(is.vector(mat)) mat <- matrix(mat, nrow=1)
  res <- Map(bartlett.test, split(mat, row(mat)), list(groups))

  gs <- matrix(groups, ncol=length(groups), nrow=nrow(mat), byrow=TRUE)
  gs[is.na(mat)] <- NA

  nGroups  <- apply(gs, 1, function(x) sum(table(x) > 1))
  nSamples <- apply(gs, 1, function(x) sum(table(x)[table(x)>1]))

  nPerGroup <- apply(gs, 1, function(x) table(factor(x, levels=unique(groups))))
  nPerGroup <- t(nPerGroup)
  nPerGroup[nPerGroup < 2] <- NA

  vtot <- t(apply(mat, 1, tapply, factor(groups, levels=unique(groups)), var, na.rm=TRUE))
  vtot <- rowSums(vtot * (nPerGroup-1), na.rm=TRUE) / (nSamples-nGroups)

  data.frame(var.tot=vtot,
             obs.tot=nSamples,
             obs.groups=nGroups,
             ksq.statistic=unname(sapply(res, getElement, "statistic")),
             p.value=unname(sapply(res, getElement, "p.value")),
             df=unname(sapply(res, getElement, "parameter")),
             stringsAsFactors=FALSE
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

test_that("edge cases give equal results", {
  x <- rnorm(4, sd=0.000001); g <- c("a","b","a","b")
  expect_equal(base_bartlett(x, g), bartlett(x, g))
  x <- c(1,1,1,1); g <- c("a","b","a","b")
  expect_equal(base_bartlett(x, g), bartlett(x, g))
})

test_that("missing values in x are dealt with", {
  x <- t(iris[,1:4])
  g <- iris$Species
  x[1,c(1:10,51:60,101:110)] <- NA
  x[2,1:50] <- NA
  x[3,1:98] <- NA
  x[4,c(3:50,53:100,103:150)] <- NA
  expect_equal(base_bartlett(x, g), suppressWarnings(bartlett(x, g)))
})

test_that("missing values in groups are dealt with", {
  x <- rnorm(9)
  g <- rep(letters[1:3], each=3)
  expect_equal(base_bartlett(x, g), suppressWarnings(bartlett(x, g)))
  gg <- ifelse(g=="a", NA, g)
  expect_equal(base_bartlett(x, gg), suppressWarnings(bartlett(x, gg)))
  gg <- g; gg[c(1,4,7)] <- NA
  expect_equal(base_bartlett(x, gg), suppressWarnings(bartlett(x, gg)))
})

test_that("missing values in both are dealth with", {
  x <- t(iris[,1:4]); x[sample(length(x), 100)] <- NA
  g <- iris$Species; g[sample(length(g), 100)] <- NA
  expect_equal(base_bartlett(x, g), suppressWarnings(bartlett(x, g)))
})

test_that("groups with one element remaining are dropped correctly", {
  x <- rnorm(12)
  g <- rep(letters[1:4], each=3)
  gg <- g; gg[1:2] <- NA
  expect_equal(base_bartlett(x[-c(1:3)], gg[-c(1:3)]), suppressWarnings(bartlett(x, gg)))
  gg <- g; gg[c(1:2,4:5,7,10)] <- NA
  expect_equal(base_bartlett(x[-c(1:6)], gg[-c(1:6)]), suppressWarnings(bartlett(x, gg)))
})

################################################################################
################################ TEST WARNINGS #################################
################################################################################

test_that("unsolvable situations produce warnings", {
  er <- "1 rows had less than 2 groups with enough observations"
  expect_warning(bartlett(NA_integer_, "a"), er)
  expect_warning(bartlett(rep(NA_integer_, 4), c("a", "a", "b", "b")), er)
  expect_warning(bartlett(1:2, c("a", "a")), er)
  expect_warning(bartlett(1:3, c("a", "a", "b")), er)
  expect_warning(bartlett(1:3, c("a", "b", "c")), er)
  expect_warning(bartlett(1:3, c("a", NA, NA)), er)
  er <- "1 columns skipped due to missing group information"
  expect_warning(bartlett(1:3, c("a", "b", NA)), er)
  er <- "2 columns skipped due to missing group information"
  expect_warning(bartlett(1:6, c("a", "a", "b", "b", NA, NA)), er)
  er <- "1 rows had groups with less than 2 observations"
  expect_warning(bartlett(1:5, c("a", "a", "b", "b", "c")), er)
  expect_warning(bartlett(c(1:5,NA), c("a", "a", "b", "b", "c", "c")), er)
})

################################################################################
############################### TEST PARAMETERS ################################
################################################################################

###################################### X #######################################

test_that("x cannot be missing", {
  er <- 'argument "x" is missing, with no default'
  expect_error(bartlett(), er)
})

test_that("invalid x values produce errors", {
  er <- "'x' must be a numeric matrix or a vector"
  expect_error(bartlett(iris), er)
  expect_error(bartlett(matrix(LETTERS[1:10], ncol=2)), er)
  expect_error(bartlett(list(1:10)), er)
  expect_error(bartlett(NA), er)
  expect_error(bartlett(NULL), er)
  expect_error(bartlett(matrix(nrow=0, ncol=5)), er)
  expect_error(bartlett(matrix(nrow=5, ncol=0)), er)
})

test_that("data.frame values pass when they are all numeric", {
  expect_error(bartlett(t(iris[,-5]), iris$Species), NA)
})

test_that("zero rows matrix pass", {
  expect_error(bartlett(matrix(NA_integer_, nrow=0, ncol=2), c("a", "a")), NA)
})

#################################### GROUPS ####################################

test_that("groups cannot be missing", {
  er <- 'argument "groups" is missing, with no default'
  expect_error(bartlett(1:5), er)
})

test_that("invalid groups values produce errors", {
  er <- "'groups' must be a character vector with length ncol\\(x\\)"
  expect_error(bartlett(1:5, NULL), er)
  expect_error(bartlett(1:5, rep("a",4)), er)
  expect_error(bartlett(1:5, list(rep("a", 5))), er)
  expect_error(bartlett(matrix(1:8, nrow=2), matrix(rep(c("a", "b"), 4), ncol=4)), er)
})

test_that("group values convertable to character pass", {
  expect_error(bartlett(1:4, factor(rep(c("a", "b"), 2))), NA)
  expect_error(bartlett(1:4, as.list(rep(c("a", "b"), 2))), NA)
  expect_error(bartlett(1:4, as.list(rep(c(1, 2), 2))), NA)
  expect_error(bartlett(1:4, matrix(rep(c(1, 2), 2), nrow=1)), NA)
  expect_error(bartlett(1:4, matrix(rep(c(1, 2), 2), ncol=1)), NA)
})

