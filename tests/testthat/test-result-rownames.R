context("Result rownames")

test_that("when no row-names in the input - numbers are added", {
  X <- matrix(rnorm(100), nrow=10)
  Y <- matrix(rnorm(100), nrow=10)
  grp <- sample(c(1,0), 10, replace=TRUE)
  rnames <- as.character(1:nrow(X))
  expect_equal(rownames(ttest_onegroup(x=X)), rnames)
  expect_equal(rownames(ttest_equalvar(x=X, y=Y)), rnames)
  expect_equal(rownames(ttest_welch(x=X, y=Y)), rnames)
  expect_equal(rownames(ttest_paired(x=X, y=Y)), rnames)
  expect_equal(rownames(oneway_equalvar(x=X, groups=grp)), rnames)
  expect_equal(rownames(bartlett(x=X, groups=grp)), rnames)
  expect_equal(rownames(ievora(x=X, groups=grp)), rnames)
})

test_that("when X doesn't have rownames - names from Y or groups are not used", {
  X <- matrix(rnorm(100), nrow=10)
  Y <- matrix(rnorm(100), nrow=10)
  grp <- sample(c(1,0), 10, replace=TRUE)
  rownames(Y) <- LETTERS[1:10]
  names(grp) <- LETTERS[1:10]
  rnames <- as.character(1:nrow(X))
  expect_equal(rownames(ttest_onegroup(x=X)), rnames)
  expect_equal(rownames(ttest_equalvar(x=X, y=Y)), rnames)
  expect_equal(rownames(ttest_welch(x=X, y=Y)), rnames)
  expect_equal(rownames(ttest_paired(x=X, y=Y)), rnames)
  expect_equal(rownames(oneway_equalvar(x=X, groups=grp)), rnames)
  expect_equal(rownames(bartlett(x=X, groups=grp)), rnames)
  expect_equal(rownames(ievora(x=X, groups=grp)), rnames)
})

test_that("when row-names are specified - they are preserved", {
  # matrix case
  X <- matrix(rnorm(100), nrow=10)
  rownames(X) <- LETTERS[1:10]
  Y <- matrix(rnorm(100), nrow=10)
  grp <- sample(c(1,0), 10, replace=TRUE)
  rnames <- rownames(X)
  expect_equal(rownames(ttest_onegroup(x=X)), rnames)
  expect_equal(rownames(ttest_equalvar(x=X, y=Y)), rnames)
  expect_equal(rownames(ttest_welch(x=X, y=Y)), rnames)
  expect_equal(rownames(ttest_paired(x=X, y=Y)), rnames)
  expect_equal(rownames(oneway_equalvar(x=X, groups=grp)), rnames)
  expect_equal(rownames(bartlett(x=X, groups=grp)), rnames)
  expect_equal(rownames(ievora(x=X, groups=grp)), rnames)
  # data.frame case
  X <- as.data.frame(X)
  Y <- as.data.frame(Y)
  rnames <- rownames(X)
  expect_equal(rownames(ttest_onegroup(x=X)), rnames)
  expect_equal(rownames(ttest_equalvar(x=X, y=Y)), rnames)
  expect_equal(rownames(ttest_welch(x=X, y=Y)), rnames)
  expect_equal(rownames(ttest_paired(x=X, y=Y)), rnames)
  expect_equal(rownames(oneway_equalvar(x=X, groups=grp)), rnames)
  expect_equal(rownames(bartlett(x=X, groups=grp)), rnames)
  expect_equal(rownames(ievora(x=X, groups=grp)), rnames)
})

test_that("when row-names are duplicated - they are modified to be unique", {
  X <- matrix(rnorm(100), nrow=10)
  Y <- matrix(rnorm(100), nrow=10)
  grp <- sample(c(1,0), 10, replace=TRUE)
  rownames(X) <- c(rep("A",5), rep("B", 5))
  rnames <- make.unique(rownames(X))
  expect_equal(rownames(ttest_onegroup(x=X)), rnames)
  expect_equal(rownames(ttest_equalvar(x=X, y=Y)), rnames)
  expect_equal(rownames(ttest_welch(x=X, y=Y)), rnames)
  expect_equal(rownames(ttest_paired(x=X, y=Y)), rnames)
  expect_equal(rownames(oneway_equalvar(x=X, groups=grp)), rnames)
  expect_equal(rownames(bartlett(x=X, groups=grp)), rnames)
  expect_equal(rownames(ievora(x=X, groups=grp)), rnames)
})

