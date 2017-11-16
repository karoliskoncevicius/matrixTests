context("Errors in main parameters")

################################################################################
################################### MISSING ####################################
################################################################################

test_that("x cannot be missing", {
  er <- 'argument "x" is missing, with no default'
  expect_error(ttest_equalvar(), er)
  expect_error(ttest_welch(), er)
  expect_error(ttest_onegroup(), er)
  expect_error(ttest_paired(), er)
  expect_error(oneway_equalvar(), er)
  expect_error(bartlett(), er)
  expect_error(ievora(), er)
})

test_that("y cannot be missing", {
  er <- 'argument "y" is missing, with no default'
  expect_error(ttest_equalvar(x=NA), er)
  expect_error(ttest_welch(x=NA), er)
  expect_error(ttest_paired(x=NA), er)
})

test_that("groups cannot be missing", {
  er <- 'argument "groups" is missing, with no default'
  expect_error(oneway_equalvar(x=NA), er)
  expect_error(bartlett(x=NA), er)
  expect_error(ievora(x=NA), er)
})

################################################################################
################################# NON NUMERIC ##################################
################################################################################

test_that("x cannot be a character", {
  matX <- matrix(c("1", "2"), nrow=1)
  er <- '"x" must be a numeric matrix or vector'
  expect_error(ttest_onegroup(x=matX), er)
  expect_error(ttest_equalvar(x=matX, y=0), er)
  expect_error(ttest_welch(x=matX, y=0), er)
  expect_error(ttest_paired(x=matX, y=0), er)
  expect_error(oneway_equalvar(x=matX, groups="a"), er)
  expect_error(bartlett(x=matX, groups="a"), er)
  expect_error(ievora(x=matX, groups="a"), er)
})

test_that("y cannot be a character", {
  matX <- matrix(c("1", "2"), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(ttest_equalvar(x=0, y=matX), er)
  expect_error(ttest_welch(x=0, y=matX), er)
  expect_error(ttest_paired(x=0, y=matX), er)
})


test_that("x cannot be partially numeric", {
  er <- '"x" must be a numeric matrix or vector'
  expect_error(ttest_onegroup(x=iris), er)
  expect_error(ttest_equalvar(x=iris, y=0), er)
  expect_error(ttest_welch(x=iris, y=0), er)
  expect_error(ttest_paired(x=iris, y=0), er)
  expect_error(oneway_equalvar(x=iris, groups="a"), er)
  expect_error(bartlett(x=iris, groups="a"), er)
  expect_error(ievora(x=iris, groups="a"), er)
})

test_that("y cannot be partially numeric", {
  er <- '"y" must be a numeric matrix or vector'
  expect_error(ttest_equalvar(x=0, y=iris), er)
  expect_error(ttest_welch(x=0, y=iris), er)
  expect_error(ttest_paired(x=0, y=iris), er)
})



test_that("x cannot be complex", {
  matX <- matrix(complex(c(1,2), c(3,4)), nrow=1)
  er <- '"x" must be a numeric matrix or vector'
  expect_error(ttest_onegroup(x=matX), er)
  expect_error(ttest_equalvar(x=matX, y=0), er)
  expect_error(ttest_welch(x=matX, y=0), er)
  expect_error(ttest_paired(x=matX, y=0), er)
  expect_error(oneway_equalvar(x=matX, groups="a"), er)
  expect_error(bartlett(x=matX, groups="a"), er)
  expect_error(ievora(x=matX, groups="a"), er)
})

test_that("y cannot be complex", {
  matX <- matrix(complex(c(1,2), c(3,4)), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(ttest_equalvar(x=0, y=matX), er)
  expect_error(ttest_welch(x=0, y=matX), er)
  expect_error(ttest_paired(x=0, y=matX), er)
})


test_that("x cannot be logical", {
  matX <- matrix(c(TRUE, FALSE), nrow=1)
  er <- '"x" must be a numeric matrix or vector'
  expect_error(ttest_onegroup(x=matX), er)
  expect_error(ttest_equalvar(x=matX, y=0), er)
  expect_error(ttest_welch(x=matX, y=0), er)
  expect_error(ttest_paired(x=matX, y=0), er)
  expect_error(oneway_equalvar(x=matX, groups="a"), er)
  expect_error(bartlett(x=matX, groups="a"), er)
  expect_error(ievora(x=matX, groups="a"), er)
})

test_that("y cannot be logical", {
  matX <- matrix(c(TRUE, FALSE), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(ttest_equalvar(x=0, y=matX), er)
  expect_error(ttest_welch(x=0, y=matX), er)
  expect_error(ttest_paired(x=0, y=matX), er)
})


test_that("x cannot be NULL", {
  er <- '"x" must be a numeric matrix or vector'
  expect_error(ttest_onegroup(x=NULL), er)
  expect_error(ttest_equalvar(x=NULL, y=0), er)
  expect_error(ttest_welch(x=NULL, y=0), er)
  expect_error(ttest_paired(x=NULL, y=0), er)
  expect_error(oneway_equalvar(x=NULL, groups="a"), er)
  expect_error(bartlett(x=NULL, groups="a"), er)
  expect_error(ievora(x=NULL, groups="a"), er)
})

test_that("y cannot be NULL", {
  matX <- matrix(c(TRUE, FALSE), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(ttest_equalvar(x=0, y=NULL), er)
  expect_error(ttest_welch(x=0, y=NULL), er)
  expect_error(ttest_paired(x=0, y=NULL), er)
})


test_that("x cannot be in a list", {
  er <- '"x" must be a numeric matrix or vector'
  expect_error(ttest_onegroup(x=list(1:5)), er)
  expect_error(ttest_equalvar(x=list(1:5), y=0), er)
  expect_error(ttest_welch(x=list(1:5), y=0), er)
  expect_error(ttest_paired(x=list(1:5), y=0), er)
  expect_error(oneway_equalvar(x=list(1:5), groups="a"), er)
  expect_error(bartlett(x=list(1:5), groups="a"), er)
  expect_error(ievora(x=list(1:5), groups="a"), er)
})

test_that("y cannot be in a list", {
  matX <- matrix(c(TRUE, FALSE), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(ttest_equalvar(x=0, y=list(1:5)), er)
  expect_error(ttest_welch(x=0, y=list(1:5)), er)
  expect_error(ttest_paired(x=0, y=list(1:5)), er)
})

test_that("x cannot be a list", {
  er <- '"x" must be a numeric matrix or vector'
  expect_error(ttest_onegroup(x=as.list(1:5)), er)
  expect_error(ttest_equalvar(x=as.list(1:5), y=0), er)
  expect_error(ttest_welch(x=as.list(1:5), y=0), er)
  expect_error(ttest_paired(x=as.list(1:5), y=0), er)
  expect_error(oneway_equalvar(x=as.list(1:5), groups="a"), er)
  expect_error(bartlett(x=as.list(1:5), groups="a"), er)
  expect_error(ievora(x=as.list(1:5), groups="a"), er)
})

test_that("y cannot be a list", {
  matX <- matrix(c(TRUE, FALSE), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(ttest_equalvar(x=0, y=as.list(1:5)), er)
  expect_error(ttest_welch(x=0, y=as.list(1:5)), er)
  expect_error(ttest_paired(x=0, y=as.list(1:5)), er)
})

################################################################################
########################## WRONG GROUP SPECIFICATION ###########################
################################################################################

test_that("groups cannot be NULL", {
  matX <- matrix(1:12, ncol=3)
  er <- '"groups" must be a vector with length ncol\\(x\\)'
  expect_error(oneway_equalvar(x=matX, groups=NULL), er)
  expect_error(bartlett(x=matX, groups=NULL), er)
  expect_error(ievora(x=matX, groups=NULL), er)
})

test_that("groups cannot be a list", {
  matX <- matrix(1:12, ncol=3)
  er <- '"groups" must be a vector with length ncol\\(x\\)'
  expect_error(oneway_equalvar(x=matX, groups=list(1:3)), er)
  expect_error(bartlett(x=matX, groups=list(1:3)), er)
  expect_error(ievora(x=matX, groups=list(1:3)), er)
})

test_that("groups cannot be a matrix", {
  matX <- matrix(1:12, ncol=4)
  grp  <- cbind(c("A", "A"), c("B", "B"))
  er <- '"groups" must be a vector with length ncol\\(x\\)'
  expect_error(oneway_equalvar(x=matX, groups=grp), er)
  expect_error(bartlett(x=matX, groups=grp), er)
  expect_error(ievora(x=matX, groups=grp), er)
})

################################################################################
############################## DIMENSION MISMATCH ##############################
################################################################################

test_that("x and y has same number of rows", {
  matX <- matrix(1:9, nrow=3)
  matY <- matrix(1:16, nrow=4)
  er <- '"x" and "y" must have the same number of rows'
  expect_error(ttest_equalvar(x=matX, y=matY), er)
  expect_error(ttest_welch(x=matX, y=matY), er)
  expect_error(ttest_paired(x=matX, y=matY), er)
})

test_that("x and y has same number of columns", {
  matX <- matrix(1:9, nrow=3)
  matY <- matrix(1:12, nrow=3)
  er <- '"x" and "y" must have the same number of columns'
  expect_error(ttest_paired(x=matX, y=matY), er)
})

test_that("group length matches number of columns", {
  matX <- matrix(1:12, nrow=3)
  er <- '"groups" must be a vector with length ncol\\(x\\)'
  expect_error(oneway_equalvar(x=matX, groups=1:3), er)
  expect_error(bartlett(x=matX, groups=1:3), er)
  expect_error(ievora(x=matX, groups=1:3), er)
})

################################################################################
########################### SPECIAL GROUP EXCEPTIONS ###########################
################################################################################

test_that("groups have required number of groups", {
  matX <- matrix(1:12, nrow=3)
  er <- '"groups" must have no more than 2 unique elements'
  expect_error(ievora(x=matX, groups=c("A","B","C","C")), er)
  expect_error(ievora(x=c(1,2,3,4,NA), groups=c(0,0,1,1,2)), er)
})

