context("Errors in main parameters")

################################################################################
################################### MISSING ####################################
################################################################################

test_that("x cannot be missing", {
  er <- 'argument "x" is missing, with no default'
  expect_error(row.t.equalvar(), er)
  expect_error(row.t.welch(), er)
  expect_error(row.t.onesample(), er)
  expect_error(row.t.paired(), er)
  expect_error(row.oneway.equalvar(), er)
  expect_error(row.oneway.welch(), er)
  expect_error(row.kruskalwallis(), er)
  expect_error(row.bartlett(), er)
  expect_error(row.cor.pearson(), er)
  expect_error(row.ievora(), er)
})

test_that("y cannot be missing", {
  er <- 'argument "y" is missing, with no default'
  expect_error(row.t.equalvar(x=NA), er)
  expect_error(row.t.welch(x=NA), er)
  expect_error(row.t.paired(x=NA), er)
  expect_error(row.cor.pearson(x=NA), er)
})

test_that("groups cannot be missing", {
  er <- 'argument "groups" is missing, with no default'
  expect_error(row.oneway.equalvar(x=NA), er)
  expect_error(row.oneway.welch(x=NA), er)
  expect_error(row.kruskalwallis(x=NA), er)
  expect_error(row.bartlett(x=NA), er)
  expect_error(row.ievora(x=NA), er)
})

################################################################################
################################# NON NUMERIC ##################################
################################################################################

test_that("x cannot be a character", {
  matX <- matrix(c("1", "2"), nrow=1)
  er <- '"x" must be a numeric matrix or vector'
  expect_error(row.t.onesample(x=matX), er)
  expect_error(row.t.equalvar(x=matX, y=0), er)
  expect_error(row.t.welch(x=matX, y=0), er)
  expect_error(row.t.paired(x=matX, y=0), er)
  expect_error(row.oneway.equalvar(x=matX, groups="a"), er)
  expect_error(row.oneway.welch(x=matX, groups="a"), er)
  expect_error(row.kruskalwallis(x=matX, groups="a"), er)
  expect_error(row.bartlett(x=matX, groups="a"), er)
  expect_error(row.cor.pearson(x=matX, y=0), er)
  expect_error(row.ievora(x=matX, groups="a"), er)
})

test_that("y cannot be a character", {
  matX <- matrix(c("1", "2"), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(row.t.equalvar(x=0, y=matX), er)
  expect_error(row.t.welch(x=0, y=matX), er)
  expect_error(row.t.paired(x=0, y=matX), er)
  expect_error(row.cor.pearson(x=0, y=matX), er)
})


test_that("x cannot be partially numeric", {
  er <- '"x" must be a numeric matrix or vector'
  expect_error(row.t.onesample(x=iris), er)
  expect_error(row.t.equalvar(x=iris, y=0), er)
  expect_error(row.t.welch(x=iris, y=0), er)
  expect_error(row.t.paired(x=iris, y=0), er)
  expect_error(row.oneway.equalvar(x=iris, groups="a"), er)
  expect_error(row.oneway.welch(x=iris, groups="a"), er)
  expect_error(row.kruskalwallis(x=iris, groups="a"), er)
  expect_error(row.bartlett(x=iris, groups="a"), er)
  expect_error(row.cor.pearson(x=iris, y=0), er)
  expect_error(row.ievora(x=iris, groups="a"), er)
})

test_that("y cannot be partially numeric", {
  er <- '"y" must be a numeric matrix or vector'
  expect_error(row.t.equalvar(x=0, y=iris), er)
  expect_error(row.t.welch(x=0, y=iris), er)
  expect_error(row.t.paired(x=0, y=iris), er)
  expect_error(row.cor.pearson(x=0, y=iris), er)
})


test_that("x cannot be complex", {
  matX <- matrix(complex(c(1,2), c(3,4)), nrow=1)
  er <- '"x" must be a numeric matrix or vector'
  expect_error(row.t.onesample(x=matX), er)
  expect_error(row.t.equalvar(x=matX, y=0), er)
  expect_error(row.t.welch(x=matX, y=0), er)
  expect_error(row.t.paired(x=matX, y=0), er)
  expect_error(row.oneway.equalvar(x=matX, groups="a"), er)
  expect_error(row.oneway.welch(x=matX, groups="a"), er)
  expect_error(row.kruskalwallis(x=matX, groups="a"), er)
  expect_error(row.bartlett(x=matX, groups="a"), er)
  expect_error(row.cor.pearson(x=matX, y=0), er)
  expect_error(row.ievora(x=matX, groups="a"), er)
})

test_that("y cannot be complex", {
  matX <- matrix(complex(c(1,2), c(3,4)), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(row.t.equalvar(x=0, y=matX), er)
  expect_error(row.t.welch(x=0, y=matX), er)
  expect_error(row.t.paired(x=0, y=matX), er)
  expect_error(row.cor.pearson(x=0, y=matX), er)
})


test_that("x cannot be logical", {
  matX <- matrix(c(TRUE, FALSE), nrow=1)
  er <- '"x" must be a numeric matrix or vector'
  expect_error(row.t.onesample(x=matX), er)
  expect_error(row.t.equalvar(x=matX, y=0), er)
  expect_error(row.t.welch(x=matX, y=0), er)
  expect_error(row.t.paired(x=matX, y=0), er)
  expect_error(row.oneway.equalvar(x=matX, groups="a"), er)
  expect_error(row.oneway.welch(x=matX, groups="a"), er)
  expect_error(row.kruskalwallis(x=matX, groups="a"), er)
  expect_error(row.bartlett(x=matX, groups="a"), er)
  expect_error(row.cor.pearson(x=matX, y=0), er)
  expect_error(row.ievora(x=matX, groups="a"), er)
})

test_that("y cannot be logical", {
  matX <- matrix(c(TRUE, FALSE), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(row.t.equalvar(x=0, y=matX), er)
  expect_error(row.t.welch(x=0, y=matX), er)
  expect_error(row.t.paired(x=0, y=matX), er)
  expect_error(row.cor.pearson(x=0, y=matX), er)
})


test_that("x cannot be NULL", {
  er <- '"x" must be a numeric matrix or vector'
  expect_error(row.t.onesample(x=NULL), er)
  expect_error(row.t.equalvar(x=NULL, y=0), er)
  expect_error(row.t.welch(x=NULL, y=0), er)
  expect_error(row.t.paired(x=NULL, y=0), er)
  expect_error(row.oneway.equalvar(x=NULL, groups="a"), er)
  expect_error(row.oneway.welch(x=NULL, groups="a"), er)
  expect_error(row.kruskalwallis(x=NULL, groups="a"), er)
  expect_error(row.bartlett(x=NULL, groups="a"), er)
  expect_error(row.cor.pearson(x=NULL, y=0), er)
  expect_error(row.ievora(x=NULL, groups="a"), er)
})

test_that("y cannot be NULL", {
  matX <- matrix(c(TRUE, FALSE), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(row.t.equalvar(x=0, y=NULL), er)
  expect_error(row.t.welch(x=0, y=NULL), er)
  expect_error(row.t.paired(x=0, y=NULL), er)
  expect_error(row.cor.pearson(x=0, y=NULL), er)
})


test_that("x cannot be in a list", {
  er <- '"x" must be a numeric matrix or vector'
  expect_error(row.t.onesample(x=list(1:5)), er)
  expect_error(row.t.equalvar(x=list(1:5), y=0), er)
  expect_error(row.t.welch(x=list(1:5), y=0), er)
  expect_error(row.t.paired(x=list(1:5), y=0), er)
  expect_error(row.oneway.equalvar(x=list(1:5), groups="a"), er)
  expect_error(row.oneway.welch(x=list(1:5), groups="a"), er)
  expect_error(row.kruskalwallis(x=list(1:5), groups="a"), er)
  expect_error(row.bartlett(x=list(1:5), groups="a"), er)
  expect_error(row.cor.pearson(x=list(1:5), y=0), er)
  expect_error(row.ievora(x=list(1:5), groups="a"), er)
})

test_that("y cannot be in a list", {
  matX <- matrix(c(TRUE, FALSE), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(row.t.equalvar(x=0, y=list(1:5)), er)
  expect_error(row.t.welch(x=0, y=list(1:5)), er)
  expect_error(row.t.paired(x=0, y=list(1:5)), er)
  expect_error(row.cor.pearson(x=0, y=list(1:5)), er)
})

test_that("x cannot be a list", {
  er <- '"x" must be a numeric matrix or vector'
  expect_error(row.t.onesample(x=as.list(1:5)), er)
  expect_error(row.t.equalvar(x=as.list(1:5), y=0), er)
  expect_error(row.t.welch(x=as.list(1:5), y=0), er)
  expect_error(row.t.paired(x=as.list(1:5), y=0), er)
  expect_error(row.oneway.equalvar(x=as.list(1:5), groups="a"), er)
  expect_error(row.oneway.welch(x=as.list(1:5), groups="a"), er)
  expect_error(row.kruskalwallis(x=as.list(1:5), groups="a"), er)
  expect_error(row.bartlett(x=as.list(1:5), groups="a"), er)
  expect_error(row.cor.pearson(x=as.list(1:5), y=0), er)
  expect_error(row.ievora(x=as.list(1:5), groups="a"), er)
})

test_that("y cannot be a list", {
  matX <- matrix(c(TRUE, FALSE), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(row.t.equalvar(x=0, y=as.list(1:5)), er)
  expect_error(row.t.welch(x=0, y=as.list(1:5)), er)
  expect_error(row.t.paired(x=0, y=as.list(1:5)), er)
  expect_error(row.cor.pearson(x=0, y=as.list(1:5)), er)
})

################################################################################
########################## WRONG GROUP SPECIFICATION ###########################
################################################################################

test_that("groups cannot be NULL", {
  matX <- matrix(1:12, ncol=3)
  er <- '"groups" must be a vector with length ncol\\(x\\)'
  expect_error(row.oneway.equalvar(x=matX, groups=NULL), er)
  expect_error(row.oneway.welch(x=matX, groups=NULL), er)
  expect_error(row.kruskalwallis(x=matX, groups=NULL), er)
  expect_error(row.bartlett(x=matX, groups=NULL), er)
  expect_error(row.ievora(x=matX, groups=NULL), er)
})

test_that("groups cannot be a list", {
  matX <- matrix(1:12, ncol=3)
  er <- '"groups" must be a vector with length ncol\\(x\\)'
  expect_error(row.oneway.equalvar(x=matX, groups=list(1:3)), er)
  expect_error(row.oneway.welch(x=matX, groups=list(1:3)), er)
  expect_error(row.kruskalwallis(x=matX, groups=list(1:3)), er)
  expect_error(row.bartlett(x=matX, groups=list(1:3)), er)
  expect_error(row.ievora(x=matX, groups=list(1:3)), er)
})

test_that("groups cannot be a matrix", {
  matX <- matrix(1:12, ncol=4)
  grp  <- cbind(c("A", "A"), c("B", "B"))
  er <- '"groups" must be a vector with length ncol\\(x\\)'
  expect_error(row.oneway.equalvar(x=matX, groups=grp), er)
  expect_error(row.oneway.welch(x=matX, groups=grp), er)
  expect_error(row.kruskalwallis(x=matX, groups=grp), er)
  expect_error(row.bartlett(x=matX, groups=grp), er)
  expect_error(row.ievora(x=matX, groups=grp), er)
})

################################################################################
############################## DIMENSION MISMATCH ##############################
################################################################################

test_that("x and y has same number of rows", {
  matX <- matrix(1:9, nrow=3)
  matY <- matrix(1:16, nrow=4)
  er <- '"x" and "y" must have the same number of rows'
  expect_error(row.t.equalvar(x=matX, y=matY), er)
  expect_error(row.t.welch(x=matX, y=matY), er)
  expect_error(row.t.paired(x=matX, y=matY), er)
  expect_error(row.cor.pearson(x=matX, y=matY), er)
})

test_that("x and y has same number of columns", {
  matX <- matrix(1:9, nrow=3)
  matY <- matrix(1:12, nrow=3)
  er <- '"x" and "y" must have the same number of columns'
  expect_error(row.t.paired(x=matX, y=matY), er)
  expect_error(row.cor.pearson(x=matX, y=matY), er)
})

test_that("group length matches number of columns", {
  matX <- matrix(1:12, nrow=3)
  er <- '"groups" must be a vector with length ncol\\(x\\)'
  expect_error(row.oneway.equalvar(x=matX, groups=1:3), er)
  expect_error(row.oneway.welch(x=matX, groups=1:3), er)
  expect_error(row.kruskalwallis(x=matX, groups=1:3), er)
  expect_error(row.bartlett(x=matX, groups=1:3), er)
  expect_error(row.ievora(x=matX, groups=1:3), er)
})

################################################################################
########################### SPECIAL GROUP EXCEPTIONS ###########################
################################################################################

test_that("groups have required number of groups", {
  matX <- matrix(1:12, nrow=3)
  er <- '"groups" must have no more than 2 unique elements'
  expect_error(row.ievora(x=matX, groups=c("A","B","C","C")), er)
  expect_error(row.ievora(x=c(1,2,3,4,NA), groups=c(0,0,1,1,2)), er)
})

