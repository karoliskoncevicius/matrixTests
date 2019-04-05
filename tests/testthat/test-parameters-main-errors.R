context("Errors in main parameters")

################################################################################
################################### MISSING ####################################
################################################################################

test_that("x cannot be missing", {
  er <- 'argument "x" is missing, with no default'
  expect_error(row_t_equalvar(), er)
  expect_error(row_t_welch(), er)
  expect_error(row_t_onesample(), er)
  expect_error(row_t_paired(), er)
  expect_error(row_oneway_equalvar(), er)
  expect_error(row_oneway_welch(), er)
  expect_error(row_kruskalwallis(), er)
  expect_error(row_bartlett(), er)
  expect_error(row_cor_pearson(), er)
  expect_error(row_ievora(), er)
  expect_error(row_jarquebera(), er)
  expect_error(row_flignerkilleen(), er)
})

test_that("y cannot be missing", {
  er <- 'argument "y" is missing, with no default'
  expect_error(row_t_equalvar(x=NA), er)
  expect_error(row_t_welch(x=NA), er)
  expect_error(row_t_paired(x=NA), er)
  expect_error(row_cor_pearson(x=NA), er)
})

test_that("groups cannot be missing", {
  er <- 'argument "g" is missing, with no default'
  expect_error(row_oneway_equalvar(x=NA), er)
  expect_error(row_oneway_welch(x=NA), er)
  expect_error(row_kruskalwallis(x=NA), er)
  expect_error(row_bartlett(x=NA), er)
  expect_error(row_flignerkilleen(x=NA), er)
})

test_that("binary cannot be missing", {
  er <- 'argument "b" is missing, with no default'
  expect_error(row_ievora(x=NA), er)
})

################################################################################
################################# NON NUMERIC ##################################
################################################################################

test_that("x cannot be a character", {
  matX <- matrix(c("1", "2"), nrow=1)
  er <- '"x" must be a numeric matrix or vector'
  expect_error(row_t_onesample(x=matX), er)
  expect_error(row_t_equalvar(x=matX, y=0), er)
  expect_error(row_t_welch(x=matX, y=0), er)
  expect_error(row_t_paired(x=matX, y=0), er)
  expect_error(row_oneway_equalvar(x=matX, g="a"), er)
  expect_error(row_oneway_welch(x=matX, g="a"), er)
  expect_error(row_kruskalwallis(x=matX, g="a"), er)
  expect_error(row_bartlett(x=matX, g="a"), er)
  expect_error(row_cor_pearson(x=matX, y=0), er)
  expect_error(row_ievora(x=matX, b="a"), er)
  expect_error(row_jarquebera(x=matX), er)
  expect_error(row_flignerkilleen(x=matX, g="a"), er)
})

test_that("y cannot be a character", {
  matX <- matrix(c("1", "2"), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(row_t_equalvar(x=0, y=matX), er)
  expect_error(row_t_welch(x=0, y=matX), er)
  expect_error(row_t_paired(x=0, y=matX), er)
  expect_error(row_cor_pearson(x=0, y=matX), er)
})


test_that("x cannot be partially numeric", {
  er <- '"x" must be a numeric matrix or vector'
  expect_error(row_t_onesample(x=iris), er)
  expect_error(row_t_equalvar(x=iris, y=0), er)
  expect_error(row_t_welch(x=iris, y=0), er)
  expect_error(row_t_paired(x=iris, y=0), er)
  expect_error(row_oneway_equalvar(x=iris, g="a"), er)
  expect_error(row_oneway_welch(x=iris, g="a"), er)
  expect_error(row_kruskalwallis(x=iris, g="a"), er)
  expect_error(row_bartlett(x=iris, g="a"), er)
  expect_error(row_cor_pearson(x=iris, y=0), er)
  expect_error(row_ievora(x=iris, b="a"), er)
  expect_error(row_jarquebera(x=iris), er)
  expect_error(row_flignerkilleen(x=iris, g="a"), er)
})

test_that("y cannot be partially numeric", {
  er <- '"y" must be a numeric matrix or vector'
  expect_error(row_t_equalvar(x=0, y=iris), er)
  expect_error(row_t_welch(x=0, y=iris), er)
  expect_error(row_t_paired(x=0, y=iris), er)
  expect_error(row_cor_pearson(x=0, y=iris), er)
})


test_that("x cannot be complex", {
  matX <- matrix(complex(c(1,2), c(3,4)), nrow=1)
  er <- '"x" must be a numeric matrix or vector'
  expect_error(row_t_onesample(x=matX), er)
  expect_error(row_t_equalvar(x=matX, y=0), er)
  expect_error(row_t_welch(x=matX, y=0), er)
  expect_error(row_t_paired(x=matX, y=0), er)
  expect_error(row_oneway_equalvar(x=matX, g="a"), er)
  expect_error(row_oneway_welch(x=matX, g="a"), er)
  expect_error(row_kruskalwallis(x=matX, g="a"), er)
  expect_error(row_bartlett(x=matX, g="a"), er)
  expect_error(row_cor_pearson(x=matX, y=0), er)
  expect_error(row_ievora(x=matX, b="a"), er)
  expect_error(row_jarquebera(x=matX), er)
  expect_error(row_flignerkilleen(x=matX, g="a"), er)
})

test_that("y cannot be complex", {
  matX <- matrix(complex(c(1,2), c(3,4)), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(row_t_equalvar(x=0, y=matX), er)
  expect_error(row_t_welch(x=0, y=matX), er)
  expect_error(row_t_paired(x=0, y=matX), er)
  expect_error(row_cor_pearson(x=0, y=matX), er)
})


test_that("x cannot be logical", {
  matX <- matrix(c(TRUE, FALSE), nrow=1)
  er <- '"x" must be a numeric matrix or vector'
  expect_error(row_t_onesample(x=matX), er)
  expect_error(row_t_equalvar(x=matX, y=0), er)
  expect_error(row_t_welch(x=matX, y=0), er)
  expect_error(row_t_paired(x=matX, y=0), er)
  expect_error(row_oneway_equalvar(x=matX, g="a"), er)
  expect_error(row_oneway_welch(x=matX, g="a"), er)
  expect_error(row_kruskalwallis(x=matX, g="a"), er)
  expect_error(row_bartlett(x=matX, g="a"), er)
  expect_error(row_cor_pearson(x=matX, y=0), er)
  expect_error(row_ievora(x=matX, b="a"), er)
  expect_error(row_jarquebera(x=matX), er)
  expect_error(row_flignerkilleen(x=matX, g="a"), er)
})

test_that("y cannot be logical", {
  matX <- matrix(c(TRUE, FALSE), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(row_t_equalvar(x=0, y=matX), er)
  expect_error(row_t_welch(x=0, y=matX), er)
  expect_error(row_t_paired(x=0, y=matX), er)
  expect_error(row_cor_pearson(x=0, y=matX), er)
})


test_that("x cannot be NULL", {
  er <- '"x" must be a numeric matrix or vector'
  expect_error(row_t_onesample(x=NULL), er)
  expect_error(row_t_equalvar(x=NULL, y=0), er)
  expect_error(row_t_welch(x=NULL, y=0), er)
  expect_error(row_t_paired(x=NULL, y=0), er)
  expect_error(row_oneway_equalvar(x=NULL, g="a"), er)
  expect_error(row_oneway_welch(x=NULL, g="a"), er)
  expect_error(row_kruskalwallis(x=NULL, g="a"), er)
  expect_error(row_bartlett(x=NULL, g="a"), er)
  expect_error(row_cor_pearson(x=NULL, y=0), er)
  expect_error(row_ievora(x=NULL, b="a"), er)
  expect_error(row_jarquebera(x=NULL), er)
  expect_error(row_flignerkilleen(x=NULL, g="a"), er)
})

test_that("y cannot be NULL", {
  matX <- matrix(c(TRUE, FALSE), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(row_t_equalvar(x=0, y=NULL), er)
  expect_error(row_t_welch(x=0, y=NULL), er)
  expect_error(row_t_paired(x=0, y=NULL), er)
  expect_error(row_cor_pearson(x=0, y=NULL), er)
})


test_that("x cannot be in a list", {
  er <- '"x" must be a numeric matrix or vector'
  expect_error(row_t_onesample(x=list(1:5)), er)
  expect_error(row_t_equalvar(x=list(1:5), y=0), er)
  expect_error(row_t_welch(x=list(1:5), y=0), er)
  expect_error(row_t_paired(x=list(1:5), y=0), er)
  expect_error(row_oneway_equalvar(x=list(1:5), g="a"), er)
  expect_error(row_oneway_welch(x=list(1:5), g="a"), er)
  expect_error(row_kruskalwallis(x=list(1:5), g="a"), er)
  expect_error(row_bartlett(x=list(1:5), g="a"), er)
  expect_error(row_cor_pearson(x=list(1:5), y=0), er)
  expect_error(row_ievora(x=list(1:5), b="a"), er)
  expect_error(row_jarquebera(x=list(1:5)), er)
  expect_error(row_flignerkilleen(x=list(1:5), g="a"), er)
})

test_that("y cannot be in a list", {
  matX <- matrix(c(TRUE, FALSE), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(row_t_equalvar(x=0, y=list(1:5)), er)
  expect_error(row_t_welch(x=0, y=list(1:5)), er)
  expect_error(row_t_paired(x=0, y=list(1:5)), er)
  expect_error(row_cor_pearson(x=0, y=list(1:5)), er)
})

test_that("x cannot be a list", {
  er <- '"x" must be a numeric matrix or vector'
  expect_error(row_t_onesample(x=as.list(1:5)), er)
  expect_error(row_t_equalvar(x=as.list(1:5), y=0), er)
  expect_error(row_t_welch(x=as.list(1:5), y=0), er)
  expect_error(row_t_paired(x=as.list(1:5), y=0), er)
  expect_error(row_oneway_equalvar(x=as.list(1:5), g="a"), er)
  expect_error(row_oneway_welch(x=as.list(1:5), g="a"), er)
  expect_error(row_kruskalwallis(x=as.list(1:5), g="a"), er)
  expect_error(row_bartlett(x=as.list(1:5), g="a"), er)
  expect_error(row_cor_pearson(x=as.list(1:5), y=0), er)
  expect_error(row_ievora(x=as.list(1:5), b="a"), er)
  expect_error(row_jarquebera(x=as.list(1:5)), er)
  expect_error(row_flignerkilleen(x=as.list(1:5), g="a"), er)
})

test_that("y cannot be a list", {
  matX <- matrix(c(TRUE, FALSE), nrow=1)
  er <- '"y" must be a numeric matrix or vector'
  expect_error(row_t_equalvar(x=0, y=as.list(1:5)), er)
  expect_error(row_t_welch(x=0, y=as.list(1:5)), er)
  expect_error(row_t_paired(x=0, y=as.list(1:5)), er)
  expect_error(row_cor_pearson(x=0, y=as.list(1:5)), er)
})

################################################################################
########################## WRONG GROUP SPECIFICATION ###########################
################################################################################

test_that("groups cannot be NULL", {
  matX <- matrix(1:12, ncol=3)
  er <- '"g" must be a vector with length ncol\\(x\\)'
  expect_error(row_oneway_equalvar(x=matX, g=NULL), er)
  expect_error(row_oneway_welch(x=matX, g=NULL), er)
  expect_error(row_kruskalwallis(x=matX, g=NULL), er)
  expect_error(row_bartlett(x=matX, g=NULL), er)
  expect_error(row_flignerkilleen(x=matX, g=NULL), er)
})

test_that("groups cannot be a list", {
  matX <- matrix(1:12, ncol=3)
  er <- '"g" must be a vector with length ncol\\(x\\)'
  expect_error(row_oneway_equalvar(x=matX, g=list(1:3)), er)
  expect_error(row_oneway_welch(x=matX, g=list(1:3)), er)
  expect_error(row_kruskalwallis(x=matX, g=list(1:3)), er)
  expect_error(row_bartlett(x=matX, g=list(1:3)), er)
  expect_error(row_flignerkilleen(x=matX, g=list(1:3)), er)
})

test_that("groups cannot be a matrix", {
  matX <- matrix(1:12, ncol=4)
  grp  <- cbind(c("A", "A"), c("B", "B"))
  er <- '"g" must be a vector with length ncol\\(x\\)'
  expect_error(row_oneway_equalvar(x=matX, g=grp), er)
  expect_error(row_oneway_welch(x=matX, g=grp), er)
  expect_error(row_kruskalwallis(x=matX, g=grp), er)
  expect_error(row_bartlett(x=matX, g=grp), er)
  expect_error(row_flignerkilleen(x=matX, g=grp), er)
})

################################################################################
########################## WRONG BINARY SPECIFICATION ##########################
################################################################################

test_that("binary cannot be NULL", {
  matX <- matrix(1:12, ncol=3)
  er <- '"b" must be a vector with length ncol\\(x\\)'
  expect_error(row_ievora(x=matX, b=NULL), er)
})

test_that("binary cannot be a list", {
  matX <- matrix(1:12, ncol=3)
  er <- '"b" must be a vector with length ncol\\(x\\)'
  expect_error(row_ievora(x=matX, b=list(1:3)), er)
})

test_that("binary cannot be a matrix", {
  matX <- matrix(1:12, ncol=4)
  grp  <- cbind(c("A", "A"), c("B", "B"))
  er <- '"b" must be a vector with length ncol\\(x\\)'
  expect_error(row_ievora(x=matX, b=grp), er)
})

test_that("binary has required number of groups", {
  matX <- matrix(1:12, nrow=3)
  er <- '"b" must have no more than 2 unique elements'
  expect_error(row_ievora(x=matX, b=c("A","B","C","C")), er)
  expect_error(row_ievora(x=c(1,2,3,4,NA), b=c(0,0,1,1,2)), er)
})

################################################################################
############################## DIMENSION MISMATCH ##############################
################################################################################

test_that("x and y has same number of rows", {
  matX <- matrix(1:9, nrow=3)
  matY <- matrix(1:16, nrow=4)
  er <- '"x" and "y" must have the same number of rows'
  expect_error(row_t_equalvar(x=matX, y=matY), er)
  expect_error(row_t_welch(x=matX, y=matY), er)
  expect_error(row_t_paired(x=matX, y=matY), er)
  expect_error(row_cor_pearson(x=matX, y=matY), er)
})

test_that("x and y has same number of columns", {
  matX <- matrix(1:9, nrow=3)
  matY <- matrix(1:12, nrow=3)
  er <- '"x" and "y" must have the same number of columns'
  expect_error(row_t_paired(x=matX, y=matY), er)
  expect_error(row_cor_pearson(x=matX, y=matY), er)
})

test_that("group length matches number of columns", {
  matX <- matrix(1:12, nrow=3)
  er <- '"g" must be a vector with length ncol\\(x\\)'
  expect_error(row_oneway_equalvar(x=matX, g=1:3), er)
  expect_error(row_oneway_welch(x=matX, g=1:3), er)
  expect_error(row_kruskalwallis(x=matX, g=1:3), er)
  expect_error(row_bartlett(x=matX, g=1:3), er)
  expect_error(row_flignerkilleen(x=matX, g=1:3), er)
})

test_that("binary length matches number of columns", {
  matX <- matrix(1:12, nrow=3)
  er <- '"b" must be a vector with length ncol\\(x\\)'
  expect_error(row_ievora(x=matX, b=1:3), er)
})

