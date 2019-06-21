context("row and col versions give equal outputs")

################################################################################
##################### STANDARD MATRICES GIVE EQUAL RESULTS #####################
################################################################################

test_that("row eq col on matrix xs and ys", {
  X <- matrix(rnorm(17*6), ncol=6)
  Y <- matrix(rnorm(17*6), ncol=6)
  g <- c("a", "a", "a", "b", "b", "b")
  expect_equal(row_t_onesample(X), col_t_onesample(t(X)))
  expect_equal(row_t_equalvar(X, Y), col_t_equalvar(t(X), t(Y)))
  expect_equal(row_t_welch(X, Y), col_t_welch(t(X), t(Y)))
  expect_equal(row_t_paired(X, Y), col_t_paired(t(X), t(Y)))
  expect_equal(row_oneway_equalvar(X, g), col_oneway_equalvar(t(X), g))
  expect_equal(row_oneway_welch(X, g), col_oneway_welch(t(X), g))
  expect_equal(row_kruskalwallis(X, g), col_kruskalwallis(t(X), g))
  expect_equal(row_bartlett(X, g), col_bartlett(t(X), g))
  expect_equal(row_flignerkilleen(X, g), col_flignerkilleen(t(X), g))
  expect_equal(row_cor_pearson(X, Y), col_cor_pearson(t(X), t(Y)))
  expect_equal(row_jarquebera(X), col_jarquebera(t(X)))
  expect_equal(row_wilcoxon_onesample(X), col_wilcoxon_onesample(t(X)))
  expect_equal(row_wilcoxon_twosample(X, Y), col_wilcoxon_twosample(t(X), t(Y)))
  expect_equal(row_wilcoxon_paired(X, Y), col_wilcoxon_paired(t(X), t(Y)))
  expect_equal(row_ievora(X, g), col_ievora(t(X), g))
})

################################################################################
############################# NUMERIC DATA FRAMES ##############################
################################################################################

test_that("row eq col on data.frame xs and ys", {
  X <- as.data.frame(matrix(rnorm(17*6), ncol=6))
  Y <- as.data.frame(matrix(rnorm(17*6), ncol=6))
  g <- c("a", "a", "a", "b", "b", "b")
  expect_equal(row_t_onesample(X), col_t_onesample(t(X)))
  expect_equal(row_t_equalvar(X, Y), col_t_equalvar(t(X), t(Y)))
  expect_equal(row_t_welch(X, Y), col_t_welch(t(X), t(Y)))
  expect_equal(row_t_paired(X, Y), col_t_paired(t(X), t(Y)))
  expect_equal(row_oneway_equalvar(X, g), col_oneway_equalvar(t(X), g))
  expect_equal(row_oneway_welch(X, g), col_oneway_welch(t(X), g))
  expect_equal(row_kruskalwallis(X, g), col_kruskalwallis(t(X), g))
  expect_equal(row_bartlett(X, g), col_bartlett(t(X), g))
  expect_equal(row_flignerkilleen(X, g), col_flignerkilleen(t(X), g))
  expect_equal(row_cor_pearson(X, Y), col_cor_pearson(t(X), t(Y)))
  expect_equal(row_jarquebera(X), col_jarquebera(t(X)))
  expect_equal(row_wilcoxon_onesample(X), col_wilcoxon_onesample(t(X)))
  expect_equal(row_wilcoxon_twosample(X, Y), col_wilcoxon_twosample(t(X), t(Y)))
  expect_equal(row_wilcoxon_paired(X, Y), col_wilcoxon_paired(t(X), t(Y)))
  expect_equal(row_ievora(X, g), col_ievora(t(X), g))
})

################################################################################
############################### NUMERIC VECTORS ################################
################################################################################

test_that("row eq col on vector xs and ys", {
  X <- rnorm(6)
  Y <- rnorm(6)
  g <- c("a", "a", "a", "b", "b", "b")
  expect_equal(row_t_onesample(X), col_t_onesample(X))
  expect_equal(row_t_equalvar(X, Y), col_t_equalvar(X, Y))
  expect_equal(row_t_welch(X, Y), col_t_welch(X, Y))
  expect_equal(row_t_paired(X, Y), col_t_paired(X, Y))
  expect_equal(row_oneway_equalvar(X, g), col_oneway_equalvar(X, g))
  expect_equal(row_oneway_welch(X, g), col_oneway_welch(X, g))
  expect_equal(row_kruskalwallis(X, g), col_kruskalwallis(X, g))
  expect_equal(row_bartlett(X, g), col_bartlett(X, g))
  expect_equal(row_flignerkilleen(X, g), col_flignerkilleen(X, g))
  expect_equal(row_cor_pearson(X, Y), col_cor_pearson(X, Y))
  expect_equal(row_jarquebera(X), col_jarquebera(X))
  expect_equal(row_wilcoxon_onesample(X), col_wilcoxon_onesample(X))
  expect_equal(row_wilcoxon_twosample(X, Y), col_wilcoxon_twosample(X, Y))
  expect_equal(row_wilcoxon_paired(X, Y), col_wilcoxon_paired(X, Y))
  expect_equal(row_ievora(X, g), col_ievora(X, g))
})

################################################################################
################################ EMPTY MATRICES ################################
################################################################################

test_that("row eq col on 0 dimension matrix xs and ys", {
  X <- matrix(NA_integer_, nrow=0, ncol=0)
  Y <- matrix(NA_integer_, nrow=0, ncol=0)
  g <- character()
  expect_equal(row_t_onesample(X), col_t_onesample(X))
  expect_equal(row_t_equalvar(X, Y), col_t_equalvar(X, Y))
  expect_equal(row_t_welch(X, Y), col_t_welch(X, Y))
  expect_equal(row_t_paired(X, Y), col_t_paired(X, Y))
  expect_equal(row_oneway_equalvar(X, g), col_oneway_equalvar(X, g))
  expect_equal(row_oneway_welch(X, g), col_oneway_welch(X, g))
  expect_equal(row_kruskalwallis(X, g), col_kruskalwallis(X, g))
  expect_equal(row_bartlett(X, g), col_bartlett(X, g))
  expect_equal(row_flignerkilleen(X, g), col_flignerkilleen(X, g))
  expect_equal(row_cor_pearson(X, Y), col_cor_pearson(X, Y))
  expect_equal(row_jarquebera(X), col_jarquebera(X))
  expect_equal(row_wilcoxon_onesample(X), col_wilcoxon_onesample(X))
  expect_equal(row_wilcoxon_twosample(X, Y), col_wilcoxon_twosample(X, Y))
  expect_equal(row_wilcoxon_paired(X, Y), col_wilcoxon_paired(X, Y))
  expect_equal(row_ievora(X, g), col_ievora(X, g))
})

################################################################################
############################## X MATRIX Y VECTOR ###############################
################################################################################

test_that("row eq col when x is a matrix and y is a vector", {
  X <- matrix(rnorm(17*6), ncol=6)
  Y <- rnorm(6)
  g <- c("a", "a", "a", "b", "b", "b")
  expect_equal(row_t_equalvar(X, Y), col_t_equalvar(t(X), Y))
  expect_equal(row_t_welch(X, Y), col_t_welch(t(X), Y))
  expect_equal(row_t_paired(X, Y), col_t_paired(t(X), Y))
  expect_equal(row_cor_pearson(X, Y), col_cor_pearson(t(X), Y))
  expect_equal(row_wilcoxon_twosample(X, Y), col_wilcoxon_twosample(t(X), Y))
  expect_equal(row_wilcoxon_paired(X, Y), col_wilcoxon_paired(t(X), Y))
})

