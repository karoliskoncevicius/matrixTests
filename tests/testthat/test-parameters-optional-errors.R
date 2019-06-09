context("Errors in optional parameters")

################################################################################
######################## VALUES MUST HAVE CORRECT TYPE #########################
################################################################################

test_that("alternative must be a character", {
  er <- '"alternative" must be a character vector with length 1 or nrow\\(x\\)'
  mat <- matrix(1:12, nrow=3)
  # NULL
  expect_error(row_t_onesample(x=mat, alternative=NULL), er)
  expect_error(row_t_equalvar(x=mat, y=mat, alternative=NULL), er)
  expect_error(row_t_welch(x=mat, y=mat, alternative=NULL), er)
  expect_error(row_t_paired(x=mat, y=mat, alternative=NULL), er)
  expect_error(row_cor_pearson(x=mat, y=mat, alternative=NULL), er)
  expect_error(row_wilcoxon_onesample(x=mat, alternative=NULL), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, alternative=NULL), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, alternative=NULL), er)
  # NA
  expect_error(row_t_onesample(x=mat, alternative=NA), er)
  expect_error(row_t_equalvar(x=mat, y=mat, alternative=NA), er)
  expect_error(row_t_welch(x=mat, y=mat, alternative=NA), er)
  expect_error(row_cor_pearson(x=mat, y=mat, alternative=NA), er)
  expect_error(row_wilcoxon_onesample(x=mat, alternative=NA), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, alternative=NA), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, alternative=NA), er)
  # numeric
  expect_error(row_t_onesample(x=mat, alternative=1), er)
  expect_error(row_t_equalvar(x=mat, y=mat, alternative=1), er)
  expect_error(row_t_welch(x=mat, y=mat, alternative=1), er)
  expect_error(row_t_paired(x=mat, y=mat, alternative=1), er)
  expect_error(row_cor_pearson(x=mat, y=mat, alternative=1), er)
  expect_error(row_wilcoxon_onesample(x=mat, alternative=1), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, alternative=1), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, alternative=1), er)
  # complex
  expect_error(row_t_onesample(x=mat, alternative=complex(1)), er)
  expect_error(row_t_equalvar(x=mat, y=mat, alternative=complex(1)), er)
  expect_error(row_t_welch(x=mat, y=mat, alternative=complex(1)), er)
  expect_error(row_t_paired(x=mat, y=mat, alternative=complex(1)), er)
  expect_error(row_cor_pearson(x=mat, y=mat, alternative=complex(1)), er)
  expect_error(row_wilcoxon_onesample(x=mat, alternative=complex(1)), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, alternative=complex(1)), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, alternative=complex(1)), er)
  # in list
  expect_error(row_t_onesample(x=mat, alternative=list("less")), er)
  expect_error(row_t_equalvar(x=mat, y=mat, alternative=list("less")), er)
  expect_error(row_t_welch(x=mat, y=mat, alternative=list("less")), er)
  expect_error(row_t_paired(x=mat, y=mat, alternative=list("less")), er)
  expect_error(row_cor_pearson(x=mat, y=mat, alternative=list("less")), er)
  expect_error(row_wilcoxon_onesample(x=mat, alternative=list("less")), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, alternative=list("less")), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, alternative=list("less")), er)
  # data frame
  expect_error(row_t_onesample(x=mat, alternative=data.frame("less")), er)
  expect_error(row_t_equalvar(x=mat, y=mat, alternative=data.frame("less")), er)
  expect_error(row_t_welch(x=mat, y=mat, alternative=data.frame("less")), er)
  expect_error(row_t_paired(x=mat, y=mat, alternative=data.frame("less")), er)
  expect_error(row_cor_pearson(x=mat, y=mat, alternative=data.frame("less")), er)
  expect_error(row_wilcoxon_onesample(x=mat, alternative=data.frame("less")), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, alternative=data.frame("less")), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, alternative=data.frame("less")), er)
})

test_that("mu must be numeric", {
  er <- '"mu" must be a numeric vector with length 1 or nrow\\(x\\)'
  mat <- matrix(1:12, nrow=3)
  # NULL
  expect_error(row_t_onesample(x=mat, mu=NULL), er)
  expect_error(row_t_equalvar(x=mat, y=mat, mu=NULL), er)
  expect_error(row_t_welch(x=mat, y=mat, mu=NULL), er)
  expect_error(row_t_paired(x=mat, y=mat, mu=NULL), er)
  expect_error(row_wilcoxon_onesample(x=mat, mu=NULL), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, mu=NULL), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, mu=NULL), er)
  # NA
  expect_error(row_t_onesample(x=mat, mu=NA), er)
  expect_error(row_t_equalvar(x=mat, y=mat, mu=NA), er)
  expect_error(row_t_welch(x=mat, y=mat, mu=NA), er)
  expect_error(row_t_paired(x=mat, y=mat, mu=NA), er)
  expect_error(row_wilcoxon_onesample(x=mat, mu=NA), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, mu=NA), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, mu=NA), er)
  # character
  expect_error(row_t_onesample(x=mat, mu="1"), er)
  expect_error(row_t_equalvar(x=mat, y=mat, mu="1"), er)
  expect_error(row_t_welch(x=mat, y=mat, mu="1"), er)
  expect_error(row_t_paired(x=mat, y=mat, mu="1"), er)
  expect_error(row_wilcoxon_onesample(x=mat, mu="1"), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, mu="1"), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, mu="1"), er)
  # complex
  expect_error(row_t_onesample(x=mat, mu=complex(1)), er)
  expect_error(row_t_equalvar(x=mat, y=mat, mu=complex(1)), er)
  expect_error(row_t_welch(x=mat, y=mat, mu=complex(1)), er)
  expect_error(row_t_paired(x=mat, y=mat, mu=complex(1)), er)
  expect_error(row_wilcoxon_onesample(x=mat, mu=complex(1)), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, mu=complex(1)), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, mu=complex(1)), er)
  # in list
  expect_error(row_t_onesample(x=mat, mu=list(1)), er)
  expect_error(row_t_equalvar(x=mat, y=mat, mu=list(1)), er)
  expect_error(row_t_welch(x=mat, y=mat, mu=list(1)), er)
  expect_error(row_t_paired(x=mat, y=mat, mu=list(1)), er)
  expect_error(row_wilcoxon_onesample(x=mat, mu=list(1)), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, mu=list(1)), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, mu=list(1)), er)
  # data frame
  expect_error(row_t_onesample(x=mat, mu=data.frame(1)), er)
  expect_error(row_t_equalvar(x=mat, y=mat, mu=data.frame(1)), er)
  expect_error(row_t_welch(x=mat, y=mat, mu=data.frame(1)), er)
  expect_error(row_t_paired(x=mat, y=mat, mu=data.frame(1)), er)
  expect_error(row_wilcoxon_onesample(x=mat, mu=data.frame(1)), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, mu=data.frame(1)), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, mu=data.frame(1)), er)
})


test_that("conf.level must be numeric", {
  er <- '"conf.level" must be a numeric vector with length 1 or nrow\\(x\\)'
  mat <- matrix(1:12, nrow=3)
  # NULL
  expect_error(row_t_onesample(x=mat, conf.level=NULL), er)
  expect_error(row_t_equalvar(x=mat, y=mat, conf.level=NULL), er)
  expect_error(row_t_welch(x=mat, y=mat, conf.level=NULL), er)
  expect_error(row_t_paired(x=mat, y=mat, conf.level=NULL), er)
  expect_error(row_cor_pearson(x=mat, y=mat, conf.level=NULL), er)
  # NA
  expect_error(row_t_onesample(x=mat, conf.level=NA), er)
  expect_error(row_t_equalvar(x=mat, y=mat, conf.level=NA), er)
  expect_error(row_t_welch(x=mat, y=mat, conf.level=NA), er)
  expect_error(row_t_paired(x=mat, y=mat, conf.level=NA), er)
  expect_error(row_cor_pearson(x=mat, y=mat, conf.level=NA), er)
  # character
  expect_error(row_t_onesample(x=mat, conf.level="1"), er)
  expect_error(row_t_equalvar(x=mat, y=mat, conf.level="1"), er)
  expect_error(row_t_welch(x=mat, y=mat, conf.level="1"), er)
  expect_error(row_t_paired(x=mat, y=mat, conf.level="1"), er)
  expect_error(row_cor_pearson(x=mat, y=mat, conf.level="1"), er)
  # complex
  expect_error(row_t_onesample(x=mat, conf.level=complex(1)), er)
  expect_error(row_t_equalvar(x=mat, y=mat, conf.level=complex(1)), er)
  expect_error(row_t_welch(x=mat, y=mat, conf.level=complex(1)), er)
  expect_error(row_t_paired(x=mat, y=mat, conf.level=complex(1)), er)
  expect_error(row_cor_pearson(x=mat, y=mat, conf.level=complex(1)), er)
  # in list
  expect_error(row_t_onesample(x=mat, conf.level=list(1)), er)
  expect_error(row_t_equalvar(x=mat, y=mat, conf.level=list(1)), er)
  expect_error(row_t_welch(x=mat, y=mat, conf.level=list(1)), er)
  expect_error(row_t_paired(x=mat, y=mat, conf.level=list(1)), er)
  expect_error(row_cor_pearson(x=mat, y=mat, conf.level=list(1)), er)
  # data frame
  expect_error(row_t_onesample(x=mat, conf.level=data.frame(1)), er)
  expect_error(row_t_equalvar(x=mat, y=mat, conf.level=data.frame(1)), er)
  expect_error(row_t_welch(x=mat, y=mat, conf.level=data.frame(1)), er)
  expect_error(row_t_paired(x=mat, y=mat, conf.level=data.frame(1)), er)
  expect_error(row_cor_pearson(x=mat, y=mat, conf.level=data.frame(1)), er)
})

test_that("exact indicator must be logical", {
  er1 <- '"exact" must be a logical vector with length 1 or nrow\\(x\\)'
  mat <- matrix(1:12, nrow=3)
  # NULL
  expect_error(row_wilcoxon_onesample(x=mat, exact=NULL), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, exact=NULL), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, exact=NULL), er1)
  # NA (non logical)
  expect_error(row_wilcoxon_onesample(x=mat, exact=NA_character_), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, exact=NA_character_), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, exact=NA_character_), er1)
  # character
  expect_error(row_wilcoxon_onesample(x=mat, exact="TRUE"), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, exact="TRUE"), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, exact="TRUE"), er1)
  # numeric
  expect_error(row_wilcoxon_onesample(x=mat, exact=1), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, exact=1), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, exact=1), er1)
  # complex
  expect_error(row_wilcoxon_onesample(x=mat, exact=complex(1)), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, exact=complex(1)), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, exact=complex(1)), er1)
  # in list
  expect_error(row_wilcoxon_onesample(x=mat, exact=list(TRUE)), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, exact=list(TRUE)), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, exact=list(TRUE)), er1)
  # data frame
  expect_error(row_wilcoxon_onesample(x=mat, exact=data.frame(TRUE)), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, exact=data.frame(TRUE)), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, exact=data.frame(TRUE)), er1)
})


test_that("correct indicator must be logical or NA", {
  er1 <- '"correct" must be a logical vector with length 1 or nrow\\(x\\)'
  mat <- matrix(1:12, nrow=3)
  # NULL
  expect_error(row_wilcoxon_onesample(x=mat, correct=NULL), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, correct=NULL), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, correct=NULL), er1)
  # NA
  expect_error(row_wilcoxon_onesample(x=mat, correct=NA_character_), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, correct=NA_character_), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, correct=NA_character_), er1)
  # character
  expect_error(row_wilcoxon_onesample(x=mat, correct="TRUE"), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, correct="TRUE"), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, correct="TRUE"), er1)
  # numeric (not 1 or 0)
  expect_error(row_wilcoxon_onesample(x=mat, correct=1), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, correct=1), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, correct=1), er1)
  # complex (not 1 or 0)
  expect_error(row_wilcoxon_onesample(x=mat, correct=complex(1)), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, correct=complex(1)), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, correct=complex(1)), er1)
  # in list
  expect_error(row_wilcoxon_onesample(x=mat, correct=list(TRUE)), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, correct=list(TRUE)), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, correct=list(TRUE)), er1)
  # data frame
  expect_error(row_wilcoxon_onesample(x=mat, correct=data.frame(TRUE)), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, correct=data.frame(TRUE)), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, correct=data.frame(TRUE)), er1)
})

test_that("p-value cutoffs must be numeric", {
  er1 <- '"cutT" must be a numeric vector with length 1'
  er2 <- '"cutBfdr" must be a numeric vector with length 1'
  mat <- matrix(1:12, nrow=3)
  grp <- c(1,1,0,0)
  # NULL
  expect_error(row_ievora(x=mat, b=grp, cutT=NULL), er1)
  expect_error(row_ievora(x=mat, b=grp, cutBfdr=NULL), er2)
  # NA
  expect_error(row_ievora(x=mat, b=grp, cutT=NA), er1)
  expect_error(row_ievora(x=mat, b=grp, cutBfdr=NA), er2)
  # character
  expect_error(row_ievora(x=mat, b=grp, cutT="0.05"), er1)
  expect_error(row_ievora(x=mat, b=grp, cutBfdr="0.05"), er2)
  # complex
  expect_error(row_ievora(x=mat, b=grp, cutT=complex(1)), er1)
  expect_error(row_ievora(x=mat, b=grp, cutBfdr=complex(1)), er2)
  # in list
  expect_error(row_ievora(x=mat, b=grp, cutT=list(1)), er1)
  expect_error(row_ievora(x=mat, b=grp, cutBfdr=list(1)), er2)
  # data frame
  expect_error(row_ievora(x=mat, b=grp, cutT=data.frame(1)), er1)
  expect_error(row_ievora(x=mat, b=grp, cutBfdr=data.frame(1)), er2)
})

################################################################################
##################### VALUES MUST HAVE CORRECT DIMENSIONS ######################
################################################################################

test_that("alternative has correct dimensions", {
  er <- '"alternative" must be a character vector with length 1 or nrow\\(x\\)'
  mat <- matrix(1:12, ncol=3)
  # too short
  expect_error(row_t_onesample(x=mat, alternative=c("less", "greater")), er)
  expect_error(row_t_equalvar(x=mat, y=mat, alternative=c("less", "greater")), er)
  expect_error(row_t_welch(x=mat, y=mat, alternative=c("less", "greater")), er)
  expect_error(row_t_paired(x=mat, y=mat, alternative=c("less", "greater")), er)
  expect_error(row_cor_pearson(x=mat, y=mat, alternative=c("less", "greater")), er)
  expect_error(row_wilcoxon_onesample(x=mat, alternative=c("less", "greater")), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, alternative=c("less", "greater")), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, alternative=c("less", "greater")), er)
  # too long
  expect_error(row_t_onesample(x=mat, alternative=rep("less", 5)), er)
  expect_error(row_t_equalvar(x=mat, y=mat, alternative=rep("less", 5)), er)
  expect_error(row_t_welch(x=mat, y=mat, alternative=rep("less", 5)), er)
  expect_error(row_t_paired(x=mat, y=mat, alternative=rep("less", 5)), er)
  expect_error(row_cor_pearson(x=mat, y=mat, alternative=rep("less", 5)), er)
  expect_error(row_wilcoxon_onesample(x=mat, alternative=rep("less", 5)), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, alternative=rep("less", 5)), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, alternative=rep("less", 5)), er)
  # matrix format
  alt <- matrix(rep("less", 4), ncol=2)
  expect_error(row_t_onesample(x=mat, alternative=alt), er)
  expect_error(row_t_equalvar(x=mat, y=mat, alternative=alt), er)
  expect_error(row_t_welch(x=mat, y=mat, alternative=alt), er)
  expect_error(row_t_paired(x=mat, y=mat, alternative=alt), er)
  expect_error(row_cor_pearson(x=mat, y=mat, alternative=alt), er)
  expect_error(row_wilcoxon_onesample(x=mat, alternative=alt), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, alternative=alt), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, alternative=alt), er)
})

test_that("mu has correct dimensions", {
  er <- '"mu" must be a numeric vector with length 1 or nrow\\(x\\)'
  mat <- matrix(1:12, ncol=3)
  # too short
  expect_error(row_t_onesample(x=mat, mu=c(1,2)), er)
  expect_error(row_t_equalvar(x=mat, y=mat, mu=c(1,2)), er)
  expect_error(row_t_welch(x=mat, y=mat, mu=c(1,2)), er)
  expect_error(row_t_paired(x=mat, y=mat, mu=c(1,2)), er)
  expect_error(row_wilcoxon_onesample(x=mat, mu=c(1,2)), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, mu=c(1,2)), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, mu=c(1,2)), er)
  # too long
  expect_error(row_t_onesample(x=mat, mu=rep(0, 5)), er)
  expect_error(row_t_equalvar(x=mat, y=mat, mu=rep(0, 5)), er)
  expect_error(row_t_welch(x=mat, y=mat, mu=rep(0, 5)), er)
  expect_error(row_t_paired(x=mat, y=mat, mu=rep(0, 5)), er)
  expect_error(row_wilcoxon_onesample(x=mat, mu=rep(0, 5)), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, mu=rep(0, 5)), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, mu=rep(0, 5)), er)
  # matrix format
  mus <- matrix(rep(1, 4), ncol=2)
  expect_error(row_t_onesample(x=mat, mu=mus), er)
  expect_error(row_t_equalvar(x=mat, y=mat, mu=mus), er)
  expect_error(row_t_welch(x=mat, y=mat, mu=mus), er)
  expect_error(row_t_paired(x=mat, y=mat, mu=mus), er)
  expect_error(row_wilcoxon_onesample(x=mat, mu=mus), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, mu=mus), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, mu=mus), er)
})

test_that("conf.level has correct dimensions", {
  er <- '"conf.level" must be a numeric vector with length 1 or nrow\\(x\\)'
  mat <- matrix(1:12, ncol=3)
  # too short
  expect_error(row_t_onesample(x=mat, conf.level=c(1,2)), er)
  expect_error(row_t_equalvar(x=mat, y=mat, conf.level=c(1,2)), er)
  expect_error(row_t_welch(x=mat, y=mat, conf.level=c(1,2)), er)
  expect_error(row_t_paired(x=mat, y=mat, conf.level=c(1,2)), er)
  expect_error(row_cor_pearson(x=mat, y=mat, conf.level=c(1,2)), er)
  # too long
  expect_error(row_t_onesample(x=mat, conf.level=rep(0, 5)), er)
  expect_error(row_t_equalvar(x=mat, y=mat, conf.level=rep(0, 5)), er)
  expect_error(row_t_welch(x=mat, y=mat, conf.level=rep(0, 5)), er)
  expect_error(row_t_paired(x=mat, y=mat, conf.level=rep(0, 5)), er)
  expect_error(row_cor_pearson(x=mat, y=mat, conf.level=rep(0, 5)), er)
  # matrix format
  cfs <- matrix(rep(1, 4), ncol=2)
  expect_error(row_t_onesample(x=mat, conf.level=cfs), er)
  expect_error(row_t_equalvar(x=mat, y=mat, conf.level=cfs), er)
  expect_error(row_t_welch(x=mat, y=mat, conf.level=cfs), er)
  expect_error(row_t_paired(x=mat, y=mat, conf.level=cfs), er)
  expect_error(row_cor_pearson(x=mat, y=mat, conf.level=cfs), er)
})

test_that("exact indicator has correct dimensions", {
  er1 <- '"exact" must be a logical vector with length 1 or nrow\\(x\\)'
  mat <- matrix(1:12, ncol=3)
  # too short
  expect_error(row_wilcoxon_onesample(x=mat, exact=c(TRUE, FALSE)), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, exact=c(TRUE, FALSE)), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, exact=c(TRUE, FALSE)), er1)
  # too long
  expect_error(row_wilcoxon_onesample(x=mat, exact=rep(TRUE,6)), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, exact=rep(TRUE,6)), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, exact=rep(TRUE,6)), er1)
  # matrix format
  exs <- matrix(rep(TRUE, 4), ncol=2)
  expect_error(row_wilcoxon_onesample(x=mat, exact=exs), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, exact=exs), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, exact=exs), er1)
})

test_that("correct indicator has correct dimensions", {
  er1 <- '"correct" must be a logical vector with length 1 or nrow\\(x\\)'
  mat <- matrix(1:12, ncol=3)
  # too short
  expect_error(row_wilcoxon_onesample(x=mat, correct=c(TRUE, FALSE)), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, correct=c(TRUE, FALSE)), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, correct=c(TRUE, FALSE)), er1)
  # too long
  expect_error(row_wilcoxon_onesample(x=mat, correct=rep(TRUE,6)), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, correct=rep(TRUE,6)), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, correct=rep(TRUE,6)), er1)
  # matrix format
  exs <- matrix(rep(TRUE, 4), ncol=2)
  expect_error(row_wilcoxon_onesample(x=mat, correct=exs), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, correct=exs), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, correct=exs), er1)
})

test_that("p-value cutoffs have the right dimensions", {
  er1 <- '"cutT" must be a numeric vector with length 1'
  er2 <- '"cutBfdr" must be a numeric vector with length 1'
  mat <- matrix(1:12, ncol=3)
  grp <- c(0,0,1)
  # too long
  expect_error(row_ievora(x=mat, b=grp, cutT=c(1,2)), er1)
  expect_error(row_ievora(x=mat, b=grp, cutBfdr=c(1,2)), er2)
  # matrix format
  cts <- matrix(rep(1, 4), ncol=2)
  expect_error(row_ievora(x=mat, b=grp, cutT=cts), er1)
  expect_error(row_ievora(x=mat, b=grp, cutBfdr=cts), er2)
})

################################################################################
######################## VALUES MUST BE IN CORRECT SET #########################
################################################################################

test_that("alternative is in: less, greater, two-sided)", {
  er <- 'all "alternative" values must be in: two\\.sided, less, greater'
  mat <- matrix(1:12, nrow=3)
  # one value
  expect_error(row_t_onesample(x=mat, alternative="ga"), er)
  expect_error(row_t_equalvar(x=mat, y=mat, alternative="ga"), er)
  expect_error(row_t_welch(x=mat, y=mat, alternative="ga"), er)
  expect_error(row_t_paired(x=mat, y=mat, alternative="ga"), er)
  expect_error(row_cor_pearson(x=mat, y=mat, alternative="ga"), er)
  expect_error(row_wilcoxon_onesample(x=mat, alternative="ga"), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, alternative="ga"), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, alternative="ga"), er)
  # for each row and one incorrect
  expect_error(row_t_onesample(x=mat, alternative=c("t","l","c")), er)
  expect_error(row_t_equalvar(x=mat, y=mat, alternative=c("t","l","c")), er)
  expect_error(row_t_welch(x=mat, y=mat, alternative=c("t","l","c")), er)
  expect_error(row_t_paired(x=mat, y=mat, alternative=c("t","l","c")), er)
  expect_error(row_cor_pearson(x=mat, y=mat, alternative=c("t","l","c")), er)
  expect_error(row_wilcoxon_onesample(x=mat, alternative=c("t","l","c")), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, alternative=c("t","l","c")), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, alternative=c("t","l","c")), er)
})

test_that("conf.level is in: 0-1)", {
  er <- 'all "conf.level" values must be between: 0 and 1'
  mat <- matrix(1:12, nrow=3)
  # slightly below
  expect_error(row_t_onesample(x=mat, conf.level=-0.001), er)
  expect_error(row_t_equalvar(x=mat, y=mat, conf.level=-0.001), er)
  expect_error(row_t_welch(x=mat, y=mat, conf.level=-0.001), er)
  expect_error(row_t_paired(x=mat, y=mat, conf.level=-0.001), er)
  expect_error(row_cor_pearson(x=mat, y=mat, conf.level=-0.001), er)
  # slightly above
  expect_error(row_t_onesample(x=mat, conf.level=1.001), er)
  expect_error(row_t_equalvar(x=mat, y=mat, conf.level=1.001), er)
  expect_error(row_t_welch(x=mat, y=mat, conf.level=1.001), er)
  expect_error(row_t_paired(x=mat, y=mat, conf.level=1.001), er)
  expect_error(row_cor_pearson(x=mat, y=mat, conf.level=1.001), er)
  # special values
  expect_error(row_t_onesample(x=mat, conf.level=NA_integer_), er)
  expect_error(row_t_equalvar(x=mat, y=mat, conf.level=NaN), er)
  expect_error(row_t_welch(x=mat, y=mat, conf.level=Inf), er)
  expect_error(row_t_paired(x=mat, y=mat, conf.level=-Inf), er)
  expect_error(row_cor_pearson(x=mat, y=mat, conf.level=-Inf), er)
})

test_that("when mu is in: [-Inf,Inf])", {
  er <- 'all "mu" values must be between: -Inf and Inf'
  mat <- matrix(1:12, nrow=3)
  # NA
  expect_error(row_t_onesample(x=mat, mu=NA_integer_), er)
  expect_error(row_t_equalvar(x=mat, y=mat, mu=NA_integer_), er)
  expect_error(row_t_welch(x=mat, y=mat, mu=NA_integer_), er)
  expect_error(row_t_paired(x=mat, y=mat, mu=NA_integer_), er)
  # NaN
  expect_error(row_t_onesample(x=mat, mu=NaN), er)
  expect_error(row_t_equalvar(x=mat, y=mat, mu=NaN), er)
  expect_error(row_t_welch(x=mat, y=mat, mu=NaN), er)
  expect_error(row_t_paired(x=mat, y=mat, mu=NaN), er)
})

test_that("when mu is in: (-Inf,Inf))", {
  er <- 'all "mu" values must be greater than -Inf and lower than Inf'
  mat <- matrix(1:12, nrow=3)
  # NA
  expect_error(row_wilcoxon_onesample(x=mat, mu=NA_integer_), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, mu=NA_integer_), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, mu=NA_integer_), er)
  # NaN
  expect_error(row_wilcoxon_onesample(x=mat, mu=NaN), er)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, mu=NaN), er)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, mu=NaN), er)
})

test_that("correct indicator cannot be NA", {
  er1 <- 'all "correct" values must be in: TRUE, FALSE'
  mat <- matrix(1:12, nrow=3)
  expect_error(row_wilcoxon_onesample(x=mat, correct=NA), er1)
  expect_error(row_wilcoxon_twosample(x=mat, y=mat, correct=NA), er1)
  expect_error(row_wilcoxon_paired(x=mat, y=mat, correct=NA), er1)
})

test_that("p-value cut-offs must be in 0-1", {
  er1 <- 'all "cutT" values must be between: 0 and 1'
  er2 <- 'all "cutBfdr" values must be between: 0 and 1'
  mat <- matrix(1:12, nrow=3)
  grp <- c(0,0,1,1)
  # slightly below
  expect_error(row_ievora(x=mat, b=grp, cutT=-0.001), er1)
  expect_error(row_ievora(x=mat, b=grp, cutBfdr=-0.001), er2)
  # slightly above
  expect_error(row_ievora(x=mat, b=grp, cutT=1.001), er1)
  expect_error(row_ievora(x=mat, b=grp, cutBfdr=1.001), er2)
  # special values
  expect_error(row_ievora(x=mat, b=grp, cutT=NA_integer_), er1)
  expect_error(row_ievora(x=mat, b=grp, cutBfdr=NaN), er2)
  expect_error(row_ievora(x=mat, b=grp, cutT=Inf), er1)
  expect_error(row_ievora(x=mat, b=grp, cutBfdr=-Inf), er2)
})

