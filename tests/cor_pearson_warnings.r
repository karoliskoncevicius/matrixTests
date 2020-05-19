library(matrixTests)
source("utils/capture.r")

#--- less than 3 complete observations -----------------------------------------

wrn <- 'row_cor_pearson: 1 of the rows had less than 3 complete observations.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue", "conf.low", "conf.high")

# 0 observations
res <- capture(row_cor_pearson(numeric(), numeric()))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))

# 1 observation
res <- capture(row_cor_pearson(1, 1))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))

# 2 observations
res <- capture(row_cor_pearson(1:2, 3:4))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))

# 2 matched observations after NAs
res <- capture(row_cor_pearson(c(NA,1:3), c(3,NA,5,4)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))


#--- exactly 3 complete observations -------------------------------------------

wrn <- 'row_cor_pearson: 1 of the rows had exactly 3 complete observations: no confidence intervals produced.\nFirst occurrence at row 1'
nacolumns <- c("conf.low", "conf.high")

# exactly 3 observations
res <- capture(row_cor_pearson(c(1,1,2), c(1,2,3)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))

# exactly 3 observations after removing NA values
res <- capture(row_cor_pearson(c(NA,1:4), c(3,NA,1,3,2)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))


#--- x has zero variance -------------------------------------------------------

wrn <- 'row_cor_pearson: 1 of the rows had zero standard deviation in x.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue", "conf.low", "conf.high")

# all values are the same
res <- capture(row_cor_pearson(c(1,1,1,1), 1:4))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))

# all values are the same after removing NAs
res <- capture(row_cor_pearson(c(1,1,1,1,2), c(1:4,NA)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))


#--- y has zero variance -------------------------------------------------------

wrn <- 'row_cor_pearson: 1 of the rows had zero standard deviation in y.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue", "conf.low", "conf.high")

# all values are the same
res <- capture(row_cor_pearson(1:4, c(1,1,1,1)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))

# all values are the same after removing NAs
res <- capture(row_cor_pearson(c(1:4,NA), c(1,1,1,1,2)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))


#--- both have 0 variance ------------------------------------------------------

wrnx <- 'row_cor_pearson: 1 of the rows had zero standard deviation in x.\nFirst occurrence at row 1'
wrny <- 'row_cor_pearson: 1 of the rows had zero standard deviation in y.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue", "conf.low", "conf.high")

# all values are constant in both groups
res <- capture(row_cor_pearson(c(1,1,1,1), c(1,1,1,1)))
stopifnot(all.equal(res$warning, c(wrnx, wrny)))
stopifnot(all(is.na(res$value[,nacolumns])))

# all values are constant in both groups after removing NAs
res <- capture(row_cor_pearson(c(NA,1,1,1,1,1,2), c(1,2,2,2,2,2,NA)))
stopifnot(all.equal(res$warning, c(wrnx, wrny)))
stopifnot(all(is.na(res$value[,nacolumns])))


#--- correlation is perfect ----------------------------------------------------

wrn <- 'row_cor_pearson: 1 of the rows had essentially perfect fit: results might be unreliable for small sample sizes.\nFirst occurrence at row 1'

# perfect positive correlation
res <- capture(row_cor_pearson(c(1:4), c(1:4)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$statistic, Inf))

# perfect negative correlation
res <- capture(row_cor_pearson(c(1:4), -c(1:4)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$statistic, -Inf))

# perfect correlation with NAs
res <- capture(row_cor_pearson(-c(1:4,NA), -c(1:4,2)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$statistic, Inf))

