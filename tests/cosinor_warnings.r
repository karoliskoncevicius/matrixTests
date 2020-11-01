library(matrixTests)
source("utils/capture.r")

#--- removing na timepoints ----------------------------------------------------

wrn <- '2 columns dropped due to missing time information'

# 2 NAs
x <- 1:10
t <- c(1,2,3,4,NA,NA,7,8,9,10)
res <- capture(row_cosinor(x, t))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs, 8))


#--- less than 3 observations --------------------------------------------------

wrn <- 'row_cosinor: 1 of the rows had less than 3 complete observations: no p-values produced, amplitude and acrophase will be unreliable.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# 1 observations
res <- capture(row_cosinor(1, 1))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 1))

# 2 observations
res <- capture(row_cosinor(1:2, 1:2))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 2))

# 2 observations with NAs
res <- capture(row_cosinor(c(1:2,NA,NA), 1:4))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 2))


#--- exactly 3 observations ----------------------------------------------------

wrn <- 'row_cosinor: 1 of the rows had exactly 3 complete observations: no p-values produced.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# 3 observations
res <- capture(row_cosinor(1:3, 1:3))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 3))

# with NAs present
res <- capture(row_cosinor(c(1:3,NA), 1:4))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 3))


#--- less than 3 unique timepoints ---------------------------------------------

wrn <- 'row_cosinor: 1 of the rows had less than 3 unique timepoints within the specified period: amplitude and acrophase will be unreliable.\nFirst occurrence at row 1'

# TODO: check what to do about turning statistics to NA.
#       right now only the case with 1 time point turns them to NA.

# one distinct point, duplicated multiple times
res <- capture(row_cosinor(1:4, rep(1,4)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs, 4))

# two distinct points, duplicated multiple times
res <- capture(row_cosinor(1:4, c(1,2,1,2)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs, 4))

# with NAs present
res <- capture(res <- row_cosinor(c(1,2,3,4,NA,NA), c(1,2,1,2,3,4)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs, 4))

# four points, but period is such that they are not unique
res <- capture(row_cosinor(1:4, c(1,2,3,4), 2))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs, 4))


#--- constant values -----------------------------------------------------------

wrn <- 'row_cosinor: 1 of the rows had essentially constant values.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# all the values are 0
res <- capture(row_cosinor(c(0,0,0,0), c(1,2,3,4)))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs, 4))

# all the values are constant
res <- capture(row_cosinor(c(1,1,1,1), c(1,2,3,4)))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs, 4))

# all the vallues are constant with NAs
res <- capture(row_cosinor(c(1,1,1,1,4), c(2,2,2,2,NA)))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$warning[2], wrn))
stopifnot(all.equal(res$value$obs, 4))


#--- perfect fit ---------------------------------------------------------------

wrn <- 'row_cosinor: 1 of the rows had essentially perfect fit.\nFirst occurrence at row 1'

# perfect sine
res <- capture(row_cosinor(sin(2*pi*1:24/24), 1:24))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$statistic, Inf))

# only three distinct points present
res <- capture(row_cosinor(c(1:3,3), c(1:3,3), 24))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$statistic, Inf))

