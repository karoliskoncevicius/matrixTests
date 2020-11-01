library(matrixTests)
source("utils/capture.r")

#--- x has less than 2 observations --------------------------------------------

wrn <- 'row_t_onesample: 1 of the rows had less than 2 "x" observations.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue", "conf.low", "conf.high")

# 0 observations
res <- capture(row_t_onesample(numeric()))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 0))

# 1 observation
res <- capture(row_t_onesample(1))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 1))

# 1 observations with NAs
res <- capture(row_t_onesample(c(NA,1,NA)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 1))


#--- all values are constant ---------------------------------------------------

wrn <- 'row_t_onesample: 1 of the rows had essentially constant values.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue", "conf.low", "conf.high")

# 2 equal observations
res <- capture(row_t_onesample(c(1,1)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 2))
stopifnot(all.equal(res$value$var, 0))

# 3 equal observations and some NAs
res <- capture(row_t_onesample(c(NA,2,2,2,NA)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs, 3))
stopifnot(all.equal(res$value$var, 0))

