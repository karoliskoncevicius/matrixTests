library(matrixTests)
source("utils/capture.r")

#--- less than 2 paired observations -------------------------------------------

wrn <- 'row_t_paired: 1 of the rows had less than 2 paired observations.\nFirst occurrence at row 1'
nacolumns <- c("stderr", "df", "statistic", "pvalue", "conf.low", "conf.high")

# 0 observations
res <- capture(row_t_paired(NA_integer_, NA_integer_))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.pair, 0))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.y, 0))

# 0 observations in x and 2 in y
res <- capture(row_t_paired(c(NA_integer_, NA_integer_), 1:2))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.pair, 0))
stopifnot(all.equal(res$value$obs.x, 0))
stopifnot(all.equal(res$value$obs.y, 2))

# 1 observation
res <- capture(row_t_paired(1, 2))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.pair, 1))
stopifnot(all.equal(res$value$obs.x, 1))
stopifnot(all.equal(res$value$obs.y, 1))

# 2 observations in x and 1 in y
res <- capture(row_t_paired(1:2, c(1, NA_integer_)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.pair, 1))
stopifnot(all.equal(res$value$obs.x, 2))
stopifnot(all.equal(res$value$obs.y, 1))

# 1 paired observation after removing NAs
res <- capture(row_t_paired(c(1,2,3,NA,NA), c(NA,NA,2,1,4)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.pair, 1))
stopifnot(all.equal(res$value$obs.x, 3))
stopifnot(all.equal(res$value$obs.y, 3))


#--- constant paired values ----------------------------------------------------

wrn <- 'row_t_paired: 1 of the rows had essentially constant values.\nFirst occurrence at row 1'
nacolumns <- c("stderr", "df", "statistic", "pvalue", "conf.low", "conf.high")

# all values are the same
res <- capture(row_t_paired(c(1,1,1,1), c(1,1,1,1)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$var.x, 0))
stopifnot(all.equal(res$value$var.y, 0))
stopifnot(all.equal(res$value$var.diff, 0))

# all values are equal within each group
res <- capture(row_t_paired(c(1,1,1,1), c(2,2,2,2)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$var.x, 0))
stopifnot(all.equal(res$value$var.y, 0))
stopifnot(all.equal(res$value$var.diff, 0))

# all values are equal within each group after NAs
res <- capture(row_t_paired(c(NA,1,1,1,1,3), c(4,2,2,2,2,NA)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$var.diff, 0))

# values are equal after taking differene
res <- capture(row_t_paired(c(3,2,0), c(2,1,-1)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$var.diff, 0))

# values are equal after taking differene with NAs
res <- capture(row_t_paired(c(100,3,2,0,NA), c(NA,2,1,-1,-50)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$var.diff, 0))

