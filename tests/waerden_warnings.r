library(matrixTests)
source("utils/capture.r")

#--- removing NA groups --------------------------------------------------------

wrn <- '2 columns dropped due to missing group information'

# 2 NAs
x <- 1:10
g <- c(1,1,1,1,NA,NA,2,2,2,2)
res <- capture(row_waerden(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs.tot, 8))
stopifnot(all.equal(res$value$obs.groups, 2))


#--- less than 2 observations --------------------------------------------------

wrn <- 'row_waerden: 1 of the rows had less than 2 total observations.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# one value one group
res <- capture(row_waerden(1, "a"))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.groups, 1))
stopifnot(all.equal(res$value$obs.tot, 1))

# one value one group with NAs
res <- capture(row_waerden(c(1,NA,NA,NA), c(1,1,2,2)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.groups, 1))
stopifnot(all.equal(res$value$obs.tot, 1))


#--- less than 2 groups --------------------------------------------------------

wrn <- 'row_waerden: 1 of the rows had less than 2 groups with enough observations.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# single group with 10 observations
res <- capture(row_waerden(1:10, rep(1,10)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.groups, 1))
stopifnot(all.equal(res$value$obs.tot, 10))

# 3 groups but others have only NA values
x <- c(1:8, NA, NA)
g <- c(rep("a", 8), "b", "c")
res <- capture(row_waerden(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.groups, 1))
stopifnot(all.equal(res$value$obs.tot, 8))


#--- single unique value -------------------------------------------------------

wrn <- 'row_waerden: 1 of the rows had essentially constant values.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# two groups
res <- capture(row_waerden(rep(1, 3), c("a","a","b")))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.groups, 2))
stopifnot(all.equal(res$value$obs.tot, 3))

# two groups with NA values
x <- c(rep(1,8),NA,NA)
g <- c(rep("a", 4), rep("b", 4), "c", "d")
res <- capture(row_waerden(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.groups, 2))
stopifnot(all.equal(res$value$obs.tot, 8))

