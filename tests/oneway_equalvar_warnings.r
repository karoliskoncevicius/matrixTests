library(matrixTests)
source("utils/capture.r")

#--- removing na groups --------------------------------------------------------

wrn <- '2 columns dropped due to missing group information'

# 2 NAs
x <- 1:10
g <- c(1,1,1,1,NA,NA,2,2,2,2)
res <- capture(row_oneway_equalvar(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs.tot, 8))
stopifnot(all.equal(res$value$obs.groups, 2))


#--- one observation per group -------------------------------------------------

wrn <- 'row_oneway_equalvar: 1 of the rows had one observation per group.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# two groups with a single obervation each
x <- 1:2
g <- 1:2
res <- capture(row_oneway_equalvar(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 2))
stopifnot(all.equal(res$value$obs.groups, 2))

# 10 groups 10 observations
x <- 1:10
g <- 1:10
res <- capture(row_oneway_equalvar(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 10))
stopifnot(all.equal(res$value$obs.groups, 10))

# 2 groups 3 observations but one is NA
x <- c(rnorm(2), NA)
g <- c("a", "b", "b")
res <- capture(row_oneway_equalvar(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 2))
stopifnot(all.equal(res$value$obs.groups, 2))


#--- less than 2 groups --------------------------------------------------------

wrn <- 'row_oneway_equalvar: 1 of the rows had less than 2 groups with enough observations.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# single group with 10 observations
res <- capture(row_oneway_equalvar(1:10, rep(1,10)))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.groups, 1))
stopifnot(all.equal(res$value$obs.tot, 10))

# 3 groups but others have only NA values
x <- c(rnorm(8), NA, NA)
g <- c(rep("a", 8), "b", "c")
res <- capture(row_oneway_equalvar(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.groups, 1))
stopifnot(all.equal(res$value$obs.tot, 8))


#--- single unique value -------------------------------------------------------

wrn <- 'row_oneway_equalvar: 1 of the rows had essentially constant values.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# two groups
res <- capture(row_oneway_equalvar(rep(1, 3), c("a","a","b")))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.groups, 2))
stopifnot(all.equal(res$value$obs.tot, 3))

# two groups with NA values
x <- c(rep(1,8),NA,NA)
g <- c(rep("a", 4), rep("b", 4), "c", "d")
res <- capture(row_oneway_equalvar(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.groups, 2))
stopifnot(all.equal(res$value$obs.tot, 8))


#--- single unique values within each group ------------------------------------

wrn <- 'row_oneway_equalvar: 1 of the rows had zero within group variance: result might be unreliable.\nFirst occurrence at row 1'

# two groups - constant values within group
x <- c(1,1,0,0)
g <- c(1,1,2,2)
res <- capture(row_oneway_equalvar(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs.tot, 4))
stopifnot(all.equal(res$value$obs.groups, 2))

# three groups - constant values within each group + NAs
x <- c(3,3,NA,0,0,NA,-1,-1,NA)
g <- c(1,1,1,2,2,2,3,3,3)
res <- capture(row_oneway_equalvar(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs.tot, 6))
stopifnot(all.equal(res$value$obs.groups, 3))

