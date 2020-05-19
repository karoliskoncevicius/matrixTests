library(matrixTests)
source("utils/capture.r")

#--- removing na groups --------------------------------------------------------

wrn <- '2 columns dropped due to missing group information'

# 2 NAs
x <- 1:10
g <- c(1,1,1,1,NA,NA,2,2,2,2)
res <- capture(row_flignerkilleen(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all.equal(res$value$obs.tot, 8))
stopifnot(all.equal(res$value$obs.groups, 2))


#--- less than 2 groups --------------------------------------------------------

wrn <- 'row_flignerkilleen: 1 of the rows had less than 2 groups with enough observations.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# all values in one group
x <- rnorm(10)
g <- rep("a", 10)
res <- capture(row_flignerkilleen(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 10))
stopifnot(all.equal(res$value$obs.groups, 1))

# two groups but one has only NAs
x <- c(rnorm(3), NA, NA, NA)
g <- rep(c("a","b"), each=3)
res <- capture(row_flignerkilleen(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 3))
stopifnot(all.equal(res$value$obs.groups, 1))


#--- one observation per group -------------------------------------------------

wrn <- 'row_flignerkilleen: 1 of the rows had one observation per group.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# two groups with a single obervation each
x <- 1:2
g <- 1:2
res <- capture(row_flignerkilleen(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 2))
stopifnot(all.equal(res$value$obs.groups, 2))

# 10 groups 10 observations
x <- 1:10
g <- 1:10
res <- capture(row_flignerkilleen(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 10))
stopifnot(all.equal(res$value$obs.groups, 10))

# 2 groups 3 observations but one is NA
x <- c(rnorm(2), NA)
g <- c("a", "b", "b")
res <- capture(row_flignerkilleen(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 2))
stopifnot(all.equal(res$value$obs.groups, 2))


#--- all groups have zero variance ---------------------------------------------

wrn <- 'row_flignerkilleen: 1 of the rows had zero variance in all of the groups.\nFirst occurrence at row 1'
nacolumns <- c("statistic", "pvalue")

# two groups - all values are constant
x <- c(1,1,1,1)
g <- c("a","a","b","b")
res <- capture(row_flignerkilleen(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 4))
stopifnot(all.equal(res$value$obs.groups, 2))

# two groups - constant values within each group
x <- c(3,3,0,0)
g <- c("a","a","b","b")
res <- capture(row_flignerkilleen(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 4))
stopifnot(all.equal(res$value$obs.groups, 2))

# two groups - constant values plus NA within each group
x <- c(3,3,NA,0,0,NA)
g <- c("a","a","a","b","b","b")
res <- capture(row_flignerkilleen(x, g))
stopifnot(all.equal(res$warning, wrn))
stopifnot(all(is.na(res$value[,nacolumns])))
stopifnot(all.equal(res$value$obs.tot, 4))
stopifnot(all.equal(res$value$obs.groups, 2))

