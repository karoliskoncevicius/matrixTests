# Matrix Tests #

A package dedicated to running statistical hypothesis tests on rows of matrices.

## Example ##

```r

X <- matrix(rnorm(10000000), ncol=10)
Y <- matrix(rnorm(10000000), ncol=10)

# The usual way
res1 <- vector(nrow(X), mode="list")
for(i in 1:nrow(X)) {
  res[[i]] <- t.test(X[i,], Y[i,])
}

# matrixTest way
res2 <- ttest_welch(X, Y)

```

## Goals ##

1. Fast execution via vectorization.
2. Output that is detailed and easy to use.
3. Result compatibility with tests that are implemented in R.

## Available Tests ##

* `ttest_single`   - one sample t.test.
* `ttest_welch`    - t.test with Welch adjustment.
* `ttest_equalvar` - t.test for two groups with equal variance.
* `ttest_paired`   - paired t.test.

## Planned ##

* one-way ANOVA
* test for Pearson's correlation coefficient
* linear regression


