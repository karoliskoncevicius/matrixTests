# Matrix Tests #

A package dedicated to running statistical hypothesis tests on rows of matrices.

## Example ##

```r

X <- matrix(rnorm(10000000), ncol=10)
Y <- matrix(rnorm(10000000), ncol=10)

# The usual way
res1 <- vector(nrow(X), mode="list")
for(i in 1:nrow(X)) {
  res1[[i]] <- t.test(X[i,], Y[i,])
}

# matrixTest way
res2 <- ttest_welch(X, Y)

```

## Goals ##

1. Fast execution via vectorization.
2. Output that is detailed and easy to use.
3. Result compatibility with tests that are implemented in R.

## Available Tests ##

* `ttest_onegroup`  - t.test for a single group.
* `ttest_welch`     - t.test with Welch adjustment.
* `ttest_equalvar`  - t.test for two groups with equal variance.
* `ttest_paired`    - t.test for paired observations.
* `oneway_equalvar` - oneway ANOVA for groups with equal variance.
* `kruskalwallis`   - Kruskal-Wallis test.
* `bartlett`        - Bartlett's test.

## Test-Based Algorithms ##

* `ievora` - detect differential variability.

## Planned ##

* oneway test with Welch correction
* test for Pearson's correlation coefficient
* linear regression

## Installation ##

Using the `devtools` library:

```r
library(devtools)
install_github("KKPMW/matrixTests")
```

## See Also ##

1. *Computing thousands of test statistics simultaneously in R*,
Holger Schwender, Tina Müller. Statistical Computing & Graphics.
Volume 18, No 1, June 2007.
2. `lmFit` in the **limma** package.
3. `rowttests()` in the **genefilter** package.
4. `mt.teststat()` in the **multtest** package.
5. `row.T.test()` in the **HybridMTest** package.
6. `rowTtest()` in the **viper** package.


