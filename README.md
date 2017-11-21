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

|             Description             | matrixTests       | R equivalent
------------------------------------------------------------------------------------------------
| t.test for a single group.          | `ttest_onegroup`  | t.test(x)
| t.test with Welch adjustment.       | `ttest_welch`     | t.test(x,y)
| t.test with equal variance.         | `ttest_equalvar`  | t.test(x,y,var.equal=TRUE)
| t.test for paired observations.     | `ttest_paired`    | t.test(x,y,paired=TRUE)
| Pearson correlation test.           | `cor_pearson`     | cor.test(x,y)
| oneway ANOVA with equal variance.   | `oneway_equalvar` | oneway.test(x,g,var.equal=TRUE)
| oneway ANOVA with Welch adjustment. | `oneway_welch`    | oneway.test(x,g)
| Kruskal-Wallis test.                | `kruskalwallis`   | kruskal.test(x,g)
| Bartlett's test.                    | `bartlett`        | bartlett.test(x,g)

## Test-Based Algorithms ##

* `ievora` - detect differential variability.

## Planned ##

* test for Spearman and Kendall correlations
* Fisher's exact test

## Installation ##

Using the `devtools` library:

```r
library(devtools)
install_github("KKPMW/matrixTests")
```

## Dependencies ##

1. `matrixStats` package.

## See Also ##

1. *Computing thousands of test statistics simultaneously in R*,
Holger Schwender, Tina Müller. Statistical Computing & Graphics.
Volume 18, No 1, June 2007.
2. `lmFit` in the **limma** package.
3. `rowttests()` in the **genefilter** package.
4. `mt.teststat()` in the **multtest** package.
5. `row.T.test()` in the **HybridMTest** package.
6. `rowTtest()` in the **viper** package.
7. `ttests()` in the **Rfast** package.


