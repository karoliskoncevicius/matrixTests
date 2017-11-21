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

|             Name                   |      matrixTests        |       R equivalent
|------------------------------------|-------------------------|-------------------------------------
| Single group t.test                | `ttest_onegroup(x)`     | `t.test(x)`
| Welch t.test                       | `ttest_welch(x, y)`     | `t.test(x, y)`
| Equal variance t.test              | `ttest_equalvar(x, y)`  | `t.test(x, y, var.equal=TRUE)`
| Paired t.test                      | `ttest_paired(x, y)`    | `t.test(x, y, paired=TRUE)`
| Pearson's correlation test         | `cor_pearson(x, y)`     | `cor.test(x, y)`
| Welch oneway ANOVA                 | `oneway_welch(x, g)`    | `oneway.test(x, g)`
| Equal variance oneway ANOVA        | `oneway_equalvar(x, g)` | `oneway.test(x, g, var.equal=TRUE)`
| Kruskal-Wallis test                | `kruskalwallis(x, g)`   | `kruskal.test(x, g)`
| Bartlett's test                    | `bartlett(x, g)`        | `bartlett.test(x, g)`

## Test-Based Algorithms ##

|             Description             |      matrixTests       |       R equivalent
|-------------------------------------|------------------------|-------------------------------------
| EVORA                               | `ievora(x, g)`         | ---

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


