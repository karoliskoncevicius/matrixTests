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
res2 <- row.t.welch(X, Y)

```

## Goals ##

1. Fast execution via vectorization.
2. Output that is detailed and easy to use.
3. Result compatibility with tests that are implemented in R.

## Available Tests ##

|             Name                   |      matrixTests            |       R equivalent
|------------------------------------|-----------------------------|-------------------------------------
| Single sample t.test               | `row.t.onesample(x)`        | `t.test(x)`
| Welch t.test                       | `row.t.welch(x, y)`         | `t.test(x, y)`
| Equal variance t.test              | `row.t.equalvar(x, y)`      | `t.test(x, y, var.equal=TRUE)`
| Paired t.test                      | `row.t.paired(x, y)`        | `t.test(x, y, paired=TRUE)`
| Pearson's correlation test         | `row.cor.pearson(x, y)`     | `cor.test(x, y)`
| Welch oneway ANOVA                 | `row.oneway.welch(x, g)`    | `oneway.test(x, g)`
| Equal variance oneway ANOVA        | `row.oneway.equalvar(x, g)` | `oneway.test(x, g, var.equal=TRUE)`
| Kruskal-Wallis test                | `row.kruskalwallis(x, g)`   | `kruskal.test(x, g)`
| Bartlett's test                    | `row.bartlett(x, g)`        | `bartlett.test(x, g)`

## Test-Based Algorithms ##

|             Description             |      matrixTests       |       R equivalent
|-------------------------------------|------------------------|-------------------------------------
| EVORA                               | `row.ievora(x, g)`     | ---

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

## Design Decisions ##

The following design decision were taken (in no particular order):

#### Function Names ####

Function names are composed of 3 elements separated by dots where needed:

`[row/col].testname.variant`

The variant part can be dropped when not applicable or when only a single
variant for that test is implemented so far.

A few examples: `row.oneway.equalvar`, `row.bartlett`

#### Compatibility with R ####

Results of tests should be as compatible as possible with the ones implemented
in base R. Allowed exceptions are cases where R's implementation is incorrect
or limiting.

A good example is `oneway.test()` which works only if all groups have more than
2 observations even when *var.equal* is set to *TRUE*. The strict requirement
for the test to run technically is for at least one group to have more than 1
observation. Therefore in such cases `row.oneway.equalvar()` works even if base
R version throws an error.

For another such example consider `bartlett.test()`. This function works without
any warnings when supplied with constant data and returns NA values:

`bartlett.test(rep(1,4), c("a","a","b","b"))`.

The typical behaviour in such situations for base R tests is to throw an error:

`t.test(rep(1,4) ~ c("a","a","b",b"))`.

Functions in this package try to be consistent with each other and be as
informative as possible. Therefore in such cases `row.bartlett()` will throw a
warning even if base R function does not.

#### Warnings and Errors ####

Errors are produced only when the input parameters are not correctly specified.

Warnings are shown in situations when there was something wrong with doing a
test itself given the specified input parameters.

Such a decision for warnings was taken because users will typically perform
multiple tests (one for each row). The function cannot fail when one or few
of those tests cannot be completed. So even when R base tests throw an error
the functions in this package will instead produce an informative warning for
the row that failed and if needed will set all it's return values related to the
test to NA so that the user will not be able to use them by mistake.

Note that in these cases only test-related values like test statistic, p-value
and confidence interval are set to NA. Other returned values: number of
observations, means, variances and similar will still be returned as usual.

As an example of such behaviour consider the case when base t-test with Welch
correction fails because it has not enough observations:

`t.test(c(1,2), 3)`

Function in this package proceeds, but throws a warning and takes care to set
the failed outputs to NA:

`row.t.welch(c(1,2), 3)`

This allows us to continue working in cases where typically we have enough
observations per group but some rows might not have enough due to NA values.

```r
mat1 <- rbind(c(1,2), c(3,NA))
mat2 <- rbind(c(2,3), c(0,4))
`row.t.welch(mat1, mat2)`
```

#### NA and NaN values ####

**NA** and **NaN** values from the input matrices are silently removed and each
row is treated like a vector that has no NA/NaN values.

When **NA** or **NaN** values are present in the parameter specifying the groups
the corresponding values from the input matrices are dropped before doing the
tests. For example if the specified group variable is `c(NA,"a", "a", "b", "b")`
then the entire first column from the input matrix corresponding to that group
will be removed.

Other parameters might allow or not allow NA values depending on context. For
example you cannot specify **NA** as wanted confidence level when doing a test
because not knowing your confidence level makes no sense.

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


