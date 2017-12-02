# Matrix Tests #

A package dedicated to running statistical hypothesis tests on rows of matrices.

## Example ##

T-test on each row of 2 matrices each with a million rows.

```r
X <- matrix(rnorm(10000000), ncol=10)
Y <- matrix(rnorm(10000000), ncol=10)
```

#### The usual way ####

```r
res1 <- vector(nrow(X), mode="list")
for(i in 1:nrow(X)) {
  res1[[i]] <- t.test(X[i,], Y[i,])
}
# RUN TIME: 2 minutes 13 seconds
```

Output for first 2 rows:

```
  res1[1:2]
[[1]]

        Welch Two Sample t-test

data:  X[i, ] and Y[i, ]
t = -0.42194, df = 17.989, p-value = 0.6781
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -1.2311709  0.8193669
sample estimates:
 mean of x  mean of y
0.06435162 0.27025362

[[2]]

        Welch Two Sample t-test

data:  X[i, ] and Y[i, ]
t = 0.18962, df = 15.213, p-value = 0.8521
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.9089581  1.0867183
sample estimates:
 mean of x  mean of y
-0.1360916 -0.2249717
```

#### matrixTest way ####

```r
res2 <- row.t.welch(X, Y)
# RUN TIME: 2.3 seconds
```

Output for first 2 rows:

```
> res2[1:2,]
       mean.x     mean.y   mean.diff     var.x    var.y obs.x obs.y obs.tot statistic.t   p.value     ci.low   ci.high    stderr       df mean.null conf.level alternative
1  0.06435162  0.2702536 -0.20590200 1.2207363 1.160574    10    10      20  -0.4219418 0.6780672 -1.2311709 0.8193669 0.4879867 17.98852         0       0.95   two.sided
2 -0.13609158 -0.2249717  0.08888009 0.6282957 1.568692    10    10      20   0.1896228 0.8521116 -0.9089581 1.0867183 0.4687204 15.21276         0       0.95   two.sided
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

## Test-Based Procedures ##

|             Description             |      matrixTests       |       R equivalent
|-------------------------------------|------------------------|-------------------------------------
| EVORA                               | `row.ievora(x, g)`     | ---

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

### Function Names ###

Function names are composed of 3 elements separated by dots where needed:

`[row/col].testname.variant`

The variant part can be dropped when not applicable or when only a single
variant for that test is implemented so far.

A few examples: `row.oneway.equalvar`, `row.bartlett`

In order to make the function names shorter the word *test* is not included.

### Single Test Per Function ###

Functions should provide a single type of test.

This means that some of the tests that in base R are implemented under a single
function will be split into different functions. For example in base R function
`t.test()` has parameters that can specify the type of t.test to use:
equal variance or Welch adjusted, paired or non paired.

In this package those types of choices are separated into separate functions
for the following reasons:

1. Output structure returned by a function is not dependant on the input values.
2. Users is made to choose the test explicitly with no hidden defaults.
3. This convention makes it easier to add more test types later.

### Input Values ###

The functions try to make sense of the provided inputs when possible.
All the cases when the inputs are incorrectly specified should throw an error.

Edge cases should be handled gracefully. For example when input is a numeric
`matrix` with 0 rows - the result is a 0 row `data.frame`.

Below is a short list of implemented input rules:

1. Main parameters are numeric matrices.
2. Vectors are transformed into a 1-row matrix.
3. Data frames are automatically transformed into matrices if all of their
   columns are numeric.
4. When two matrices are required - number of their rows should match.
5. Group specifications and additional parameters typically can have only one
   value that will be applied to all the rows.
6. In some cases additional parameters can have a separate value for each row.
   Those cases are specified in the help documentation.


### Outputs ###

Outputs are contained in a `data.frame` with each row providing the result of
the test performed on the corresponding row of the input matrix.

#### Output categories ####

Columns hold the relevant results that can be divided into 3 main categories:

1. Descriptive statistics related to the test. Typically ordered by increasing
complexity. For example: 1) number of observations, 2) means, 3) variances.

2. Outputs related to the test itself. Typically ordered by increasing complexity.
For example: 1) degrees of freedom, 2) test statistic, 3) p-value, 4) confidence interval.

3. Input parameters that were chosen for the row. Ordered by their appearance in
the function call. For example: 1) alternative hypothesis type, 2) mean of null hypothesis, 3) confidence level.

#### Column names ####

Column names of the output are written in consistent fashion and typically have
two parts: type of result and specification, separated by a dot. Some examples:

* obs.x - number of x observations
* obs.tot - total number of observations
* mean.x - mean of x
* mean.diff - mean of x and y difference
* ci.low - lower confidence interval
* statistic.t - t statistic of the test

Exception from this rule are values that have a dot or dash in their name.
Like *p.value*.

#### Row names ####

Row names are transfered from the main input matrix. If the row names of the
matrix were not unique - the are made unique using `make.unique()`. In case
input matrix had no row names the numbers `1:nrow(x)` are used.


### Compatibility with R ###

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

`bartlett.test(rep(1,4), c("a","a","b","b"))`

The typical behaviour in such situations for base R tests is to throw an error:

`t.test(rep(1,4) ~ c("a","a","b",b"))`

Functions in this package try to be consistent with each other and be as
informative as possible. Therefore in such cases `row.bartlett()` will throw a
warning even if base R function does not.

### Warnings and Errors ###

Errors are produced only when the input parameters are not correctly specified.

Warnings are shown in situations when there was something wrong with doing a
test itself given the specified input parameters.

Such a decision for warnings was taken because users will typically perform
multiple tests (one for each row). The function cannot fail when one or few
of those tests cannot be completed. So even when R base tests throw an error
the functions in this package will instead produce an informative warning for
the row that failed and if needed will set all it's return values related to the
test to `NA` so that the user will not be able to use them by mistake.

Note that in these cases only test-related values like test statistic, p-value
and confidence interval are set to NA. Other returned values: number of
observations, means, variances and similar will still be returned as usual.

As an example of such behaviour consider the case when base t-test with Welch
correction fails because it has not enough observations:

```r
t.test(c(1,2), 3)
```
```
Error in t.test.default(c(1, 2), 3) : not enough 'y' observations
```

Function in this package proceeds, but throws a warning and takes care to set
the failed outputs to NA:

```r
row.t.welch(c(1,2), 3)
```
```
  mean.x mean.y mean.diff var.x var.y obs.x obs.y obs.tot statistic.t p.value ci.low ci.high stderr  df mean.null conf.level alternative
1    1.5      3      -1.5   0.5   NaN     2     1       3          NA      NA     NA      NA    NaN NaN         0       0.95   two.sided
Warning message:
In showWarning(w2, "had less than 2 \"y\" observations") :
  1 of the rows had less than 2 "y" observations. First occurrence at row 1
```

This allows the function continue working in cases where typically we have enough
observations per group but some rows might not have enough due to NA values.

```r
mat1 <- rbind(c(1,2), c(3,NA))
mat2 <- rbind(c(2,3), c(0,4))
row.t.welch(mat1, mat2)
```
```
  mean.x mean.y mean.diff var.x var.y obs.x obs.y obs.tot statistic.t   p.value    ci.low  ci.high    stderr  df mean.null conf.level alternative
1    1.5    2.5        -1   0.5   0.5     2     2       4   -1.414214 0.2928932 -4.042435 2.042435 0.7071068   2         0       0.95   two.sided
2    3.0    2.0         1   NaN   8.0     1     2       3          NA        NA        NA       NA       NaN NaN         0       0.95   two.sided
Warning message:
In showWarning(w1, "had less than 2 \"x\" observations") :
  1 of the rows had less than 2 "x" observations. First occurrence at row 2
```

### NA and NaN values ###

`NA` and `NaN` values from the input matrices are silently removed and each
row is treated like a vector that has no `NA`/`NaN` values.

When `NA` or `NaN` values are present in the parameter specifying the groups
the corresponding values from the input matrices are dropped before doing the
tests. For example if the specified group variable has a `NA`:

```r
x <- rnorm(5)
g <- c(NA,"a", "a", "b", "b")
row.oneway.welch(x=x, g=g)
```
```
  obs.tot obs.groups df.treatment df.residuals statistic.F   p.value
1       4          2            1     1.488032    5.461371 0.1863233
```

then the entire first column from the input matrix x corresponding to that group
will be removed. And the result will be equivalent to:

```r
row.oneway.welch(x=x[-1], g=g[-1])
```
```
  obs.tot obs.groups df.treatment df.residuals statistic.F   p.value
1       4          2            1     1.488032    5.461371 0.1863233
```

Other parameters might allow or not allow `NA` values depending on context. For
example you cannot specify `NA` as wanted confidence level when doing a test
because not knowing your confidence level makes little sense.

## Notes ##

All the tests are implemented in R. So when running a test on a single row there
should be no increase in execution speed compared with the base R versions. In
most cases a slight decrease is expected due to the more detailed output.

For now the column-wise versions of the tests simply transposes the input
matrix and calls the equivalent row-wise test.

Candidates of tests that will be implemented next:

1. Shapiro-Wilks test for normality
2. Spearman and Kendall correlation tests
3. Test for proportions
4. Fisher's exact test

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


