[![CRAN version](http://www.r-pkg.org/badges/version/matrixTests)](https://cran.r-project.org/package=matrixTests)
[![Build Status](https://travis-ci.com/KKPMW/matrixTests.svg?branch=master)](https://travis-ci.com/KKPMW/matrixTests)
[![codecov](https://codecov.io/gh/KKPMW/matrixTests/branch/master/graph/badge.svg)](https://codecov.io/gh/KKPMW/matrixTests)
[![Monthly Downloads](https://cranlogs.r-pkg.org/badges/matrixTests)](https://cranlogs.r-pkg.org/badges/matrixTests)

# Matrix Tests #

A package dedicated to running multiple statistical hypothesis tests on rows and columns of matrices.

## Goals ##

1. Fast execution via vectorization.
2. Handling of edge cases (NA values, 0 row inputs).
3. Output that is detailed and easy to use.
4. Result compatibility with tests that are implemented in R.

## Examples ##

**1) Running one way ANOVA on every column of iris data using Species as groups**

```r
col_oneway_equalvar(iris[,-5], iris$Species)
```
```
             obs.tot obs.groups  sumsq.between  sumsq.within  meansq.between  meansq.within df.between df.within  statistic       pvalue
Sepal.Length     150          3       63.21213       38.9562       31.606067     0.26500816          2       147  119.26450 1.669669e-31
Sepal.Width      150          3       11.34493       16.9620        5.672467     0.11538776          2       147   49.16004 4.492017e-17
Petal.Length     150          3      437.10280       27.2226      218.551400     0.18518776          2       147 1180.16118 2.856777e-91
Petal.Width      150          3       80.41333        6.1566       40.206667     0.04188163          2       147  960.00715 4.169446e-85
```

**2) t-test on each row of 2 matrices each with a million rows (matrixTests version vs simple t.test())**

```r
X <- matrix(rnorm(10000000), ncol=10)
Y <- matrix(rnorm(10000000), ncol=10)
```

**matrixTest row_t_welch() way** &#9200; 2.4 seconds

```r
res1 <- row_t_welch(X, Y)
```
```
> res1[1:2,]
  obs.x obs.y obs.tot      mean.x     mean.y  mean.diff    var.x     var.y stderr          df  statistic    pvalue   conf.low conf.high alternative mean.null conf.level
1    10    10      20 -0.06643757 -0.2985907  0.2321531 1.627547 0.9140158 0.5041392 16.68493  0.4604941 0.6511065 -0.8330197 1.2973259   two.sided         0       0.95
2    10    10      20 -0.02447724  0.4805317 -0.5050090 1.424720 1.2936936 0.5213841 17.95828 -0.9685930 0.3456133 -1.6005787 0.5905608   two.sided         0       0.95
```

**The usual t.test() way** &#9200; 2 minutes 16 seconds

```r
res2 <- vector(nrow(X), mode="list")

for(i in 1:nrow(X)) {
  res2[[i]] <- t.test(X[i,], Y[i,])
}
```

```
  res2[1:2]
[[1]]

        Welch Two Sample t-test

data:  X[i, ] and Y[i, ]
t = 0.46049, df = 16.685, p-value = 0.6511
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.8330197  1.2973259
sample estimates:
  mean of x   mean of y
-0.06643757 -0.29859071


[[2]]

        Welch Two Sample t-test

data:  X[i, ] and Y[i, ]
t = -0.96859, df = 17.958, p-value = 0.3456
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -1.6005787  0.5905608
sample estimates:
  mean of x   mean of y
-0.02447724  0.48053173
```

## Available Tests ##

|             Name                   |      matrixTests            |       R equivalent
|------------------------------------|-----------------------------|-------------------------------------
| Single sample t.test               | `row_t_onesample(x)`        | `t.test(x)`
| Welch t.test                       | `row_t_welch(x, y)`         | `t.test(x, y)`
| Equal variance t.test              | `row_t_equalvar(x, y)`      | `t.test(x, y, var.equal=TRUE)`
| Paired t.test                      | `row_t_paired(x, y)`        | `t.test(x, y, paired=TRUE)`
| Pearson's correlation test         | `row_cor_pearson(x, y)`     | `cor.test(x, y)`
| Welch oneway ANOVA                 | `row_oneway_welch(x, g)`    | `oneway.test(x, g)`
| Equal variance oneway ANOVA        | `row_oneway_equalvar(x, g)` | `oneway.test(x, g, var.equal=TRUE)`
| Kruskal-Wallis test                | `row_kruskalwallis(x, g)`   | `kruskal.test(x, g)`
| Bartlett's test                    | `row_bartlett(x, g)`        | `bartlett.test(x, g)`
| Fligner-Killeen test               | `row_flignerkilleen(x)`     | `fligner.test(x)`
| Jarque-Bera test                   | `row_jarquebera(x)`         | `moments::jarque.test(x)`

## Test-Based Procedures ##

|             Description            |      matrixTests            |       R equivalent
|------------------------------------|-----------------------------|-----------------------------------------
| EVORA                              | `row_ievora(x, b)`          | ---

## Installation ##

From **CRAN**:

```r
install.packages("matrixTests")
```

Using the `devtools` library:

```r
library(devtools)
install_github("KKPMW/matrixTests")
```

To install the **developement version** (stable updates not yet on **CRAN**):

```r
library(devtools)
install_github("KKPMW/matrixTests", ref="dev")
```


## Further Information ##

For more information please refer to the [Wiki](https://github.com/KKPMW/matrixTests/wiki) page:

1. [Design Decisions](https://github.com/KKPMW/matrixTests/wiki/Design-Decisions)
2. [Future Plans](https://github.com/KKPMW/matrixTests/wiki/Future-Plans)
3. [Possible Shortcomings](https://github.com/KKPMW/matrixTests/wiki/Possible-Shortcomings)

## See Also ##

### Literature ###

**Computing thousands of test statistics simultaneously in R**,
*Holger Schwender, Tina Müller*. Statistical Computing & Graphics. Volume 18, No 1, June 2007.

### Packages ###

CRAN:

1. `ttests()` in the [**Rfast**](https://CRAN.R-project.org/package=Rfast) package.
2. `row.ttest.stat()` in the [**metaMA**](https://CRAN.R-project.org/package=metaMA) package.
3. `MultiTtest()` in the [**ClassComparison**](https://CRAN.R-project.org/package=ClassComparison) package.
4. `bartlettTests()` in the [**heplots**](https://CRAN.R-project.org/package=heplots) package.

BioConductor:

1. `lmFit()` in the [**limma**](https://bioconductor.org/packages/release/bioc/html/limma.html) package.
2. `rowttests()` in the [**genefilter**](https://bioconductor.org/packages/release/bioc/html/genefilter.html) package.
3. `mt.teststat()` in the [**multtest**](https://www.bioconductor.org/packages/release/bioc/html/multtest.html) package.
4. `row.T.test()` in the [**HybridMTest**](https://www.bioconductor.org/packages/release/bioc/html/HybridMTest.html) package.
5. `rowTtest()` in the [**viper**](https://bioconductor.org/packages/release/bioc/html/viper.html) package.
6. `lmPerGene()` in the [**GSEAlm**](https://www.bioconductor.org/packages/release/bioc/html/GSEAlm.html) package.

GitHub:

1. `rowWilcoxonTests()` in the [**sanssouci**](https://github.com/pneuvial/sanssouci) package.
