[![CRAN version](http://www.r-pkg.org/badges/version/matrixTests)](https://cran.r-project.org/package=matrixTests)
[![Build Status](https://travis-ci.com/karoliskoncevicius/matrixTests.svg?branch=master)](https://travis-ci.com/karoliskoncevicius/matrixTests)
[![codecov](https://codecov.io/gh/karoliskoncevicius/matrixTests/branch/master/graph/badge.svg)](https://codecov.io/gh/karoliskoncevicius/matrixTests)
[![dependencies](https://tinyverse.netlify.com/badge/matrixTests)](https://CRAN.R-project.org/package=matrixTests)
[![Monthly Downloads](https://cranlogs.r-pkg.org/badges/matrixTests)](https://cranlogs.r-pkg.org/badges/matrixTests)

# Matrix Tests #

A package dedicated to running multiple statistical hypothesis tests on rows and columns of matrices.

![illustration](http://karolis.koncevicius.lt/data/matrixtests/illustration.png)

## Goals ##

1. Fast execution via vectorization.
2. Convenient and detailed output format.
3. Compatibility with tests implemented in base R.
4. Careful handling of missing values and edge cases.

## Examples ##

#### 1. Bartlett's test on columns ####

Bartlett's test on every column of iris dataset using Species as groups:

```r
col_bartlett(iris[,-5], iris$Species)
```
```
             obs.tot obs.groups var.pooled df statistic                pvalue
Sepal.Length     150          3 0.26500816  2 16.005702 0.0003345076070163084
Sepal.Width      150          3 0.11538776  2  2.091075 0.3515028004158132768
Petal.Length     150          3 0.18518776  2 55.422503 0.0000000000009229038
Petal.Width      150          3 0.04188163  2 39.213114 0.0000000030547839322
```

#### 2. Welch t-test on rows ####

Welch t-test performed on each row of 2 large (million row) matrices:

```r
X <- matrix(rnorm(10000000), ncol = 10)
Y <- matrix(rnorm(10000000), ncol = 10)

row_t_welch(X, Y)  # running time: 2.4 seconds
```

Confidence interval computations can be turned-off for further increase in speed:

```r
row_t_welch(X, Y, conf.level = NA)  # running time: 1 second
```

## Available Tests ##

|           Test                   |           Function                       |
|----------------------------------|------------------------------------------|
| **Location tests (1 group)**     |                                          |
| Single sample Student's t.test   | `row_t_onesample(x)`                     |
| Single sample Wilcoxon's test    | `row_wilcoxon_onesample(x)`              |
|                                  |                                          |
| **Location tests (2 groups)**    |                                          |
| Equal variance Student's t.test  | `row_t_equalvar(x, y)`                   |
| Welch adjusted Student's t.test  | `row_t_welch(x, y)`                      |
| Two sample Wilcoxon's test       | `row_wilcoxon_twosample(x, y)`           |
|                                  |                                          |
| **Location tests (paired)**      |                                          |
| Paired Student's t.test          | `row_t_paired(x, y)`                     |
| Paired Wilcoxon's test           | `row_wilcoxon_paired(x, y)`              |
|                                  |                                          |
| **Location tests (2+ groups)**   |                                          |
| Equal variance oneway anova      | `row_oneway_equalvar(x, g)`              |
| Welch's oneway anova             | `row_oneway_welch(x, g)`                 |
| Kruskal-Wallis test              | `row_kruskalwallis(x, g)`                |
| van der Waerden's test           | `row_waerden(x, g)`                      |
|                                  |                                          |
| **Scale tests (2 groups)**       |                                          |
| F variance test                  | `row_f_var(x, y)`                        |
|                                  |                                          |
| **Scale tests (2+ groups)**      |                                          |
| Bartlett's test                  | `row_bartlett(x, g)`                     |
| Fligner-Killeen test             | `row_flignerkilleen(x, g)`               |
| Levene's test                    | `row_levene(x, g)`                       |
| Brown-Forsythe test              | `row_brownforsythe(x, g)`                |
|                                  |                                          |
| **Association tests**            |                                          |
| Pearson's correlation test       | `row_cor_pearson(x, y)`                  |
|                                  |                                          |
| **Periodicity tests**            |                                          |
| Cosinor                          | `row_cosinor(x, t, period)`              |
|                                  |                                          |
| **Distribution tests**           |                                          |
| Kolmogorov-Smirnov test          | `row_kolmogorovsmirnov_twosample(x, y)`  |
|                                  |                                          |
| **Normality tests**              |                                          |
| Jarque-Bera test                 | `row_jarquebera(x)`                      |
| Anderson-Darling test            | `row_andersondarling(x)`                 |


## Further Information ##

For more information please refer to the [Wiki](https://github.com/karoliskoncevicius/matrixTests/wiki) page:

1. [Installation Instructions](https://github.com/karoliskoncevicius/matrixTests/wiki/Installation)
2. [Design Decisions](https://github.com/karoliskoncevicius/matrixTests/wiki/Design-Decisions)
3. [Speed Benchmarks](https://github.com/karoliskoncevicius/matrixTests/wiki/Benchmarks)
4. [Bug Fixes and Improvements to Base R](https://github.com/karoliskoncevicius/matrixTests/wiki/Bug-Fixes-and-Improvements-to-Base-R)


## See Also ##

### Literature ###

**Computing thousands of test statistics simultaneously in R**, *Holger Schwender, Tina Müller*.\
Statistical Computing & Graphics. Volume 18, No 1, June 2007.

### Packages ###

CRAN:

1. `ttests()` in the [**Rfast**](https://CRAN.R-project.org/package=Rfast) package.
2. `row.ttest.stat()` in the [**metaMA**](https://CRAN.R-project.org/package=metaMA) package.
3. `MultiTtest()` in the [**ClassComparison**](https://CRAN.R-project.org/package=ClassComparison) package.
4. `bartlettTests()` in the [**heplots**](https://CRAN.R-project.org/package=heplots) package.
5. `harmonic.regression()` in the [**HarmonicRegression**](https://CRAN.R-project.org/package=HarmonicRegression) package.

BioConductor:

1. `lmFit()` in the [**limma**](https://bioconductor.org/packages/release/bioc/html/limma.html) package.
2. `rowttests()` in the [**genefilter**](https://bioconductor.org/packages/release/bioc/html/genefilter.html) package.
3. `mt.teststat()` in the [**multtest**](https://www.bioconductor.org/packages/release/bioc/html/multtest.html) package.
4. `row.T.test()` in the [**HybridMTest**](https://www.bioconductor.org/packages/release/bioc/html/HybridMTest.html) package.
5. `rowTtest()` in the [**viper**](https://bioconductor.org/packages/release/bioc/html/viper.html) package.
6. `lmPerGene()` in the [**GSEAlm**](https://www.bioconductor.org/packages/release/bioc/html/GSEAlm.html) package.

GitHub:

1. `rowWilcoxonTests()` in the [**sanssouci**](https://github.com/pneuvial/sanssouci) package.
2. `matrix.t.test()` in the [**pi0**](https://github.com/gitlongor/pi0) package.
3. `wilcoxauc()` in the [**presto**](https://github.com/immunogenomics/presto) package.
