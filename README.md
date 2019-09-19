[![CRAN version](http://www.r-pkg.org/badges/version/matrixTests)](https://cran.r-project.org/package=matrixTests)
[![Build Status](https://travis-ci.com/KKPMW/matrixTests.svg?branch=master)](https://travis-ci.com/KKPMW/matrixTests)
[![dependencies](https://tinyverse.netlify.com/badge/matrixTests)](https://CRAN.R-project.org/package=matrixTests)
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

**1) Welch t-test on each row of 2 large (million row) matrices**

```r
X <- matrix(rnorm(10000000), ncol=10)
Y <- matrix(rnorm(10000000), ncol=10)

row_t_welch(X, Y)  # running time: 2.4 seconds
```
```
  obs.x obs.y obs.tot      mean.x     mean.y  mean.diff    var.x     var.y stderr          df  statistic    pvalue   conf.low conf.high alternative mean.null conf.level
1    10    10      20 -0.06643757 -0.2985907  0.2321531 1.627547 0.9140158 0.5041392 16.68493  0.4604941 0.6511065 -0.8330197 1.2973259   two.sided         0       0.95
2    10    10      20 -0.02447724  0.4805317 -0.5050090 1.424720 1.2936936 0.5213841 17.95828 -0.9685930 0.3456133 -1.6005787 0.5905608   two.sided         0       0.95
```

**2) One way ANOVA on every column of iris data using Species as groups**

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

## Available Tests ##

|             Name                             |      matrixTests               |       R equivalent
|----------------------------------------------|--------------------------------|-------------------------------------
| **Location tests (1 group)**                 |                                |
| &nbsp; &nbsp; Single sample t.test           | `row_t_onesample(x)`           | `t.test(x)`
| &nbsp; &nbsp; Single sample Wilcoxon test    | `row_wilcoxon_onesample(x)`    | `wilcox.test(x)`
| **Location tests (2 groups)**                |                                |
| &nbsp; &nbsp; Equal variance t.test          | `row_t_equalvar(x, y)`         | `t.test(x, y, var.equal=TRUE)`
| &nbsp; &nbsp; Welch t.test                   | `row_t_welch(x, y)`            | `t.test(x, y)`
| &nbsp; &nbsp; Two sample Wilcoxon test       | `row_wilcoxon_twosample(x, y)` | `wilcox.test(x, y)`
| **Location tests (2+ groups)**               |                                |
| &nbsp; &nbsp; Equal variance oneway ANOVA    | `row_oneway_equalvar(x, g)`    | `oneway.test(x ~ g, var.equal=TRUE)`
| &nbsp; &nbsp; Welch oneway ANOVA             | `row_oneway_welch(x, g)`       | `oneway.test(x ~ g)`
| &nbsp; &nbsp; Kruskal-Wallis test            | `row_kruskalwallis(x, g)`      | `kruskal.test(x, g)`
| **Location tests (paired)**                  |                                |
| &nbsp; &nbsp; Paired t.test                  | `row_t_paired(x, y)`           | `t.test(x, y, paired=TRUE)`
| &nbsp; &nbsp; Paired Wilcoxon test           | `row_wilcoxon_paired(x, y)`    | `wilcox.test(x, y, paired=TRUE)`
| **Scale tests (2 groups)**                   |                                |
| &nbsp; &nbsp; F variance test                | `row_f_var(x, y)`              | `var.test(x, y)`
| **Scale tests (2+ groups)**                  |                                |
| &nbsp; &nbsp; Bartlett's test                | `row_bartlett(x, g)`           | `bartlett.test(x, g)`
| &nbsp; &nbsp; Fligner-Killeen test           | `row_flignerkilleen(x, g)`     | `fligner.test(x, g)`
| &nbsp; &nbsp; Levene's test                  | `row_levene(x, g)`             | `car:leveneTest(x, g, center="mean")`
| &nbsp; &nbsp; Brown-Forsythe test            | `row_brownforsythe(x, g)`      | `car:leveneTest(x, g, center="median")`
| **Assosiations tests**                       |                                |
| &nbsp; &nbsp; Pearson's correlation test     | `row_cor_pearson(x, y)`        | `cor.test(x, y)`
| **Distribution tests**                       |                                |
| &nbsp; &nbsp; Jarque-Bera test               | `row_jarquebera(x)`            | `moments::jarque.test(x)`
| **Test based procedures**                    |                                |
| &nbsp; &nbsp; EVORA                          | `row_ievora(x, b)`             | ---

## Further Information ##

For more information please refer to the [Wiki](https://github.com/KKPMW/matrixTests/wiki) page:

1. [Installation Instructions](https://github.com/KKPMW/matrixTests/wiki/Installation)
2. [Design Decisions](https://github.com/KKPMW/matrixTests/wiki/Design-Decisions)
3. [Speed Benchmarks](https://github.com/KKPMW/matrixTests/wiki/Benchmarks)
4. [Future Plans](https://github.com/KKPMW/matrixTests/wiki/Future-Plans)
5. [Possible Shortcomings](https://github.com/KKPMW/matrixTests/wiki/Possible-Shortcomings)

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
