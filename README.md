
# rankresponse

<!-- badges: start -->
[![R-CMD-check](https://github.com/cmatKhan/rankresponse/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cmatKhan/rankresponse/actions/workflows/R-CMD-check.yaml)
[![test-codecov](https://github.com/cmatKhan/rankresponse/actions/workflows/test-codecov.yaml/badge.svg)](https://github.com/cmatKhan/rankresponse/actions/workflows/test-codecov.yaml)
<!-- badges: end -->

The purpose of `rankresponse` is to provide a method of evaluating the 
concordance between data derived from a transcription factor binding signal 
assay (eg Calling Cards or ChIP) and expression data derived from perturbing 
the same transcription factor (eg knock-out or overexpression). The procedure 
is simple:  

1. expression and binding data are (inner) joined on the gene feature

2. the expression data effect (eg log2FoldChange of perturbed / control) are 
binarized by the effect magnitude and p-value. Default thresholds are 0 for the
log2FoldChange and 0.05 for the p-value.

3. The table is next ranked from least to greatest by the binding data p-value. 
For example:

| gene         | binding_signal     | responsive |
|--------------|--------------------|------------|
| a            | .0000000001        | TRUE       |
| f            | .000000001         | TRUE       |
| e            | .00000001          | FALSE      |
| z            | .0000001           | TRUE       |
| v            | .000001            | FALSE      |
| w            | .00001             | FALSE      |
| b            | .0001              | TRUE       |
| d            | .001               | FALSE      |
| g            | .01                | FALSE      |
| h            | .1                 | FALSE      |

4. Given a bin resolution (default 5) the cumulative responsiveness, determined 
by the binarized expression data, is calculated. From the example above, 
the summary would be:

| bin         | responsive_ratio    |
|-------------|---------------------|
| 1           | 0.60 (3 TRUE of 5)  |
| 2           | 0.40 (4 TRUE of 10) |


## Installation

To install the current release:

``` r
install.packages('rankresponse')
```

To install the development version:

``` r
# install.packages("devtools")
devtools::install_github("cmatKhan/rankresponse")
```


