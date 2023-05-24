PhantasusLite – a tool designed to integrate the work with RNA-Seq count
matrices into a single and fast R-pipeline. This R-package supports
loading the data from GEO by GSE. It provides url access to the remote
repository with the archs4, archs4_zoo and dee2 HDF5-files for getting
the count matrix. Finally phantasusLite allows to get an ExpressionSet
with the expression matrix for future differential expression analysis.

To install the package you should:

1.  Install R-package rhdf5client from the github
2.  Install R-package phantasusLite from the github

``` r
library(devtools)
install_github("assaron/rhdf5client")
install_github("rsablina/phantasusLite")
```

To run the package enter the code sample below. ExpressionSet from the
GEO doesn’t contain the expression matrix – exprs(es) is empty. Function
loadCountsFromHSDS returns an ExpressionSet with the expression matrix –
the second exprs(es) contains an expression matrix.

The remote repository URL is
‘<https://ctlab.itmo.ru/hsds/?domain=/counts>’.

``` r
library(GEOquery)
library(rhdf5client)
library(phantasusLite)
ess <- getGEO("GSE85653")
url <- 'https://ctlab.itmo.ru/hsds/?domain=/counts'
es <- ess[[1]]
head(exprs(es))
es <- loadCountsFromHSDS(es, url)
head(exprs(es))
```
