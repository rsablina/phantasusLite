``` r
library(GEOquery)
```

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Setting options('download.file.method.GEOquery'='auto')

    ## Setting options('GEOquery.inmemory.gpl'=FALSE)

``` r
library(rhdf5client)
```

    ## Loading required package: DelayedArray

    ## Loading required package: stats4

    ## Loading required package: Matrix

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## Warning: package 'matrixStats' was built under R version 4.2.3

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     rowMedians

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:Matrix':
    ## 
    ##     expand, unname

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:base':
    ## 
    ##     apply, rowsum, scale, sweep

``` r
library(phantasusLite)
ess <- getGEO("GSE164173")
```

    ## Found 1 file(s)

    ## GSE164173_series_matrix.txt.gz

``` r
src <- HSDSSource('https://ctlab.itmo.ru/hsds/')
ess <- ess[[1]]
ess <- loadCountsFromHSDS(ess, src, '/counts')
head(exprs(ess))
```

    ##               GSM4998773 GSM4998774 GSM4998775 GSM4998776 GSM4998777 GSM4998778
    ## 0610007P14Rik        672        531        556        611        789        892
    ## 0610009B22Rik         59         61         96        135        111        155
    ## 0610009L18Rik          1         30         36         12         21         37
    ## 0610009O20Rik        147        183        191        140        211        208
    ## 0610010F05Rik        145        115         96        109         98        194
    ## 0610010K14Rik        566        408        488        406        663        518
    ##               GSM4998779 GSM4998780
    ## 0610007P14Rik        743        938
    ## 0610009B22Rik        168        298
    ## 0610009L18Rik         11         27
    ## 0610009O20Rik        168        238
    ## 0610010F05Rik         80         47
    ## 0610010K14Rik        526        650
