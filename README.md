```         
library(GEOquery)
ess <- getGEO("GSE164173")
src <- HSDSSource('https://ctlab.itmo.ru/hsds/')
ess <- ess[[1]]
ess <- loadCountsFromHSDS(es, src, '/counts')
```
