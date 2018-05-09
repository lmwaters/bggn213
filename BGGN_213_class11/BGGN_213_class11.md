BGGN 213 class 11
================

PDB statistics
==============

``` r
pdbstat <- read.csv("Data Export Summary.csv", row.names=1)
percent <- (pdbstat$Total/sum(pdbstat$Total))*100
names(percent) <- row.names(pdbstat)
percent
```

    ##               X-Ray                 NMR Electron Microscopy 
    ##         89.51673340          8.71321614          1.51239392 
    ##               Other        Multi Method 
    ##          0.16986775          0.08778879
