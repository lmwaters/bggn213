---
title: "Untitled"
output: 
  html_document: 
    keep_md: yes
---

map.colors <- function (value,high.low,palette) {
  proportion <- ((value-high.low[1])/(high.low[2]-high.low[1]))
  index <- round ((length(palette)-1)*proportion)+1
  return (palette[index])
}



```r
map.colors2 <- function(x, high.low = range(x) , palette = cm.colors(100)) {

 #determine where in hig.low range value x lies
   percent <- ((x - high.low[1]) / (high.low[2] - high.low[1]))
 
 #where in palette vector is this percent
  index <- round( (length(palette)-1) *percent)+1

   return(palette[index])
 }
```
1st function  

```r
add <- function(x, y=1) {
 # Sum the input x and y
 x + y
}
```
2nd function

```r
rescale <- function(x) {
 rng <-range(x)
 (x - rng[1]) / (rng[2] - rng[1])
}
```
test 2nd function

```r
rescale(c(1,2,NA,3,10))
```

```
## [1] NA NA NA NA NA
```

```r
rescale(1:10)
```

```
##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
##  [8] 0.7777778 0.8888889 1.0000000
```

```r
rescale2<- function(x) {
 rng <-range(x, na.rm=TRUE)
 (x - rng[1]) / (rng[2] - rng[1])
}
```
test rescale2

```r
rescale2(c(1,2,NA,3,10))
```

```
## [1] 0.0000000 0.1111111        NA 0.2222222 1.0000000
```


```r
rescale3 <- function(x, na.rm=TRUE, plot=TRUE) {
 if(na.rm) {
 rng <-range(x, na.rm=na.rm)
 } else {
 rng <-range(x)
 }
 print("Hello")
 answer <- (x - rng[1]) / (rng[2] - rng[1])
 
 print("is it me you are looking for?")
 if(plot) {
 plot(answer, typ="b", lwd=4)
 }
 print("I can see it in ...")
  return(answer)
}
```
test rescale3

```r
rescale3(c(1,2,NA,3,10))
```

```
## [1] "Hello"
## [1] "is it me you are looking for?"
```

![](mapcolor_markdown_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```
## [1] "I can see it in ..."
```

```
## [1] 0.0000000 0.1111111        NA 0.2222222 1.0000000
```
## Section 2b class 6

```r
library(bio3d)
s1 <- read.pdb("4AKE") # kinase with drug
```

```
##   Note: Accessing on-line PDB file
```

```r
s1
```

```
## 
##  Call:  read.pdb(file = "4AKE")
## 
##    Total Models#: 1
##      Total Atoms#: 3459,  XYZs#: 10377  Chains#: 2  (values: A B)
## 
##      Protein Atoms#: 3312  (residues/Calpha atoms#: 428)
##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
## 
##      Non-protein/nucleic Atoms#: 147  (residues: 147)
##      Non-protein/nucleic resid values: [ HOH (147) ]
## 
##    Protein sequence:
##       MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVT
##       DELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRI
##       VGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIG
##       YYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILGMRIILLGAPGA...<cut>...KILG
## 
## + attr: atom, xyz, seqres, helix, sheet,
##         calpha, remark, call
```


```r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s1.chainA
```

```
## 
##  Call:  trim.pdb(pdb = s1, chain = "A", elety = "CA")
## 
##    Total Models#: 1
##      Total Atoms#: 214,  XYZs#: 642  Chains#: 1  (values: A)
## 
##      Protein Atoms#: 214  (residues/Calpha atoms#: 214)
##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
## 
##      Non-protein/nucleic Atoms#: 0  (residues: 0)
##      Non-protein/nucleic resid values: [ none ]
## 
##    Protein sequence:
##       MRIILLGAPGAGKGTQAQFIMEKYGIPQISTGDMLRAAVKSGSELGKQAKDIMDAGKLVT
##       DELVIALVKERIAQEDCRNGFLLDGFPRTIPQADAMKEAGINVDYVLEFDVPDELIVDRI
##       VGRRVHAPSGRVYHVKFNPPKVEGKDDVTGEELTTRKDDQEETVRKRLVEYHQMTAPLIG
##       YYSKEAEAGNTKYAKVDGTKPVAEVRADLEKILG
## 
## + attr: atom, helix, sheet, seqres, xyz,
##         calpha, call
```


```r
library(bio3d)
s1 <- read.pdb("4AKE") # kinase with drug
```

```
##   Note: Accessing on-line PDB file
```

```
## Warning in get.pdb(file, path = tempdir(), verbose = FALSE): C:\Users\laura
## \AppData\Local\Temp\RtmpyuaOzX/4AKE.pdb exists. Skipping download
```

```r
s2 <- read.pdb("1AKE") # kinase no drug
```

```
##   Note: Accessing on-line PDB file
##    PDB has ALT records, taking A only, rm.alt=TRUE
```

```r
s3 <- read.pdb("1E4Y") # kinase with drug
```

```
##   Note: Accessing on-line PDB file
```

```r
s1.chainA <- trim.pdb(s1, chain="A", elety="CA")
s2.chainA <- trim.pdb(s2, chain="A", elety="CA")
s3.chainA <- trim.pdb(s1, chain="A", elety="CA")
s1.b <- s1.chainA$atom$b
s2.b <- s2.chainA$atom$b
s3.b <- s3.chainA$atom$b
```


```r
hc <- hclust( dist( rbind(s1.b, s2.b, s3.b) ) )
plot(hc)
```

![](mapcolor_markdown_files/figure-html/unnamed-chunk-13-1.png)<!-- -->


