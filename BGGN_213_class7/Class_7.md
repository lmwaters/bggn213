---
title: "Class 7"
output: 
  html_document: 
    keep_md: yes
---


#Section 1 FUnctions
#can source any file of code with source()

```r
source("http://tinyurl.com/rescale-R")
```
#viewing environment

```r
ls()
```

```
##  [1] "both_na"         "both_na2"        "both_na3"       
##  [4] "df1"             "df2"             "df3"            
##  [7] "gene_intersect"  "gene_intersect2" "gene_intersect3"
## [10] "gene_intersect4" "rescale"         "rescale2"
```
#check rescale function

```r
rescale(1:10)
```

```
##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
##  [8] 0.7777778 0.8888889 1.0000000
```
#break function

```r
rescale(1:10, "string")
```

```r
rescale2(c(1:10, "string"))
```
#Function for finding missing values
#write a both_na() function

```r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
is.na(x)
```

```
## [1] FALSE FALSE  TRUE FALSE  TRUE
```

```r
which(is.na(x))
```

```
## [1] 3 5
```

```r
sum( is.na(x))
```

```
## [1] 2
```

```r
is.na(x)&is.na(y)
```

```
## [1] FALSE FALSE  TRUE FALSE FALSE
```

```r
sum(is.na(x)&is.na(y))
```

```
## [1] 1
```
both_na <- function(x, y) {
  ## Check for NA elements in both input vectors 
  sum( is.na(x) & is.na(y) )

```r
both_na(x,y)
```

```
## [1] 1
```

```r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
both_na(x, y2)
```
function(x, y) {
  ## Check for NA elements in both input vectors and don't allow re-cycling 
  if(length(x) != length(y)) {
    stop("Input x and y should be vectors of the same length", call.=FALSE)
  }
  sum( is.na(x) & is.na(y) )

```r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
both_na2(x, y1)
```

```
## [1] 2
```


```r
both_na3 <- function(x, y) {
  ## Print some info on where NA's are as well as the number of them 
  if(length(x) != length(y)) {
    stop("Input x and y should be vectors of the same length", call.=FALSE)
  }
  na.in.both <- ( is.na(x) & is.na(y) )
  na.number  <- sum(na.in.both)
  na.which   <- which(na.in.both)

  message("Found ", na.number, " NA's at position(s):", 
          paste(na.which, collapse=", ") ) 
  
  return( list(number=na.number, which=na.which) )
}
both_na3(x,y)
```
#Section 2

```r
x <- df1$IDs
y <- df2$IDs
x
```

```
## [1] "gene1" "gene2" "gene3"
```

```r
y
```

```
## [1] "gene2" "gene4" "gene3" "gene5"
```
#intersect() annd %in%

```r
intersect(x,y)
```

```
## [1] "gene2" "gene3"
```

```r
x%in%y
```

```
## [1] FALSE  TRUE  TRUE
```


```r
# Putting together
cbind( x[ x %in% y ], y[ y %in% x ] )
```

```
##      [,1]    [,2]   
## [1,] "gene2" "gene2"
## [2,] "gene3" "gene3"
```

```r
gene_intersect <- function(x, y) { 
   cbind( x[ x %in% y ], y[ y %in% x ] )
}
gene_intersect(x,y)
```

```
##      [,1]    [,2]   
## [1,] "gene2" "gene2"
## [2,] "gene3" "gene3"
```
gene_intersect2 <- function(df1, df2) { 
   cbind( df1[ df1$IDs %in% df2$IDs, ], 
          df2[ df2$IDs %in% df1$IDs, "exp"] )
}

```r
gene_intersect2(df1,df2)
```

```
##     IDs exp df2[df2$IDs %in% df1$IDs, "exp"]
## 2 gene2   1                               -2
## 3 gene3   1                                1
```
gene_intersect3 <- function(df1, df2, gene.colname="IDs") { 
   cbind( df1[ df1[,gene.colname] %in% df2[,gene.colname], ], 
          exp2=df2[ df2[,gene.colname] %in% df1[,gene.colname], "exp"] )
}

```r
gene_intersect3(df1,df2)
```

```
##     IDs exp exp2
## 2 gene2   1   -2
## 3 gene3   1    1
```
gene_intersect4 <- function(df1, df2, gene.colname="IDs") { 

  df1.name <- df1[,gene.colname]
  df2.name <- df2[,gene.colname]

  df1.inds <- df1.name %in% df2.name
  df2.inds <- df2.name %in% df1.name

   cbind( df1[ df1.inds, ], 
          exp2=df2[ df2.inds, "exp"] )
}

```r
gene_intersect4(df1,df2)
```

```
##     IDs exp exp2
## 2 gene2   1   -2
## 3 gene3   1    1
```


```r
merge(df1, df2, by="IDs")
```

```
##     IDs exp.x exp.y
## 1 gene2     1    -2
## 2 gene3     1     1
```
