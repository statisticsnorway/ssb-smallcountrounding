# Small count rounding by various methods

Run
[`makeroundtabs()`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/makeroundtabs.md)
or
[`RoundViaDummy2()`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/RoundViaDummy2.md)

## Usage

``` r
Round2(data, control, ..., method = c("roundtabs", "viadummy", "viadummyAll"))
```

## Arguments

- data:

  Input data.frame. Same as input parameter `A` to
  [`makeroundtabs()`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/makeroundtabs.md)

- control:

  As input to
  [`makeroundtabs()`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/makeroundtabs.md)

- ...:

  Other paramameters as input to
  [`makeroundtabs()`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/makeroundtabs.md)

- method:

  One of "roundtabs" (default), "viadummy", "viadummyAll"
  (`allTerms=TRUE`)

## Value

Output from
[`makeroundtabs`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/makeroundtabs.md)
or
[`RoundViaDummy2`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/RoundViaDummy2.md)

## See also

`RoundKostra` (SSB internal),
[`RoundViaDummy`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/RoundViaDummy.md),
[`Lists2formula`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/Lists2formula.md),
[`MakeControl`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/MakeControl.md),
[`FindMaxDiff`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/FindMaxDiff.md)

## Examples

``` r
if (FALSE) { # \dontrun{
z <- SmallCountData("sosialFiktiv")
d <- list(c("region","mnd") , c("hovedint","mnd2") , c("fylke","hovedint","mnd2") , 
          c("kostragr","hovedint","mnd"))
con <- MakeControl(d,z)
sor <- names(z)[c(4,5,3,2,1)]

roundedA <-  Round2(data=z,b=3,d=d,micro=FALSE,sort=sor,control=con, nin="ant",nout="Rndtall",
                    minit=2,maxit=10,maxdiff=5,seed=123,method="roundtabs")
roundedB <-  Round2(data=z,b=3,d=d,micro=FALSE,sort=sor,control=con, nin="ant",nout="Rndtall",
                    minit=2,maxit=10,maxdiff=5,seed=123,method="viadummy")
#10 rows of rounded data
roundedA$Ar[1:10,]  #"roundtabs"
cbind(z,roundedB$yInner)[1:10,] #"viadummy"

# recalculate maxdiff  nMaxdiff
dA <- FindMaxDiff(roundedA$Ar,con,"ant","m")
dB <- FindMaxDiff(z,con,roundedB$yInner[,1],roundedB$yInner[,2])

# Formula from d and control
Lists2formula(d,con,z)

# Formula from another d
d2 <-list(sor)
Lists2formula(d2,con,z)
Lists2formula(d2,con) # Without knowing data
Lists2formula(d2,data=z) # Without control
Lists2formula(d2,data=z) # Without control and data
} # }
```
