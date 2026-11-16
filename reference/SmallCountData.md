# Function that returns a dataset

Function that returns a dataset

## Usage

``` r
SmallCountData(dataset, path = NULL)
```

## Arguments

- dataset:

  Name of data set within the SmallCountRounding package

- path:

  When non-NULL the data set is read from "path/dataset.RData"

## Value

The dataset

## Note

Except for `"europe6"`, `"eHrc"`, `"eDimList"` and `"exPSD"`, the
function returns the same datasets as
[`SSBtoolsData`](https://statisticsnorway.github.io/ssb-ssbtools/reference/SSBtoolsData.html).

## See also

[`SSBtoolsData`](https://statisticsnorway.github.io/ssb-ssbtools/reference/SSBtoolsData.html),
[`Hrc2DimList`](https://statisticsnorway.github.io/ssb-ssbtools/reference/DimList2Hrc.html)

## Examples

``` r
 SmallCountData("z1")
#>    region hovedint ant
#> 1       A    annet  11
#> 2       B    annet   7
#> 3       C    annet   5
#> 4       D    annet  13
#> 5       E    annet   9
#> 6       F    annet  12
#> 7       G    annet   6
#> 8       H    annet   9
#> 9       A   arbeid  11
#> 10      B   arbeid   1
#> 11      C   arbeid   8
#> 12      D   arbeid   2
#> 13      E   arbeid  14
#> 14      F   arbeid   9
#> 15      G   arbeid   4
#> 16      H   arbeid   3
#> 17      A soshjelp  55
#> 18      B soshjelp  29
#> 19      C soshjelp  35
#> 20      D soshjelp  17
#> 21      E soshjelp  63
#> 22      F soshjelp  24
#> 23      G soshjelp  22
#> 24      H soshjelp  38
#> 25      A    trygd  36
#> 26      B    trygd  18
#> 27      C    trygd  25
#> 28      D    trygd  13
#> 29      E    trygd  52
#> 30      F    trygd  22
#> 31      G    trygd   8
#> 32      H    trygd  15
 SmallCountData("e6")
#>        geo    eu year freq
#> 1  Iceland nonEU 2018    2
#> 2 Portugal    EU 2018    3
#> 3    Spain    EU 2018    7
#> 4  Iceland nonEU 2019    1
#> 5 Portugal    EU 2019    5
#> 6    Spain    EU 2019    6
 SmallCountData("eHrc")      #  TauArgus coded hierarchies
#> $geo
#> [1] "EU"        "@Portugal" "@Spain"    "nonEU"     "@Iceland" 
#> 
#> $year
#> [1] "2018" "2019"
#> 
 SmallCountData("eDimList")  #  sdcTable coded hierarchies
#> $geo
#>   levels    codes
#> 1      @    Total
#> 2     @@       EU
#> 3    @@@ Portugal
#> 4    @@@    Spain
#> 5     @@    nonEU
#> 6    @@@  Iceland
#> 
#> $year
#>   levels codes
#> 1      @ Total
#> 2     @@  2018
#> 3     @@  2019
#> 
 SmallCountData("exPSD")     #  Example data in presentation at Privacy in statistical databases
#>    rows cols freq
#> 1  row1 col1    6
#> 2  row2 col1    1
#> 3  row3 col1    0
#> 4  row1 col2    0
#> 5  row2 col2    2
#> 6  row3 col2    1
#> 7  row1 col3    1
#> 8  row2 col3    3
#> 9  row3 col3    1
#> 10 row1 col4    3
#> 11 row2 col4    1
#> 12 row3 col4    0
#> 13 row1 col5    4
#> 14 row2 col5    2
#> 15 row3 col5    2
```
