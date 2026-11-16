# Two-way table from PLSrounding output

Output from
[`PLSrounding`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/PLSrounding.md)
is presented as two-way table(s) in cases where this is possible. A
requirement is that the number of main dimensional variables is two.

## Usage

``` r
PLS2way(obj, variable = c("rounded", "original", "difference", "code"))
```

## Arguments

- obj:

  Output object from
  [`PLSrounding`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/PLSrounding.md)

- variable:

  One of `"rounded"` (default), `"original"`, `"difference"` or
  `"code"`.

## Value

A data frame

## Details

When parameter `"variable"` is `"code"`, output is coded as `"#"`
(publish), `"."` (inner) and `"&"` (both).

## Examples

``` r
# Making tables from PLSrounding examples 
z <- SmallCountData("e6")
a <- PLSrounding(z, "freq", formula = ~eu * year + geo)
PLS2way(a, "original")
#>          X2018 X2019 Total
#> Iceland      2     1     3
#> Portugal     3     5     8
#> Spain        7     6    13
#> EU          10    11    21
#> Total       12    12    24
#> nonEU        2     1     3
PLS2way(a, "difference")
#>          X2018 X2019 Total
#> Iceland      1    -1     0
#> Portugal     0     0     0
#> Spain        0     0     0
#> EU           0     0     0
#> Total        1    -1     0
#> nonEU        1    -1     0
PLS2way(a, "code")
#>          X2018 X2019 Total
#> Iceland      .     .     #
#> Portugal     .     .     #
#> Spain        .     .     #
#> EU           #     #     #
#> Total        #     #     #
#> nonEU        #     #     #
PLS2way(PLSrounding(z, "freq", formula = ~eu * year + geo * year), "code")
#>          X2018 X2019 Total
#> Iceland      &     &     #
#> Portugal     &     &     #
#> Spain        &     &     #
#> EU           #     #     #
#> Total        #     #     #
#> nonEU        #     #     #
eHrc2 <- list(geo = c("EU", "@Portugal", "@Spain", "Iceland"), year = c("2018", "2019"))
PLS2way(PLSrounding(z, "freq", hierarchies = eHrc2))
#>          X2018 X2019 Total
#> Iceland      3     0     3
#> Portugal     3     5     8
#> Spain        7     6    13
#> EU          10    11    21
#> Total       13    11    24
```
