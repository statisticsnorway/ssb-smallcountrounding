# PLS inspired rounding

Small count rounding of necessary inner cells are performed so that all
small frequencies of cross-classifications to be published (publishable
cells) are rounded. The publishable cells can be defined from a model
formula, hierarchies or automatically from data.

## Usage

``` r
PLSrounding(
  data,
  freqVar = NULL,
  roundBase = 3,
  hierarchies = NULL,
  formula = NULL,
  dimVar = NULL,
  maxRound = roundBase - 1,
  printInc = nrow(data) > 1000,
  output = NULL,
  extend0 = FALSE,
  preAggregate = NA,
  zeroCandidates = FALSE,
  forceInner = FALSE,
  allSmall = FALSE,
  step = 0,
  aggregatePackage = "base",
  aggregateNA = TRUE,
  aggregateBaseOrder = FALSE,
  rowGroupsPackage = aggregatePackage,
  ...,
  action_unused_dots = getOption("SmallCountRounding.action_unused_dots", "warn"),
  allowed_unused_dots = getOption("SmallCountRounding.allowed_unused_dots", character(0))
)

PLSroundingInner(..., output = "inner")

PLSroundingPublish(..., output = "publish")
```

## Arguments

- data:

  Input data (inner cells), typically a data frame, tibble, or
  data.table. If `data` is not a classic data frame, it will be coerced
  to one internally unless `preAggregate` is `TRUE` and
  `aggregatePackage` is `"data.table"`.

- freqVar:

  Variable holding counts (inner cells frequencies). When `NULL`
  (default), microdata is assumed.

- roundBase:

  Rounding base

- hierarchies:

  List of hierarchies

- formula:

  Model formula defining publishable cells

- dimVar:

  The main dimensional variables and additional aggregating variables.
  This parameter can be useful when hierarchies and formula are
  unspecified.

- maxRound:

  Inner cells contributing to original publishable cells equal to or
  less than maxRound will be rounded

- printInc:

  Printing iteration information to console when TRUE

- output:

  Possible non-NULL values are `"input"`, `"inner"` and `"publish"`.
  Then a single data frame is returned.

- extend0:

  When `extend0` is set to `TRUE`, the data is automatically extended.
  This is relevant when `zeroCandidates = TRUE` (see
  [`RoundViaDummy`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/RoundViaDummy.md)).
  Additionally, `extend0` can be specified as a list, representing the
  `varGroups` parameter in the
  [`Extend0`](https://statisticsnorway.github.io/ssb-ssbtools/reference/Extend0.html)
  function. Can also be set to `"all"` which means that input codes in
  hierarchies are considered in addition to those in data.

- preAggregate:

  When `TRUE`, the data will be aggregated within the function to an
  appropriate level. This is defined by the dimensional variables
  according to `dimVar`, `hierarchies` or `formula`. When `FALSE`, no
  aggregation is performed. When `NA` (default), the function will
  automatically decide whether to aggregate: aggregation is applied
  unless `freqVar` is present and the data contain no duplicated rows
  with respect to the dimensional variables. Exception: if a non-`NULL`
  `x` (the model matrix) is supplied via `...`, `NA` is treated as
  `FALSE`.

- zeroCandidates:

  When `TRUE`, inner cells with zero frequencies are treated as
  candidates to be rounded up. See
  [`RoundViaDummy()`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/RoundViaDummy.md)
  for details. Note: When using microdata input (`freqVar = NULL`),
  `extend0` must be specified for `zeroCandidates` to have any effect.

- forceInner:

  When `TRUE`, all inner cells are forced to be rounded. See
  [`RoundViaDummy()`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/RoundViaDummy.md)
  for details.

- allSmall:

  When `TRUE`, all small inner cells (`<= maxRound`) are rounded. See
  [`RoundViaDummy()`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/RoundViaDummy.md)
  for details.

- step:

  A parameter controlling forward and backward steps within the
  algorithm. The `step` parameter can significantly impact performance
  on large datasets. For example, using `step = list(100)` may be a
  useful approach. See the parameter description and details in the
  underlying function
  [`RoundViaDummy()`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/RoundViaDummy.md).

- aggregatePackage:

  Package used to preAggregate. Parameter `pkg` to
  [`aggregate_by_pkg`](https://statisticsnorway.github.io/ssb-ssbtools/reference/aggregate_by_pkg.html).

- aggregateNA:

  Whether to include NAs in the grouping variables while preAggregate.
  Parameter `include_na` to
  [`aggregate_by_pkg`](https://statisticsnorway.github.io/ssb-ssbtools/reference/aggregate_by_pkg.html).

- aggregateBaseOrder:

  Parameter `base_order` to
  [`aggregate_by_pkg`](https://statisticsnorway.github.io/ssb-ssbtools/reference/aggregate_by_pkg.html),
  used when preAggregate. The default is set to `FALSE` to avoid
  unnecessary sorting operations. When `TRUE`, an attempt is made to
  return the same result with `data.table` as with base R. This cannot
  be guaranteed due to potential variations in sorting behavior across
  different systems.

- rowGroupsPackage:

  Parameter `pkg` to
  [`RowGroups`](https://statisticsnorway.github.io/ssb-ssbtools/reference/RowGroups.html).
  The parameter is input to
  [`Formula2ModelMatrix`](https://statisticsnorway.github.io/ssb-ssbtools/reference/FormulaSums.html)
  via
  [`ModelMatrix`](https://statisticsnorway.github.io/ssb-ssbtools/reference/ModelMatrix.html).

- ...:

  Further parameters sent to `RoundViaDummy`

- action_unused_dots:

  Character string controlling how unused arguments in `...` are
  handled. Internally uses
  [`ellipsis::check_dots_used()`](https://rlang.r-lib.org/reference/check_dots_used.html)
  with a custom action. One of "warn", "abort", "inform", or "none". The
  value "none" disables the check entirely. The default is taken from
  `getOption("SmallCountRounding.action_unused_dots")`, falling back to
  "warn" if the option is not set. Users can change the default globally
  with e.g. `options(SmallCountRounding.action_unused_dots = "abort")`.

- allowed_unused_dots:

  Character vector of argument names ignored by the unused-argument
  check. May be useful when this function is wrapped by another
  function, or in other cases where a correctly spelled argument is
  nevertheless not registered as used. The default is taken from
  `getOption("SmallCountRounding.allowed_unused_dots")`, falling back to
  `character(0)` if the option is not set. Users can change the default
  globally with e.g.
  `options(SmallCountRounding.allowed_unused_dots = c("plotColor", "lineType"))`.

## Value

Output is a four-element list with class attribute "PLSrounded", which
ensures informative printing and enables the use of
[`FormulaSelection`](https://statisticsnorway.github.io/ssb-ssbtools/reference/FormulaSelection.html)
on this object.

- inner:

  Data frame corresponding to input data with the main dimensional
  variables and with cell frequencies (original, rounded, difference).

- publish:

  Data frame of publishable data with the main dimensional variables and
  with cell frequencies (original, rounded, difference).

- metrics:

  A named character vector of various statistics calculated from the two
  output data frames ("`inner_`" used to distinguish). See examples
  below and the function
  [`HDutility`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/HDutility.md).

- freqTable:

  Matrix of frequencies of cell frequencies and absolute differences.
  For example, row "`rounded`" and column "`inn.4+`" is the number of
  rounded inner cell frequencies greater than or equal to `4`.

## Details

This function is a user-friendly wrapper for `RoundViaDummy` with data
frame output and with computed summary of the results. See
[`RoundViaDummy`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/RoundViaDummy.md)
for more details.

## References

Langsrud, Ø. and Heldal, J. (2018): “An Algorithm for Small Count
Rounding of Tabular Data”. Presented at: *Privacy in statistical
databases*, Valencia, Spain. September 26-28, 2018.
<https://www.researchgate.net/publication/327768398_An_Algorithm_for_Small_Count_Rounding_of_Tabular_Data>

## See also

[`RoundViaDummy`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/RoundViaDummy.md),
[`PLS2way`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/PLS2way.md),
[`ModelMatrix`](https://statisticsnorway.github.io/ssb-ssbtools/reference/ModelMatrix.html)

## Examples

``` r
# Small example data set
z <- SmallCountData("e6")
print(z)
#>        geo    eu year freq
#> 1  Iceland nonEU 2018    2
#> 2 Portugal    EU 2018    3
#> 3    Spain    EU 2018    7
#> 4  Iceland nonEU 2019    1
#> 5 Portugal    EU 2019    5
#> 6    Spain    EU 2019    6

# Publishable cells by formula interface
a <- PLSrounding(z, "freq", roundBase = 5,  formula = ~geo + eu + year)
print(a)
#> 
#> PLSrounding summary:  
#> 
#>        maxdiff      HDutility    meanAbsDiff rootMeanSquare 
#>              3          0.938           1.25         1.6583 
#> 
#> Frequencies of cell frequencies and absolute differences:  
#> 
#>          inn.0 inn.1-4 inn.5 inn.6+ inn.all pub.0 pub.1-4 pub.5 pub.6+ pub.all
#> original     .       3     1      2       6     .       2     .      6       8
#> rounded      1       1     2      2       6     .       .     2      6       8
#> absDiff      4       2     .      .       6     3       5     .      .       8
#> 
print(a$inner)
#>        geo year original rounded difference
#> 1  Iceland 2018        2       5          3
#> 2 Portugal 2018        3       3          0
#> 3    Spain 2018        7       7          0
#> 4  Iceland 2019        1       0         -1
#> 5 Portugal 2019        5       5          0
#> 6    Spain 2019        6       6          0
print(a$publish)
#>        geo  year original rounded difference
#> 1    Total Total       24      26          2
#> 2  Iceland Total        3       5          2
#> 3 Portugal Total        8       8          0
#> 4    Spain Total       13      13          0
#> 5       EU Total       21      21          0
#> 6    nonEU Total        3       5          2
#> 7    Total  2018       12      15          3
#> 8    Total  2019       12      11         -1
print(a$metrics)
#>            roundBase             maxRound              maxdiff 
#>            5.0000000            4.0000000            3.0000000 
#>      inner_HDutility            HDutility    inner_meanAbsDiff 
#>            0.8131709            0.9380433            0.6666667 
#>          meanAbsDiff inner_rootMeanSquare       rootMeanSquare 
#>            1.2500000            1.2909944            1.6583124 
print(a$freqTable)
#>          inn.0 inn.1-4 inn.5 inn.6+ inn.all pub.0 pub.1-4 pub.5 pub.6+ pub.all
#> original     0       3     1      2       6     0       2     0      6       8
#> rounded      1       1     2      2       6     0       0     2      6       8
#> absDiff      4       2     0      0       6     3       5     0      0       8

# Using FormulaSelection()
FormulaSelection(a$publish, ~eu + year)
#>     geo  year original rounded difference
#> 1 Total Total       24      26          2
#> 5    EU Total       21      21          0
#> 6 nonEU Total        3       5          2
#> 7 Total  2018       12      15          3
#> 8 Total  2019       12      11         -1
FormulaSelection(a, ~eu + year) # same as above
#>     geo  year original rounded difference
#> 1 Total Total       24      26          2
#> 5    EU Total       21      21          0
#> 6 nonEU Total        3       5          2
#> 7 Total  2018       12      15          3
#> 8 Total  2019       12      11         -1
FormulaSelection(a)             # just a$publish
#>        geo  year original rounded difference
#> 1    Total Total       24      26          2
#> 2  Iceland Total        3       5          2
#> 3 Portugal Total        8       8          0
#> 4    Spain Total       13      13          0
#> 5       EU Total       21      21          0
#> 6    nonEU Total        3       5          2
#> 7    Total  2018       12      15          3
#> 8    Total  2019       12      11         -1

# Recalculation of maxdiff, HDutility, meanAbsDiff and rootMeanSquare
max(abs(a$publish[, "difference"]))
#> [1] 3
HDutility(a$publish[, "original"], a$publish[, "rounded"])
#> [1] 0.9380433
mean(abs(a$publish[, "difference"]))
#> [1] 1.25
sqrt(mean((a$publish[, "difference"])^2))
#> [1] 1.658312

# Five lines below produce equivalent results 
# Ordering of rows can be different
PLSrounding(z, "freq", dimVar = c("geo", "eu", "year"))
#> 
#> PLSrounding summary:  
#> 
#>        maxdiff      HDutility    meanAbsDiff rootMeanSquare 
#>              1         0.9117         0.3333         0.5774 
#> 
#> Frequencies of cell frequencies and absolute differences:  
#> 
#>          inn.0 inn.1-2 inn.3 inn.4+ inn.all pub.0 pub.1-2 pub.3 pub.4+ pub.all
#> original     .       2     1      3       6     .       4     3     11      18
#> rounded      1       .     2      3       6     2       .     5     11      18
#> absDiff      4       2     .      .       6    12       6     .      .      18
#> 
PLSrounding(z, "freq", formula = ~eu * year + geo * year)
#> 
#> PLSrounding summary:  
#> 
#>        maxdiff      HDutility    meanAbsDiff rootMeanSquare 
#>              1         0.9117         0.3333         0.5774 
#> 
#> Frequencies of cell frequencies and absolute differences:  
#> 
#>          inn.0 inn.1-2 inn.3 inn.4+ inn.all pub.0 pub.1-2 pub.3 pub.4+ pub.all
#> original     .       2     1      3       6     .       4     3     11      18
#> rounded      1       .     2      3       6     2       .     5     11      18
#> absDiff      4       2     .      .       6    12       6     .      .      18
#> 
PLSrounding(z[, -2], "freq", hierarchies = SmallCountData("eHrc"))
#> 
#> PLSrounding summary:  
#> 
#>        maxdiff      HDutility    meanAbsDiff rootMeanSquare 
#>              1         0.9117         0.3333         0.5774 
#> 
#> Frequencies of cell frequencies and absolute differences:  
#> 
#>          inn.0 inn.1-2 inn.3 inn.4+ inn.all pub.0 pub.1-2 pub.3 pub.4+ pub.all
#> original     .       2     1      3       6     .       4     3     11      18
#> rounded      1       .     2      3       6     2       .     5     11      18
#> absDiff      4       2     .      .       6    12       6     .      .      18
#> 
PLSrounding(z[, -2], "freq", hierarchies = SmallCountData("eDimList"))
#> 
#> PLSrounding summary:  
#> 
#>        maxdiff      HDutility    meanAbsDiff rootMeanSquare 
#>              1         0.9117         0.3333         0.5774 
#> 
#> Frequencies of cell frequencies and absolute differences:  
#> 
#>          inn.0 inn.1-2 inn.3 inn.4+ inn.all pub.0 pub.1-2 pub.3 pub.4+ pub.all
#> original     .       2     1      3       6     .       4     3     11      18
#> rounded      1       .     2      3       6     2       .     5     11      18
#> absDiff      4       2     .      .       6    12       6     .      .      18
#> 
PLSrounding(z[, -2], "freq", hierarchies = SmallCountData("eDimList"), formula = ~geo * year)
#> 
#> PLSrounding summary:  
#> 
#>        maxdiff      HDutility    meanAbsDiff rootMeanSquare 
#>              1         0.9117         0.3333         0.5774 
#> 
#> Frequencies of cell frequencies and absolute differences:  
#> 
#>          inn.0 inn.1-2 inn.3 inn.4+ inn.all pub.0 pub.1-2 pub.3 pub.4+ pub.all
#> original     .       2     1      3       6     .       4     3     11      18
#> rounded      1       .     2      3       6     2       .     5     11      18
#> absDiff      4       2     .      .       6    12       6     .      .      18
#> 

# Define publishable cells differently by making use of formula interface
PLSrounding(z, "freq", formula = ~eu * year + geo)
#> 
#> PLSrounding summary:  
#> 
#>        maxdiff      HDutility    meanAbsDiff rootMeanSquare 
#>              1          0.931         0.3333         0.5774 
#> 
#> Frequencies of cell frequencies and absolute differences:  
#> 
#>          inn.0 inn.1-2 inn.3 inn.4+ inn.all pub.0 pub.1-2 pub.3 pub.4+ pub.all
#> original     .       2     1      3       6     .       2     2      8      12
#> rounded      1       .     2      3       6     1       .     3      8      12
#> absDiff      4       2     .      .       6     8       4     .      .      12
#> 

# Define publishable cells differently by making use of hierarchy interface
eHrc2 <- list(geo = c("EU", "@Portugal", "@Spain", "Iceland"), year = c("2018", "2019"))
PLSrounding(z, "freq", hierarchies = eHrc2)
#> 
#> PLSrounding summary:  
#> 
#>        maxdiff      HDutility    meanAbsDiff rootMeanSquare 
#>              1         0.9357         0.2667         0.5164 
#> 
#> Frequencies of cell frequencies and absolute differences:  
#> 
#>          inn.0 inn.1-2 inn.3 inn.4+ inn.all pub.0 pub.1-2 pub.3 pub.4+ pub.all
#> original     .       2     1      3       6     .       2     2     11      15
#> rounded      1       .     2      3       6     1       .     3     11      15
#> absDiff      4       2     .      .       6    11       4     .      .      15
#> 

# Also possible to combine hierarchies and formula
PLSrounding(z, "freq", hierarchies = SmallCountData("eDimList"), formula = ~geo + year)
#> 
#> PLSrounding summary:  
#> 
#>        maxdiff      HDutility    meanAbsDiff rootMeanSquare 
#>              0              1              0              0 
#> 
#> Frequencies of cell frequencies and absolute differences:  
#> 
#>          inn.0 inn.1-2 inn.3 inn.4+ inn.all pub.0 pub.1-2 pub.3 pub.4+ pub.all
#> original     .       2     1      3       6     .       .     2      7       9
#> rounded      .       2     1      3       6     .       .     2      7       9
#> absDiff      6       .     .      .       6     9       .     .      .       9
#> 

# Single data frame output and total codes 
PLSroundingInner(z, "freq", roundBase = 5, formula = ~geo + eu + year)
#>        geo year original rounded difference
#> 1  Iceland 2018        2       5          3
#> 2 Portugal 2018        3       3          0
#> 3    Spain 2018        7       7          0
#> 4  Iceland 2019        1       0         -1
#> 5 Portugal 2019        5       5          0
#> 6    Spain 2019        6       6          0
PLSroundingPublish(z, roundBase = 5, formula = ~geo + eu + year)
#>        geo  year original rounded difference
#> 1    Total Total        6       5         -1
#> 2  Iceland Total        2       0         -2
#> 3 Portugal Total        2       0         -2
#> 4    Spain Total        2       5          3
#> 5       EU Total        4       5          1
#> 6    nonEU Total        2       0         -2
#> 7    Total  2018        3       5          2
#> 8    Total  2019        3       0         -3
PLSroundingPublish(z, roundBase = 5, formula = ~geo + eu + year,
                   total = c(geo = "Europe", year = "All_years"))
#>        geo      year original rounded difference
#> 1   Europe All_years        6       5         -1
#> 2  Iceland All_years        2       0         -2
#> 3 Portugal All_years        2       0         -2
#> 4    Spain All_years        2       5          3
#> 5       EU All_years        4       5          1
#> 6    nonEU All_years        2       0         -2
#> 7   Europe      2018        3       5          2
#> 8   Europe      2019        3       0         -3

# Microdata input
PLSroundingInner(rbind(z, z), roundBase = 5, formula = ~geo + eu + year)
#>        geo year original rounded difference
#> 1 Portugal 2018        2       0         -2
#> 2    Spain 2018        2       5          3
#> 3  Iceland 2018        2       0         -2
#> 4 Portugal 2019        2       0         -2
#> 5    Spain 2019        2       0         -2
#> 6  Iceland 2019        2       5          3

# Zero perturbed due to both  extend0 = TRUE and zeroCandidates = TRUE 
set.seed(12345)
PLSroundingInner(z[sample.int(5, 12, replace = TRUE), 1:3], 
                 formula = ~geo + eu + year, roundBase = 5, 
                 extend0 = TRUE, zeroCandidates = TRUE, printInc = TRUE)
#> [preAggregate 12*3->5*4]
#> [extend0 5*4->6*4]
#> [-**..:=]
#>        geo year original rounded difference
#> 1 Portugal 2018        4       0         -4
#> 2    Spain 2018        3       0         -3
#> 3  Iceland 2018        2       5          3
#> 4 Portugal 2019        1       0         -1
#> 5  Iceland 2019        2       0         -2
#> 6    Spain 2019        0       5          5

# Parameter avoidHierarchical (see RoundViaDummy and ModelMatrix) 
PLSroundingPublish(z, roundBase = 5, formula = ~geo + eu + year, avoidHierarchical = TRUE)
#>        geo    eu  year original rounded difference
#> 1    Total Total Total        6       5         -1
#> 2  Iceland Total Total        2       0         -2
#> 3 Portugal Total Total        2       0         -2
#> 4    Spain Total Total        2       5          3
#> 5    Total    EU Total        4       5          1
#> 6    Total nonEU Total        2       0         -2
#> 7    Total Total  2018        3       5          2
#> 8    Total Total  2019        3       0         -3


# To illustrate hierarchical_extend0 
#    (parameter to underlying function, SSBtools::Extend0fromModelMatrixInput)
PLSroundingInner(z[-c(2:3), ], roundBase = 5, formula = ~geo + eu + year, 
   avoidHierarchical = TRUE, zeroCandidates = TRUE, extend0 = TRUE)
#>         geo    eu year original rounded difference
#> 1   Iceland nonEU 2018        1       0         -1
#> 2  Portugal    EU 2019        1       0         -1
#> 3     Spain    EU 2019        1       0         -1
#> 4   Iceland nonEU 2019        1       0         -1
#> 5  Portugal nonEU 2018        0       0          0
#> 6     Spain nonEU 2018        0       0          0
#> 7   Iceland    EU 2018        0       0          0
#> 8  Portugal    EU 2018        0       0          0
#> 9     Spain    EU 2018        0       0          0
#> 10 Portugal nonEU 2019        0       0          0
#> 11    Spain nonEU 2019        0       0          0
#> 12  Iceland    EU 2019        0       5          5
PLSroundingInner(z[-c(2:3), ], roundBase = 5, formula = ~geo + eu + year, 
   avoidHierarchical = TRUE, zeroCandidates = TRUE, extend0 = TRUE, 
   hierarchical_extend0 = TRUE)
#>        geo    eu year original rounded difference
#> 1  Iceland nonEU 2018        1       0         -1
#> 2 Portugal    EU 2019        1       0         -1
#> 3    Spain    EU 2019        1       0         -1
#> 4  Iceland nonEU 2019        1       5          4
#> 5 Portugal    EU 2018        0       0          0
#> 6    Spain    EU 2018        0       0          0

# Package sdcHierarchies can be used to create hierarchies. 
# The small example code below works if this package is available. 
if (require(sdcHierarchies)) {
  z2 <- cbind(geo = c("11", "21", "22"), z[, 3:4], stringsAsFactors = FALSE)
  h2 <- list(
    geo = hier_compute(inp = unique(z2$geo), dim_spec = c(1, 1), root = "Tot", as = "df"),
    year = hier_convert(hier_create(root = "Total", nodes = c("2018", "2019")), as = "df"))
  PLSrounding(z2, "freq", hierarchies = h2)
}
#> Loading required package: sdcHierarchies
#> Loading required package: shinythemes
#> Package 'sdcHierarchies' 0.23.2 has been loaded.
#> 
#> PLSrounding summary:  
#> 
#>        maxdiff      HDutility    meanAbsDiff rootMeanSquare 
#>              1         0.9117         0.3333         0.5774 
#> 
#> Frequencies of cell frequencies and absolute differences:  
#> 
#>          inn.0 inn.1-2 inn.3 inn.4+ inn.all pub.0 pub.1-2 pub.3 pub.4+ pub.all
#> original     .       2     1      3       6     .       4     3     11      18
#> rounded      1       .     2      3       6     2       .     5     11      18
#> absDiff      4       2     .      .       6    12       6     .      .      18
#> 

# Use PLS2way to produce tables as in Langsrud and Heldal (2018) and to demonstrate 
# parameters maxRound, zeroCandidates and identifyNew (see RoundViaDummy).   
# Parameter rndSeed used to ensure same output as in reference.
exPSD <- SmallCountData("exPSD")
a <- PLSrounding(exPSD, "freq", 5, formula = ~rows + cols, rndSeed=124)
PLS2way(a, "original")  # Table 1
#>       col1 col2 col3 col4 col5 Total
#> row1     6    0    1    3    4    14
#> row2     1    2    3    1    2     9
#> row3     0    1    1    0    2     4
#> Total    7    3    5    4    8    27
PLS2way(a)  # Table 2
#>       col1 col2 col3 col4 col5 Total
#> row1     6    0    5    0    4    15
#> row2     1    0    0    5    2     8
#> row3     0    5    0    0    0     5
#> Total    7    5    5    5    6    28
a <- PLSrounding(exPSD, "freq", 5, formula = ~rows + cols, identifyNew = FALSE, rndSeed=124)
PLS2way(a)  # Table 3
#>       col1 col2 col3 col4 col5 Total
#> row1     6    0    1    0    4    11
#> row2     1    0    3    5    2    11
#> row3     0    5    0    0    0     5
#> Total    7    5    4    5    6    27
a <- PLSrounding(exPSD, "freq", 5, formula = ~rows + cols, maxRound = 7)
PLS2way(a)  # Values in col1 rounded
#>       col1 col2 col3 col4 col5 Total
#> row1     5    0    0    5    0    10
#> row2     0    0    5    0    5    10
#> row3     0    5    0    0    0     5
#> Total    5    5    5    5    5    25
a <- PLSrounding(exPSD, "freq", 5, formula = ~rows + cols, zeroCandidates = TRUE)
PLS2way(a)  # (row3, col4): original is 0 and rounded is 5
#>       col1 col2 col3 col4 col5 Total
#> row1     6    0    5    0    4    15
#> row2     1    5    0    0    2     8
#> row3     0    0    0    5    0     5
#> Total    7    5    5    5    6    28

# Using formula followed by FormulaSelection 
output <- PLSrounding(data = SmallCountData("example1"), 
                      formula = ~age * geo * year + eu * year, 
                      freqVar = "freq", 
                      roundBase = 5)
FormulaSelection(output, ~(age + eu) * year)
#>      age   geo  year original rounded difference
#> 1  Total Total Total       59      59          0
#> 2    old Total Total       38      37         -1
#> 3  young Total Total       21      22          1
#> 7  Total Total  2014       20      21          1
#> 8  Total Total  2015       18      16         -2
#> 9  Total Total  2016       21      22          1
#> 10 Total    EU Total       46      49          3
#> 11 Total nonEU Total       13      10         -3
#> 18   old Total  2014       13      16          3
#> 19   old Total  2015       13      11         -2
#> 20   old Total  2016       12      10         -2
#> 21 young Total  2014        7       5         -2
#> 22 young Total  2015        5       5          0
#> 23 young Total  2016        9      12          3
#> 33 Total    EU  2014       15      16          1
#> 34 Total nonEU  2014        5       5          0
#> 35 Total    EU  2015       15      16          1
#> 36 Total nonEU  2015        3       0         -3
#> 37 Total    EU  2016       16      17          1
#> 38 Total nonEU  2016        5       5          0

# Example similar to the one in the documentation of tables_by_formulas,
# but using PLSroundingPublish with roundBase = 4.
tables_by_formulas(SSBtoolsData("magnitude1"),
                   table_fun = PLSroundingPublish, 
                   table_formulas = list(table_1 = ~region * sector2, 
                                         table_2 = ~region1:sector4 - 1, 
                                         table_3 = ~region + sector4 - 1), 
                   substitute_vars = list(region = c("geo", "eu"), region1 = "eu"), 
                   collapse_vars = list(sector = c("sector2", "sector4")), 
                   roundBase = 4) 
#>      region        sector original rounded difference table_1 table_2 table_3
#> 1     Total         Total       20      21          1    TRUE   FALSE   FALSE
#> 2   Iceland         Total        4       4          0    TRUE   FALSE    TRUE
#> 3  Portugal         Total        8       8          0    TRUE   FALSE    TRUE
#> 4     Spain         Total        8       9          1    TRUE   FALSE    TRUE
#> 5        EU         Total       16      17          1    TRUE   FALSE    TRUE
#> 6     nonEU         Total        4       4          0    TRUE   FALSE    TRUE
#> 7     Total       private       16      17          1    TRUE   FALSE   FALSE
#> 8     Total        public        4       4          0    TRUE   FALSE   FALSE
#> 9     Total   Agriculture        4       4          0   FALSE   FALSE    TRUE
#> 10    Total Entertainment        6       5         -1   FALSE   FALSE    TRUE
#> 11    Total  Governmental        4       4          0   FALSE   FALSE    TRUE
#> 12    Total      Industry        6       8          2   FALSE   FALSE    TRUE
#> 13  Iceland       private        4       4          0    TRUE   FALSE   FALSE
#> 14 Portugal       private        6       4         -2    TRUE   FALSE   FALSE
#> 15 Portugal        public        2       4          2    TRUE   FALSE   FALSE
#> 16    Spain       private        6       9          3    TRUE   FALSE   FALSE
#> 17    Spain        public        2       0         -2    TRUE   FALSE   FALSE
#> 18       EU       private       12      13          1    TRUE   FALSE   FALSE
#> 19       EU        public        4       4          0    TRUE   FALSE   FALSE
#> 20    nonEU       private        4       4          0    TRUE   FALSE   FALSE
#> 21       EU   Agriculture        4       4          0   FALSE    TRUE   FALSE
#> 22       EU Entertainment        5       5          0   FALSE    TRUE   FALSE
#> 23       EU  Governmental        4       4          0   FALSE    TRUE   FALSE
#> 24       EU      Industry        3       4          1   FALSE    TRUE   FALSE
#> 25    nonEU Entertainment        1       0         -1   FALSE    TRUE   FALSE
#> 26    nonEU      Industry        3       4          1   FALSE    TRUE   FALSE
```
