# Small Count Rounding of Tabular Data

Small count rounding via a dummy matrix and by an algorithm inspired by
PLS

## Usage

``` r
RoundViaDummy(
  data,
  freqVar,
  formula = NULL,
  roundBase = 3,
  singleRandom = FALSE,
  crossTable = TRUE,
  total = "Total",
  maxIterRows = 1000,
  maxIter = 1e+07,
  x = NULL,
  hierarchies = NULL,
  xReturn = FALSE,
  maxRound = roundBase - 1,
  zeroCandidates = FALSE,
  forceInner = FALSE,
  identifyNew = TRUE,
  step = 0,
  preRounded = NULL,
  leverageCheck = FALSE,
  easyCheck = TRUE,
  printInc = TRUE,
  rndSeed = 123,
  dimVar = NULL,
  plsWeights = NULL,
  preDifference = NULL,
  allSmall = FALSE,
  ...
)
```

## Arguments

- data:

  Input data as a data frame (inner cells)

- freqVar:

  Variable holding counts (name or number)

- formula:

  Model formula defining publishable cells. Will be used to calculate
  `x` (via
  [`ModelMatrix`](https://statisticsnorway.github.io/ssb-ssbtools/reference/ModelMatrix.html)).
  When NULL, x must be supplied.

- roundBase:

  Rounding base

- singleRandom:

  Single random draw when TRUE (instead of algorithm)

- crossTable:

  When TRUE, cross table in output and caculations via FormulaSums()

- total:

  String(s) used to name totals, passed to
  [`SSBtools::ModelMatrix()`](https://statisticsnorway.github.io/ssb-ssbtools/reference/ModelMatrix.html).
  Named vectors or lists are supported.

- maxIterRows:

  See details

- maxIter:

  Maximum number of iterations

- x:

  Dummy matrix defining publishable cells

- hierarchies:

  List of hierarchies, which can be converted by
  [`AutoHierarchies`](https://statisticsnorway.github.io/ssb-ssbtools/reference/AutoHierarchies.html).
  Thus, a single string as hierarchy input is assumed to be a total
  code. Exceptions are `"rowFactor"` or `""`, which correspond to only
  using the categories in the data.

- xReturn:

  Dummy matrix in output when TRUE (as input parameter `x`)

- maxRound:

  Inner cells contributing to original publishable cells equal to or
  less than maxRound will be rounded.

- zeroCandidates:

  When TRUE, inner cells in input with zero count (and multiple of
  roundBase when maxRound is in use) contributing to publishable cells
  will be included as candidates to obtain roundBase value. With vector
  input, the rule is specified individually for each cell. This can be
  specified as a vector, a variable in data or a function generating it
  (see details).

- forceInner:

  When TRUE, all inner cells will be rounded. Use vector input to force
  individual cells to be rounded. This can be specified as a vector, a
  variable in data or a function generating it (see details). Can be
  combined with parameter zeroCandidates to allow zeros and roundBase
  multiples to be rounded up.

- identifyNew:

  When `TRUE`, new cells may be identified after initial rounding to
  ensure all rounded publishable cells equal to or less than `maxRound`
  to be `roundBase` multiples. Use `NA` for the a less conservative
  behavior (old behavior). Then it is ensured that no nonzero rounded
  publishable cells are smaller than `roundBase`. When `maxRound` is
  default, there is no difference between `TRUE` and `NA`.

- step:

  When `step>1`, the original forward part of the algorithm is replaced
  by a kind of stepwise. After `step` steps forward, backward steps may
  be performed. The `step` parameter is also used for backward-forward
  iteration at the end of the algorithm; `step` backward steps may be
  performed. For greater control, the `step` parameter can be specified
  as a vector. Additionally, it can be provided as a list to trigger a
  final re-run iteration. See details.

- preRounded:

  A vector or a variable in data that contains a mixture of missing
  values and predetermined values of rounded inner cells. Can also be
  specified as a function generating it (see details).

- leverageCheck:

  When TRUE, all inner cells that depends linearly on the published
  cells and with small frequencies (`<=maxRound`) will be rounded. The
  computation of leverages can be very time and memory consuming. The
  function
  [`Reduce0exact`](https://statisticsnorway.github.io/ssb-ssbtools/reference/Reduce0exact.html)
  is called. The default leverage limit is `0.999999`. Another limit can
  be sent as input instead of `TRUE`. Checking is performed before and
  after (since new zeros) rounding. Extra iterations are performed when
  needed.

- easyCheck:

  A light version of the above leverage checking. Checking is performed
  after rounding. Extra iterations are performed when needed.
  `Reduce0exact` is called with `reduceByLeverage=FALSE` and
  `reduceByColSums=TRUE`.

- printInc:

  Printing iteration information to console when TRUE

- rndSeed:

  If non-NULL, a random generator seed to be used locally within the
  function without affecting the random value stream in R.

- dimVar:

  The main dimensional variables and additional aggregating variables.
  This parameter can be useful when hierarchies and formula are
  unspecified.

- plsWeights:

  A vector of weights for each cell to be published or a function
  generating it (see details). For use in the algorithm criterion.

- preDifference:

  A data.frame with differences already obtained from rounding another
  subset of data. There must be columns that match `crossTable`.
  Differences must be in the last column.

- allSmall:

  When TRUE, all small inner cells (`<= maxRound`) are rounded. This
  parameter is a simplified alternative to specifying `forceInner` (see
  details).

- ...:

  Further parameters sent to
  [`ModelMatrix`](https://statisticsnorway.github.io/ssb-ssbtools/reference/ModelMatrix.html).
  In particular, one can specify `removeEmpty=TRUE` to omit empty
  combinations. The parameter `inputInOutput` can be used to specify
  whether to include codes from input. The parameter `avoidHierarchical`
  ([`Formula2ModelMatrix`](https://statisticsnorway.github.io/ssb-ssbtools/reference/FormulaSums.html))
  can be combined with formula input.

## Value

A list where the two first elements are two column matrices. The first
matrix consists of inner cells and the second of cells to be published.
In each matrix the first and the second column contains, respectively,
original and rounded values. By default the cross table is the third
element of the output list.

## Details

Small count rounding of necessary inner cells are performed so that all
small frequencies of cross-classifications to be published (publishable
cells) are rounded. This is equivalent to changing micro data since
frequencies of unique combinations are changed. Thus, additivity and
consistency are guaranteed. The matrix multiplication formula is:
`yPublish` `=` `t(x)` `%*%` `yInner`, where `x` is the dummy matrix.

Parameters `zeroCandidates`, `forceInner`, `preRounded` and `plsWeights`
can be specified as functions. The supplied functions take the following
arguments: `data`, `yPublish`, `yInner`, `crossTable`, `x`, `roundBase`,
`maxRound`, and `...`, where the first two are numeric vectors of
original counts. When `allSmall` is `TRUE`, `forceInner` is set to
`function(yInner, maxRound, ...)` `yInner <= maxRound`.

Details about the `step` parameter:

- `step` as a numeric vector is converted to three parameters by

  - `step1 <- step[1]`

  - `step2 <- ifelse(length(step)>=2, step[2], round(step/2))`

  - `step3 <- ifelse(length(step)>=3, step[3], step[1])`

  After `step1` steps forward, up to `step2` backward steps may be
  performed. At the end of the algorithm; up to `step3` backward steps
  may be executed repeatedly.

- `step` when provided as a list (of numeric vectors), is adjusted to a
  length of 3 using `rep_len(step, 3)`.

  - `step[[1]]` is used in the main iterations.

  - `step[[2]]`, when non-`NULL`, is used in a final re-run iteration.

  - `step[[3]]` is used in extra iterations caused by `easyCheck` or
    `leverageCheck`.

  Setting `step = list(0)` will result in standard behavior, with the
  exception that an extra re-run iteration is performed. The most
  detailed setting is achieved by setting `step` to a length-3 list
  where each element has length 3.

## Note

Iterations are needed since after initial rounding of identified cells,
new cells are identified. If cases of a high number of identified cells
the algorithm can be too memory consuming (unless singleRandom=TRUE). To
avoid problems, not more than maxIterRows cells are rounded in each
iteration. The iteration limit (maxIter) is by default set to be high
since a low number of maxIterRows may need a high number of iterations.

## See also

See the user-friendly wrapper
[`PLSrounding`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/PLSrounding.md)
and see `Round2` for rounding by other algorithm

## Examples

``` r
# See similar and related examples in PLSrounding documentation
RoundViaDummy(SmallCountData("e6"), "freq")
#> [-**.:=]
#> $yInner
#>      original rounded
#> [1,]        2       3
#> [2,]        3       3
#> [3,]        7       7
#> [4,]        1       0
#> [5,]        5       5
#> [6,]        6       6
#> 
#> $yPublish
#>                original rounded
#> Total:Total          24      24
#> Total:2018           12      13
#> Total:2019           12      11
#> EU:Total             21      21
#> EU:2018              10      10
#> EU:2019              11      11
#> nonEU:Total           3       3
#> nonEU:2018            2       3
#> nonEU:2019            1       0
#> Iceland:Total         3       3
#> Iceland:2018          2       3
#> Iceland:2019          1       0
#> Portugal:Total        8       8
#> Portugal:2018         3       3
#> Portugal:2019         5       5
#> Spain:Total          13      13
#> Spain:2018            7       7
#> Spain:2019            6       6
#> 
#> $crossTable
#>         geo  year
#> 1     Total Total
#> 2     Total  2018
#> 3     Total  2019
#> 4        EU Total
#> 5        EU  2018
#> 6        EU  2019
#> 7     nonEU Total
#> 8     nonEU  2018
#> 9     nonEU  2019
#> 10  Iceland Total
#> 11  Iceland  2018
#> 12  Iceland  2019
#> 13 Portugal Total
#> 14 Portugal  2018
#> 15 Portugal  2019
#> 16    Spain Total
#> 17    Spain  2018
#> 18    Spain  2019
#> 
RoundViaDummy(SmallCountData("e6"), "freq", formula = ~eu * year + geo)
#> [-**.:=]
#> $yInner
#>      original rounded
#> [1,]        2       3
#> [2,]        3       3
#> [3,]        7       7
#> [4,]        1       0
#> [5,]        5       5
#> [6,]        6       6
#> 
#> $yPublish
#>                original rounded
#> Total-Total          24      24
#> Total-EU             21      21
#> Total-nonEU           3       3
#> 2018-Total           12      13
#> 2019-Total           12      11
#> Total-Iceland         3       3
#> Total-Portugal        8       8
#> Total-Spain          13      13
#> 2018-EU              10      10
#> 2019-EU              11      11
#> 2018-nonEU            2       3
#> 2019-nonEU            1       0
#> 
#> $crossTable
#>     year      geo
#> 1  Total    Total
#> 2  Total       EU
#> 3  Total    nonEU
#> 4   2018    Total
#> 5   2019    Total
#> 6  Total  Iceland
#> 7  Total Portugal
#> 8  Total    Spain
#> 9   2018       EU
#> 10  2019       EU
#> 11  2018    nonEU
#> 12  2019    nonEU
#> 
RoundViaDummy(SmallCountData("e6"), "freq", hierarchies = 
   list(geo = c("EU", "@Portugal", "@Spain", "Iceland"), year = c("2018", "2019")))
#> [-**.:=]
#> $yInner
#>      original rounded
#> [1,]        2       3
#> [2,]        3       3
#> [3,]        7       7
#> [4,]        1       0
#> [5,]        5       5
#> [6,]        6       6
#> 
#> $yPublish
#>                original rounded
#> Total:Total          24      24
#> Total:2018           12      13
#> Total:2019           12      11
#> EU:Total             21      21
#> EU:2018              10      10
#> EU:2019              11      11
#> Iceland:Total         3       3
#> Iceland:2018          2       3
#> Iceland:2019          1       0
#> Portugal:Total        8       8
#> Portugal:2018         3       3
#> Portugal:2019         5       5
#> Spain:Total          13      13
#> Spain:2018            7       7
#> Spain:2019            6       6
#> 
#> $crossTable
#>         geo  year
#> 1     Total Total
#> 2     Total  2018
#> 3     Total  2019
#> 4        EU Total
#> 5        EU  2018
#> 6        EU  2019
#> 7   Iceland Total
#> 8   Iceland  2018
#> 9   Iceland  2019
#> 10 Portugal Total
#> 11 Portugal  2018
#> 12 Portugal  2019
#> 13    Spain Total
#> 14    Spain  2018
#> 15    Spain  2019
#> 

RoundViaDummy(SmallCountData('z2'), 
              'ant', ~region + hovedint + fylke*hovedint + kostragr*hovedint, 10)
#> [-**..:=]
#> $yInner
#>       original rounded
#>  [1,]       11      11
#>  [2,]        7      10
#>  [3,]        5       5
#>  [4,]       13      13
#>  [5,]        9       9
#>  [6,]       12      12
#>  [7,]        6       6
#>  [8,]        9       9
#>  [9,]        3       3
#> [10,]        9       9
#> [11,]        4       4
#> [12,]       11      11
#> [13,]        1       0
#> [14,]        8       8
#> [15,]        2       2
#> [16,]       14      14
#> [17,]        9       9
#> [18,]        4      10
#> [19,]        3       0
#> [20,]        0       0
#> [21,]        0       0
#> [22,]        2       0
#> [23,]       55      55
#> [24,]       29      29
#> [25,]       35      35
#> [26,]       17      17
#> [27,]       63      63
#> [28,]       24      24
#> [29,]       22      22
#> [30,]       38      38
#> [31,]        9       9
#> [32,]       32      32
#> [33,]       18      18
#> [34,]       36      36
#> [35,]       18      18
#> [36,]       25      25
#> [37,]       13      13
#> [38,]       52      52
#> [39,]       22      22
#> [40,]        8       8
#> [41,]       15      15
#> [42,]        2       2
#> [43,]       20      20
#> [44,]       11      11
#> 
#> $yPublish
#>                original rounded
#> Total-Total         706     709
#> A-Total             113     113
#> B-Total              55      57
#> C-Total              73      73
#> D-Total              45      45
#> E-Total             138     138
#> F-Total              67      67
#> G-Total              40      46
#> H-Total              65      62
#> I-Total              14      14
#> J-Total              61      61
#> K-Total              35      33
#> Total-annet          88      91
#> Total-arbeid         54      54
#> Total-soshjelp      342     342
#> Total-trygd         222     222
#> 1-Total             127     127
#> 4-Total              55      57
#> 5-Total             118     118
#> 6-Total             205     205
#> 8-Total             105     108
#> 10-Total             96      94
#> 300-Total           596     601
#> 400-Total           110     108
#> 1-annet              14      14
#> 4-annet               7      10
#> 5-annet              18      18
#> 6-annet              21      21
#> 8-annet              15      15
#> 10-annet             13      13
#> 1-arbeid             11      11
#> 4-arbeid              1       0
#> 5-arbeid             10      10
#> 6-arbeid             23      23
#> 8-arbeid              7      10
#> 10-arbeid             2       0
#> 1-soshjelp           64      64
#> 4-soshjelp           29      29
#> 5-soshjelp           52      52
#> 6-soshjelp           87      87
#> 8-soshjelp           60      60
#> 10-soshjelp          50      50
#> 1-trygd              38      38
#> 4-trygd              18      18
#> 5-trygd              38      38
#> 6-trygd              74      74
#> 8-trygd              23      23
#> 10-trygd             31      31
#> 300-annet            72      75
#> 400-annet            16      16
#> 300-arbeid           52      54
#> 400-arbeid            2       0
#> 300-soshjelp        283     283
#> 400-soshjelp         59      59
#> 300-trygd           189     189
#> 400-trygd            33      33
#> 
#> $crossTable
#>    region hovedint
#> 1   Total    Total
#> 2       A    Total
#> 3       B    Total
#> 4       C    Total
#> 5       D    Total
#> 6       E    Total
#> 7       F    Total
#> 8       G    Total
#> 9       H    Total
#> 10      I    Total
#> 11      J    Total
#> 12      K    Total
#> 13  Total    annet
#> 14  Total   arbeid
#> 15  Total soshjelp
#> 16  Total    trygd
#> 17      1    Total
#> 18      4    Total
#> 19      5    Total
#> 20      6    Total
#> 21      8    Total
#> 22     10    Total
#> 23    300    Total
#> 24    400    Total
#> 25      1    annet
#> 26      4    annet
#> 27      5    annet
#> 28      6    annet
#> 29      8    annet
#> 30     10    annet
#> 31      1   arbeid
#> 32      4   arbeid
#> 33      5   arbeid
#> 34      6   arbeid
#> 35      8   arbeid
#> 36     10   arbeid
#> 37      1 soshjelp
#> 38      4 soshjelp
#> 39      5 soshjelp
#> 40      6 soshjelp
#> 41      8 soshjelp
#> 42     10 soshjelp
#> 43      1    trygd
#> 44      4    trygd
#> 45      5    trygd
#> 46      6    trygd
#> 47      8    trygd
#> 48     10    trygd
#> 49    300    annet
#> 50    400    annet
#> 51    300   arbeid
#> 52    400   arbeid
#> 53    300 soshjelp
#> 54    400 soshjelp
#> 55    300    trygd
#> 56    400    trygd
#> 
mf <- ~region*mnd + hovedint*mnd + fylke*hovedint*mnd + kostragr*hovedint*mnd
a <- RoundViaDummy(SmallCountData('z3'), 'ant', mf, 5)
#> [-**........=-**.:=-**.:=]
b <- RoundViaDummy(SmallCountData('sosialFiktiv'), 'ant', mf, 4)
#> [-**.........:=]
print(cor(b[[2]]),digits=12) # Correlation between original and rounded
#>                original        rounded
#> original 1.000000000000 0.999999987033
#> rounded  0.999999987033 1.000000000000

# Demonstrate parameter leverageCheck 
# The 42nd inner cell must be rounded since it can be revealed from the published cells.
mf2 <- ~region + hovedint + fylke * hovedint + kostragr * hovedint
RoundViaDummy(SmallCountData("z2"), "ant", mf2, leverageCheck = FALSE)$yInner[42, ]
#> [-**.:=]
#> original  rounded 
#>        2        2 
RoundViaDummy(SmallCountData("z2"), "ant", mf2, leverageCheck = TRUE)$yInner[42, ]
#> [{-x-=H-}(check:1)-**..:={-zx-=H-}]
#> original  rounded 
#>        2        3 

if (FALSE) { # \dontrun{
# Demonstrate parameters maxRound, zeroCandidates and forceInner 
# by tabulating the inner cells that have been changed.
z4 <- SmallCountData("sosialFiktiv")
for (forceInner in c("FALSE", "z4$ant < 10")) 
  for (zeroCandidates in c(FALSE, TRUE)) 
    for (maxRound in c(2, 5)) {
      set.seed(123)
      a <- RoundViaDummy(z4, "ant", formula = mf, maxRound = maxRound, 
                         zeroCandidates = zeroCandidates, 
                         forceInner = eval(parse(text = forceInner)))
      change <- a$yInner[, "original"] != a$yInner[, "rounded"]
      cat("\n\n---------------------------------------------------\n")
      cat("      maxRound:", maxRound, "\n")
      cat("zeroCandidates:", zeroCandidates, "\n")
      cat("    forceInner:", forceInner, "\n\n")
      print(table(original = a$yInner[change, "original"], rounded = a$yInner[change, "rounded"]))
      cat("---------------------------------------------------\n")
    }
} # }
```
