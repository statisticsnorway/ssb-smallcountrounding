# PLSrounding on portions of data at a time

The
[`PLSrounding`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/PLSrounding.md)
runs are coordinated by using preliminary differences as input for the
next run (parameter `preDifference`)

## Usage

``` r
PLSroundingLoop(
  data,
  loopId,
  ...,
  zeroCandidates = FALSE,
  forceInner = FALSE,
  preRounded = NULL,
  plsWeights = NULL,
  printInc = TRUE,
  preDifference = TRUE,
  preOutput = NULL,
  rndSeed = 123,
  action_unused_dots = "warn"
)
```

## Arguments

- data:

  Input data as a data frame (inner cells)

- loopId:

  Variable holding id for loops

- ...:

  `PLSrounding` parameters

- zeroCandidates:

  `PLSrounding` parameter (see details)

- forceInner:

  `PLSrounding` parameter (see details)

- preRounded:

  `PLSrounding` parameter (see details)

- plsWeights:

  `PLSrounding` parameter (see details)

- printInc:

  Printing iteration information to console when TRUE

- preDifference:

  When TRUE, the `preDifference` parameter to `PLSrounding` is used.
  Each time with the differences obtained so far.

- preOutput:

  preOutput The function can continue from output from a previous run

- rndSeed:

  If non-NULL, a random generator seed to be set locally at the
  beginning of `PLSroundingLoop` without affecting the random value
  stream in R. Within `PLSroundingLoop`, `PLSrounding` is called with
  `rndSeed = NULL`.

- action_unused_dots:

  `PLSrounding` parameter.

## Value

As output from
[`PLSrounding`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/PLSrounding.md)

## Details

Note that in this function `zeroCandidates`, `forceInner`, `preRounded`
and `plsWeights` cannot be supplied as vectors. They may be specified as
functions or as variables in the input data.

## Examples

``` r
mf2 <- ~region + fylke * hovedint
z2 <- SmallCountData("z2")
a <- PLSroundingLoop(z2, loopId = "kostragr", freqVar = "ant", formula = mf2)
#>    1: []
#>    2: [-]
a
#> 
#> PLSrounding summary:  
#> 
#>        maxdiff      HDutility    meanAbsDiff rootMeanSquare 
#>              1         0.9873         0.2391          0.489 
#> 
#> Frequencies of cell frequencies and absolute differences:  
#> 
#>          inn.0 inn.1-2 inn.3 inn.4+ inn.all pub.0 pub.1-2 pub.3 pub.4+ pub.all
#> original     2       4     2     36      44     .       2     .     44      46
#> rounded      3       1     4     36      44     1       .     1     44      46
#> absDiff     41       3     .      .      44    35      11     .      .      46
#> 
```
