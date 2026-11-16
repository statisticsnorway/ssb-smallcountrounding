# FormulaSelection method for PLSrounded

FormulaSelection method for PLSrounded

## Usage

``` r
# S3 method for class 'PLSrounded'
FormulaSelection(x, formula = NULL, intercept = NA, logical = FALSE)
```

## Arguments

- x:

  PLSrounded object

- formula:

  `formula` parameter to
  [`FormulaSelection`](https://statisticsnorway.github.io/ssb-ssbtools/reference/FormulaSelection.html).
  When `NULL` (default), the publish data frame is returned without any
  limitation.

- intercept:

  `intercept` parameter to `FormulaSelection`.

- logical:

  `logical` parameter to `FormulaSelection`.

## Value

Limited version of the publish data frame
