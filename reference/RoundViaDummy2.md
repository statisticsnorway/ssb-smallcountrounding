# RoundViaDummy2

RoundViaDummy with input as makeroundtabs

## Usage

``` r
RoundViaDummy2(
  data,
  b,
  d,
  nin,
  micro = FALSE,
  control = NULL,
  allTerms = FALSE,
  singleRandom = FALSE
)
```

## Arguments

- data:

  data

- b:

  b

- d:

  d

- nin:

  nin

- micro:

  micro

- control:

  control

- allTerms:

  Use all interaction terms in formula instead of using control

- singleRandom:

  Single random draw when TRUE (instead of algorithm)

## Value

Output from RoundViaDummy extended with "yControl" "maxdiff" "nMaxdiff"
"formula"
