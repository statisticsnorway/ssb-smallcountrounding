# Hellinger Distance (Utility)

Hellinger distance (`HD`) and a related utility measure (`HDutility`)
described in the reference below. The utility measure is made to be
bounded between 0 and 1.

## Usage

``` r
HD(f, g)

HDutility(f, g)
```

## Arguments

- f:

  Vector of original counts

- g:

  Vector of perturbed counts

## Value

Hellinger distance or related utility measure

## Details

HD is defined as "`sqrt(sum((sqrt(f) - sqrt(g))^2)/2)`" and HDutility is
defined as "`1 - HD(f, g)/sqrt(sum(f))`".

## References

Shlomo, N., Antal, L., & Elliot, M. (2015). Measuring Disclosure Risk
and Data Utility for Flexible Table Generators, Journal of Official
Statistics, 31(2), 305-324.
[doi:10.1515/jos-2015-0019](https://doi.org/10.1515/jos-2015-0019)

## Examples

``` r
f <- 1:6
g <- c(0, 3, 3, 3, 6, 6)
print(c(
  HD = HD(f, g), 
  HDutility = HDutility(f, g), 
  maxdiff = max(abs(g - f)), 
  meanAbsDiff = mean(abs(g - f)), 
  rootMeanSquare = sqrt(mean((g - f)^2))
))
#>             HD      HDutility        maxdiff    meanAbsDiff rootMeanSquare 
#>      0.7805018      0.8296805      1.0000000      0.6666667      0.8164966 
```
