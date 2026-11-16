# function makeroundtabs

This function creates a set of consistently rounded frequency count
tables or hypercubes by means of a version of small count rounding.

## Usage

``` r
makeroundtabs(
  A,
  b = 3,
  d,
  micro = "TRUE",
  sort,
  control,
  nin = "n",
  nout = "n",
  minit = 3,
  maxit = 3,
  maxdiff = 5,
  seed
)
```

## Arguments

- A:

  A data frame representing a micro dataset or a frequency count
  hypercube. The (first) columns define the variables. If A is a
  hypercube the last column contains the number of units in each cell.
  If A is a micro dataset it is reduced to hypercube by the function
  aggrtab.

- b:

  Rounding base. Counts in A less than b tat are contributing to counts
  less than b in the marginal cubes D are selected from A. The selected
  dataframe is called B

- d:

  A list d{\[\[j\]\]} whose elements are vectors of variable names from
  A defining marginal tables/cubes D of A that we are interested in.

- micro:

  Logical. TRUE if A is a micro dataset (default). FALSE if A i a
  frequency count hypercube.

- sort:

  An ordered list of variables in hypercubes in D meant for priority
  sorting of the reduced hypercube B before rounding. Not all variables
  in D should be included.

- control:

  A list of marginals of the hypercubes in D where deviations of
  aggregated rounded counts are checked against original counts.

- nin:

  Name of count variable if A is a hypercube. Default name: "n".

- nout:

  Name of the frequency count variable in the output tables.

- minit:

  Minimum number of searches to be carried out.

- maxit:

  Maximum number of searches to be carried out.

- maxdiff:

  If maximum difference in "control" is no larger than maxit, the stop
  search.

- seed:

  Input seed for first systematic random search.

## Value

Ar: The rounded version of A

Br: The rounded version of B

D: The original hypercube of interest.

Dr: The rounded version of D. The final table of interest.

maxdiff: The largest absolute difference between cells D and Dr among
cells in the control list.

nmaxdiff: The number of occurences if Maxdiff

## See also

Dependencies:
[`aggrtab`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/aggrtab.md),
[`redcube`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/redcube.md),
[`roundcube`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/roundcube.md)

## Author

Johan Heldal, January 2018
