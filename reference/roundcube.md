# function roundcube

This function rounds small counts in a set of hypercubes D produced by
the function redcube and searches for a solution with smallest possible
deviations from the original hypercube at some aggregated levels.

## Usage

``` r
roundcube(rc, sort, control, minit, maxit, maxdiff, seed)
```

## Arguments

- rc:

  The list of outpts from redcube

- sort:

  An ordered list of variables in hypercubes in D meant for priority
  sorting of the reduced hypercube B before rounding. Not all variables
  in D should be included.

- control:

  A list of marginals of the hypercubes in D where deviations of
  aggregated rounded counts are checked against original counts.

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

Ar: The rounded version of A Br: The rounded version of B D: The
original hypercube of interest. Dr: The rounded version of D. The final
table of interest. maxdiff: The largest absolute difference between
cells D and Dr among cells in the control list. nmaxdiff: The number of
occurences if Maxdiff

## Author

Johan Heldal, January 2018
