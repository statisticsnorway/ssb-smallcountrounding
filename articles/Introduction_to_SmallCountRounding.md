# Introduction to 'SmallCountRounding'

  
  

## Introductory example

### Example without code

First some data to be rounded.  

|       | col1 | col2 | col3 | col4 | col5 | Total |
|:------|-----:|-----:|-----:|-----:|-----:|------:|
| row1  |    6 |    0 |    1 |    3 |    4 |    14 |
| row2  |    1 |    2 |    3 |    1 |    2 |     9 |
| row3  |    0 |    1 |    1 |    0 |    2 |     4 |
| Total |    7 |    3 |    5 |    4 |    8 |    27 |

**Table 1**: Original data in tabular form with row and column totals
{.table .table .table-bordered
style="font-size: 16px; width: auto !important; "}

  
Some of the inner cells are rounded. Thereafter new totals are computed.
The underlying algorithm tries to keep the values of these totals close
to the original ones.  

|       | col1 | col2 | col3 | col4 | col5 | Total |
|:------|-----:|-----:|-----:|-----:|-----:|------:|
| row1  |    6 |    0 |    0 |    5 |    5 |    16 |
| row2  |    0 |    0 |    5 |    0 |    0 |     5 |
| row3  |    0 |    0 |    0 |    0 |    5 |     5 |
| Total |    6 |    0 |    5 |    5 |   10 |    26 |

**Table 2**: All small inner cell values (1-4) are rounded using 5 as
rounding base. {.table .table .table-bordered
style="font-size: 16px; width: auto !important; "}

  
When the inner cells are not going to be published, the number of cells
to be rounded can be limited.  

|       | col1 | col2 | col3 | col4 | col5 | Total |
|:------|-----:|-----:|-----:|-----:|-----:|------:|
| row1  |    6 |    0 |    0 |    5 |    4 |    15 |
| row2  |    1 |    0 |    5 |    0 |    2 |     8 |
| row3  |    0 |    5 |    0 |    0 |    0 |     5 |
| Total |    7 |    5 |    5 |    5 |    6 |    28 |

**Table 3**: Assuming only row and column totals to be published,
necessary small inner cell values (1-4) are rounded using 5 as rounding
base. {.table .table .table-bordered
style="font-size: 16px; width: auto !important; "}

------------------------------------------------------------------------

### The example dataset in Table 1

``` r

library(SmallCountRounding)
z <- SmallCountData("exPSD")
z
```

       rows cols freq
    1  row1 col1    6
    2  row2 col1    1
    3  row3 col1    0
    4  row1 col2    0
    5  row2 col2    2
    6  row3 col2    1
    7  row1 col3    1
    8  row2 col3    3
    9  row3 col3    1
    10 row1 col4    3
    11 row2 col4    1
    12 row3 col4    0
    13 row1 col5    4
    14 row2 col5    2
    15 row3 col5    2

### Rounding all small cells (Table 2)

To avoid any small values in the range 1-4 we can use 5 as rounding
base.

``` r

a <- PLSrounding(z, freqVar = "freq", roundBase = 5, dimVar = c("rows", "cols"))
```

The result is given in Table 2 and can bee seen in the output elements
below.

``` r

a$inner
```

       rows cols original rounded difference
    1  row1 col1        6       6          0
    2  row2 col1        1       0         -1
    3  row3 col1        0       0          0
    4  row1 col2        0       0          0
    5  row2 col2        2       0         -2
    6  row3 col2        1       0         -1
    7  row1 col3        1       0         -1
    8  row2 col3        3       5          2
    9  row3 col3        1       0         -1
    10 row1 col4        3       5          2
    11 row2 col4        1       0         -1
    12 row3 col4        0       0          0
    13 row1 col5        4       5          1
    14 row2 col5        2       0         -2
    15 row3 col5        2       5          3

``` r

a$publish
```

        rows  cols original rounded difference
    1  Total Total       27      26         -1
    2  Total  col1        7       6         -1
    3  Total  col2        3       0         -3
    4  Total  col3        5       5          0
    5  Total  col4        4       5          1
    6  Total  col5        8      10          2
    7   row1 Total       14      16          2
    8   row1  col1        6       6          0
    9   row1  col2        0       0          0
    10  row1  col3        1       0         -1
    11  row1  col4        3       5          2
    12  row1  col5        4       5          1
    13  row2 Total        9       5         -4
    14  row2  col1        1       0         -1
    15  row2  col2        2       0         -2
    16  row2  col3        3       5          2
    17  row2  col4        1       0         -1
    18  row2  col5        2       0         -2
    19  row3 Total        4       5          1
    20  row3  col1        0       0          0
    21  row3  col2        1       0         -1
    22  row3  col3        1       0         -1
    23  row3  col4        0       0          0
    24  row3  col5        2       5          3

The output element `publish` contains the original and rounded versions
of the all the 24 values in Table 2. The corresponding element `inner`
contains only the 15 inner cells and is similar to the input data. The
values in publish are `additive`. That is, marginal cells (Totals) can
be computed straightforwardly from `inner` for both original and rounded
counts.

### Rounding necessary small inner cell (Table 3)

Assuming only row and column totals to be published, the publishable
cells can be defined by the formula `~rows+cols`. Rounding can now be
performed by:

``` r

b <- PLSrounding(z, "freq", 5, formula = ~rows + cols)
```

The result is given in Table 3 and can bee seen in the output elements
below.

``` r

b$inner
```

       rows cols original rounded difference
    1  row1 col1        6       6          0
    2  row2 col1        1       1          0
    3  row3 col1        0       0          0
    4  row1 col2        0       0          0
    5  row2 col2        2       0         -2
    6  row3 col2        1       5          4
    7  row1 col3        1       0         -1
    8  row2 col3        3       5          2
    9  row3 col3        1       0         -1
    10 row1 col4        3       5          2
    11 row2 col4        1       0         -1
    12 row3 col4        0       0          0
    13 row1 col5        4       4          0
    14 row2 col5        2       2          0
    15 row3 col5        2       0         -2

``` r

b$publish
```

       rows  cols original rounded difference
    1 Total Total       27      28          1
    2  row1 Total       14      15          1
    3  row2 Total        9       8         -1
    4  row3 Total        4       5          1
    5 Total  col1        7       7          0
    6 Total  col2        3       5          2
    7 Total  col3        5       5          0
    8 Total  col4        4       5          1
    9 Total  col5        8       6         -2

### Unique output obtained by local random generator seed

The underlying algorithm is sequential. Within a loop, the next cell to
be given the rounding base value is selected according to a criterion.
Random draw is used when draw criterion. To ensure unique output, a
fixed random generator seed is used locally within the function without
affecting the random value stream in R. See the documentation of
`rndSeed`, a parameter to `RoundViaDummy`.

### The output object

The result of printing the output from `PLSrounding` is (`a` and `b` as
above):

``` r

a
```


    PLSrounding summary:  

           maxdiff      HDutility    meanAbsDiff rootMeanSquare 
                 4          0.744         1.3333         1.6833 

    Frequencies of cell frequencies and absolute differences:  

             inn.0 inn.1-4 inn.5 inn.6+ inn.all pub.0 pub.1-4 pub.5 pub.6+ pub.all
    original     3      11     .      1      15     3      14     1      6      24
    rounded     10       .     4      1      15    11       .     8      5      24
    absDiff      4      11     .      .      15     5      19     .      .      24

``` r

b
```


    PLSrounding summary:  

           maxdiff      HDutility    meanAbsDiff rootMeanSquare 
                 2          0.941              1         1.2019 

    Frequencies of cell frequencies and absolute differences:  

             inn.0 inn.1-4 inn.5 inn.6+ inn.all pub.0 pub.1-4 pub.5 pub.6+ pub.all
    original     3      11     .      1      15     .       3     1      5       9
    rounded      8       3     3      1      15     .       .     4      5       9
    absDiff      7       8     .      .      15     2       7     .      .       9

First some utility measures are printet. For example `maxdiff` is the
maximum difference between an original and rounded cells within
`publish`. Thereafter a table of frequencies of cell frequencies and
absolute differences are printed. Summary of `inner` and `publish` are
shown in the left and right parts of the table, respectively. For
example, row `rounded` and column `inn.6+` is the number of rounded
inner cell frequencies greater than or equal to 6. The last row
(`absDiff`) is based on the differences without signs.

It is possible to compute manually the printed utility measures by:

``` r

f <- b$publish$original
g <- b$publish$rounded
print(c(
  maxdiff        =  max(abs(g - f)), 
  HDutility      =  HDutility(f, g), 
  meanAbsDiff    =  mean(abs(g - f)), 
  rootMeanSquare =  sqrt(mean((g - f)^2))
))
```

           maxdiff      HDutility    meanAbsDiff rootMeanSquare 
          2.000000       0.940951       1.000000       1.201850 

These measures are also found in the output element `metrics` together
with the same measures based on `inner`. See
[`?HDutility`](https://statisticsnorway.github.io/ssb-smallcountrounding/reference/HDutility.md)
for more information about the utility measure based on the Hellinger
distance.

Apart from printing, output is a usual list and `summary` works as
usual.

``` r

summary(b)
```

              Length Class      Mode   
    inner      5     data.frame list   
    publish    5     data.frame list   
    metrics    9     -none-     numeric
    freqTable 30     -none-     numeric

The output element `freqTable` is the table seen when the output object
is printed (frequencies of cell frequencies and absolute differences).

  
  

## Hierarchical data

### Example without code

Below is a small data set to be used as input.

| geo      | eu    | year | freq |
|:---------|:------|:-----|-----:|
| Iceland  | nonEU | 2018 |    2 |
| Portugal | EU    | 2018 |    3 |
| Spain    | EU    | 2018 |    7 |
| Iceland  | nonEU | 2019 |    1 |
| Portugal | EU    | 2019 |    5 |
| Spain    | EU    | 2019 |    6 |

**Table 4**: Input data {.table .table .table-bordered
style="font-size: 14px; width: auto !important; "}

  

The variables `geo` and `eu` is hierarchical related. This data set can
be processed in several ways. In some cases, the entire table will be
input and in other cases the `eu` column can be omitted. Then, the
hierarchical information is sent as input in another way. One
possibility is the table below, where the hierarchy is coded as in the r
package sdcTable.

| levels | codes    |
|:-------|:---------|
| @      | Total    |
| @@     | EU       |
| @@@    | Portugal |
| @@@    | Spain    |
| @@     | nonEU    |
| @@@    | Iceland  |

**Table 5**: Hierarchy, `geo` {.table .table .table-bordered
style="font-size: 14px; width: auto !important; "}

  

Another possibility is TauArgus coding. More general coding is also
possible. See `?AutoHierarchies` for more information.

Below is output in the case were all possible combinations (including
the inner cells) are to be published. Also in this example we use 5 as a
rounding base. As can be seen below, this output can be generated in
several ways. The inner cells are colored according to the rounding.

| geo      | year  | original | rounded | difference |
|:---------|:------|---------:|--------:|-----------:|
| Total    | Total |       24 |      23 |         -1 |
| Total    | 2018  |       12 |      12 |          0 |
| Total    | 2019  |       12 |      11 |         -1 |
| EU       | Total |       21 |      23 |          2 |
| EU       | 2018  |       10 |      12 |          2 |
| EU       | 2019  |       11 |      11 |          0 |
| nonEU    | Total |        3 |       0 |         -3 |
| nonEU    | 2018  |        2 |       0 |         -2 |
| nonEU    | 2019  |        1 |       0 |         -1 |
| Iceland  | Total |        3 |       0 |         -3 |
| Iceland  | 2018  |        2 |       0 |         -2 |
| Iceland  | 2019  |        1 |       0 |         -1 |
| Portugal | Total |        8 |      10 |          2 |
| Portugal | 2018  |        3 |       5 |          2 |
| Portugal | 2019  |        5 |       5 |          0 |
| Spain    | Total |       13 |      13 |          0 |
| Spain    | 2018  |        7 |       7 |          0 |
| Spain    | 2019  |        6 |       6 |          0 |

**Table 6**: Ouput data (publish) {.table .table .table-bordered
style="font-size: 14px; width: auto !important; "}

### Data (Table 4) and hierarchies (Table 5)

``` r

e6  <-  SmallCountData("e6")             # As Table 4 
eDimList <- SmallCountData("eDimList")
eDimList
```

    $geo
      levels    codes
    1      @    Total
    2     @@       EU
    3    @@@ Portugal
    4    @@@    Spain
    5     @@    nonEU
    6    @@@  Iceland

    $year
      levels codes
    1      @ Total
    2     @@  2018
    3     @@  2019

As seen above, a hierarchy is specified for both variables.
`eDimList$geo` is given in Table 5 and `eDimList$year` is a plain
hierarchy with total code.

### Four ways to produce Table 6

The four lines below produce the same results with element `publish` as
in Table 6. Ordering of rows can be different.

``` r

#  PLSrounding(e6, "freq", 5) # This is no longer possible. See ver. 1.1.0 news # a)
PLSrounding(e6, "freq", 5, dimVar = c("geo", "eu", "year"))                     # b) 
PLSrounding(e6, "freq", 5, formula = ~eu * year + geo * year)                   # c)
PLSrounding(e6[, -2], "freq", 5, hierarchies = eDimList)                        # d)
PLSrounding(e6[, -2], "freq", 5, hierarchies = eDimList, formula = ~geo * year) # e)
```

- In b) and d), the function uses hierarchies for the calculations.
  In b) the hierarchies are found automatically from the input data. In
  a), `dimVar` is assumed to be all variables except `freq`.
- In c) the cross-classifications are found from the formula. In
  addition, hierarchical relations in the input data are analysed so
  that `geo` and `eu` are combined into the same output column.  
- In e), how to cross the hierarchies are defined by a formula.

### Remarks and other parameters

A difference occur when all combinations are not contained in input
data. Then c) above will limit output to combinations available in
input. In the other cases zeroes will be added. The extra zeroes can be
avoided by using `removeEmpty=TRUE`. Note also the parameter
`inputInOutput` which can be used to specify whether to include codes
from input. Below is an example with incomplete input data using both
these parameters.

``` r

out <- PLSrounding(e6[-1, ], "freq", 5, removeEmpty = TRUE, inputInOutput = c(FALSE,TRUE),
                   dimVar = c("geo", "eu", "year"))
out
```


    PLSrounding summary:  

           maxdiff      HDutility    meanAbsDiff rootMeanSquare 
                 1         0.8925            0.5         0.7071 

    Frequencies of cell frequencies and absolute differences:  

             inn.0 inn.1-4 inn.5 inn.6+ inn.all pub.0 pub.1-4 pub.5 pub.6+ pub.all
    original     .       2     1      2       5     .       2     .      6       8
    rounded      1       1     1      2       5     2       .     .      6       8
    absDiff      4       1     .      .       5     4       4     .      .       8

``` r

out$inner
```

           geo year original rounded difference
    2 Portugal 2018        3       3          0
    3    Spain 2018        7       7          0
    4  Iceland 2019        1       0         -1
    5 Portugal 2019        5       5          0
    6    Spain 2019        6       6          0

``` r

out$publish
```

        geo  year original rounded difference
    1 Total Total       22      21         -1
    2 Total  2018       10      10          0
    3 Total  2019       12      11         -1
    4    EU Total       21      21          0
    5    EU  2018       10      10          0
    6    EU  2019       11      11          0
    7 nonEU Total        1       0         -1
    8 nonEU  2019        1       0         -1

In this case only a single inner cell needed to be rounded (Iceland,
2019). The original small value of (Portugal, 2018) could be retained.
