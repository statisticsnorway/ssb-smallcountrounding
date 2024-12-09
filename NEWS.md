
## SmallCountRounding	1.1.0
* Added a check to ensure that at least one of `dimVar`, `hierarchies`, or `formula` is specified.
  - This is a breaking change that may affect previous code.
  - Previously, if all were unspecified, `dimVar` was automatically generated from the remaining columns.
  - While this behavior was correctly implemented, it often stemmed from user input errors and could lead to unexpected behavior or crashes.
  - This change now requires explicit input, making the function more robust and reducing the risk of user errors.
* Improved support for `tibble` and `data.table` input (parameter `data`).
  - Input is now explicitly coerced to a data frame using `as.data.frame()` where necessary to ensure consistent behavior.
  - When `preAggregate` is `TRUE` and `aggregatePackage` is `"data.table"`, the use of `as.data.frame()` is skipped to avoid unnecessary back-and-forth conversion of `data.table` objects, preserving efficiency. 
  - Applies to `PLSrounding()` and its wrappers.
* Note new hierarchy possibilities due to the new version of 
     [the SSBtools package](https://CRAN.R-project.org/package=SSBtools) (version 1.6.0).
  -  Output from functions like `get_klass()` in the 
      [klassR package](https://cran.r-project.org/package=klassR) 
      or `hier_create()` in the 
      [sdcHierarchies package](https://cran.r-project.org/package=sdcHierarchies) 
      can now be used directly as input. Example of usage:
     ```r
      a <- get_klass(classification = "24")
      b <- hier_create(root = "Total", nodes = LETTERS[1:5])
      mydata <- data.frame(tree = sample(a$code[nchar(a$code) > 1], 200, replace = TRUE), 
                           letter = LETTERS[1:5])
      PLSroundingPublish(mydata, roundBase = 5, hierarchies = list(tree = a, letter = b)) 
     ```
  - New possibilities for working with both formulas and hierarchies are now available through the `map_hierarchies_to_data()` function. 
  - Improved functionality for combining formulas with the `Formula2ModelMatrix()` parameter `avoidHierarchical = TRUE`, 
    thanks to the new `total_collapse()` function which can be applied to output.

## SmallCountRounding	1.0.8
* `FormulaSelection()` now works with the output from `PLSrounding()`.
  - The output dataset corresponding to a restricted part of the input formula can now be easily retrieved.
  - See the examples in the documentation.
* `extend0`  is new parameter to `PLSrounding()`, enabling data to be automatically extended by zero frequency rows.
  - This is relevant when `zeroCandidates = TRUE`.
  - The old parameter in `PLSroundingFits()` has been renamed from `extend0` to `extend0Fits`. Code that used the old parameter will now behave differently.
  - Note that `extend0` and `extend0Fits` can now be specified in more advanced ways beyond just TRUE/FALSE.
* Improvements to the `step` parameter, which can be passed to `PLSrounding()` and is documented in the underlying function `RoundViaDummy()`:
  - A bug that could cause a hang when using `step` has been fixed.
  - The `step` parameter can now be specified as a vector for greater control.
  - Additionally, it can be provided as a list to trigger a final re-run iteration.
  - The `step` parameter can significantly impact performance on large datasets. For example, using `step = list(100)` may be a useful approach.
* Due to updates in [the SSBtools package](https://CRAN.R-project.org/package=SSBtools) (version 1.5.4), 
  it is now meaningful to include NA’s in the grouping variables. 
  - Note the parameter `NAomit` to `SSBtools::Formula2ModelMatrix()`: 
    * When `TRUE`, NAs in the grouping variables are omitted in output and not included as a separate category.
    * This parameter can be input to `PLSrounding()` and its wrappers.
  - `aggregateNA` is new parameter to `PLSrounding()`:
    * Whether to include NAs in the grouping variables while preAggregate.
    * Needs to be `TRUE` (default) to utilize the above `NAomit` parameter.
* Due to updates in [the SSBtools package](https://CRAN.R-project.org/package=SSBtools) (version 1.5.4), 
  where [data.table](https://cran.r-project.org/package=data.table) is now listed under *Suggests*, 
  some functionality can be speeded up. 
  - Set the new parameter `aggregatePackage` to  `"data.table"` to utilize this possibility.
    * `aggregatePackage` is parameter to `PLSrounding()` and its wrappers.
    * Also note the related new parameters `aggregateBaseOrder`. 



## SmallCountRounding	1.0.5
* Minor updates with no changes in functionality
  - Changed package license to MIT, in accordance with the policy at Statistics Norway.
  - Some technical changes in documentation to comply with standards.


## SmallCountRounding	1.0.3
* Workaround for old `R` versions where the `isFALSE` function is not defined.


## SmallCountRounding	1.0.2
* Improved behavior of the `identifyNew` parameter when the `maxRound` parameter is used. 
  - New description of the `identifyNew` parameter: 
                     When `TRUE`, new cells may be identified after initial rounding to ensure all rounded publishable 
                     cells equal to or less than `maxRound` to be `roundBase` multiples. Use `NA` for the a less conservative 
                     behavior (old behavior). Then it is ensured that no nonzero rounded publishable cells are smaller 
                     than `roundBase`. When `maxRound` is default, there is no difference between `TRUE` and `NA`.


## SmallCountRounding	1.0.0
* New function, `PLSroundingLoop`: PLSrounding on portions of data at a time.
  - The \code{\link{PLSrounding}} runs are coordinated by using preliminary differences as input for the next run (parameter `preDifference`)
* Parameters `zeroCandidates`, `forceInner`,  `preRounded` and `plsWeights` can now be specified as functions.
  - These cannot be supplied as vectors in `PLSroundingLoop`.
* New parameter, `allSmall`.
  - When TRUE, all small inner cells (`<= maxRound`) are rounded. A simplified alternative to specifying `forceInner`.
* Adaption needed after Matrix ver. 1.4-2 (not a user-visible change)

## SmallCountRounding	0.9.0
* New function, `PLSroundingFits`, for post-processing to expected frequencies
  - Expected inner cell frequencies are generated by iterative proportional fitting
* `plsWeights`	is new parameter to `RoundViaDummy` (and `PLSrounding`)
  - A vector of weights for each cell to be published or a function generating it. For use in the algorithm criterion. 


## SmallCountRounding	0.8.0
* Now, microdata input is allowed. This is due to
  - Allowing empty `freqVar` in input.
  - The new parameter `preAggregate`: When `TRUE`, the data will be aggregated beforehand within the function by the dimensional variables.
* It is possible to avoid handling hierarchical variables when using the formula interface.
  - This is a consequence of parameter `avoidHierarchical` to `Formula2ModelMatrix` in the SSBtools package.


## SmallCountRounding	0.7.0
* Now, a random generator seed is used locally within the function without affecting the random value stream in R.
  - Handled by `rndSeed`, a new parameter to `RoundViaDummy` (and `PLSrounding`). 
  - By default, `rndSeed = 123`. This means that repeated runs with equal input will result in equal output. 
  - To get back the old behaviour of the function, set `rndSeed` to `NULL`.
* Possible to return a single data frame: `"inner"` or `"publish"`.
  - Handled by `output`, a new parameter to `PLSrounding`. 
  - New wrapper functions, `PLSroundingInner` and  `PLSroundingPublish`.
* `dimVar`	is new parameter to `RoundViaDummy` and `PLSrounding`
  - The main dimensional variables and additional aggregating variables. This parameter can be  useful when hierarchies and formula are unspecified. 


## SmallCountRounding	0.6.0
* `preRounded`	is new parameter to `RoundViaDummy` (and `PLSrounding`)
  - A mixture of missing values and predetermined values of rounded inner cells


## SmallCountRounding	0.5.0

* Formula combined with hierarchies is now possible
  - This is a consequence of the function `HierarchiesAndFormula2ModelMatrix` in the SSBtools package
* `leverageCheck` and `easyCheck` are new parameters to `RoundViaDummy`
  - This provides protection against possible disclosure of small numbers by linear relationships (difference attack)
  - The function `Reduce0exact` in the SSBtools package is utilised
* `printInc` is new parameter to `PLSrounding` and `RoundViaDummy`
  - Whether to print iteration information to the console
* Possible to specify `removeEmpty=TRUE` to omit empty combinations
  -  Parameter to `Hierarchies2ModelMatrix` and `HierarchiesAndFormula2ModelMatrix` in the SSBtools package
* The parameter `inputInOutput` is also mentioned in the RoundViaDummy documentation 
  -  Parameter to same functions as above
  -  Can be used to specify whether to include codes from input
* A vignette, "Introduction to 'SmallCountRounding'", is included  

  
## SmallCountRounding	0.4.0

* Last version before any news
