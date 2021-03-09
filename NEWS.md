

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