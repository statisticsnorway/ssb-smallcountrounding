
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