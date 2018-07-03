# R package SmallCountRounding

### Small Count Rounding of Tabular Data


```r
library(devtools)                                     # Load package containing install_github
install_github("statisticsnorway/SSBtools")           # Install SSBtools from GitHub
install_github("statisticsnorway/SmallCountRounding") # Install SmallCountRounding from GitHub
library(SmallCountRounding)                           # Load SmallCountRounding 
?HierarchyCompute                                     # Help documentation of function RoundViaDummy
```

-----------

### The function RoundViaDummy

Small count rounding of necessary inner cells (Heldal, 2017) are performed so that all small frequencies of cross-classifications to be published (publishable cells) are rounded. This is equivalent to changing micro data since frequencies of unique combinations are changed. Thus, additivity and consistency are guaranteed.
This is performed by an algorithm inspired by partial least squares regression (Langsrud and Heldal, 2018).


### References

Heldal, J.: The European Census Hub 2011 Hypercubes - Norwegian SDC Experiences. In: Work Session on Statistical Data Confidentiality (2017), Skopje, The former Yugoslav Republic of Macedonia, September 20-22 , 2017.

Langsrud, Ø. and  Johan Heldal, J.: An Algorithm for Small Count Rounding of Tabular Data Øyvind Langsrud. 
To be presented at: Privacy in statistical databases, Valencia, Spain. September 26-28, 2018.