# R package SmallCountRounding 
[![Mentioned in Awesome Official Statistics ](https://awesome.re/mentioned-badge.svg)](http://www.awesomeofficialstatistics.org)


## Small Count Rounding of Tabular Data 


See the package vignette: 
[Introduction to ‘SmallCountRounding’](https://cran.r-project.org/web/packages/SmallCountRounding/vignettes/Introduction_to_SmallCountRounding.html)

-----------

#### Installation from CRAN 
(Recommended, unless you want to test the newest changes.)

```r
install.packages("SmallCountRounding")
```

#### Installation from GitHub

(If needed, install the `devtools` package from CRAN first.)


```r
library(devtools)                                     # Load package containing install_github
install_github("statisticsnorway/SSBtools")           # Install SSBtools from GitHub 
install_github("statisticsnorway/SmallCountRounding") # Install SmallCountRounding from GitHub
```

####  Load the package 

```r
library(SmallCountRounding)                           # Load SmallCountRounding 
?PLSrounding                                          # Help documentation of function PLSrounding
```


-----------

#### A statistical disclosure control tool to protect frequency tables in cases where small values are sensitive

The main function, PLSrounding(), performs small count rounding of necessary inner cells (Heldal, 2017)
so that all small frequencies of cross-classifications to be published (publishable cells) are rounded. This is equivalent to changing micro data since frequencies of unique combinations are changed. Thus, additivity and consistency are guaranteed.
This is performed by an algorithm inspired by partial least squares regression (Langsrud and Heldal, 2018).


#### References

Heldal, J.: The European Census Hub 2011 Hypercubes - Norwegian SDC Experiences. In: Work Session on Statistical Data Confidentiality (2017), Skopje, The former Yugoslav Republic of Macedonia, September 20-22 , 2017.

Langsrud, Ø. and Heldal, J.: An Algorithm for Small Count Rounding of Tabular Data. 
To be presented at: Privacy in statistical databases, Valencia, Spain. September 26-28, 2018.
 [Full-text](https://www.researchgate.net/publication/327768398_An_Algorithm_for_Small_Count_Rounding_of_Tabular_Data)
 and
[presentation](https://www.researchgate.net/publication/327916165_An_Algorithm_for_Small_Count_Rounding_of_Tabular_Data_-_Presentation) 
 available at Researchgate.
 
 -----------
 Official version on CRAN: [https://cran.r-project.org/package=SmallCountRounding](https://cran.r-project.org/package=SmallCountRounding)
