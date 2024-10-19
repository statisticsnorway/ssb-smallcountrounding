
printInc <- FALSE


tre <- function(x) {
  c(nrow(x), sum(x$ipFit > 0), sum(x$difference == 5))
}


test_that("PLSroundingFits with extend0/extend0Fits", {
  my_km2 <- SSBtools::SSBtoolsData("my_km2")
  
  # Default automatic extension (extend0 = TRUE)
  a <- PLSroundingFits(my_km2, "freq", roundBase = 5, zeroCandidates = TRUE,
                       printInc = printInc,
          formula = ~(Sex + Age) * Municipality * Square1000m + Square250m)[["inner"]]
  
  expect_equal(tre(a), c(60, 27, 0)) 
  
  b <- PLSroundingFits(my_km2, "freq", roundBase = 5, zeroCandidates = TRUE,
                       printInc = printInc,
                       formula = ~(Sex + Age) * Municipality * Square1000m + Square250m,
                  extend0Fits = list(c("Sex", "Age"), 
          c("Municipality", "Square1000m", "Square250m")))[["inner"]]
  
  expect_equal(tre(b), c(36, 21, 0)) 
  
  
  d <- PLSroundingFits(my_km2, "freq", roundBase = 5, zeroCandidates = TRUE, 
                       printInc = printInc,
                       formula = ~(Sex + Age) * Municipality * Square1000m + Square250m, 
                       extend0 = list(c("Sex", "Age"), 
                                      c("Municipality", "Square1000m", "Square250m")))[["inner"]]
  
  expect_equal(tre(d), c(60, 27, 1)) 
  
  
  e <- PLSroundingFits(my_km2, "freq", roundBase = 5, zeroCandidates = TRUE, 
                       printInc = printInc,
                       formula = ~(Sex + Age) * Municipality * Square1000m + Square250m, 
                       extend0 = TRUE)[["inner"]]
  
  expect_equal(tre(e), c(60, 27, 3)) 
  
  
})
