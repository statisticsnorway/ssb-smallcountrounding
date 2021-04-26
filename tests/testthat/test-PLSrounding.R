
printInc = FALSE

test_that("PLSrounding works", {
  # Small example data set
  z <- SmallCountData("e6")

  printInc = FALSE
    
  a <- PLSrounding(z, "freq", printInc = printInc)
  expect_equivalent(a$metrics, PLSrounding(z, "freq", formula = ~eu * year + geo * year, printInc = printInc)$metrics)
  expect_equivalent(a$metrics, PLSrounding(z[, -2], "freq", hierarchies = SmallCountData("eHrc"), printInc = printInc)$metrics)
  expect_equivalent(a$metrics, PLSrounding(z[, -2], "freq", hierarchies = SmallCountData("eDimList"), printInc = printInc)$metrics)
  expect_equivalent(a$metrics, PLSrounding(z[, -2], "freq", hierarchies = SmallCountData("eDimList"), formula = ~geo * year, printInc = printInc)$metrics)
  
  expect_equivalent(PLSrounding(z[, -2], "freq", hierarchies = SmallCountData("eDimList"), formula = ~geo + year, printInc = printInc)$metrics["maxdiff"], 0)
  
  expect_equivalent(PLSroundingInner(z, "freq", printInc = printInc, dimVar = c("geo", "eu"), roundBase = 5)$difference, c(3, 0, 0, -1, 0, 0))
  expect_equivalent(PLSroundingInner(z, "freq", printInc = printInc, dimVar = c("geo", "eu"), roundBase = 5, rndSeed = 1234)$difference, c(-2, 0, 0, 4, 0, 0))
  set.seed(12345)
  expect_equivalent(PLSroundingInner(z, "freq", printInc = printInc, dimVar = c("geo", "eu"), roundBase = 5, rndSeed = NULL)$difference, c(-2, 0, 0, 4, 0, 0))
  expect_equivalent(PLSroundingPublish(z, "freq", printInc = printInc, dimVar = c("geo", "eu"), roundBase = 5)$difference, c(2, 0, 2, 2, 0, 0))
  

  mf2 <- ~region + hovedint + fylke * hovedint + kostragr * hovedint
  z2 = SmallCountData("z2")
  a <- PLSrounding( z2, "ant", formula = mf2, xReturn = TRUE, printInc = printInc)
  expect_equivalent(t(as.matrix(a$x)) %*% as.matrix(a$inner[, c("original", "rounded")]), as.matrix(a$publish[, c("original", "rounded")]))
  expect_equivalent(sum(a$publish[, "rounded"] == 2), 0)
  expect_true(a$inner[42, "rounded"] == 2)
  a <- PLSrounding( z2, "ant", formula = mf2, leverageCheck = TRUE, printInc = printInc)
  expect_false(a$inner[42, "rounded"] == 2)
  a <- PLSrounding( z2, "ant", formula = mf2, leverageCheck = 0.9999999, printInc = printInc)
  expect_false(a$inner[42, "rounded"] == 2)
  a <- PLSrounding( z2, "ant", formula = mf2, leverageCheck = 1.1, printInc = printInc)
  expect_true(a$inner[42, "rounded"] == 2)
  
  z <- z2[-c(1,3,7,11,13,17), ]
  a0 <- PLSrounding( z, "ant", printInc = printInc, removeEmpty=FALSE)
  a1 <- PLSrounding( z, "ant", printInc = printInc, removeEmpty=TRUE)
  a2 <- PLSrounding( z, "ant", printInc = printInc, formula = ~region * hovedint + fylke * hovedint + kostragr * hovedint)
  expect_equivalent(a1$freqTable,a2$freqTable)
  expect_false(a0$freqTable[3,10]==a1$freqTable[3,10])
  
  
  mf3 <- ~region*mnd + region*hovedint + fylke*hovedint*mnd + kostragr*hovedint*mnd
  z = SmallCountData("z3")
  a <- PLSrounding(z, "ant", 50, formula = mf3, easyCheck = FALSE, printInc = printInc)
  z$ant2 <- a$inner$rounded
  b0 <- PLSrounding(z, "ant2", 50, formula = mf3, easyCheck = FALSE, printInc = printInc)
  b1 <- PLSrounding(z, "ant2", 50, formula = mf3, printInc = printInc)
  expect_true(b0$metrics["maxdiff"]==0)
  expect_false(b1$metrics["maxdiff"]==0)
  b2 <- PLSrounding(z, "ant2", 50, formula = mf3, leverageCheck = TRUE, printInc = printInc)
  expect_identical(b1,b2)
  
  z <- z[z$ant>0, ]
  dL <- FindDimLists(z[,-c(3,6,7)])
  a0 <- PLSrounding( z, "ant", hierarchies= dL, formula = ~region*hovedint*mnd-region:hovedint:mnd, printInc = printInc, removeEmpty=FALSE)
  a1 <- PLSrounding( z, "ant", hierarchies= dL, formula = ~region*hovedint*mnd-region:hovedint:mnd, printInc = printInc, removeEmpty=TRUE)
  expect_false(a0$freqTable[1,6]==0)
  expect_true(a1$freqTable[1,6]==0)
  
  exPSD <- SmallCountData("exPSD")
  a <- PLSrounding(exPSD, "freq", 5, formula = ~rows + cols, printInc = printInc)
  expect_equivalent(a$publish$rounded, c(28, 15, 8, 5, 7, 5, 5, 5, 6))
  
  a <- PLSrounding(exPSD, "freq", 5, formula = ~rows + cols, identifyNew = FALSE, printInc = printInc)
  expect_equivalent(a$publish$rounded, c(27, 16, 6, 5, 7, 5, 4, 5, 6))
  
  a <- PLSrounding(exPSD, "freq", 5, formula = ~rows + cols, maxRound = 7, printInc = printInc)
  expect_equivalent(a$inner$rounded, c(5, 0, 0, 0, 0, 5, 0, 5, 0, 5, 0, 0, 4, 2, 0))
  
  a <- PLSrounding(exPSD, "freq", 5, formula = ~rows + cols, zeroCandidates = TRUE, printInc = printInc)
  expect_equivalent(a$inner$rounded, c(6, 1, 0, 0, 5, 0, 5, 0, 0, 0, 0, 5, 4, 2, 0))
  
})



test_that("preAggregate works", {
  printInc <- FALSE
  
  formula <- ~eu * year + geo * year
  zM <- SSBtools::MakeMicro(SmallCountData("e6"), "freq")[, -4]
  zF <- aggregate(list(freq = zM[[1]]), zM, length)
  
  aM <- PLSrounding(zM, formula = formula, printInc = printInc)
  aF <- PLSrounding(zF, "freq", formula = formula, printInc = printInc)
  
  expect_equal(nrow(aM[[1]]), nrow(zF))
  expect_equal(diff(range(diff(sort(SSBtools::Match(aM[[1]], aF[[1]]))))), 0)
  expect_identical(aM[[2]], aF[[2]][names(aM[[2]])])
  expect_identical(aM[3:4], aF[3:4])
  
  aM <- PLSrounding(zM, printInc = printInc)
  expect_equal(diff(range(diff(sort(SSBtools::Match(aM[[1]], aF[[1]]))))), 0)
  expect_equal(diff(range(diff(sort(SSBtools::Match(aM[[2]], aF[[2]]))))), 0)
  expect_identical(aM[3:4], aF[3:4])
  
  aM <- PLSrounding(zM, hierarchies = SmallCountData("eDimList"), printInc = printInc)
  expect_equal(diff(range(diff(sort(SSBtools::Match(aM[[1]], aF[[1]]))))), 0)
  expect_equal(diff(range(diff(sort(SSBtools::Match(aM[[2]], aF[[2]]))))), 0)
  expect_identical(aM[3:4], aF[3:4])
  
  aM <- PLSrounding(zM, formula = formula, printInc = printInc,  preAggregate = FALSE)
  expect_equal(nrow(aM[[1]]), nrow(zM))
  expect_equal(diff(range(diff(sort(SSBtools::Match(aM[[2]], aF[[2]]))))), 0)

})

test_that("Parameter preRounded", {
  mf2 <- ~region + hovedint + fylke * hovedint + kostragr * hovedint
  z2 <- SmallCountData("z2")
  z2$freq <- c(2, 2, 0, 2, 0, 1, 2, 2, 2, 0, 2, 0, 0, 0, 2, 0, 0, 1, 0, 2, 0, 0, 2, 0, 0, 0, 1, 0, 2, 1, 2, 1, 2, 0, 0, 1, 1, 2, 0, 1, 2, 1, 0, 0)
  
  a <- PLSrounding(z2, "ant", formula = mf2, easyCheck = TRUE, leverageCheck = FALSE)
  z2$antR <- a$inner$rounded
  expect_warning(b <- PLSrounding(z2, "ant", formula = mf2, easyCheck = TRUE, leverageCheck = TRUE, preRounded = "antR"))
  expect_identical(a, b)
  
  a <- PLSrounding(z2, "freq", formula = mf2, easyCheck = FALSE, leverageCheck = FALSE)
  z2$freqR <- a$inner$rounded
  expect_warning(b <- PLSrounding(z2, "freq", formula = mf2, easyCheck = TRUE, leverageCheck = FALSE, preRounded = "freqR"))
  expect_identical(a, b)
  z2$freqR[a$inner$difference < 0] <- NA
  b <- PLSrounding(z2, "freq", formula = mf2, easyCheck = FALSE, leverageCheck = FALSE, preRounded = "freqR")
  expect_identical(a, b)
  b <- PLSrounding(z2, "freq", formula = mf2, easyCheck = FALSE, leverageCheck = FALSE, preRounded =  z2$freqR, maxIterRows = 3)
  expect_identical(a, b)
  b <- PLSrounding(z2, "freq", formula = mf2, easyCheck = FALSE, leverageCheck = FALSE, preRounded = "freqR", maxIterRows = 1)
  expect_identical(a, b)
  
  z2$freqR <- NA
  z2$freqR[1] <- 18
  b <- PLSrounding(z2, "freq", formula = mf2, easyCheck = FALSE, leverageCheck = FALSE, preRounded = "freqR")
  expect_true(sum(b$inner$difference) < 2)
  b <- PLSrounding(z2, "freq", formula = mf2, easyCheck = FALSE, leverageCheck = FALSE, preRounded = "freqR", maxIterRows = 2)
  expect_true(sum(b$inner$difference) < 2)
  b <- PLSrounding(z2, "freq", formula = mf2, easyCheck = FALSE, leverageCheck = FALSE, preRounded =  z2$freqR, maxIterRows = 1)
  expect_true(sum(b$inner$difference) < 2)
  
  z2$freqR[1] <- 1000
  b <- PLSrounding(z2, "freq", formula = mf2, easyCheck = FALSE, leverageCheck = FALSE, preRounded = "freqR")
  expect_true(sum(b$inner$difference > 0) == 1)
  
  b <- PLSrounding(z2, "freq", formula = mf2, easyCheck = FALSE, leverageCheck = FALSE, preRounded = "freqR", forceInner = TRUE)
  expect_true(sum(b$inner$rounded) == 1000)
  
  z2$freq[1] <- 1000
  z2$freqR[1] <- 100
  b <- PLSrounding(z2, "freq", formula = mf2, easyCheck = FALSE, leverageCheck = FALSE, preRounded = "freqR")
  expect_true(sum(b$inner$difference < 0) == 1)
  
  b <- PLSrounding(z2, "freq", formula = mf2, easyCheck = FALSE, leverageCheck = FALSE, preRounded = "freqR", forceInner = TRUE)
  expect_equivalent(unique(b$inner$rounded), c(100, 3, 0))
  
  b <- PLSrounding(z2, "freq", formula = mf2, easyCheck = FALSE, leverageCheck = FALSE, preRounded = "freqR", forceInner = TRUE, zeroCandidates = TRUE)
  expect_equivalent(unique(b$inner$rounded), c(100, 3))
})



PLStest = function(..., seed, Version){
  set.seed(seed)
  capture.output({ a <- PLSrounding(..., Version = Version, rndSeed = NULL)})
  set.seed(seed)
  b <-PLSrounding(..., printInc = printInc, rndSeed = NULL)
  expect_identical(a,b)
}


test_that("Same as Version_0.3.0", {
  skip_on_cran()
  seed = 123
  mf <- ~region*mnd + hovedint*mnd + fylke*hovedint*mnd + kostragr*hovedint*mnd
  PLStest(SmallCountData('z3'), 'ant', 3, formula = mf, seed= seed, Version = "0.3.0")
  # PLSrounding(SmallCountData('z3'), 'ant', 5, formula = mf, seed= seed, Version = "0.3.0", maxIterRows = 30)
  PLStest(SmallCountData('z3'), 'ant', 7, formula = mf, seed= seed, Version = "0.3.0", singleRandom = TRUE)
  mf <- ~region*mnd + hovedint*mnd + fylke*hovedint*mnd
  PLStest(SmallCountData('z3'), 'ant', 10, formula = mf, seed= seed, Version = "0.3.0")
  PLStest(SmallCountData('z3'), 'ant', 5, seed= seed, Version = "0.3.0")
})



test_that("Same as Version_0.3.0 many tests", {
  skip("Too many tests") 
  seed = 123
  mf <- ~region*mnd + hovedint*mnd + fylke*hovedint*mnd + kostragr*hovedint*mnd
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 3, formula = mf, seed= seed, Version = "0.3.0")
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 4, formula = mf, seed= seed, Version = "0.3.0")
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 5, formula = mf, seed= seed, Version = "0.3.0")
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 7, formula = mf, seed= seed, Version = "0.3.0")
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 10, formula = mf, seed= seed, Version = "0.3.0")
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 20, formula = mf, seed= seed, Version = "0.3.0",  maxIterRows = 10000)
  
  mf <- ~region*mnd + hovedint*mnd + fylke*hovedint*mnd
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 5, formula = mf, seed= seed, Version = "0.3.0")
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 10, formula = mf, seed= seed, Version = "0.3.0")
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 20, formula = mf, seed= seed, Version = "0.3.0",  maxIterRows = 10000)
  
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 5, seed= seed, Version = "0.3.0")
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 20, seed= seed, Version = "0.3.0")
  
  mf <- ~region*mnd + hovedint*mnd + fylke*hovedint*mnd + kostragr*hovedint*mnd
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 3, formula = mf, seed= seed, Version = "0.3.0", singleRandom = TRUE)
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 4, formula = mf, seed= seed, Version = "0.3.0", singleRandom = TRUE)
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 5, formula = mf, seed= seed, Version = "0.3.0", singleRandom = TRUE)
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 7, formula = mf, seed= seed, Version = "0.3.0", singleRandom = TRUE)
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 10, formula = mf, seed= seed, Version = "0.3.0", singleRandom = TRUE)
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 20, formula = mf, seed= seed, Version = "0.3.0", singleRandom = TRUE)
  
  mf <- ~region*mnd + hovedint*mnd + fylke*hovedint*mnd
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 5, formula = mf, seed= seed, Version = "0.3.0", singleRandom = TRUE)
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 10, formula = mf, seed= seed, Version = "0.3.0", singleRandom = TRUE)
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 20, formula = mf, seed= seed, Version = "0.3.0", singleRandom = TRUE)
  
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 5, seed= seed, Version = "0.3.0", singleRandom = TRUE)
  PLStest(SmallCountData('sosialFiktiv'), 'ant', 20, seed= seed, Version = "0.3.0", singleRandom = TRUE)
})
