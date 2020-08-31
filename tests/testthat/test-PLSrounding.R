test_that("PLSrounding works", {
  # Small example data set
  z <- SmallCountData("e6")
  
  a <- PLSrounding(z, "freq")
  expect_equivalent(a$metrics, PLSrounding(z, "freq", formula = ~eu * year + geo * year)$metrics)
  expect_equivalent(a$metrics, PLSrounding(z[, -2], "freq", hierarchies = SmallCountData("eHrc"))$metrics)
  expect_equivalent(a$metrics, PLSrounding(z[, -2], "freq", hierarchies = SmallCountData("eDimList"))$metrics)
  expect_equivalent(a$metrics, PLSrounding(z[, -2], "freq", hierarchies = SmallCountData("eDimList"), formula = ~geo * year)$metrics)
  
  expect_equivalent(PLSrounding(z[, -2], "freq", hierarchies = SmallCountData("eDimList"), formula = ~geo + year)$metrics["maxdiff"], 0)

  mf2 <- ~region + hovedint + fylke * hovedint + kostragr * hovedint
  a <- PLSrounding(SmallCountData("z2"), "ant", formula = mf2, xReturn = TRUE)
  expect_equivalent(t(as.matrix(a$x)) %*% as.matrix(a$inner[, c("original", "rounded")]), as.matrix(a$publish[, c("original", "rounded")]))
  expect_equivalent(sum(a$publish[, "rounded"] == 2), 0)
  expect_true(a$inner[42, "rounded"] == 2)
  a <- PLSrounding(SmallCountData("z2"), "ant", formula = mf2, leverageCheck = TRUE)
  expect_false(a$inner[42, "rounded"] == 2)
  a <- PLSrounding(SmallCountData("z2"), "ant", formula = mf2, leverageCheck = 0.9999999)
  expect_false(a$inner[42, "rounded"] == 2)
  a <- PLSrounding(SmallCountData("z2"), "ant", formula = mf2, leverageCheck = 1.1)
  expect_true(a$inner[42, "rounded"] == 2)
  
  
  exPSD <- SmallCountData("exPSD")
  set.seed(12345)
  a <- PLSrounding(exPSD, "freq", 5, formula = ~rows + cols)
  expect_equivalent(a$publish$rounded, c(28, 15, 8, 5, 7, 5, 5, 5, 6))
  
  set.seed(12345)
  a <- PLSrounding(exPSD, "freq", 5, formula = ~rows + cols, identifyNew = FALSE)
  expect_equivalent(a$publish$rounded, c(27, 16, 6, 5, 7, 5, 4, 5, 6))
  
  set.seed(12345)
  a <- PLSrounding(exPSD, "freq", 5, formula = ~rows + cols, maxRound = 7)
  expect_equivalent(a$inner$rounded, c(5, 0, 0, 0, 0, 5, 0, 5, 0, 5, 0, 0, 4, 2, 0))
  
  set.seed(12345)
  a <- PLSrounding(exPSD, "freq", 5, formula = ~rows + cols, zeroCandidates = TRUE)
  expect_equivalent(a$inner$rounded, c(6, 1, 0, 0, 5, 0, 5, 0, 0, 0, 0, 5, 4, 2, 0))
  
})



PLStest = function(..., seed, Version){
  set.seed(seed)
  a = PLSrounding(..., Version = Version)
  set.seed(seed)
  b = PLSrounding(...)
  expect_identical(a,b)
}


test_that("Same as Version_0.3.0", {
  seed = 123
  mf <- ~region*mnd + hovedint*mnd + fylke*hovedint*mnd + kostragr*hovedint*mnd
  PLStest(SmallCountData('z3'), 'ant', 3, formula = mf, seed= seed, Version = "0.3.0")
  PLSrounding(SmallCountData('z3'), 'ant', 5, formula = mf, seed= seed, Version = "0.3.0", maxIterRows = 30)
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
