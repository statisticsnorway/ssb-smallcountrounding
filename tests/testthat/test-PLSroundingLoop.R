test_that("PLSroundingLoop  works", {
  printInc <- FALSE
  set.seed(12)
  mf <- ~mnd * (region + fylke * hovedint)
  z3 <- SmallCountData("z3")
  a <- PLSroundingLoop(z3, loopId = "kostragr", freqVar = "ant", formula = mf, seed = NULL, printInc = printInc)
  set.seed(12)
  b1 <- PLSroundingLoop(z3[z3$kostragr < 250, ], loopId = "kostragr", freqVar = "ant", formula = mf, seed = NULL, printInc = printInc)
  b <- PLSroundingLoop(z3[z3$kostragr > 250, ], loopId = "kostragr", freqVar = "ant", formula = mf, seed = NULL, preOutput = b1, printInc = printInc)
  expect_equal(a, b)
})
