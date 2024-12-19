test_that("basic index computation works", {
  mat = exp(cbind(1:10, 2:11, 3:12))
  nd = data.frame(time = 1:10)
  ind1 = index(mat, timevar = "time", newdata = nd)
  expect_equal(ind1$index, exp(1:nrow(ind1)-1))
  expect_equal(ind1$ci_2.5, exp(1:nrow(ind1)-1))
  expect_equal(ind1$se, rep(0,nrow(ind1)))
  expect_equal(ind1$baseline, 1/3*rep(exp(1) + exp(2) + exp(3),nrow(ind1)))

  # Weights and repeat observations
  ind2 = index(rbind(mat, mat[8:10, ]), timevar = "time", newdata = rbind(nd, nd[8:10,, drop = FALSE]), weights = c(rep(1, 7), rep(.5, 6)))
  expect_equal(ind1$index, ind2$index) # TODO: Test whole data frame.

  # Shuffle columns
  ind2b = ind1 = index(mat[, sample(ncol(mat))], timevar = "time", newdata = nd)
  expect_identical(ind1$index, ind2b$index) # TODO: Test whole data frame.

  # Shuffle rows
  ord = c(4,  6,  8,  7,  2,  5,  1,  9,  3, 10)
  ind3 = index(mat[ord, ], timevar = "time", newdata = nd[ord,, drop=FALSE])
  expect_identical(ind1[,-3], ind3[sort(ord, index.return = TRUE)$ix, -3]) # TODO: Remove -3

  # Type = delta
  ind4 = index(mat, timevar = "time", newdata = nd, type = "delta")
  expect_identical(diff(ind4$index), rep(0, nrow(ind4) - 1))

  # Type = raw
  ind5 = index(mat, timevar = "time", newdata = nd, type = "raw")
  expect_equal(ind5$index, rowMeans(mat))
})

test_that("Group indices work.", {
  mat = cbind(matrix())
})

test_that("Zero abundances ok.", {
  mat = matrix(1, ncol = 6, nrow = 10)
  nd = data.frame(time = 1:10)
  mat[10,1:3] = 0
  ind1 = index(mat, timevar = "time", newdata = nd)
  expect_identical(ind1$index, c(rep(1,9), .5))
  expect_identical(ind1$se_log, c(rep(0,9), NaN))

  mat[10,] = 0
  ind2 = index(mat, timevar = "time", newdata = nd)
  expect_identical(ind2$index, c(rep(1,9),0))

  mat[1, 1:3] = 0
  ind3 = index(mat, timevar = "time", newdata = nd)
  expect_identical(ind3$index, c(NaN, rep(Inf, 8), NaN))
  expect_identical(ind3$ci_2.5, c(NA, rep(1, 8), NA))
  expect_identical(ind3$baseline, rep(0.5, 10))

  mat[1, ] = 0
  ind4 = index(mat, timevar = "time", newdata = nd)
  expect_identical(ind4$index, c(NaN, rep(Inf, 8), NaN))
  expect_identical(ind4$ci_2.5, c(NA, rep(Inf, 8), NA))
  expect_identical(ind4$baseline, rep(0, 10))

})

test_that("Single column works", {
  mat = matrix(1, nrow = 10, ncol = 1)
  nd = data.frame(time = 1:10)
  ind = index(mat, timevar = "time", newdata = nd)
  expect_identical(ind$index, rep(1, 10))
  expect_identical(ind$se_log, NA + numeric(10))
})

test_that("Unbalanced newdata is ok.", {
  #TODO
})

test_that("Test of gam.", {
  gam_fit = mgcv::gam(count ~ s(yr, lat), data = cuckoo, family = quasipoisson())
  #brms_fit = brms::brm(count ~ s(yr, lat), data = cuckoo, family = "negbinomial")
  #readRDS(paste0(test_path(), "/testdata/brms_test.rda"))
  nd = expand.grid(yr = 2000:2020, lat = c(55,60,65))
  preds = predict(gam_fit, newdata = nd, type = "response")
  indT = c(preds[1:21] / preds[1], preds[22:42] / preds[22], preds[43:63] / preds[43])
  names(indT) = NULL
  ind = index(gam_fit, timevar = "yr", byvar = "lat", newdata = nd)
  expect_equal(indT, ind$index)
  indT = c(preds / sum(preds[c(1,22,43)]))
  names(indT) = NULL
  ind = index(gam_fit, timevar = "yr", byvar = "lat", newdata = nd, type = "global")
  expect_equal(indT, ind$index)

})


test_that("Test of brms.", {
  # TODO
})


