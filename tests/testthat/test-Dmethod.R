test_that("Dmethod (VLS) runs on theophylline data and returns expected structure", {
  d <- as.data.frame(Theoph)
  Expr <- conc ~ Dose * exp(ai2 + ai3 - ai1) *
    (exp(-Time * exp(ai3)) - exp(-Time * exp(ai2))) /
    (exp(ai2) - exp(ai3))
  start  <- c(ai1 = -3.22, ai2 = 0.47, ai3 = -2.45)
  random <- c("ai1 ~ B1 + bi1", "ai2 ~ B2 + bi2", "ai3 ~ B3 + bi3")

  DVLS <- Dmethod(d, Expr, group = "Subject",
                  random = random, start = start, method = "VLS",
                  verbose = 0)

  expect_s3_class(DVLS, "Dmethod")
  expect_true(all(c("Dhat", "Sigma2", "Beta", "method", "R",
                    "Ypred", "residuals") %in% names(DVLS)))
  expect_equal(DVLS$method, "VLS")
  expect_true(is.matrix(DVLS$Dhat))
  expect_equal(dim(DVLS$Dhat), c(3, 3))
  expect_true(DVLS$Sigma2 > 0)
  expect_true(!is.null(attr(DVLS, "internal")))
  expect_true(all(c("data", "Expr", "group", "random", "start") %in%
                    names(attr(DVLS, "internal"))))
})

test_that("Tstat computes a finite scalar from a Dmethod object", {
  d <- as.data.frame(Theoph)
  Expr <- conc ~ Dose * exp(ai2 + ai3 - ai1) *
    (exp(-Time * exp(ai3)) - exp(-Time * exp(ai2))) /
    (exp(ai2) - exp(ai3))
  start  <- c(ai1 = -3.22, ai2 = 0.47, ai3 = -2.45)
  random <- c("ai1 ~ B1 + bi1", "ai2 ~ B2 + bi2", "ai3 ~ B3 + bi3")

  DVLS <- Dmethod(d, Expr, group = "Subject",
                  random = random, start = start, method = "VLS",
                  verbose = 0)
  Tobs <- Tstat(DVLS)

  expect_true(is.numeric(Tobs))
  expect_length(Tobs, 1)
  expect_true(is.finite(Tobs))
  expect_true(Tobs >= 0)
})

test_that("bootstrap_se runs on a small theophylline fit and returns SE tables", {
  d <- as.data.frame(Theoph)
  Expr <- conc ~ Dose * exp(ai2 + ai3 - ai1) *
    (exp(-Time * exp(ai3)) - exp(-Time * exp(ai2))) /
    (exp(ai2) - exp(ai3))
  start  <- c(ai1 = -3.22, ai2 = 0.47, ai3 = -2.45)
  random <- c("ai1 ~ B1 + bi1", "ai2 ~ B2 + bi2", "ai3 ~ B3 + bi3")

  DVLS <- Dmethod(d, Expr, group = "Subject",
                  random = random, start = start, method = "VLS",
                  verbose = 0)

  boot <- bootstrap_se(DVLS, nboot = 5, type = "case", seed = 1, verbose = 0)

  expect_s3_class(boot, "Dboot")
  expect_true(all(c("Estimate", "SE") %in% colnames(boot$Beta)))
  expect_true(all(c("Estimate", "SE") %in% colnames(boot$Dhat)))
  expect_true(all(c("Estimate", "SE") %in% colnames(boot$Sigma2)))
})

test_that("Dhypothesis_test runs the joint test and returns a valid Dtest object", {
  d <- as.data.frame(Theoph)
  Expr <- conc ~ Dose * exp(ai2 + ai3 - ai1) *
    (exp(-Time * exp(ai3)) - exp(-Time * exp(ai2))) /
    (exp(ai2) - exp(ai3))
  start  <- c(ai1 = -3.22, ai2 = 0.47, ai3 = -2.45)
  random <- c("ai1 ~ B1 + bi1", "ai2 ~ B2 + bi2", "ai3 ~ B3 + bi3")

  DVLS <- Dmethod(d, Expr, group = "Subject",
                  random = random, start = start, method = "VLS",
                  verbose = 0)

  H_all <- Dhypothesis_test(d, Expr, group = "Subject",
                            random = random, start = start,
                            Dhatt = DVLS, nperm = 10, seed = 1,
                            verbose = 0)

  expect_s3_class(H_all, "Dtest")
  expect_true(all(c("Decision", "pvalue", "Tobs", "Tperm",
                    "Dhatt", "bi_out", "plot") %in% names(H_all)))
  expect_true(H_all$Decision %in% c("Reject H0", "Do not reject H0"))
  expect_true(is.numeric(attr(H_all, "internal")$pvalue_num))
  expect_true(attr(H_all, "internal")$pvalue_num >= 0 && attr(H_all, "internal")$pvalue_num <= 1)
  expect_true(is.finite(H_all$Tobs))
  expect_length(H_all$Tperm, 10)
  expect_null(H_all$bi_out)
})

test_that("Dhypothesis_test runs a subset test via bi_out", {
  d <- as.data.frame(Theoph)
  Expr <- conc ~ Dose * exp(ai2 + ai3 - ai1) *
    (exp(-Time * exp(ai3)) - exp(-Time * exp(ai2))) /
    (exp(ai2) - exp(ai3))
  start  <- c(ai1 = -3.22, ai2 = 0.47, ai3 = -2.45)
  random <- c("ai1 ~ B1 + bi1", "ai2 ~ B2 + bi2", "ai3 ~ B3 + bi3")

  DVLS <- Dmethod(d, Expr, group = "Subject",
                  random = random, start = start, method = "VLS",
                  verbose = 0)

  H_bi3 <- Dhypothesis_test(d, Expr, group = "Subject",
                            random = random, start = start,
                            Dhatt = DVLS, bi_out = "bi3",
                            nperm = 10, seed = 1, verbose = 0)

  expect_s3_class(H_bi3, "Dtest")
  expect_equal(H_bi3$bi_out, "bi3")
  expect_true(is.numeric(attr(H_bi3, "internal")$pvalue_num))
  expect_true(attr(H_bi3, "internal")$pvalue_num >= 0 && attr(H_bi3, "internal")$pvalue_num <= 1)
})

test_that("Dhypothesis_test computes Dhatt internally when not supplied", {
  d <- as.data.frame(Theoph)
  Expr <- conc ~ Dose * exp(ai2 + ai3 - ai1) *
    (exp(-Time * exp(ai3)) - exp(-Time * exp(ai2))) /
    (exp(ai2) - exp(ai3))
  start  <- c(ai1 = -3.22, ai2 = 0.47, ai3 = -2.45)
  random <- c("ai1 ~ B1 + bi1", "ai2 ~ B2 + bi2", "ai3 ~ B3 + bi3")

  H <- Dhypothesis_test(d, Expr, group = "Subject",
                        random = random, start = start,
                        nperm = 10, seed = 1, verbose = 0)

  expect_s3_class(H, "Dtest")
  expect_s3_class(H$Dhatt, "Dmethod")
  expect_true(is.numeric(attr(H, "internal")$pvalue_num))
})
