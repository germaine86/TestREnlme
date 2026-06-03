```r
################################################################################
## Theophylline example: full random-effects model
##
## This script reproduces the Theophylline analysis in which all three
## Similar approach was used for skill acqisition data
## pharmacokinetic parameters are allowed to have corresponding random effects.
##
## Model:
##   C_i  = exp(ai1)  : clearance-related parameter
##   Ka_i = exp(ai2)  : absorption rate
##   Ke_i = exp(ai3)  : elimination rate
##
## Random-effects structure:
##   ai1 = B1 + bi1
##   ai2 = B2 + bi2
##   ai3 = B3 + bi3
##
## The script estimates variance components using VLS, MM, and MMF, and then
## performs permutation tests for:
##   1. all random effects jointly;
##   2. bi1 and bi3 jointly, conditional on retaining bi2;
##   3. bi3 only, conditional on retaining bi1 and bi2.
################################################################################

## Load data
data("Theoph", package = "datasets")
d <- as.data.frame(Theoph)

## Nonlinear one-compartment pharmacokinetic model
Expr <- conc ~ Dose * exp(ai2 + ai3 - ai1) *
  (exp(-Time * exp(ai3)) - exp(-Time * exp(ai2))) /
  (exp(ai2) - exp(ai3))

## Starting values for fixed-effect parameters
start <- c(ai1 = -3.22, ai2 = 0.47, ai3 = -2.45)

## Full random-effects structure
random <- c(
  "ai1 ~ B1 + bi1",
  "ai2 ~ B2 + bi2",
  "ai3 ~ B3 + bi3"
)

## Number of permutations used in the paper
nperm0 <- 200



################################################################################ 
## Load reproducibility code 
################################################################################ 
source("../Dmethod.R") 
source("../Dhypothesis_test.R")
source("../Utils.R")



################################################################################
## 1. VLS estimation and permutation tests
################################################################################

## Estimate variance components using VLS
DVLS <- Dmethod(
  data = d,
  Expr = Expr,
  group = "Subject",
  random = random,
  start = start,
  method = "VLS",
  verbose = 1
)

summary(DVLS)

## Test 1: H0: all random effects are zero
HVLS1 <- Dhypothesis_test(
  data = d,
  Expr = Expr,
  group = "Subject",
  random = random,
  start = start,
  Dhatt = DVLS,
  method = "VLS",
  nperm = nperm0,
  seed = 2,
  verbose = 1
)

## Test 2: H0: bi1 = bi3 = 0, conditional on retaining bi2
HVLS13 <- Dhypothesis_test(
  data = d,
  Expr = Expr,
  group = "Subject",
  random = random,
  start = start,
  bi_out = c("bi1", "bi3"),
  Dhatt = DVLS,
  method = "VLS",
  nperm = nperm0,
  seed = 2,
  verbose = 1
)

## Test 3: H0: bi3 = 0, conditional on retaining bi1 and bi2
HVLS3 <- Dhypothesis_test(
  data = d,
  Expr = Expr,
  group = "Subject",
  random = random,
  start = start,
  bi_out = "bi3",
  Dhatt = DVLS,
  method = "VLS",
  nperm = nperm0,
  seed = 2,
  verbose = 1
)


################################################################################
## 2. MM estimation and permutation tests
################################################################################

## MM and MMF share the same first-stage subject-specific NLS fits.
## We compute these once using MM_base() and reuse them below.
mb <- MM_base(
  data = d,
  Expr = Expr,
  group = "Subject",
  random = random,
  start = start,
  verbose = 1
)

## Estimate variance components using MM
DMM <- Dmethod(
  data = d,
  Expr = Expr,
  group = "Subject",
  random = random,
  start = start,
  method = "MM",
  MM_base_obj = mb,
  verbose = 0
)

summary(DMM)

## Test 1: H0: all random effects are zero
HMM1 <- Dhypothesis_test(
  data = d,
  Expr = Expr,
  group = "Subject",
  random = random,
  start = start,
  Dhatt = DMM,
  method = "MM",
  nperm = nperm0,
  seed = 2,
  verbose = 1
)

## Test 2: H0: bi1 = bi3 = 0, conditional on retaining bi2
HMM13 <- Dhypothesis_test(
  data = d,
  Expr = Expr,
  group = "Subject",
  random = random,
  start = start,
  bi_out = c("bi1", "bi3"),
  Dhatt = DMM,
  method = "MM",
  nperm = nperm0,
  seed = 2,
  verbose = 1
)

## Test 3: H0: bi3 = 0, conditional on retaining bi1 and bi2
HMM3 <- Dhypothesis_test(
  data = d,
  Expr = Expr,
  group = "Subject",
  random = random,
  start = start,
  bi_out = "bi3",
  Dhatt = DMM,
  method = "MM",
  nperm = nperm0,
  seed = 2,
  verbose = 1
)


################################################################################
## 3. MMF estimation and permutation tests
################################################################################

## Estimate variance components using MMF
DMMF <- Dmethod(
  data = d,
  Expr = Expr,
  group = "Subject",
  random = random,
  start = start,
  method = "MMF",
  MM_base_obj = mb,
  verbose = 0
)

summary(DMMF)

## Test 1: H0: all random effects are zero
HMMF1 <- Dhypothesis_test(
  data = d,
  Expr = Expr,
  group = "Subject",
  random = random,
  start = start,
  Dhatt = DMMF,
  method = "MMF",
  nperm = nperm0,
  seed = 2,
  verbose = 1
)

## Test 2: H0: bi1 = bi3 = 0, conditional on retaining bi2
HMMF13 <- Dhypothesis_test(
  data = d,
  Expr = Expr,
  group = "Subject",
  random = random,
  start = start,
  bi_out = c("bi1", "bi3"),
  Dhatt = DMMF,
  method = "MMF",
  nperm = nperm0,
  seed = 2,
  verbose = 1
)

## Test 3: H0: bi3 = 0, conditional on retaining bi1 and bi2
HMMF3 <- Dhypothesis_test(
  data = d,
  Expr = Expr,
  group = "Subject",
  random = random,
  start = start,
  bi_out = "bi3",
  Dhatt = DMMF,
  method = "MMF",
  nperm = nperm0,
  seed = 2,
  verbose = 1
)


################################################################################
## 4. Summary of hypothesis-test results
################################################################################

list(
  VLS = list(
    all_random_effects = HVLS1$pval,
    bi1_bi3_given_bi2 = HVLS13$pval,
    bi3_given_bi1_bi2 = HVLS3$pval
  ),
  MM = list(
    all_random_effects = HMM1$pval,
    bi1_bi3_given_bi2 = HMM13$pval,
    bi3_given_bi1_bi2 = HMM3$pval
  ),
  MMF = list(
    all_random_effects = HMMF1$pval,
    bi1_bi3_given_bi2 = HMMF13$pval,
    bi3_given_bi1_bi2 = HMMF3$pval
  )
)
```
