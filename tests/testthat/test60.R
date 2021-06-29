# Reproduce and Test Fix #60

context("Meta Regression Continuous Only")

data(afib)

my_afib <- afib

my_afib$stroke <- "Decimal"

View(my_afib)

dataprep <- data.prep(arm.data = my_afib,
                      varname.t = "treatment",
                      varname.s = "study")

str(dataprep)


test_that("Non continuous covariate produces error in meta regression", 
          {expect_error(
            nma.model(data=dataprep,
                                  outcome="events",
                                  N="sampleSize",
                                  reference="02",
                                  family="binomial",
                                  link="logit",
                                  effects="fixed",
                                  covariate="stroke",
                                  prior.beta="EXCHANGEABLE"
                      ), "Error in nma.model(data = dataprep, outcome = \"events\", N = \"sampleSize\",  : \n  Covariate should have more than one unique value"
          )})


nma.model(data=dataprep,
          outcome="events",
          N="sampleSize",
          reference="02",
          family="binomial",
          link="logit",
          effects="fixed",
          covariate="stroke",
          prior.beta="EXCHANGEABLE")
          