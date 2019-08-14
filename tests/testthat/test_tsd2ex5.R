
context("TSD2 Example 5 (Parkinson's)")

tc <- textConnection(NULL, "w")
sink(tc)

library(BUGSnet)

set.seed(BUGSnet:::tsd2ex5$bugsnet$seed[1])

dat <- data.prep(BUGSnet:::tsd2ex5$data, varname.t = "Treatment", varname.s = "Study")

#TSD2 Example 5 Fixed Effects Model
model_fe <- nma.model(dat, outcome = "y", N = "n", sd = "sd", reference = "1", family = "normal",
                      link = "identity", effects = "fixed")
results_fe <- nma.run(model_fe, n.adapt = 5000, n.burnin = 50000, n.iter = 100000)
s_fe <- summary(results_fe$samples[,2:5])
tbl_fe <- cbind(s_fe$statistics[,1:2], s_fe$quantiles[,c(3,1,5)])

#TSD2 Example 5 Random Effects Model
model_re <- nma.model(dat, outcome = "y", N = "n", sd = "sd", reference = "1", family = "normal",
                      link = "identity", effects = "random")
results_re <- nma.run(model_re, n.adapt = 5000, n.burnin = 50000, n.iter = 100000)
s_re <- summary(results_re$samples[,2:5])
tbl_re <- cbind(s_re$statistics[,1:2], s_re$quantiles[,c(3,1,5)])

#round and append tables
results <- as.data.frame(round(rbind(tbl_fe, tbl_re), 2))
rownames(results) <- NULL

sink()
close(tc)

benchmark <- BUGSnet:::tsd2ex5$bugsnet[,3:7]
rownames(benchmark) <- NULL

test_that("nma.run results match benchmark", { expect_equal(benchmark, results) })
