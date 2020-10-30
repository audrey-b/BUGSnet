## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(BUGSnet)

data(afib)
dataprep <- data.prep(arm.data = afib,
                      varname.t = "treatment",
                      varname.s = "study")

## ---- net.plot, echo=TRUE, fig.width=5, fig.height=5--------------------------
 net.plot(dataprep, node.scale = 1, 
          edge.scale=1,
          label.offset1 = 2)

## ---- results = "hide"--------------------------------------------------------
network.char <- net.tab(data = dataprep,
                        outcome = "events",
                        N = "sampleSize",
                        type.outcome = "binomial")

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(network.char$network)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(network.char$intervention)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(network.char$comparison)

## -----------------------------------------------------------------------------
fixed_effects_model <- nma.model(data=dataprep,
                                  outcome="events",
                                  N="sampleSize",
                                  reference="02",
                                  family="binomial",
                                  link="logit",
                                  effects="fixed",
                                  covariate="stroke",
                                  prior.beta="EXCHANGEABLE")

random_effects_model <- nma.model(data=dataprep,
                                  outcome="events",
                                  N="sampleSize",
                                  reference="02",
                                  family="binomial",
                                  link="logit",
                                  effects="random",
                                  covariate="stroke",
                                  prior.beta="EXCHANGEABLE")


## ---- results = "hide"--------------------------------------------------------
fixed_effects_results <- nma.run(fixed_effects_model,
                           n.adapt=1000,
                           n.burnin=1000,
                           n.iter=5000)

random_effects_results <- nma.run(random_effects_model,
                           n.adapt=1000,
                           n.burnin=1000,
                           n.iter=5000)


## ---- fig.width=7, fig.height=4, results = "hide"-----------------------------
par(mfrow = c(1,2))
nma.fit(fixed_effects_results, main = "Fixed Effects Model" )
nma.fit(random_effects_results, main= "Random Effects Model")

## ---- echo=TRUE, fig.width=7, fig.height=4------------------------------------
nma.regplot(random_effects_results)

## ---- echo=TRUE, fig.width=7, fig.height=4, dpi = 95--------------------------
sucra.out <- nma.rank(random_effects_results, largerbetter=FALSE, cov.value=0.1, sucra.palette= "Set1")
sucra.out$sucraplot

## ---- echo=TRUE, fig.width=15, fig.height=10----------------------------------
league.out <- nma.league(random_effects_results, 
                             central.tdcy = "median",
                             order = as.vector(t(dataprep$treatments)),
                             cov.value=0.1, 
                             log.scale = FALSE)
league.out$heatplot

## ---- echo=TRUE, fig.width=10, fig.height=4-----------------------------------
nma.forest(random_effects_results, 
           comparator="02", 
           central.tdcy = "median",
           cov.value=0.1,
           x.trans="log")

