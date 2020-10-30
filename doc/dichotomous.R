## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")


## -----------------------------------------------------------------------------
library(BUGSnet)

data(thrombolytic)
age <- rnorm(nrow(thrombolytic), 60, 10)
dich.slr <- data.prep(arm.data = cbind(thrombolytic, age),
                     varname.t = "treatment",
                     varname.s = "study")

## ---- net.plot, echo=TRUE, fig.width=6, fig.height=6--------------------------
 net.plot(dich.slr, node.scale = 1.5, edge.scale=0.5)

## ---- results = "hide"--------------------------------------------------------
network.char <- net.tab(data = dich.slr,
                        outcome = "events",
                        N = "sampleSize", 
                        type.outcome = "binomial",
                        time = NULL)

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(network.char$network)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(network.char$intervention)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(network.char$comparison)

## ---- echo=FALSE--------------------------------------------------------------
#pairwise.dat <- by.comparison(data.nma=dich.slr, outcome="events", type.outcome="binomial", N="sampleSize")
 
#pairwise.dat.with.c <- pairwise.dat %>% 
#   filter(grepl("SK", comparison)) %>%
#   filter((trt.e == "tPA"| trt.c == "tPA"))
 
#  tPA_vs_SK <- pairwise(data.nma = dich.slr,
#                           name.trt1 = "SK", 
#                           name.trt2 = "tPA", 
#                           outcome = "events",
#                           N = "sampleSize")
 
# pairwise.all(data.nma=dich.slr, 
#             outcome="events",
#             N = "sampleSize",
#             method.tau="DL",
#             sm= "RR") 

## -----------------------------------------------------------------------------
fixed_effects_model <- nma.model(data=dich.slr,
                     outcome="events",
                     N="sampleSize",
                     reference="SK",
                     family="binomial",
                     link="log",
                     effects="fixed")

random_effects_model <- nma.model(data=dich.slr,
                     outcome="events",
                     N="sampleSize",
                     reference="SK",
                     family="binomial",
                     link="log",
                     effects="random")


## ---- results = "hide"--------------------------------------------------------
set.seed(20190829)
fixed_effects_results <- nma.run(fixed_effects_model,
                           n.adapt=1000,
                           n.burnin=1000,
                           n.iter=10000)

random_effects_results <- nma.run(random_effects_model,
                           n.adapt=1000,
                           n.burnin=1000,
                           n.iter=10000)


## ---- fig.width=7, fig.height=4, results = "hide"-----------------------------
par(mfrow = c(1,2))
nma.fit(fixed_effects_results, main = "Fixed Effects Model" )
nma.fit(random_effects_results, main= "Random Effects Model")


## ---- echo=TRUE, fig.width=7, fig.height=4, dpi = 95--------------------------
sucra.out <- nma.rank(fixed_effects_results, largerbetter=FALSE, sucra.palette= "Set1")
sucra.out$sucraplot

## ---- echo=TRUE, fig.width=7, fig.height=4------------------------------------
league.out <- nma.league(fixed_effects_results,  
                         central.tdcy="median",
                         order = sucra.out$order)
league.out$heatplot

## ---- echo=TRUE, fig.width=7, fig.height=4------------------------------------
nma.forest(fixed_effects_results,
           central.tdcy="median",
           comparator = "SK")

## ---- results = "hide", fig.show = 'hide',  fig.width=8, fig.height = 8-------
fe_inconsistency_model <- nma.model(data=dich.slr,
                     outcome="events",
                     N="sampleSize",
                     reference="SK",
                     family="binomial",
                     link="log",
                     type = "inconsistency",
                     effects="fixed")

fe_inconsistency_results <- nma.run(fe_inconsistency_model,
                           n.adapt=1000,
                           n.burnin=1000,
                           n.iter=10000)

## ---- fig.height=4, fig.width=7-----------------------------------------------
par(mfrow = c(1,2))
fe_model_fit <- nma.fit(fixed_effects_results)
inconsist_model_fit <- nma.fit(fe_inconsistency_results)

## ---- fig.height=6, fig.width = 6---------------------------------------------
nma.compare(fe_model_fit, inconsist_model_fit)

