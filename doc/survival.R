## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include = FALSE---------------------------------------------------------
library(BUGSnet)

## ---- echo = FALSE------------------------------------------------------------
data(diabetes.sim)
knitr::kable(diabetes.sim)

## -----------------------------------------------------------------------------
library(BUGSnet)
data(diabetes.sim)
rate2.slr <- data.prep(arm.data = diabetes.sim,
                     varname.t = "Treatment",
                     varname.s = "Study")

## ---- net.plot, echo=TRUE, fig.width=6, fig.height=5, fig.align="centre"------
par(mar=c(1,1,1,1)) #reduce margins for nicer plotting
net.plot(rate2.slr, node.scale = 3, edge.scale=1.5)

## ---- echo=TRUE, fig.width=6, fig.height=5, fig.align="centre"----------------
par(mar=c(1,1,1,1)) #reduce margins for nicer plotting
net.plot(rate2.slr, node.scale = 3, edge.scale=1.5, flag = "ARB")

## ---- results = "hide"--------------------------------------------------------
network.char <- net.tab(data = rate2.slr,
                        outcome = "diabetes",
                        N = "n",
                        type.outcome = "rate2",
                        time = "followup")

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(network.char$network)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(network.char$intervention)

## ---- echo=FALSE--------------------------------------------------------------
knitr::kable(network.char$comparison)

## ---- fig.width=8, fig.height = 5, fig.align="centre"-------------------------
data.plot(data = rate2.slr,
              covariate = "age", 
              half.length = "age_SD",
              by = "treatment",
              avg.hline=TRUE, #add overall average line?
              text.size = 12) 

## -----------------------------------------------------------------------------
random_effects_model <- nma.model(data=rate2.slr,
                                 outcome="diabetes",
                                 N="n",
                                 reference="Diuretic",
                                 family="binomial",
                                 link="cloglog",
                                 time = "followup",
                                 effects= "random")

fixed_effects_model <- nma.model(data=rate2.slr,
                                 outcome="diabetes",
                                 N="n",
                                 reference="Diuretic",
                                 family="binomial",
                                 link="cloglog",
                                 time = "followup",
                                 effects="fixed")

## ---- results = "hide"--------------------------------------------------------
set.seed(20190828)
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
fe_model_fit <- nma.fit(fixed_effects_results, main = "Fixed Effects Model" )
re_model_fit <- nma.fit(random_effects_results, main= "Random Effects Model")

## ---- results = "hide", fig.show = 'hide',  fig.width=8, fig.height = 8, fig.align="centre"----
re_inconsistency_model <- nma.model(data=rate2.slr,
                                   outcome="diabetes",
                                   N="n",
                                   reference="Diuretic",
                                   family="binomial",
                                   link="cloglog",
                                   time = "followup",
                                   type = "inconsistency", #specifies inconsistency model
                                   effects="random")

re_inconsistency_results <- nma.run(re_inconsistency_model,
                                         n.adapt=1000,
                                         n.burnin=1000,
                                         n.iter=10000)


## ---- fig.width=7, fig.height=4, results = "hide"-----------------------------
par(mfrow = c(1,2))
re_model_fit <- nma.fit(random_effects_results, main = "Consistency Model" )
inconsist_model_fit <- nma.fit(re_inconsistency_results, main= "Inconsistency Model")

## ---- fig.height=6, fig.width = 6---------------------------------------------
nma.compare(re_model_fit, inconsist_model_fit)

## ---- echo=TRUE, fig.width=7, fig.height=4, dpi = 95--------------------------
sucra.out <- nma.rank(random_effects_results, largerbetter=FALSE, sucra.palette= "Set1")
sucra.out$rankogram

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(sucra.out$ranktable)

## ---- echo=FALSE, fig.width=7, fig.height=4, dpi = 95, fig.align="centre"-----
sucra.out$sucraplot

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(sucra.out$sucratable)

## ---- echo=FALSE, fig.width=7, fig.height=4, dpi = 95, fig.align="centre"-----
sucra.out$sucraplot + ggplot2::theme_gray()

## ---- echo=TRUE, fig.width=7, fig.height=4, fig.align="centre"----------------
league.out <- nma.league(random_effects_results,  
                         central.tdcy="median",
                         order = sucra.out$order,
                         log.scale = FALSE,
                         low.colour = "springgreen4",
                         mid.colour = "white",
                         high.colour = "red",
                         digits = 2)
league.out$heatplot

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(league.out$table)

## ---- echo=TRUE, fig.width=7, fig.height=4, fig.align="centre"----------------
nma.forest(random_effects_results,
           central.tdcy="median",
           comparator = "Placebo",
           log.scale = FALSE)

## -----------------------------------------------------------------------------
cat(random_effects_model$bugs)

## ---- fig.height = 15, fig.width=8--------------------------------------------
nma.diag(random_effects_results, plot_prompt = FALSE)

