library(BUGSnet)

data(diabetes.sim)

rate2.slr <- data.prep(arm.data = diabetes.sim,
                       varname.t = "Treatment",
                       varname.s = "Study")

png("netplot.png", width = 30, height = 20, units = "cm", res = 600, pointsize = 14)
## ---- net.plot, echo=TRUE, fig.width=6, fig.height=5, fig.align="centre"----
#par(mar=c(1,1,1,1)) #reduce margins for nicer plotting
par(mfrow=c(1,2))
net.plot(rate2.slr, outcome = "diabetes", node.scale = 3, edge.scale=1.5)

## ---- echo=TRUE, fig.width=6, fig.height=5, fig.align="centre"-----------
#par(mar=c(1,1,1,1)) #reduce margins for nicer plotting
net.plot(rate2.slr, outcome = "diabetes", node.scale = 3, edge.scale=1.5, flag = "ARB")
graphics.off()

## ---- results = "hide"---------------------------------------------------
network.char <- net.tab(data = rate2.slr,
                        outcome = "diabetes",
                        N = "n",
                        type.outcome = "rate2",
                        time = "followup")

## ---- echo = FALSE-------------------------------------------------------
network.char$network

## ---- echo=FALSE---------------------------------------------------------
network.char$intervention

## ---- echo=FALSE---------------------------------------------------------
network.char$comparison

## ---- fig.width=8, fig.height = 5, fig.align="centre"--------------------
png("dataplot.png", width = 18, height = 14, units = "cm", res = 600, pointsize = 14)
data.plot(data = rate2.slr,
          covariate = "age", 
          half.length = "age_SD",
          by = "treatment",
          fill.str = "age_type", #which value is reported (e.g mean/median)
          avg.hline=TRUE,
          text.size=60) #add overall average line?
graphics.off()

## ------------------------------------------------------------------------
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

## ---- results = "hide"---------------------------------------------------
fixed_effects_results <- nma.run(fixed_effects_model,
                                 n.adapt=10000,
                                 n.burnin=50000,
                                 n.iter=100000)

random_effects_results <- nma.run(random_effects_model,
                                  n.adapt=10000,
                                  n.burnin=50000,
                                  n.iter=100000)


## ---- fig.width=7, fig.height=4, results = "hide"------------------------
png("fitplot.png", width = 25, height = 15, units = "cm", res = 600, pointsize = 16)
par(mfrow = c(1,2))
fe_model_fit <- nma.fit(fixed_effects_results, main = "Fixed Effects Model" )
re_model_fit <- nma.fit(random_effects_results, main= "Random Effects Model")
graphics.off()

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


## ---- fig.width=7, fig.height=4, results = "hide"------------------------
par(mfrow = c(1,2))
re_model_fit <- nma.fit(random_effects_results, main = "Consistency Model" )
inconsist_model_fit <- nma.fit(re_inconsistency_results, main= "Inconsistency Model")


## ---- fig.height=6, fig.width = 6----------------------------------------
png("compareplot.png", width = 14, height = 14, units = "cm", res = 600, pointsize = 12)
nma.compare(re_model_fit, inconsist_model_fit)
graphics.off()

## ---- echo=TRUE, fig.width=7, fig.height=4, dpi = 95---------------------
sucra.out <- nma.rank(random_effects_results, largerbetter=FALSE, sucra.palette= "Set1")
sucra.out$rankogram

## ---- echo = FALSE-------------------------------------------------------
sucra.out$ranktable

## ---- echo=FALSE, fig.width=7, fig.height=4, dpi = 95, fig.align="centre"----
gg1 <- sucra.out$sucraplot
ggsave(filename = "sucraplot.png", plot = gg1, width = 12, height = 10, dpi = 300, units = "cm")

## ---- echo = FALSE-------------------------------------------------------
sucra.out$sucratable

## ---- echo=FALSE, fig.width=7, fig.height=4, dpi = 95, fig.align="centre"----
sucra.out$sucraplot + theme_gray()

## ---- echo=TRUE, fig.width=7, fig.height=4, fig.align="centre"-----------
league.out <- nma.league(random_effects_results,  
                         central.tdcy="median",
                         order = sucra.out$order,
                         log.scale = FALSE,
                         low.colour = "springgreen4",
                         mid.colour = "white",
                         high.colour = "red")
gg2 <- league.out$heatplot
ggsave(filename = "heatplot.png", plot = gg2, width = 12, height = 10, dpi = 300, units = "cm")

## ---- echo = FALSE-------------------------------------------------------
league.out$table

## ---- echo=TRUE, fig.width=7, fig.height=4, fig.align="centre"-----------
nma.forest(random_effects_results,
           central.tdcy="median",
           comparator = "Placebo",
           log.scale = FALSE)

## ------------------------------------------------------------------------
cat(random_effects_model$bugs)

## ---- fig.height = 15, fig.width=8---------------------------------------
nma.trace(random_effects_results)
