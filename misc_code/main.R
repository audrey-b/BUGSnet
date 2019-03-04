#source("misc_code/load.packages.R")
#source("R/functions.R")
pathA <-"C:/Users/audre/Documents/Compile BUGSnet/BUGSnet/"
sepA <- "/"
pathJ <- "C:\\Users\\justi\\Desktop\\Lighthouse\\nmapackage\\BUGSnet\\"
sepJ <- "\\"
path <- pathA
sep <- sepA
source(paste0(path,"misc_code",sep,"load.packages.R"))
library(devtools)
library(roxygen2)
library(gemtc)
library(knitr)
library(readxl)
load_all(path = paste0(path,"R"))


rawdata <- read_excel(paste0(path,"data","\\","rate_example.xlsx"))
rawdata <- read_excel("data/continuous_example.xlsx", 
                      col_types = c("text", "numeric", "numeric", 
                                    "text", "text", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric"))
dich.slr <- data.prep(arm.data = rawdata,
                      varname.t = "trt_name",
                      varname.s = "trial")


# Import data --------------------------------------------------------
age <- rnorm(nrow(thrombolytic$data.ab), 60, 10)

dich.slr <- data.prep(arm.data = cbind(tbl_df(thrombolytic$data.ab),age),
                      varname.t = "treatment",
                      varname.s = "study")

#ensure that treatment and study variables are both of class character
dich.slr$arm.data$treatment <- as.character(dich.slr$arm.data$treatment)
dich.slr$arm.data$study <- as.character(dich.slr$arm.data$study)

pma(data = dich.slr,
    name.trt1 = "SK", 
    name.trt2 = "Ret",
    outcome = "responders",
    N = "sampleSize")


# Network Characteristics -------------------------------------------------

net.plot(dich.slr,  flag="PLCB", node.scale = 1.5, edge.scale=0.5, label.offset1 = 4, label.offset2 = 4)

network.char <- net.tab(data = dich.slr,
                        outcome = "responders",
                        N = "sampleSize",
                        type.outcome = "continuous",
                        time = NULL)

network.char$network
network.char$intervention
network.char$comparison

#NMA

fixed_effects_model <- nma.model(data=dich.slr,
                                 outcome="responders",
                                 N="sampleSize",
                                 reference="SK",
                                 type="inconsistency",
                                 family="binomial",
                                 link="logit",
                                 effects="fixed",
                                 covariate="age",
                                 prior.beta="UNRELATED")

random_effects_model <- nma.model(data=dich.slr,
                                  outcome="responders",
                                  N="sampleSize",
                                  type="inconsistency",
                                  reference="SK",
                                  family="binomial",
                                  link="logit",
                                  effects="random")

sink("Z:/ResearchDocuments/Research/BUGSnet/code.bug")
cat(random_effects_model$model)
sink()


random_effects_results <- nma.run(random_effects_model,
                                  monitor = c("d", "dev", "r", "n","totresdev","rhat","sigma","beta"),
                                  n.adapt=1000,
                                  n.burnin=1000,
                                  n.iter=10000)

random_effects_fit <- nma.fit(random_effects_results, main = "Random Effects Model" )
random_effects_fit$DIC
random_effects_fit$pD
random_effects_fit$pmdev

sucra.out <- nma.rank(random_effects_results, largerbetter=FALSE, cov.value=NULL)
sucra.out$sucraplot
sucra.out$rankogram

nma.forest(random_effects_results, 
           comparator="SK", 
           central.tdcy = "median")
nma.league(random_effects_results, 
           central.tdcy = "median",
           log.scale = TRUE)

nma.regplot(random_effects_results, x.range=c(38,84))

# Network Plots -----------------------------------------------------------




# 
# # Patient Characteristics plots--------------------------------------------
# source("baselines.plot.R")
# 
# #png(filename="data/test_pdfs", res=1000, height=7500, width=10000)
# comp_graph_comparison(x=slr$baseline.data, 
#                       trial.id = "trial", 
#                       treat.id = "trt", 
#                       measure.val = "age_estimate",
#                       measure.sd  = "age_SD",
#                       x.lim = c(25,65), 
#                       x.label = "Average Age (Years)" )
# #graphics.off()
# 
# #png(filename="data/test_pdfs", res=1000, height=7500, width=10000)
# comp_graph_trial(data=slr$baseline.data, 
#                  trial.id="trial",
#                  treat.id="trt",
#                  measure.val = "age_estimate", 
#                  x.lim = c(25,65), 
#                  x.label = "Average Age (Years)" )
# #graphics.off()

# Pairwise Comparisons ----------------------------------------------------

# NOTE: specify correct treatment: slr$treatments
# NOTE: specify correct oucome: colnames(slr$arm.data)
# 
# source("pairwise.R")
# png("output/pairwise.plot.png", height=720, width=1920, res=200)
# pairwise.output <- pairwise(slr,
#                             name.trt1 = "Placebo",
#                             name.trt2 = "Vernakalant IV",
#                             outcome="r1",
#                             N="N",
#                             method = "MH",
#                             method.tau="DL",
#                             sm = "RR")
# graphics.off()
# 
# tmp1 <- pairwise.output$summary
# 
# pairwise.output2 <- pairwise(slr,
#                              name.trt1 = "Placebo", 
#                              name.trt2 = "Vitamin D",
#                              outcome="o2",
#                              N="n2",
#                              method = "MH",
#                              method.tau="DL",
#                              sm = "OR")
# 
# tmp2 <- pairwise.output2$summary
# 
# # NOTE: okay to ignore warning message here
# pairwise.results <- bind_rows(tmp1, tmp2)
# rm(tmp1, tmp2)
# 
# source("pairwise.all.R")
# 
# pairwise.output.all <- pairwise.all(slr, 
#                                     outcome="r1",
#                                     N="N",
#                                     sm="RR")
# 
# pairwise.output.all$se.effect %>% max
# 
# 
# # Nma with JAGS -----------------------------------------------------------
# 
# # NOTE: turned off: print("The baseline treatment was ...")
# 
# source("nma.bugs.R")
# 
# makebugs <- nma.bugs(slr,
#                      outcome="r1", 
#                      N="N",
#                      baseline.name="Placebo",
#                      family="binomial",
#                      link="logit",
#                      effects="random")
# 
# makebugs
# 
# source("nma.analysis.R")
# 
# jagsoutput <- nma.analysis(makebugs,
#                            monitor = c("d"),
#                            n.adapt=10000,
#                            n.burnin=10000,
#                            n.iter=100000)
# 
# jagsmodel <- jagsoutput$model
# jagssamples <- jagsoutput$samples
# jagstrtkey <- jagsoutput$trt.key
# 
# # # trace plots
# #png("output/traceplot%02d.png")
# #plot(jagssamples)
# #graphics.off()
# # 
# # summary(jagssamples)
# 
# # Assess Model fit --------------------------------------------------------
# source("assess.model.fit.R")
# 
# 
# 
# # League Table ------------------------------------------------------------
# 
# source("league.table.R")
# 
# league.out <- leaguetable(jagsoutput)
# 
# leaguetable(jagsoutput, central.tdcy="mean")
# 
# leaguetable(jagsoutput, central.tdcy="median", layout="long")
# 
# write.csv(league.out, "output/leaguetable.csv")
# 
# 
# 
# 
# # SUCRA -------------------------------------------------------------------
# 
# source("sucra.table.R")
# 
# sucra.out <- sucra(jagsoutput, largerbetter=FALSE)
# 
# sucra.out.table <- sucra.out$s.table
# 
# sucra.out.plot <- sucra.out$s.plot
# 
# write.csv(sucra.out.table, "output/sucra.csv")
# 
# png("output/sucra.png", height=1080, width=1920, res=250)
# plot(sucra.out.plot)
# graphics.off()
# 
# 
# 
# # NMA Forest Plot --------------------------------------------------------
# 
# # NOTE: make sure you use the correct base.trt val
# 
# source("nma.forestplot.R")
# 
forest.out <- nma.forest(fixed_effects_results, comparator="SK")
# 
# nma.forestplot(jagsoutput, central.tdcy="median", base.trt="Placebo")
# 
# #nma.forestplot(jagsoutput, base.trt="2", line.size=1.5)
# 
# #png("output/forest.png", width=900, height=500)
# #plot(forest.out)
# #graphics.off()3
# 
# # Inconsistency pairwise vs nma ------------------------------------------
# 
# source("inconsistency.PMA.R")
# 
# inconsistency1 <- inconsistency.PMA(slr=slr, 
#                                      jagsoutput=jagsoutput,
#                                      base.trt="treatment 1", 
#                                      model="random", 
#                                      central.tdcy="median",
#                                      outcome="r1",
#                                      N="N")
# 
# inconsistency1$table
# incons.plot.pma <- inconsistency1$plot
# 
# write.csv(inconsistency1$table, file = "output/inconsistencyPMA.csv")
# 
# 
# png("output/inconsistencyPMA.png",  height=1080, width=1920, res=200)
# plot(incons.plot.pma)
# graphics.off()
# 
# 
# # Inconsistency Model TSD4-------------------------------------------------
# 
# source("nma.bugs.R")
# source("nma.analysis.R")
# 
# makebugsincons <- nma.bugs(slr,
#                      outcome="r1", 
#                      N="N",
#                      baseline.name="treatment 1",
#                      family="binomial",
#                      link="logit",
#                      effects="random",
#                      type="inconsistency")
# 
# inconsistency <- nma.analysis(makebugsincons,
#                               n.adapt=10000, 
#                               n.burnin=10000, 
#                               n.iter=100000)
# 
# jags.incons.model <- inconsistency$model
# jags.incons.samples <- inconsistency$samples
# jags.incons.trtkey <- inconsistency$trt.key
# 
# png("output/traceplot-incons.png")
# plot(jags.incons.samples)
# graphics.off()
# # 
# # summary.incons <- summary(jags.incons.samples)
# 
# source("inconsistency.model.R")
# 
# incons.results <- inconsistency.model(inconsistency, jagsoutput)
# 
# png("output/inconsistencymodel.png", height=1080, width=1920, res=250)
# plot(incons.results)
# graphics.off()
