source("load.packages.R")
source("functions.R")


# Import data --------------------------------------------------------

library(readxl)

dataset.choice <- "Afib"
#choices: Afib, diabetes, crc

if(dataset.choice=="Afib"){
  data1 <- read_excel("data/data_outcome.xlsx") 
  baselines1 <- read.csv("data/data_patients.csv") 
  #rct only
  data.nma <- data1 %>% 
    filter(!(trial %in% c("Conde 2013A",
                          "Conde 2013B", 
                          "Pohjantahti-Maaroos 2017")))
  baselines.nma <- baselines1 %>%
    filter(!(trial %in% c("Conde 2013A",
                          "Conde 2013B", 
                          "Pohjantahti-Maaroos 2017")))
  
  slr <- data.slr(data.nma,
                  baselines.nma,
                  varname.t = "trt", 
                  varname.s = "trial")
} else if(dataset.choice=="diabetes"){
  library(pcnetmeta)
  data("diabetes")
  data1 <- diabetes
  data.nma <- diabetes
  slr <- data.slr(data.nma,
                  varname.t = "t.id",
                  varname.s = "s.id")
} else if(dataset.choice=="crc") {
  data1 <- read.csv("data/data_crc.csv") 
  data.nma <- data1 %>%
    # NOTE: cannot be factor for netmetxl to work... possibly work into netmetxl function
    mutate_at(vars(study_group, trial), funs(as.character(.))) %>%
    filter(!trial=="Sandler, 2003")
  slr <- data.slr(data.nma,
                  varname.t = "study_group",
                  varname.s = "trial")
}

surv.slr <- data.slr(read_excel("data/survival_example.xlsx"),
                     varname.t = "Treatment",
                     varname.s = "Study")



# Network Characteristics -------------------------------------------------

source("netmetaxl.tables.R")

# NOTE: Specify var names; colnames(slr$raw.data)
# 
# Afib: outcome = "r1", N = "N", type.outcome = "binomial"
# diabetes: outcome = "r", N = "n", type.outcome = "binomial"

netmetaxl <- netmetaxl.tables(slr=slr,
                              outcome="r1", 
                              N="N",
                              type.outcome="binomial")

netmetaxl$network
netmetaxl$intervention
netmetaxl$comparison



# Network Plots -----------------------------------------------------------

source("network.plot.R")

# NOTE:
#
# label.offset1 controls distance from center. If label.offset1 = 0 then will plot node labels in center of node
# label.offset2 controls relative distance of labels in upper/lower quadrants relative to labels on left/right quadrants

png("output/network.plot1.png", height=1080, width=1280, res=200)
network.plot(slr, node.scale=5, edge.scale=2, label.offset1 = 0.1, label.offset2 = 0.1)
graphics.off()

#png("output/network.plot2.png", width=10, height=10, units="in", res=500)
network.plot(slr, node.scale=2, edge.scale=2, flag="Placebo" , label.offset1 = 3, label.offset2 = 1)
#graphics.off()

# Patient Characteristics plots--------------------------------------------
source("baselines.plot.R")

#png(filename="data/test_pdfs", res=1000, height=7500, width=10000)
comp_graph_comparison(x=slr$baseline.data, 
                      trial.id = "trial", 
                      treat.id = "trt", 
                      measure.val = "age_estimate",
                      measure.sd  = "age_SD",
                      x.lim = c(25,65), 
                      x.label = "Average Age (Years)" )
#graphics.off()

#png(filename="data/test_pdfs", res=1000, height=7500, width=10000)
comp_graph_trial(data=slr$baseline.data, 
                 trial.id="trial",
                 treat.id="trt",
                 measure.val = "age_estimate", 
                 x.lim = c(25,65), 
                 x.label = "Average Age (Years)" )
#graphics.off()

# Pairwise Comparisons ----------------------------------------------------

# NOTE: specify correct treatment: slr$treatments
# NOTE: specify correct oucome: colnames(slr$raw.data)

source("pairwise.R")
png("output/pairwise.plot.png", height=720, width=1920, res=200)
pairwise.output <- pairwise(slr,
                            name.trt1 = "Placebo",
                            name.trt2 = "Vernakalant IV",
                            outcome="r1",
                            N="N",
                            method = "MH",
                            method.tau="DL",
                            sm = "RR")
graphics.off()

tmp1 <- pairwise.output$summary

pairwise.output2 <- pairwise(slr,
                             name.trt1 = "Placebo", 
                             name.trt2 = "Vitamin D",
                             outcome="o2",
                             N="n2",
                             method = "MH",
                             method.tau="DL",
                             sm = "OR")

tmp2 <- pairwise.output2$summary

# NOTE: okay to ignore warning message here
pairwise.results <- bind_rows(tmp1, tmp2)
rm(tmp1, tmp2)

source("pairwise.all.R")

pairwise.output.all <- pairwise.all(slr, 
                                    outcome="r1",
                                    N="N",
                                    sm="RR")

pairwise.output.all$se.effect %>% max


# Nma with JAGS -----------------------------------------------------------

# NOTE: turned off: print("The baseline treatment was ...")

source("nma.bugs.R")

makebugs <- nma.bugs(slr,
                     outcome="r1", 
                     N="N",
                     baseline.name="Placebo",
                     family="binomial",
                     link="logit",
                     effects="random")

makebugs

source("nma.analysis.R")

jagsoutput <- nma.analysis(makebugs,
                           monitor = c("d"),
                           n.adapt=10000,
                           n.burnin=10000,
                           n.iter=100000)

jagsmodel <- jagsoutput$model
jagssamples <- jagsoutput$samples
jagstrtkey <- jagsoutput$trt.key

# # trace plots
#png("output/traceplot%02d.png")
#plot(jagssamples)
#graphics.off()
# 
# summary(jagssamples)

# Assess Model fit --------------------------------------------------------
source("assess.model.fit.R")



# League Table ------------------------------------------------------------

source("league.table.R")

league.out <- leaguetable(jagsoutput)

leaguetable(jagsoutput, central.tdcy="mean")

leaguetable(jagsoutput, central.tdcy="median", layout="long")

write.csv(league.out, "output/leaguetable.csv")




# SUCRA -------------------------------------------------------------------

source("sucra.table.R")

sucra.out <- sucra(jagsoutput, largerbetter=FALSE)

sucra.out.table <- sucra.out$s.table

sucra.out.plot <- sucra.out$s.plot

write.csv(sucra.out.table, "output/sucra.csv")

png("output/sucra.png", height=1080, width=1920, res=250)
plot(sucra.out.plot)
graphics.off()



# NMA Forest Plot --------------------------------------------------------

# NOTE: make sure you use the correct base.trt val

source("nma.forestplot.R")

forest.out <- nma.forestplot(jagsoutput, base.trt="Placebo")

nma.forestplot(jagsoutput, central.tdcy="median", base.trt="Placebo")

#nma.forestplot(jagsoutput, base.trt="2", line.size=1.5)

#png("output/forest.png", width=900, height=500)
#plot(forest.out)
#graphics.off()3

# Inconsistency pairwise vs nma ------------------------------------------

source("inconsistency.PMA.R")

inconsistency1 <- inconsistency.PMA(slr=slr, 
                                     jagsoutput=jagsoutput,
                                     base.trt="treatment 1", 
                                     model="random", 
                                     central.tdcy="median",
                                     outcome="r1",
                                     N="N")

inconsistency1$table
incons.plot.pma <- inconsistency1$plot

write.csv(inconsistency1$table, file = "output/inconsistencyPMA.csv")


png("output/inconsistencyPMA.png",  height=1080, width=1920, res=200)
plot(incons.plot.pma)
graphics.off()


# Inconsistency Model TSD4-------------------------------------------------

source("nma.bugs.R")
source("nma.analysis.R")

makebugsincons <- nma.bugs(slr,
                     outcome="r1", 
                     N="N",
                     baseline.name="treatment 1",
                     family="binomial",
                     link="logit",
                     effects="random",
                     type="inconsistency")

inconsistency <- nma.analysis(makebugsincons,
                              n.adapt=10000, 
                              n.burnin=10000, 
                              n.iter=100000)

jags.incons.model <- inconsistency$model
jags.incons.samples <- inconsistency$samples
jags.incons.trtkey <- inconsistency$trt.key

png("output/traceplot-incons.png")
plot(jags.incons.samples)
graphics.off()
# 
# summary.incons <- summary(jags.incons.samples)

source("inconsistency.model.R")

incons.results <- inconsistency.model(inconsistency, jagsoutput)

png("output/inconsistencymodel.png", height=1080, width=1920, res=250)
plot(incons.results)
graphics.off()
