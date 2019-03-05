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


rawdata <- read_excel("data/continuous_example.xlsx", 
                      col_types = c("text", "numeric", "numeric", 
                                    "text", "text", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric", "numeric", 
                                    "numeric", "numeric"))


# Import data --------------------------------------------------------

cont.slr <- data.prep(arm.data = rawdata,
                      varname.t = "trt", #errors with trt_name
                      varname.s = "trial")

#ensure that treatment and study variables are both of class character
#cont.slr$arm.data$treatment <- as.character(cont.slr$arm.data$treatment)
#cont.slr$arm.data$study <- as.character(cont.slr$arm.data$study)

pairwise <- pma(data = cont.slr,
    type.outcome="continuous",
    name.trt1 = "1", 
    name.trt2 = "2",
    outcome = "y",
    sd="sd",
    N = "n",
    sm = "MD")


# Network Characteristics -------------------------------------------------

net.plot(cont.slr,  
         flag="1", 
         node.scale = 1.5, 
         edge.scale=0.5, 
         label.offset1 = 4, 
         label.offset2 = 4)

network.char <- net.tab(data = cont.slr,
                        outcome = "y",
                        N = "n",
                        type.outcome = "continuous")

network.char$network
network.char$intervention
network.char$comparison

#NMA

fixed_effects_model <- nma.model(data=cont.slr,
                                 outcome="y",
                                 N="n",
                                 sd="sd",
                                 reference="1",
                                 type="consistency",
                                 family="normal",
                                 link="identity",
                                 effects="fixed")

random_effects_model <- nma.model(data=cont.slr,
                                  outcome="y",
                                  N="n",
                                  sd="sd",
                                  type="consistency",
                                  reference="1",
                                  family="normal",
                                  link="identity",
                                  effects="random")

#sink("Z:/ResearchDocuments/Research/BUGSnet/code.bug")
cat(random_effects_model$bugs)
#sink()


random_effects_results <- nma.run(random_effects_model,
                                  monitor = c("d","sigma"),
                                  n.adapt=1000,
                                  n.burnin=1000,
                                  n.iter=10000)

random_effects_fit <- nma.fit(random_effects_results, main = "Random Effects Model" )
random_effects_fit$DIC
random_effects_fit$pD
random_effects_fit$pmdev

sucra.out <- nma.rank(random_effects_results, 
                      largerbetter=FALSE, 
                      cov.value=NULL)
sucra.out$sucraplot
sucra.out$rankogram

nma.forest(random_effects_results, 
           comparator="1", 
           central.tdcy = "median")
league.res <- nma.league(random_effects_results, 
           central.tdcy = "median")

png("test.png", width=1500, height=1000)
league.res$heatplot
graphics.off()
