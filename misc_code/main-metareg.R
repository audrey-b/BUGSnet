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


regressor <- list(coefficient='exchangeable',
                  variable='stroke',
                  control="02")
gemtc_model <- mtc.model(atrialFibrillation,
                         type="regression",
                         link="logit",
                         likelihood="binom",
                         linearModel = "random",
                         regressor=regressor)
gemtc_results <- mtc.run(gemtc_model, 
                         n.adapt = 1000, 
                         n.iter = 50000, 
                         thin = 1)

gemtc_leaguetable <- relative.effect.table(gemtc_results, covariate=0.1) %>% as.data.frame()
#write.csv(gemtc_leaguetable, "gemtc_leaguetable.csv")

conv_exp <- function(x) {
  strsplit(x, split="[ \\(\\),]")[[1]][c(1,3,5)] %>% #split into mean, lower, upper
    as.numeric %>%
    exp() %>%
    formatC(digits=2, format= "f") %>% #keep two decimals
    str_c(c(" (", ", ", ")")) %>% #put back into string
    paste0(collapse="") %>% #put back into string
    return()
}

library(stringr)
for(i in 1:length(gemtc_leaguetable)){
  for(j in (1:length(gemtc_leaguetable))[-i]){
    gemtc_leaguetable[i,j] %<>% conv_exp()
  }
}
write.csv(gemtc_leaguetable, "gemtc_leaguetable.csv")


forest(relative.effect(gemtc_results, "02"))

gemtc_ranks <- rank.probability(gemtc_results)
plot(gemtc_ranks)
plot(gemtc_ranks, beside=TRUE) 

#png("gemtcplots%02d.png")
#plot(gemtc_results)
#graphics.off()

gemtc_model$code %>% cat
gemtc_model$om.scale^(-2)
(15*gemtc_model$om.scale)^(-2)



#plotCovariateEffect(gemtc_results, 
#                   atrialFibrillation$data.ab$treatment %>% unique, 
#                  atrialFibrillation$data.ab$treatment %>% unique, 
#                 xlim=NULL, 
#                ylim=NULL, 
#               ask=dev.interactive(orNone=TRUE))


#rawdata <- atrialFibrillation$data.ab %>%
#  left_join(atrialFibrillation$studies, by="study")

#dataprep <- data.prep(arm.data = rawdata,
#                      varname.t = "treatment",
#                      varname.s = "study")

data(afib)

dataprep <- data.prep(arm.data = afib,
                      varname.t = "treatment",
                      varname.s = "study")

random_effects_model <- nma.model(data=dataprep,
                                  outcome="events",
                                  N="sampleSize",
                                  reference="02",
                                  family="binomial",
                                  link="logit",
                                  effects="random",
                                  covariate="stroke",
                                  prior.beta="EXCHANGEABLE")

bugsnet_results <- nma.run(random_effects_model,
                           n.iter=10000,
                           n.adapt=1000,
                           n.burnin=1000)

random_effects_model$bugs %>% cat

#png("bugsnetplots%02d.png", width=2000, height=2000)
#nma.trace(bugsnet_results)
#graphics.off()



random_effects_fit <- nma.fit(bugsnet_results, main = "Random Effects Model" )
random_effects_fit$DIC
random_effects_fit$pD
random_effects_fit$pmdev
random_effects_fit$leverage

sucra.out <- nma.rank(bugsnet_results, largerbetter=FALSE, cov.value=0.1)
sucra.out$sucraplot
sucra.out$rankogram

nma.forest(bugsnet_results, 
           comparator="02", 
           central.tdcy = "median",
           cov.value=0.1)

bugsnet_league <- nma.league(bugsnet_results, 
                             central.tdcy = "median",
                             order = as.vector(t(dataprep$treatments)),
                             cov.value=0.1, 
                             log.scale = FALSE)
write.csv(bugsnet_league$table, "bugsnet_leaguetable.csv")

nma.regplot(bugsnet_results, x.range=c(0.1,1))

