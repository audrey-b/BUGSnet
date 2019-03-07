library(devtools)
library(readxl)
library(magrittr)
library(dplyr)

# dietfat <- read_excel("data-raw/rate_example.xlsx")
# dietfat <- as.data.frame(dietfat)
# use_data(dietfat, overwrite=TRUE) #rate

diabetes.sim <- read_excel("data-raw/rate2_example.xlsx")
diabetes.sim <- as.data.frame(diabetes.sim)
use_data(diabetes.sim, overwrite=TRUE) #rate2
diabetes <- diabetes.sim[,1:6]
use_data(diabetes, overwrite=TRUE)

library(gemtc)
data(thrombolytic)
thrombolytic <- thrombolytic$data.ab %>%
  rename(events=responders)
use_data(thrombolytic, overwrite=TRUE) #dichotomous

data(atrialFibrillation)
afib <- atrialFibrillation$data.ab %>%
  left_join(atrialFibrillation$studies, by="study") %>%
  rename(events=responders)
use_data(afib, overwrite=TRUE) #metareg
