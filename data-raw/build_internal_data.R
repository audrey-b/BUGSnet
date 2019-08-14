
#############################################################
# Description: builds internal use R/sysdata.rda data file  #
#   which is used to implement bugsnet.test                 #
#                                                           #
#############################################################

tsd2ex5data <- read.csv(paste0(getwd(), "/data-raw/TSD2_ex5_parkinsons.csv"), 
                    stringsAsFactors = FALSE)
tsd2ex5results <- read.csv(paste0(getwd(), "/data-raw/TSD2_ex5_parkinsons_results.csv"), 
                       stringsAsFactors = FALSE)
tsd2ex5bugsnet <- read.csv(paste0(getwd(), "/data-raw/TSD2_ex5_parkinsons_bugsnet.csv"), 
                           stringsAsFactors = FALSE, check.names = F)
tsd2ex5 <- list(data = tsd2ex5data, results = tsd2ex5results, bugsnet = tsd2ex5bugsnet)
usethis::use_data(tsd2ex5, internal = TRUE)
