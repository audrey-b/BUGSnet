library(devtools)
devtools::install(build_vignette = FALSE)
library(BUGSnet)
devtools::check()


devtools::load_all("R")
document()
devtools::check(document = FALSE)
devtools::check()

devtools::build_vignettes()

#devtools::install(build_vignettes = TRUE)


devtools::install()
library(BUGSnet)

install.packages("~/Compile BUGSnet/BUGSnet_0.0.0.9000.tar.gz", repos = NULL, type = "source")
library(BUGSnet)
??BUGSnet
