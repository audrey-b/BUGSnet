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

system("R CMD Rd2pdf . --title=Package BUGSnet --output=./manual.pdf --force --no-clean --internals")

install.packages("~/BUGSnet_0.0.0.9000.tar.gz", repos = NULL, type = "source")
library(BUGSnet)
??BUGSnet


devtools::document()
devtools::build_vignettes()
devtools::build()

system("tar zxvf ../BUGSnet_0.0.0.9000.tar.gz")
system("R CMD Rd2pdf BUGSnet")

