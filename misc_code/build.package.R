#change version
#change Depends: R (>= 3.6)

#change NEWS file

#build vignettes for website
devtools::build_vignettes()



#build pdf manual for website
system("R CMD Rd2pdf . --title=Package BUGSnet --output=./manual.pdf --force --no-clean --internals")
#delete .Rd2pdf8444 folder