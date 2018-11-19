quicknames <- function(prefixes, suffixes, sep=""){
  
  c(outer(prefixes, paste0(sep, suffixes), FUN=paste0))
  
}


