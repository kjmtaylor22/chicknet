
#' @export

tukey.select <- function(x, n, reduce=T){
  # x=tukey output table,
  # n=number of factors compared
  # reduce=NULL gives all comparisons
  # reduce=F gives only useful comparisons
  # reduce=T gives only significant useful comparisons
  z <- row.names(x)
  require(dplyr)
  x1 <- lapply(z, FUN=function(z){
    y <- strsplit(z, split="-")[[1]] %>%
          strsplit(split=':') %>% unlist()
    return(y)
  })
  x2 <- lapply(x1, FUN=function(x1){
    out <- 0
    for (i in 1:n){
      for (j in 1:n){
        y <- grep(x1[i], x1[n+j])
        if (length(y)==1){out <- out+1}
      }
    }
    if (out==n-1){y <- "yes"} else {y <- "no"}
    return(y)
  })
  z <- cbind(x, data.frame(useful=unlist(x2), signif=NA))
  vvsig <- which(z[,4] < 0.001) %>% as.integer()
  vsig <- which(z[,4] >=0.001 & z[,4] < 0.01) %>% as.integer()
  sig <- which(z[,4] >= 0.01 & z[,4] < 0.05) %>% as.integer()
  mar <- which(z[,4] >= 0.05 & z[,4] < 0.1) %>% as.integer()
  non <- which(z[,4] >= 0.1) %>% as.integer()
  z$signif[vvsig] <- "***"
  z$signif[vsig] <- "**"
  z$signif[sig] <- "*"
  z$signif[mar] <- "."
  z$signif[non] <- " "
  if (is.null(reduce)){
    return(z)
  } else {
    if (reduce==TRUE){
      z <- z[z$useful=="yes" & z$signif!=" ",]
    }
    if (reduce==F){
      z <- z[z$useful=="yes",]
    }
    return(z)
  }
}
