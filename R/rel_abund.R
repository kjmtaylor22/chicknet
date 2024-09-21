ra <- function(x){ #where x is a vector
  y <- sum(x)
  z <- sapply(x, FUN=function(x){z1 <- x/y})
  return(z)
}