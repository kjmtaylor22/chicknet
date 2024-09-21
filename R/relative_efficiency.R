Er <- function(df, list, filename=Sys.Date){
  
  library(dplyr)
  library(tidyverse)
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/Er/")==FALSE){
    dir.create("./coreplots/Er/")
  } 
  if(dir.exists(paste0("./coreplots/Er/", Sys.Date()))==FALSE){
    dir.create(paste0("./coreplots/Er/", Sys.Date()))
  } 
  
  out.df <- data.frame()
  line.df <- data.frame()
  for (i in 1:length(list)){
    
    find <- sapply(df, FUN=function(x){grep("^asv",x)[1]}) %>% unlist() %>% .[!is.na(.)]
    
    pull <- df[,c(names(find), unlist(list[[i]]))]
    
    
    if (any(names(df)%in%c("kingdom","phylum","class","order","family","genus","taxonomy"))){
      pull <- cbind(pull, df[,names(df)%in%c("kingdom","phylum","class","order","family","genus","taxonomy")])
      names(pull) <- c("asv","x","y","taxonomy")
      pull <- pull[pull$taxonomy!="",]
    } else {names(pull) <- c("asv","x","y")}
    
    
    
    pull <- mutate(pull, diffLogD=log10(x)-log10(y))
    
    my <- lmodel2::lmodel2(log10(y) ~ log10(x), data=pull, nperm=50)
    coef <- my$regression.results %>% .[.$Method=="SMA",2:3] %>% unlist() %>% round(., 3)
    conf <- my$confidence.intervals %>% .[.$Method=="SMA",2:5] %>% unlist() %>% round(., 3)
    
    SE.b <- coef[2]*sqrt((1-my$r^2)/my$n) 

    out.df <- rbind(out.df, data.frame(pull, compar=paste(unlist(list[[i]]),collapse=" (x) - (y) ")) %>%
      separate(compar, c("xlab","ylab"), " - ", remove=F))
    
    line.df <- rbind(line.df, data.frame(compar=paste(unlist(list[[i]]),collapse=" (x) - (y) "), r=my$r, n=my$n, 
                                        A.y=coef[1], b.y=coef[2], r2=round(my$rsquare,2), p=signif(my$P.param, 3),
                                        A.lwr=conf[1], A.upr=conf[2], b.lwr=conf[3], b.upr=conf[4]) %>%
                       mutate(SE.b = b.y*sqrt((1-r^2)/n), 
                              SE.A = sqrt(var(my$y))*sqrt(((1-r^2)/n)*(1+(mean(my$x)^2/var(my$x))))) %>%
                       separate(compar, c("xlab","ylab"), " - ", remove=F))
    
  }
  
  if (length(list)>5){
    dims <- ceiling(sqrt(length(list)))
    w <- 900*dims+100
    h <- 900*dims+200
  } else {
    dims <- length(list)
    w <- 900*dims+100
    h <- 1100
  }
  
  
  shift <- c((1/25)*(log10(max(out.df$x, na.rm=T))-log10(min(out.df$x, na.rm=T))), 
             (1/25)*(log10(max(out.df$y, na.rm=T))-log10(min(out.df$y, na.rm=T))))
  lab.coord <- list(X=log10(c(max(out.df$x, na.rm=T), min(out.df$x, na.rm=T))),# + c(-shift[1], shift[1]), 
                    Y=log10(c(max(out.df$y, na.rm=T), min(out.df$y, na.rm=T))))# + c(-shift[2], shift[2]))
  
  g <- ggplot(line.df) + facet_wrap(vars(compar), ncol=dims) + 
    geom_abline(aes(intercept=0,slope=1), color="grey90", linetype="dashed", size=0.8) + 
    geom_point(data=out.df, aes(x,y)) + 
    geom_abline(aes(intercept=A.y,slope=b.y), size=1.2) + 
    ggthemes::theme_few() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    geom_text(aes(x=10^lab.coord$X[2], y=10^lab.coord$Y[1], group=compar, label=paste0("y = ", A.y, " + ", b.y, "x")), color="black", vjust=1, hjust=0) +
    geom_text(aes(x=10^lab.coord$X[1], y=10^lab.coord$Y[1], group=compar, label=paste0("R^2 = ", r2)), color="black", vjust=1, hjust=1) +
    geom_text(aes(x=10^lab.coord$X[1], y=10^lab.coord$Y[2], group=compar, label=paste0("p = ", p)), color="black", vjust=0, hjust=1)
    #geom_text(aes(x=0, y=Inf, group=compar, label=paste0("y = ", A.y, " + ", b.y, "x")), color="black", vjust=1, hjust=0) +
    #geom_text(aes(x=Inf, y=Inf, group=compar, label=paste0("R^2 = ", r2)), color="black", vjust=1, hjust=1) +
    #geom_text(aes(x=Inf, y=0, group=compar, label=paste0("p = ", p)), color="black", vjust=0, hjust=1)
  
  f <- ggplot(line.df, aes(color=compar)) + 
    geom_abline(aes(intercept=0,slope=1), color="grey90", linetype="dashed", size=0.8) + 
    geom_abline(aes(intercept=A.y,slope=b.y, color=compar), size=1.2) + 
    ggthemes::theme_few() + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) 
  
  jpeg(paste0("coreplots/Er/",Sys.Date(),"/", filename,".jpeg"), width=w, height=h, res=300)
  print(g)
  dev.off()
  
  #jpeg(paste0("coreplots/Er_", filename,"_cat.jpeg"), 2000, 1100, res=300)
  #print(f)
  #dev.off()
  
  return(list("points"=out.df, "line"=line.df))
  
}