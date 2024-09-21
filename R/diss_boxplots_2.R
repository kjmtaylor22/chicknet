diss.boxplots <- function(dist, meta, select, group, legend=TRUE, fixed=NULL){
  palette(c("magenta", "red", "orange", "yellow", "lawngreen", "mediumturquoise", "dodgerblue3", "navy", "blueviolet"))
  
  if (!"dplyr"%in%row.names(installed.packages())) {install.packages("dplyr")}
  library(dplyr)
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/dissimilarity/")==FALSE){
    dir.create("./coreplots/dissimilarity/")
  }
  
  if(dir.exists(paste0("./coreplots/dissimilarity/", deparse(substitute(dist))))==FALSE){
    dir.create(paste0("./coreplots/dissimilarity/", deparse(substitute(dist))))
  }
  if(dir.exists(paste0("./coreplots/dissimilarity/", deparse(substitute(dist)), "/", group))==FALSE){
    dir.create(paste0("./coreplots/dissimilarity/", deparse(substitute(dist)), "/", group))
  }
  
  if (class(dist)%in%c("dist","data.frame")){
    dist <- as.matrix(dist)
  }
    
  sep.dist <- data.frame(SampleID=row.names(dist), as.matrix(dist))
  row.names(sep.dist) <- row.names(dist)
  
  if (!is.null(select)){
    sep.dist <- sep.dist[grep(select, sep.dist$SampleID),
                   grep(select, names(sep.dist))]
  } 
    
  centroids <- sep.dist %>% left_join(meta[,c("SampleID",group)]) %>% 
    .[,-1] %>% group_by_(group) %>% summarize_all(funs(mean)) %>% as.data.frame() %>% `row.names<-` (.[,1]) %>% 
    .[,-1]  %>% t() %>% as.data.frame() %>% cbind(SampleID=row.names(.)) %>% reshape2::melt() 


  
  id <- grep("SampleID", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "SampleID"
  } else {
    names(meta)[1] <- "SampleID"
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
  }
  if (any(!row.names(sep.dist)%in%meta[,which(names(meta)=="SampleID")])){
    warning("Check that all community samples are accounted for in metadata.")
    warning("May also need to rename sample ID column in metadata.")
    return(NULL)
  }
  
  
  tw <- left_join(centroids, meta[,c("SampleID", group)])
  tw <- tw[as.character(tw$variable)==as.character(tw[,group]),]


  clr <- c("magenta", "red", "orange", "yellow", "lawngreen", "mediumturquoise", "dodgerblue3", "navy", "blueviolet")
  clr <- colorRampPalette(clr)(length(unique(tw$variable)))
  names(clr) <- as.character(unique(tw$variable))
  palette(clr)
  
  jpeg(paste("./coreplots/dissimilarity/", deparse(substitute(dist)), "/", group,"/DIST_boxplot_", group, select, ".jpg", sep=""), 2300, 2650, res=600)
  par(cex.axis=1.1, font.axis=2, las=3, xaxt="n")
  if (is.null(fixed)){
    boxplot(value ~ variable, data=tw, boxfill=palette(), boxwex=.8, 
            main=select, ylab="dissimilarity")
  } else {
    boxplot(value ~ variable, data=tw, boxfill=palette(), 
            boxwex=.8, main=select, ylab="dissimilarity", ylim=fixed)
  }
  axis(1, at=c(1:length(levels(tw$variable))), labels=FALSE)
  text(c(1:length(levels(tw$variable))), par("usr")[3], font=2, 
       labels=c(levels(tw$variable)), srt=45, adj=c(1, 1.5), xpd=T, offset=2.5)
  if (legend==TRUE){
    legend('top', horiz = TRUE, fill = palette(), legend = levels(tw$variable), bty = 'n', cex=0.5)
  }
  dev.off()
  
  atw <- aov(value ~ variable, data=tw)
  tuktw <- TukeyHSD(atw)
  tuksel <- tukey.select(tuktw$`variable`, 1)
  
  test1 <- pairwise.t.test(tw$value, tw$variable, "bonferroni", pool.sd=T)
  test <- pairwise.wilcox.test(tw$value, tw$variable, "bonferroni")
  
  sink(paste("./coreplots/dissimilarity/",  deparse(substitute(dist)), "/", group,"/DIST_boxplot-stats_", group, select, ".txt", sep=""))
  print(list("ANOVA"=summary(atw), "MW U test"=test, "t-test"=test1, "tukey"=tuksel))
  sink()  
  
  return(tw)
}