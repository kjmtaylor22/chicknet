alpha.boxplots <- function(comm, meta, select=NULL, group, legend=TRUE, fixed=NULL, type="box"){
  
  pkgs <- c("vegan", "dplyr","fossil")
  if (any(!pkgs%in%row.names(installed.packages()))) {
    install.packages(pkgs[which(!pkgs%in%row.names(installed.packages())==T)])
  }
  library(dplyr)
  library(vegan)
  library(fossil)
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/alpha/")==FALSE){
    dir.create("./coreplots/alpha/")
  }
  if(dir.exists(paste0("./coreplots/alpha/", deparse(substitute(comm))))==FALSE){
    dir.create(paste0("./coreplots/alpha/", deparse(substitute(comm))))
  }
  
  id <- grep("sample", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "SampleID"
  } else {
    names(meta)[1] <- "SampleID"
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
    id <- 1
  }
  
  
  g <- which(names(meta)%in%group)
  
  
  if (is.null(select)==TRUE){
    sep.comm <- comm
  } else {
    sep.comm <- comm[grep(select, row.names(comm)),]
  }
  
  if (any(!row.names(sep.comm)%in%meta$SampleID)){
    warning("Check that all community samples are accounted for in metadata.")
    warning("May also need to rename sample ID column in metadata.")
    return(NULL)
  }
  
  shan <- diversity(sep.comm)
  rich <- specnumber(sep.comm)
  even <- shan/log(rich)
  #chao <- apply(sep.comm, MARGIN=1, FUN=chao1)
  
  shan <- data.frame(SampleID=names(shan), shan=shan, rich=rich, even=even)#, chao1=chao)
  
  
  
  shan <- left_join(shan, meta[,c(id,g)])
  #shan[,group] <- factor(shan[,group], levels = unique(shan[,group]))
  
  names(shan) <- c("SampleID", "shan", "rich", "even", "group.id")
  
  
  shan$group.id <- factor(shan$group.id,
                          levels=sort(unique(shan$group.id)))
  shan <- shan[order(shan$group.id),]
  
  clr <- c("magenta", "red", "orange", "yellow", "lawngreen", "mediumturquoise", "dodgerblue3", "navy", "blueviolet")
  clr <- colorRampPalette(clr)(length(unique(shan$group.id)))
  names(clr) <- as.character(unique(shan$group.id))
  palette(clr)
  
  or.par <- par( no.readonly = TRUE )
  
  if (type=="box"){
    jpeg(paste("./coreplots/alpha/",deparse(substitute(comm)),"/ALPHA_boxplot_", group, select, ".jpg", sep=""), 6000, 2300, res=800)
    par(mfrow=c(1,3), cex.axis=1.1, font.axis=2, las=3)#, xaxt="n")
    
    out <- list()
    fixedx <- 0
    for (j in c("shan", "rich", "even")){#, "chao1")){
      if (is.null(fixed)){
        boxplot(eval(parse(text=j)) ~ group.id, data=shan, boxfill=palette(), 
                boxwex=.8, main=select, ylab=j, xlab="")
      } else {
        fixedx <- fixedx+1
        boxplot(eval(parse(text=j)) ~ group.id, data=shan, boxfill=palette(), 
                boxwex=.8, main=select, ylab=j, ylim=fixed[[fixedx]])
      }
      axis(1, at=c(1:length(levels(shan$group.id))), labels=FALSE)
      text(c(1:length(levels(shan$group.id))), par("usr")[3], font=2, 
           labels=c(levels(shan$group.id)), srt=45, adj=c(1, 1.5), offset=2.5)
      if (legend==TRUE){
        if (j=="shan"){
          legend(c(0,-0.5), horiz = TRUE, fill = palette(), legend = levels(shan$group.id), xpd=T, bty = 'n', cex=1)
        }
      }
      
      atw <- aov(eval(parse(text=j)) ~ group.id, data=shan)
      tuktw <- TukeyHSD(atw)
      #source("../functions/tukey_select.R")
      tuksel <- tukey.select(tuktw$`group.id`, 1)
      
      test1 <- pairwise.t.test(with(shan, eval(parse(text=j))), shan$group.id, "bonferroni", pool.sd=T)
      test <- pairwise.wilcox.test(with(shan, eval(parse(text=j))), shan$group.id, "bonferroni")
      out[[j]] <- list("ANOVA"=summary(atw), "M-W U test"=test, "T-test"=test1, "TukeyHSD"=tuksel)
    }
    dev.off()
    
    sink(paste("./coreplots/alpha/",deparse(substitute(comm)),"/ALPHA_boxplot-stats_", group, select, ".txt", sep=""))
    print(out)
    sink() 
  }
  
  if (type=="point"){
    jpeg(paste("./coreplots/alpha/",deparse(substitute(comm)),"/ALPHA_scatter_", group, select, ".jpg", sep=""), 6000, 2100, res=800)
    par(mfrow=c(1,3), cex.axis=1.1, mar=c(4,4,3,2), font.axis=2, las=3, xaxt="n")
    
    out <- list()
    fixedx <- 0
    for (j in c("shan", "rich", "even")){#, "chao1")){
      if (is.null(fixed)){
        yax <- shan[,j]
        plot(eval(parse(text=j)) ~ jitter(as.numeric(group.id),1), data=shan, bg=group.id, main=select, ylab=j, xlab=" ", pch=22, cex=1.5)
      } else {
        fixedx <- fixedx+1
        yax <- fixed[[fixedx]]
        plot(eval(parse(text=j)) ~ jitter(as.numeric(group.id),1), data=shan, bg=group.id, main=select, ylab=j, xlab=" ", ylim=yax, pch=22, cex=1.5)
      }
      
      #axis(1, at=c(1:length(levels(shan$group.id))), labels=FALSE)
      #text(c(1:length(levels(shan$group.id))), par("usr")[3], font=2, 
      #     labels=c(levels(shan$group.id)), srt=45, adj=c(1, 1.5), xpd=T, offset=2.5)
      
      m <- lm(eval(parse(text=j)) ~ as.numeric(group.id), data=shan)
      sum <- summary(m)
      abline(m$coefficients[[1]], m$coefficients[[2]])
      
      equation <- as.expression(bquote(y == .(signif(m$coefficients[[2]], 3)) ~ x + .(signif(m$coefficients[[1]], 3)) ~~~~~~~~~ R^2 == .(signif(sum$r.squared, 3))))
      
      legend(1, (max(yax)-(max(yax)-min(yax))/10), box.lty=0, bg=alpha("white", 0.75), y.intersp=0.05, cex=0.5, legend=equation)
      
      
      if (legend==TRUE){
        if (j=="rich"){
          xax <- unique(as.numeric(shan$group.id))
          
          legend(-(max(xax)-min(xax))/4, (min(yax)-(max(yax)-min(yax))/6), 
                 x.intersp=0.2, horiz = TRUE, fill = palette(), legend = levels(shan$group.id), xpd=T, bty = 'n', pt.cex=1)
        }
      }
      
    }
    dev.off()
  }
  
  return(shan)
}