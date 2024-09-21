pcoa.centroids <- function(dist, meta, factor, col.pal="default", connect=F, connect.by=NULL, filename=Sys.Date()){
  
  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/pcoa/")==FALSE){
    dir.create("./coreplots/pcoa/")
  }
  if(dir.exists(paste0("./coreplots/pcoa/", deparse(substitute(dist))))==FALSE){
    dir.create(paste0("./coreplots/pcoa/", deparse(substitute(dist))))
  }
  
  library(ggplot2)
  library(dplyr)
  
  fac <- length(unique(meta[meta$SampleID%in%row.names(as.matrix(dist)),factor]))
  
  gg_color_hue <- function(n) {
    hues = seq(20, 380, length = n+1)
    hcl(h = hues, l = 70, c = 140)[1:n]}
  gg_color_hue2 <- function(n) {
    hues = seq(240, 580, length = n+1)
    hcl(h = hues, l = 50, c = 120)[1:n]}
  gg_color_hue3 <- function(n) {
    hues = seq(240, 580, length = n+1)
    hcl(h = hues, l = 90, c = 120)[1:n]}
  
  if (length(col.pal)>1){colors <- col.pal} else {
  
    if (col.pal=="rainbow"){colors <- colorRampPalette(c("magenta", "red", "orange", "yellow", "green", "cyan", "dodgerblue", "purple4"))(fac)}
    if (col.pal=="Set1"){if (fac>8){colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(fac)} else {colors <- RColorBrewer::brewer.pal(8, "Set1")}}
    if (col.pal=="default"){colors <- gg_color_hue(fac)}
    if (col.pal=="paired"){colors <- c(rbind(gg_color_hue2(fac/2),gg_color_hue3(fac/2)))}
    
    names(colors) <- as.character(sort(unique(meta[meta$SampleID%in%row.names(dist),factor])))
  }

  seq.cent <- function(x, axis){
    y <- c(1:length(axis))
    for (i in 1:fac-1){
      y[which(x==levels(x)[i])] <- axis[which(x==levels(x)[i+1])[1]]
    }
    y[which(x==levels(x)[fac])] <- axis[which(x==levels(x)[fac])[1]]
    return(y)
  }
  

  group.cent <- function(x, by, axis){
    inters <- as.character(unique(interaction(x, by)))
    y <- c(1:length(axis))
    for (j in as.character(unique(by))){
      pull <- grep(j, inters, value=T)
        splt <- unlist(strsplit(pull[1], split="\\."))
        splt2 <- unlist(strsplit(pull[2], split="\\."))
        y[which(x==splt[1] & by==splt[2])] <- axis[which(x==splt2[1] & by==splt2[2])[1]]
        y[which(x==splt2[1] & by==splt2[2])] <- axis[which(x==splt2[1] & by==splt2[2])[1]]
    }
    return(y)
  }

  pcoa <- ape::pcoa(dist)
  
  eig <- pcoa$values$Eigenvalues[1:length(pcoa$vectors[1,])] #only positive eigenvalues
  rescale <- eig/sum(eig)
  varexp <- round(rescale/sum(rescale)*100, 1)
  
  meta[,factor] <- as.factor(meta[,factor])
  
  if (is.null(connect.by)){
    centroids <- data.frame(SampleID=row.names(pcoa$vectors), pcoa$vectors) %>%
      left_join(meta[,c("SampleID", factor)]) %>%
      reshape2::melt() %>%
      subset(subset=variable%in%c("Axis.1", "Axis.2", "Axis.3")) %>%
      group_by_(factor, "variable") %>%
      mutate(centroid.1=mean(value))
    names(centroids)[which(names(centroids)==factor)] <- "factor"
    centroids <- centroids[order(centroids$SampleID),]
    samples <- reshape2::dcast(centroids, formula=SampleID + factor ~ variable, value.var = "value")
  } else {
    centroids <- data.frame(SampleID=row.names(pcoa$vectors), pcoa$vectors) %>%
      left_join(meta[,c("SampleID", factor, connect.by)]) %>%
      reshape2::melt() %>%
      subset(subset=variable%in%c("Axis.1", "Axis.2", "Axis.3")) %>%
      group_by_(factor, connect.by, "variable") %>%
      mutate(centroid.1=mean(value))
    names(centroids)[which(names(centroids)==factor)] <- "factor"
    names(centroids)[which(names(centroids)==connect.by)] <- "connect"
    centroids <- centroids[order(centroids$SampleID),]
    samples <- reshape2::dcast(centroids, formula=SampleID + factor + connect ~ variable, value.var = "value")
  }
  
  samples <- data.frame(samples, 
                        centroid.1=centroids$centroid.1[centroids$variable=="Axis.1"],
                        centroid.2=centroids$centroid.1[centroids$variable=="Axis.2"],
                        centroid.3=centroids$centroid.1[centroids$variable=="Axis.3"])
  if (connect==T){
    if (!is.null(connect.by)){

      samples <- dplyr::mutate(samples, FactorCent.1=group.cent(factor, connect,  centroid.1), 
                               FactorCent.2=group.cent(factor, connect,  centroid.2), 
                               FactorCent.3=group.cent(factor, connect,  centroid.3))
    } else {
      samples <- dplyr::mutate(samples, FactorCent.1=seq.cent(factor, centroid.1), 
                               FactorCent.2=seq.cent(factor, centroid.2), 
                               FactorCent.3=seq.cent(factor, centroid.3))
    }
  }
  if (connect==T){
      con.lines1 <- expression(geom_segment(aes(x=centroid.1, xend=FactorCent.1, y=centroid.2, yend=FactorCent.2), color="black", size=1))
      con.lines2 <- expression(geom_segment(aes(x=centroid.1, xend=FactorCent.1, y=centroid.3, yend=FactorCent.3), color="black", size=1))
  } else {
      con.lines1 <- expression(geom_segment(aes(x=centroid.1, xend=centroid.1, y=centroid.2, yend=centroid.2), color="black", size=1))
      con.lines2 <- expression(geom_segment(aes(x=centroid.1, xend=centroid.1, y=centroid.3, yend=centroid.3), color="black", size=1))
  }
  
  if (col.pal=="paired"){
    leg.row <- expression(guides(fill=guide_legend(nrow=2, byrow=F)))
  } else {leg.row <- expression(guides(fill=guide_legend(nrow=ceiling(sqrt(fac)), byrow=F)))}
  
  p12 <- ggplot(samples) +
    geom_segment(aes(x=Axis.1, xend=centroid.1, y=Axis.2, yend=centroid.2, color=factor), alpha=0.5) +
    geom_point(size = 4, aes(Axis.1, Axis.2, fill=factor), shape=21, color="black", alpha=0.75) +
    geom_point(size = 5, aes(centroid.1, centroid.2, fill=factor), shape=24) +
    eval(expr=con.lines1) +
    geom_hline(yintercept=0, color="grey80", linetype="dotdash") +
    geom_vline(xintercept=0, color="grey80", linetype="dotdash") +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(), 
          panel.background = element_rect(fill = "white", color = "black", size = 1),
          plot.title = element_text(size = 25, face="bold", hjust = 0.5),
          axis.title = element_text(size = 15, face="bold"),
          axis.text = element_text(color = "black", size = 12),
          legend.key = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.text = element_text(size = 18, margin=margin(r=25, l=5, b=3)),
          legend.position="bottom",
          legend.justification = "center",
          legend.direction="horizontal") +
    xlab("PC1") + ylab("PC2") + #labs(title=samples$BodySite[1]) + 
    eval(expr=leg.row)
    
  p13 <- ggplot(samples) +
    geom_segment(aes(x=Axis.1, xend=centroid.1, y=Axis.3, yend=centroid.3, color=factor), alpha=0.5) +
    geom_point(size = 4, aes(Axis.1, Axis.3, fill=factor), shape=21, color="black", alpha=0.75) +
    geom_point(size = 5, aes(centroid.1, centroid.3, fill=factor), shape=24) +
    eval(expr=con.lines2) +
    geom_hline(yintercept=0, color="grey80", linetype="dotdash") +
    geom_vline(xintercept=0, color="grey80", linetype="dotdash") +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(), 
          panel.background = element_rect(fill = "white", color = "black", size = 1),
          plot.title = element_text(size = 25, face="bold", hjust = 0.5),
          axis.title = element_text(size = 15, face="bold"),
          axis.text = element_text(color = "black", size = 12),
          legend.key = element_rect(fill="white"),
          legend.title = element_blank(),
          legend.text = element_text(size = 18, margin=margin(r=25, l=5, b=3)),
          legend.position="bottom",
          legend.justification = "center",
          legend.direction="horizontal") +
    xlab("PC1") + ylab("PC3") + #labs(title=samples$BodySite[1]) + 
    guides(fill=guide_legend(nrow=1))
  

  tiff(paste("./coreplots/pcoa/", deparse(substitute(dist)), "/PCOA_", 
             paste(factor, connect.by, filename, sep="_"), ".tiff", sep=""), 
             width=4500, height=2700, res=600)
gridExtra::grid.arrange(p12+guides(fill=F, color=F), p13+guides(fill=F, color=F), 
                        ggpubr::get_legend(p12),
                        layout_matrix=rbind(c(1,1,1,2,2,2), c(1,1,1,2,2,2), c(1,1,1,2,2,2),
                                            c(3,3,3,3,3,3)))
dev.off()

pcoa <- left_join(data.frame(pcoa$vectors, SampleID=row.names(pcoa$vectors)), 
                    meta[,c("SampleID",factor)])

sink(paste("./coreplots/pcoa/", deparse(substitute(dist)), "/PCOA_", 
           paste(factor, connect.by, filename, sep="_"), ".txt", sep=""))
cat("Variance explained: ", "\n",
             paste("PC1 (", as.character(varexp[1]), "%)"), "\n",
             paste("PC2 (", as.character(varexp[2]), "%)"),
             "\n","\n","~~~ * ~~~","\n","\n")
print(vegan::adonis(dist ~ pcoa[,factor]))
print(pairwiseAdonis::pairwise.adonis(dist, pcoa[,factor]))
print(vegan::betadisper(as.dist(dist), pcoa[,factor], type = "centroid"))
sink()

return(centroids)
}