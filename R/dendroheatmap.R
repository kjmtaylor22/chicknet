dendro.heatmap <- function(comm, tax, meta, path="", group, subgroup, core.list=NULL,
                           hi.tax="tag", filename="", stop="no", wd=NULL, ht=NULL){
  #### Series of nested functions to produce grouped heatmap

  id <- grep("sample", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "sample"
  } else {
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
    names(meta)[id[1]] <- "sample"
    id <- id[1]
  }

  if (any(!row.names(comm)%in%meta[,which(names(meta)=="sample")])){
    warning("Check that all community samples are accounted for in metadata.")
    warning("May also need to rename sample ID column in metadata.")
    return(NULL)
  }

  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/profiles/")==FALSE){
    dir.create("./coreplots/profiles/")
  }
  if(dir.exists(paste0("./coreplots/profiles/", group, "_heatmap/"))==FALSE){
    dir.create(paste0("./coreplots/profiles/", group, "_heatmap/"))
  }

  cran <- c("dplyr", "egg", "ape", "vegan", "ggplot2", "reshape2")
  if (any(!cran%in%row.names(installed.packages()))) {
    install.packages(cran[!cran%in%row.names(installed.packages())])
  }
  biocon <- c("ggtree", "ggdendro", "dendextend")
  if (any(!biocon%in%row.names(installed.packages()))){
    #source("https://bioconductor.org/biocLite.R")
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")}
    for (i in biocon){
      #biocLite(i, suppressUpdates=T)
      BiocManager::install(i)
    }
  }

  library(dplyr)
  library(egg)
  library(ape)
  library(vegan)
  library(ggplot2)
  library(ggtree)
  library(ggdendro)
  library(dendextend)

  ra <- function(x){ #where x is a vector
    y <- sum(x)
    z <- sapply(x, FUN=function(x){z1 <- ((x)/y)*100}) ## the x+1 solves the log transformation problem later without skewing the data too much
    return(z)
  }

  comm <- apply(comm, MARGIN=1, FUN=ra) %>% t() %>% as.data.frame()

  sub.tax <- function(comm, tax, meta, group, subgroup, core.list){


    if (is.null(core.list)==FALSE){
      if (class(core.list)=="numeric"){
        otus <- names(sort(colMeans(comm), decreasing = T))[1:core.list]
      } else {

        test <- grep("^asv", core.list)
        test2 <- grep("asv", core.list)

        if (length(test) < length(core.list) & length(test2) > 0){
          if (length(test2)==length(core.list)){
            if (dim(tax)[2]>12){
              otus <- tax$tag[which(tax$tag.name%in%as.character(unique(core.list)))]
            } else {otus <- tax$tag[which(tax$tag.name%in%as.character(unique(core.list)))]}
          } else {
            otus <- tax$tag[which(tax$otu.name%in%as.character(unique(core.list)))]
          }
          otus <- otus[which(duplicated(otus)==FALSE)]
          otus <- otus[which(otus %in% names(comm))]
        }
        if (length(test)==length(test2)){
          otus <- as.character(unique(core.list))
          otus <- otus[which(otus %in% names(comm))]
        }
        if (length(test2)==0){
          for (i in 1:dim(tax)[2]){
            t <- all(core.list %in% tax[,i])
            if (t==TRUE){break}
          }
          otus <- tax$tag[which(tax[,i] %in% as.character(unique(core.list)))]
          otus <- otus[which(otus %in% names(comm))]
        }
        comm <- comm[,which(names(comm) %in% otus)]
      }
    } else {otus <- colnames(comm)}

    sub <- comm

    #if (trans=="none"){
    #  rank1 <- sub
    #}
    #if (trans=="rank"){
    #  rank1 <- rank(unlist(sub))
    #  try <- matrix(rank1, nrow=dim(sub)[1], ncol=length(otus),
    #                dimnames=list(row.names(sub), otus))
    #  rank1 <- as.data.frame(try)
    #}

    #if (trans=="log"){
    #  rank1 <- log10(sub)
    #  rank1 <- apply(rank1, MARGIN=c(1,2), FUN=function(x){if (is.infinite(x)){y <- NA}else{y <- x}}) %>% as.data.frame()
    #}


    #sub <- data.frame(ID=row.names(rank1), rank1)

    sub <- data.frame(ID=row.names(sub), sub)

    sub <- right_join(meta, sub, by=c("sample"="ID"))

    if (length(unique(sub[,group]))==1|length(unique(sub[,subgroup]))==1){

      sub <- mutate(sub, !!"group.subgroup" :=
                      paste(eval(parse(text=group)), eval(parse(text=subgroup)), sample, sep="."))

    } else {
      sub <- mutate(sub, !!"group.subgroup" :=
                      paste(eval(parse(text=group)), eval(parse(text=subgroup)), sep="."))
    }
      sub <- subset(sub, select=c("group.subgroup", otus))

      row.names(sub) <- row.names(comm)


    top <- sub
    return(top)

  }




  ## function for creating cluster dendrogram for the groups and subgroups
  group.dendro <- function(data, meta, tax, path, group, subgroup){

    columns <- colnames(data)
    rows <- row.names(data)

    otus <- names(data)[grep("^asv", names(data))]
    n <- length(otus)

    new <- data.frame(ID=row.names(data), data)

    g <- which(names(meta)==group)
    sg  <- which(names(meta)==subgroup)

    data <- left_join(new, meta[,c(id, g, sg)], by=c("ID"="sample"))
    row.names(data) <- row.names(new)
    #data <- data[,-1]
    rm(new)

    if (length(unique(data[,group]))==1|length(unique(data[,subgroup]))==1){
      data <- mutate(data, !!"group.subgroup" :=
                       paste(eval(parse(text=group)), eval(parse(text=subgroup)), ID, sep="."))
    } else {
      data <- mutate(data, !!"group.subgroup" :=
                      paste(eval(parse(text=group)), eval(parse(text=subgroup)), sep="."))
    }

    data <- data[,c(which(names(data)%in%c("group.subgroup", otus)))]
    row.names(data) <- rows

    g.sg <- names(data)[-grep("^asv", names(data))]
    ngsg <- length(unique(with(data, eval(parse(text=g.sg)))))



    avgd <- data %>%
      group_by_(g.sg) %>%
      summarize_all(funs(mean=mean)) %>%
      as.data.frame()
    row.names(avgd) <- avgd[,1]
    avgd <- subset(avgd, select=grep("^asv", names(avgd)))
    names(avgd) <- otus

    avgd <- apply(avgd, MARGIN=1, FUN=ra) %>% t() %>% as.data.frame()

    tree.comm <- t(avgd) %>% `row.names<-`(tax$den.otu[match(row.names(.), tax$tag)])
    tips <- setdiff(ape::read.tree(path)$tip.label,row.names(tree.comm))
    drop.tree <- ape::drop.tip(ape::read.tree(path), tip=tips)
    bdiv <- GUniFrac::GUniFrac(t(tree.comm), drop.tree, alpha=1)$unifracs[,,1]

    #bdiv <- vegdist(avgd, "bray")

    clust <- hclust(d=as.dist(bdiv), method="average")

    dend <- as.dendrogram(clust) %>%
      dendextend::rotate(sort(as.character(unique(data$group.subgroup))))

    ddend <- dendro_data(dend, type="rectangle")

    order <- ddend$labels$label
    zoom <- ngsg*(0.0595-0.0001*ngsg)

    d <- ggplot(segment(ddend)) +
      geom_segment(aes(x=x, y=y, xend=xend, yend=yend), size=2) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = margin(0,0,0,0, "pt")) +
      coord_cartesian(xlim=c(zoom, ngsg+1-zoom),
                      ylim=c(max(ddend$segments$yend)/10,
                             max(ddend$segments$yend)))
    d[["tip.labels"]] <- order

    return(d)
  }


  ## function for creating phylogenetic cladogram for top taxa
  otu.phylo <- function(top, tax, path, hi.tax){#tax= output from bact.tax

    g.sg <- names(top)[-grep("^asv", names(top))]
    otus <- names(top)[grep("^asv", names(top))]
    ngsg <- length(unique(with(top, eval(parse(text=g.sg)))))

    tree <- ape::read.tree(path)

    tips <- as.character(tax$den.otu[-match(otus, tax$tag)])

    prune <- match(tax$den.otu, tree$tip.label)

    if (length(prune)>0){
      if(any(is.na(prune))){
        prune <- prune[which(!is.na(prune))]
      }
      tips <- c(tips, tree$tip.label[-prune])
    }

    drop.tree <- ape::drop.tip(tree, tip=tips)

    labs <- data.frame(tips=drop.tree$tip.label)
    labs <- left_join(labs, tax, by=c("tips"="den.otu"))

    if (hi.tax%in%names(tax)){
      get <- which(names(labs)==hi.tax)
      drop.tree$tip.label <- as.character(unlist(labs[,get]))
    } else {
      drop.tree$tip.label <- as.character(labs$otu.name)
    }

    p <- ggtree::ggtree(drop.tree, branch.length = "none", ladderize=FALSE) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = margin(0,0,0,0, "pt"))

    pull <- data.frame(node=p$data$node, label=p$data$label)
    dup <- which(duplicated(pull$label[1:length(otus)])==TRUE)
    new <- drop.tip(drop.tree, tip=pull$node[dup])


    n <- length(new$tip.label)

    if (n<150){
      zoom <- n*(0.0595-0.0001*n)
    }
    if (n>150){
      zoom <- n*(0.0595-0.000045*n)
    }

    cp <- ggtree::ggtree(new, branch.length = "none", size=2) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = margin(0,0,0,0, "pt")) +
      coord_cartesian(ylim=c(zoom, n+1-zoom))

    new <- cp$data[1:n, c(4,5,7)]

    out <- new[order(new$y),]


    cp[["new.labels"]] <- out$label

    return(cp)
  }

  ## function for creating heatmap using the top n otus
  top.heat <- function(top, tax, phylo, dendro, hi.tax, filename){

    g.sg <- names(top)[-grep("^asv", names(top))]
    otus <- names(top)[grep("^asv", names(top))]

    group1 <- top %>%
      group_by_(g.sg) %>%
      summarize_all(funs(mean(.,na.rm=T))) %>%
      as.data.frame()

    row.names(group1) <- group1[,1]
    group1 <- group1[,-1]

    group1 <- group1[match(dendro$tip.labels, row.names(group1)),]

    group2 <- reshape(group1, varying=otus,
                      v.names="abund",
                      timevar="otu",
                      times=otus,
                      direction="long",
                      ids=row.names(group1))

    group2 <- left_join(group2, tax, by=c("otu"="tag"))

    if (hi.tax%in%names(tax)){
      get <- which(names(group2)==hi.tax)
      group3 <- group2[,c(2, 3, get)] %>%
        group_by_("id", hi.tax) %>%
        summarize_at("abund", funs(abund=sum(.,na.rm=T))) %>%
        as.data.frame()
      names(group3)[which(names(group3)==hi.tax)] <- "otu.name"
    } else {
      group3 <- group2[,c(2, 3, 6)] %>%
        group_by(id, otu.name) %>%
        summarize_at("abund", funs(abund=max)) %>%
        as.data.frame()

    }

    #if (trans=="none"){group3 <- group3}

   # if (trans=="rank"){group3$abund <- rank(group3$abund)}

    #if (trans=="log"){
      group3$abund <- log10(group3$abund)
      group3$abund <- sapply(group3$abund, FUN=function(x){if (is.infinite(x)){y <- NA}else{y <- x}})
    #}

    group3$id <- factor(group3$id,
                        levels=dendro$tip.labels)

    group3$otu.name <- factor(group3$otu.name,
                              levels=phylo$new.labels)
    write.csv(group3, paste0("./coreplots/profiles/", group, "_heatmap/heatmapdata_", filename, ".csv"))

    min.lim <- floor(min(group3$abund,na.rm = T))

    brks <- seq(-2, 2, 1)

    colors <- rev(c("magenta", "red", "orange", "yellow", "green", "turquoise", "blue", "navy"))

    if (min.lim < brks[1]){
      colors <- c(rep("black", abs(min.lim - brks[1])+1), colors)
      brks <- seq(min.lim, 2, 1)
    }

    lbs <- paste0(10^brks, "%")

    lims <- c(min(brks), max(brks))


    h <- ggplot(group3, aes(id, otu.name)) +
      geom_raster(aes(fill=abund)) +
      scale_fill_gradientn(colors=colors,
                            guide=guide_colorbar(title=bquote(paste("Relative Abundance (", "log"[.(10)]," scale)")), title.position = "right",
                                                 ticks.colour="black", frame.colour="black", label.position = "left",
                                                 title.theme = element_text(size=55, angle=-90, hjust=0.5)),
                            limits=lims, na.value = "white", breaks=brks, labels=lbs) +
      theme(axis.ticks=element_line(color="white"),
            axis.text.x=element_text(angle=90, size=60, hjust=1, vjust=0.5, face="bold"),
            axis.text.y=element_text(size=60, vjust=0.5), axis.title = element_blank(),
            plot.margin = margin(1,1,5,1, "pt"),
            plot.caption=element_text(size=45),
            legend.text=element_text(size=50, face = "bold"),
            legend.direction = "vertical",
            legend.key.height =unit(0.13, "npc"),
            legend.key.width =unit(3, "lines")) +
      scale_y_discrete(position="right") +
      labs(caption=paste("Features: n = ", as.character(length(levels(group3$otu.name)))))

    if (stop=="data"){return(group3)}

    return(h)
  }

  top <- sub.tax(comm, tax, meta, group, subgroup, core.list)
  if (stop=="top"){return(top)}
  den <- group.dendro(comm, meta, tax, path, group, subgroup)
  if (stop=="den"){return(den)}
  phy <- otu.phylo(top, tax, path, hi.tax) #use hi.tax if you want core OR all taxa at level other than in `sub`
  if (stop=="phy"){return(phy)}
  heat <- top.heat(top, tax, phy, den, hi.tax, filename) #use hi.tax if you want core OR all taxa
  if (stop=="heat" | stop=="data"){return(heat)}

  e <- ggplot() +
    geom_blank() +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = margin(0,0,0,0, "pt"))

  n <- length(phy$new.labels)
  n2 <- length(den$tip.labels)

  if (is.null(wd)){
    wd <- 1000+65*n2
  }
  if (is.null(ht)){
    if (n>=50){
      ht <- 56*n
    }
    if (n<50 & n>=25){
      ht <- 800+65*n
    }
    if (n<25){
      ht <- 2000
    }
  }

  hts <- c(0.4-0.003*n, 3)
  wds <- c(0.5-0.003*n2, 3)


  tiff(paste0("./coreplots/profiles/", group, "_heatmap/phyloheatmap_", filename, ".tiff"), width=wd, height=ht)
  egg::ggarrange(e, den, phy, heat, heights=hts, widths=wds)
  dev.off()

}
