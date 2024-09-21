core.stack <- function(data, list, tax=NULL, factors, group, subgroup, hi.tax="tag", threshold=NULL, margin,
                       date, core.only=FALSE, fixed=FALSE, landscape=TRUE, view.all=FALSE, facet=FALSE,
                       reverse=FALSE){

  cran <- c("dplyr", "RColorBrewer", "ggplot2")
  if (any(!cran%in%row.names(installed.packages()))) {
    install.packages(cran[!cran%in%row.names(installed.packages())])
  }
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)


  if(dir.exists("./coreplots/")==FALSE){
    dir.create("./coreplots/")
  }
  if(dir.exists("./coreplots/profiles/")==FALSE){
    dir.create("./coreplots/profiles/")
  }
  if(dir.exists(paste0("./coreplots/profiles/", group, "_", date, "/"))==FALSE){
    dir.create(paste0("./coreplots/profiles/", group, "_", date, "/"))
  }

  fmx <- as.data.frame(t(data)) %>%
    summarize_all(funs(sum=sum)) %>%
    max()

  ra <- function(x){ #where x is a vector
    y <- sum(x)
    z <- sapply(x, FUN=function(x){z1 <- x/y})
    return(z)
  }
  avg <- function(x){
    z <- mean(x, na.rm=TRUE)
    return(z)
  }
  gcd <- function(x,y) {
    r <- x%%y
    return(ifelse(r, gcd(y, r), y))
  }

  id <- grep("sample", names(factors), ignore.case = T)
  if (length(id)==1){
    names(factors)[id] <- "SampleID"
  } else {
    warning(paste0("Sample ID not found. Using factor column 1 (", names(factors)[1], ") as sample ID. "))
    names(factors)[1] <- "SampleID"
  }

  if (any(!row.names(data)%in%factors[,which(names(factors)=="SampleID")])){
    warning("Check that all community samples are accounted for in metadata.")
    warning("May also need to rename sample ID column in metadata.")
    return(NULL)
  }

  features <- colnames(data)


  new <- data.frame(SampleID=row.names(data), data)
  data <- right_join(factors[,c("SampleID", group, subgroup)], new)
  data <- mutate(data, !!"supergroup" :=
                   paste(eval(parse(text=group)), eval(parse(text=subgroup)), sep="."))

  gsg <- unique(data[,c(group, subgroup, "supergroup")])

  cat("Calculating as relative abundance (%) \n")

  data <- data[,-c(1:3)] %>%
    group_by(supergroup) %>%
    summarize_all(funs(sum=sum))

  supergroup <- data[,1]


  data <- apply(data[,-1], MARGIN=1, FUN=ra) %>% t() %>% data.frame()


  names(data) <- features
  data <- cbind(supergroup, data)

  if (is.null(threshold)==TRUE){
    t.av <- function(x){
      y <- which(x>0)
      z <- mean(x[y], na.rm=TRUE)
      return(z)
    }
    new <- apply(data, MARGIN=1, FUN=t.av)
    split <- strsplit(names(new), split=NULL)
    splt <- function(x){
      y <- grep("[A-Z]", x)
      if (y[1]==y[length(y)]){
        z <- x[y[1]]
      } else {
        if (x[y[1]]=="B"){
          z <- x[y[1]]
        } else {
          z <- paste(x[y[1]], x[y[length(y)]], sep="")
        }
      }
      return(z)
    }
    fromgroup <- lapply(split, splt) %>% as.data.frame() %>% t()
    new <- data.frame(new, fromgroup=fromgroup)
    new <- summarize(group_by(new, fromgroup), mean=mean(new))
    thresh <- new$mean
    rm(new, fromgroup, split)
  }

  data <- right_join(gsg, data, by=c("supergroup"="supergroup"))

  df <- data.frame()
  for (i in 1:length(list)){
    if (exists("thresh")){
      threshold <- thresh[i]
      threshold <- round(threshold, digits=4)
    }

    cat(paste("Calculating core microbiome for ", names(list)[i], "\n", sep=""))
    group.core <- data.frame()

    if (all(is.na(unlist(list[[i]])))){next}



    core.pool <- data.frame(count=rep(1, length(unlist(list[i]))),
                            otu=unlist(list[i]))
    core <- core.pool %>% group_by(otu) %>% summarize_at(vars(count), funs(sum))

    core <- as.character(core$otu[which(core$count==length(which(is.na(list[[i]])==FALSE)))])


    if (length(core)>0 & length(which(is.na(list[[i]])==FALSE))>=1){
      sub <- data[data[,group]==names(list[i]), names(data)%in%c(subgroup, core)]

      core <- names(sub)[which(names(sub)%in%core)]

      if (is.null(tax)==FALSE){
        core <- data.frame(v1=c(1:length(core)), v2=core)
        core <- left_join(core, tax, by=c("v2"="tag"))

        if (hi.tax%in%names(tax)){

          get <- which(names(core)==hi.tax)

          core[,hi.tax] <- as.character(core[,hi.tax])

          if (any(core[,hi.tax]=="")){
            core[which(core[,hi.tax]==""),hi.tax] <- paste("unclassified", core[which(core[,hi.tax]==""),get-1])
          }

          core <- as.character(unlist(core[,hi.tax]))
        } else {core <- as.character(core$otu.name)}

        grouped <- data.frame(features=core, t(sub[,-1]))

        grouped <- grouped %>%
          group_by(features) %>%
          summarize_all(funs(sum=sum))

        names(grouped) <- c("features", unlist(as.character(sub[,1])))

        core <- as.character(unlist(grouped$features))
        grouped <- grouped[,-1]
        grouped <- t(grouped) %>% as.data.frame()
        names(grouped) <- core

        grouped <- cbind(v1=as.character(sub[,1]), grouped)

        names(grouped) <- c(subgroup, names(grouped)[-1])
        row.names(grouped) <- grouped[,1]
      }

      if (is.null(tax)==TRUE){
        grouped <- sub
      }


      grouped <- reshape(grouped, varying=core,
                         v.names="rel.abund",
                         timevar="otu",
                         times=core,
                         direction="long")

      grouped <- data.frame(grouped, groupvar=rep(names(list)[i], dim(grouped)[1]))
      names(grouped) <- c(subgroup, names(grouped)[-1])

      grouped$otu <- as.factor(grouped$otu)
      grouped$otu <- factor(grouped$otu, levels = rev(levels(grouped$otu)))
      pull <- grouped$otu[which(grouped$rel.abund > threshold)]
      pull <- pull[which(duplicated(pull)==FALSE)]
      grouped <- grouped[which(grouped$otu%in%pull),]
      group.core <- cbind(grouped, core=rep("core", dim(grouped)[1]))
    }

    if (core.only==FALSE){

      if (reverse==F){
        j.sets <- 2:length(list[[i]])
        end <- length(list[[i]])
        start <- 1
        j.stop <- expression(j-1)
      } else {
        j.sets <- (length(list[[i]])-1):1
        end <- 1
        start <- length(list[[i]])
        j.stop <- expression(j+1)
      }

      for (j in j.sets){
        if (is.na(list[[i]][j])==FALSE){
          pull <- which(unlist(list[[i]][j:end])%in%unlist(list[[i]][start:eval(j.stop)])==FALSE)
          pull <- unlist(list[[i]][j:length(list[[i]])], use.names = FALSE)[pull]
          pool <- data.frame(count=rep(1, length(pull)), otu=pull)

          otus <- names(which(summary(as.factor(pool$otu), maxsum=5000)== ##Note: this line means that the function ignores if some age groups don't meet the occurence margin parameter
                                length(which(is.na(list[[i]][j:end])==FALSE))))
          #cat(paste0(names(list)[i], ": ", names(list[[i]][j]), "-->", names(list[[i]][length(list[[i]])])), otus, sep="\n")
          if ("(Other)"%in%otus){otus <- otus[-length(otus)]}
          if (length(otus)>0 & length(which(is.na(list[[i]])==FALSE))>=1){
            sub <- data[which(data[,group]==names(list[i])),
                        which(names(data)%in%c(subgroup, otus))]

            if (is.null(tax)==FALSE){
              otus <- otus[match(names(sub)[-1], otus)]
              otus <- data.frame(v1=c(1:length(otus)), v2=otus)
              otus <- left_join(otus, tax, by=c("v2"="tag"))


              if (hi.tax%in%names(tax)){
                get <- which(names(otus)==hi.tax)

                otus[,get] <- as.character(otus[,get])

                if (any(otus[,hi.tax]=="")){
                  otus[which(otus[,hi.tax]==""),hi.tax] <- paste("unclassified", otus[which(otus[,hi.tax]==""),get-1])
                }

                otus <- unlist(otus[,get])

              } else {otus <- as.character(otus$otu.name)}

              #if (length(otus)!=dim(t(sub[,-1]))[1]){return(list(features=otus, comm=t(sub[,-1])))}
              #if (length(otus)!=dim(t(sub[,-1]))[1]){return(summary(pool$otu))}

              #otus <- otus[!is.na(otus)] ##This prevents the summary(pool$otu) from trying to return "(Other)" as an otu by accident
              grouped <- data.frame(features=otus, t(sub[,-1]))

              grouped <- grouped %>%
                group_by(features) %>%
                summarize_all(funs(sum=sum))

              otus <- as.character(unlist(grouped$features))
              grouped <- grouped[,-1]
              grouped <- t(grouped) %>% as.data.frame()
              names(grouped) <- otus


              grouped <- cbind(v1=as.character(sub[,1]), grouped)

              names(grouped) <- c(subgroup, names(grouped)[-1])
              row.names(grouped) <- grouped[,1]
            }
            if (is.null(tax)==TRUE){
              grouped <- sub
            }
            #if (i==5 & j==4){return(pool)}
            grouped <- reshape(grouped, varying=otus,
                               v.names="rel.abund",
                               timevar="otu",
                               times=otus,
                               direction="long")

            grouped <- data.frame(grouped, groupvar=rep(names(list)[i], dim(grouped)[1]))
            names(grouped) <- c(subgroup, names(grouped)[-1])
            grouped$otu <- as.factor(grouped$otu)
            grouped$otu <- factor(grouped$otu, levels = rev(levels(grouped$otu)))
            pull <- grouped$otu[which(grouped$rel.abund > threshold)]



            if(length(pull)>0){
              pull <- pull[which(duplicated(pull)==FALSE)]
              grouped <- grouped[which(grouped$otu%in%pull),]
              grouped <- cbind(grouped, core=paste(names(list[[i]][j]), "-->", names(list[[i]][end]), sep=""))
              if (dim(grouped)[1]==0){next}
            }
            if (length(pull)==0){
              grouped <- group.core[1,]
              grouped$rel.abund[1] <- 0
              grouped$core <- paste(names(list[[i]][j]), "-->", names(list[[i]][end]), sep="")
              if (dim(grouped)[1]==0){next}
            }
            if (dim(grouped)[2]<6){next}
            group.core <- rbind(group.core, grouped)
          } else {
            grouped <- group.core[1,]
            grouped$rel.abund[1] <- 0
            grouped$core <- paste(names(list[[i]][j]), "-->", names(list[[i]][end]), sep="")
            if (dim(grouped)[1]==0|dim(grouped)[2]<6){next}
            group.core <- rbind(group.core, grouped)
          }
        } else {
          grouped <- group.core[1,]
          grouped$rel.abund[1] <- 0
          grouped$core <- paste(names(list[[i]][j]), "-->", names(list[[i]][end]), sep="")
          if (dim(grouped)[1]==0|dim(grouped)[2]<6){next}
          group.core <- rbind(group.core, grouped)
        }
      }
    }

    #if (dim(group.core)[2]<6){next}
    df <- rbind(df, group.core)

  }

  colors <- c("#330033", "#F0E442", "#0072B2", "#FF0033", "#FF9933", "#999999", "#56B4E9", "#FFCC00", "#00FF00", "#FF66FF", "#666666")
  colors <- colorRampPalette(colors)(length(unique(df$otu)))
  names(colors) <- sort(as.character(unique(df$otu)))

  lev <- match(levels(factors[,subgroup]), unique(df[,subgroup])) %>% .[which(!is.na(.))]
  df[,subgroup] <- factor(df[,subgroup], levels=unique(df[,subgroup])[lev])

  df$otu <- as.character(df$otu)

  lev <- match(levels(factors[,group]), unique(df$groupvar)) %>% .[which(!is.na(.))]
  df$groupvar <- factor(df$groupvar, levels=unique(df$groupvar)[lev])

  df$core <- as.factor(df$core)
  df$core <- factor(df$core, levels=unique(df$core))

  for (i in unique(df$groupvar)){

    group.core <- subset(df, groupvar==i)

    #group.core[,subgroup] <- factor(group.core[,subgroup], levels=unique(factors[,subgroup]))


    face <- c(floor(length(unique(group.core$core))/3), ceiling(length(unique(group.core$core))/3))
    kh <- 1

    if (face[2]>face[1]){
      if (length(unique(group.core$otu)) > 25){
        cols <- floor(length(unique(group.core$otu))/25)
      } else {cols <- 1}
      leg <- list(c(0.833, 0.22), "center")
    } else {
      cols <- 4
      leg <- list("bottom", "center")
    }

    if (dim(group.core)[1]==0){next}
    group.core <- group.core[which(group.core$core%in%names(summary(group.core$core)[which(summary(group.core$core)!=1)])),]
    #group.core$core <- factor(group.core$core)

    fexp <- expression(facet_grid(core~., drop=TRUE))
    fexp2 <- expression(facet_grid(groupvar~core, drop=FALSE))

    leg <- list("right", "right")

    if (core.only==FALSE){
      if (length(unique(group.core$otu)) > 40){
        kh <- floor(length(unique(group.core$otu))/40)
      }
      if (length(unique(group.core$otu))/length(unique(group.core$core)) > 8){
        kh <- length(unique(group.core$otu))/8
      }
    }
    cols <- 1




    if (fixed==TRUE){
      fix.ylim <- expression(ylim(0,100))
    }
    if (fixed==FALSE){
      fix.ylim <- expression(ylim(0, NA))
    }
    if (fixed=="free"){
      fix.ylim <- expression(facet_wrap(vars(core), dir="v", scales="free_y", strip.position = "right", drop=T))
    }



    group.core$rel.abund <- group.core$rel.abund*100

    #lev <- match(levels(with(meta, eval(parse(text=subgroup)))),
    #             unique(with(group.core, eval(parse(text=subgroup)))))
    #lev <- lev[which(!is.na(lev))]

    #group.core[,which(names(group.core)==subgroup)] <-
    #  factor(group.core[,which(names(group.core)==subgroup)],
    #         levels=levels(with(group.core, eval(parse(text=subgroup))))[
    #           lev[which(is.na(lev)==FALSE)]])

    w <- 500*length(unique(group.core[,1]))+1000
    r <- expression(300)
    fs <- w/300
    fs <- round(c(fs+(fs/2), fs+(fs/3.2), fs+(fs/4), fs+(fs/8), fs+(fs/16), fs*2))

    if (view.all==T){
      fexp <- expression(facet_grid(groupvar~., drop=TRUE))
      if (length(unique(df[,subgroup]))<3){#length(unique(df$groupvar))){
        fexp2 <- expression(facet_grid(eval(parse(text=subgroup))~groupvar, drop=TRUE))
      } else {fexp2 <- expression(facet_grid(groupvar~., drop=TRUE))}
      h <- 2300
    } else {
      h <- 800*length(unique(group.core$core))+1200
    }

    #group.core$core <- factor(group.core$core, levels=unique(group.core$core))

    colors2 <- colors[match(unique(group.core$otu), names(colors))]

    f <- ggplot(group.core, aes(eval(parse(text=subgroup)), rel.abund, fill=otu)) +
      eval(expr=fexp) +
      geom_bar(stat="identity") +
      eval(expr=fix.ylim) +
      theme_minimal() +
      scale_fill_manual(values=colors2) +
      theme(panel.border=element_rect(linetype="solid", fill=NA),
            legend.position=leg[[1]],
            legend.justification=leg[[2]],
            panel.grid.major.x = element_line(color="grey50", linetype="dotted"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(colour="grey90", linetype="36"),
            panel.grid.minor.y = element_blank(),
            strip.text = element_text(size = fs[3], face="bold"),
            plot.title = element_text(size = fs[1], face="bold"),
            axis.title = element_text(size = fs[2], face="bold"),
            axis.text.x = element_text(color = "black", size = fs[2], angle=45, vjust=0.9, hjust=1),
            axis.text.y = element_text(color = "black", size = fs[2], face="bold"),
            legend.title = element_text(size = fs[3], face="bold"),
            legend.text = element_text(size = fs[3]),
            plot.caption = element_text(size=fs[5]),
            legend.background=element_rect(fill="white", colour="black"))

    jpeg(file=paste("./coreplots/profiles/", group, "_", date, "/", i, date, "_", threshold, "_", as.character(margin), ".jpeg", sep=""),
         width=w, height=h, units="px", res=eval(expr=r) )
    print(f + labs(title=paste("Progression of Core Microbiome Composition in", i)) +
            xlab(subgroup) + ylab("Relative Abundance (%)") +
            guides(fill=guide_legend(keywidth = fs[2]/10, keyheight=fs[2]/(10*kh), ncol=cols)) +
            labs(caption=paste("Occurence margin: ", as.character(margin*100), "%", "\n",
                               "Abundance threshold: ", as.character(threshold*100), "%")))
    dev.off()
  }
  write.csv(df, paste("./coreplots/profiles/", group, "_", date, "/", date, "_",
                      threshold, "_", as.character(margin), ".csv", sep=""))

  if (facet==TRUE){
    if ("BG" %in% unique(df$groupvar)){
      df <- subset(df, subset=groupvar!="BG")
      df$core <- factor(df$core)
    }

    sg <- which(names(df)==subgroup)
    pull <- unique(data.frame(df$groupvar, df[,sg]))
    #pull <- names(which(summary(pull[,1])==length(unique(df[,sg]))))
    #df <- df[which(df$groupvar%in%pull),]
    df$rel.abund <- df$rel.abund*100
    if (landscape==FALSE){
      w <- 4000
      h <- 8000
      legs <- expression(theme(legend.position="right",
                               legend.justification="right"))
      gui <- expression(guides(fill=guide_legend(keywidth = 1, keyheight=0.75, ncol=1)))


    } else {
      w <- 8000
      h <- 4000
      legs <- expression(theme(legend.position="bottom",
                               legend.justification="center"))
      #gui <- expression(guides(fill=guide_legend(keywidth = 1, keyheight=0.75, ncol=7)))

      if (length(unique(df[,subgroup]))<3){
        fexp2 <- expression(facet_grid(groupvar~eval(parse(text=subgroup)), drop=TRUE))
      } else {fexp2 <- expression(facet_grid(.~groupvar, drop=TRUE))}

    }


    bff <- ggplot(df, aes(eval(parse(text=subgroup)), rel.abund, fill=otu)) +
      geom_bar(stat="identity") +
      eval(expr=fix.ylim) +
      eval(expr=fexp2) +
      theme_minimal() +
      scale_fill_manual(values=colors) +
      theme(panel.border=element_rect(linetype="solid", fill=NA),
            legend.position=leg[[1]],
            legend.justification=leg[[2]],
            panel.grid.major.x = element_line(color="grey50", linetype="dotted"),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(colour="grey90", linetype="36"),
            panel.grid.minor.y = element_blank(),
            strip.text = element_text(size = fs[6], face="bold"),
            plot.title = element_text(size = fs[1], face="bold"),
            axis.title = element_text(size = fs[2], face="bold"),
            axis.text.x = element_text(color = "black", size = fs[6], angle=45, hjust=1, vjust=1),
            axis.text.y = element_text(color = "black", size = fs[6], face="bold"),
            legend.title = element_text(size = fs[3], face="bold"),
            legend.text = element_text(size = fs[1]),
            plot.caption = element_text(size=fs[5]),
            legend.background=element_rect(fill="white", colour="black")) +
      eval(expr=legs)
    jpeg(file=paste("./coreplots/profiles/", group, "_", date, "/facet", date, "_", threshold, "_", as.character(margin), ".jpeg", sep=""),
         width=w, res=600, units="px", height=h) #height=h*length(list))
    print(bff + labs(title=paste("Progression of Core Microbiome Composition")) +
            xlab(subgroup) + ylab("Relative Abundance (%)") +
            guides(fill=F)+#eval(expr=gui) +
            labs(caption=paste("Occurence margin: ", as.character(margin*100), "%", "\n",
                               "Abundance threshold: ", as.character(threshold*100), "%")))
    dev.off()

    pleg <- ggpubr::as_ggplot(ggpubr::get_legend(bff + guides(fill=guide_legend(keywidth = 1, keyheight=0.75, ncol=5))))

    jpeg(paste0("./coreplots/profiles/", group, "_", date, "/facet", date, "_", threshold, "_", as.character(margin), "_legend.jpeg"),
         width=w, height=h/2, res=600)
    print(pleg)
    dev.off()
  }
  return(df)
}
