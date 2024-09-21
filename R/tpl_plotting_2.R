#' @export

tpl.plot <- function(data, set, block, group, variable, hi.tax="tag.name", filename=Sys.Date(), tpl=TRUE,
                     subsample=FALSE, tax=NULL, meta=NULL){

  # variable = numeric data from which mean and variance are calculated
  # group = most nested; the grouping variable by which to segregate samples into sets for direct comparison
  # block = middle nest; sets of samples that SHOULD NOT be directly compared
  # set = largest nest; different experiments that SHOULD NOT be directly compared
  require(dplyr)
  require(ggplot2)

  if (all(sapply(data, FUN=is.numeric))==TRUE &
      is.null(tax)==FALSE & is.null(meta)==FALSE){

    grouped <- data.frame()
    for (i in as.character(unique(meta[meta$SampleID%in%row.names(data),set]))){

      pooled <- data[row.names(data)%in%meta$SampleID[meta[,set]==i],] %>%
        .[,colSums(.)>0] %>% t() %>% data.frame(otus=row.names(.)) %>%
        left_join(tax[,c("tag",hi.tax)], by=c("otus"="tag"))  %>%
        .[,-(dim(.)[2]-1)] %>% group_by_(hi.tax) %>%
        summarize_all(funs(sum)) %>% as.data.frame()

      otus <- as.character(unlist(pooled[,1]))


      sub.pooled <- pooled[,-1] %>% .[rowSums(.)>0,]
      row.names(sub.pooled) <- otus

      sub.pooled <- data.frame(OTU=row.names(sub.pooled), sub.pooled)

      drop.empty.asv <- reshape2::melt(sub.pooled, variable.name="SampleID", value.name="abund")

      grouped <- rbind(grouped, drop.empty.asv)
    }

    if (is.null(variable)){
      variable <- "OTU"
      joined <- right_join(meta[,c("SampleID",set,block,group)], grouped)
    }

    if (is.null(group)){
      group <- "OTU"
      joined <- right_join(meta[,c("SampleID",set,block,variable)], grouped)
    }

    data <- joined
  }


  data <- data[!is.na(data[,set]),]

  data[,block] <- as.factor(data[,block])

  un.set <- as.character(unique(data[,set]))

  un.block <- as.character(levels(data[,block])[levels(data[,block])%in%unique(data[,block])])

  gg_color_hue <- function(n) {
    hues = seq(20, 380, length = n+1)
    hcl(h = hues, l = 70, c = 140)[1:n]}

  clr <- gg_color_hue(length(un.block))
  names(clr) <- un.block

  shps <- c(7,8,18,2,16,3:6)
  names(shps) <- un.block

  parse.list <- list()
  for (i in un.set){
    exp <- list()
    for (j in un.block){
      blc <- as.character(unique(data[data[,block]==j & data[,set]==i, group]))
      exp[[j]] <- blc
    }
    parse.list[[i]] <- exp
  }

  out <- data.frame()
  out2 <- data.frame()
  for (k in 1:length(parse.list)){
    test <- data.frame()
    for (j in 1:length(parse.list[[k]])){
       x <- subset(data, subset=eval(parse(text=block))==names(parse.list[[k]][j]) &
                     eval(parse(text=set))==un.set[k])
       if (dim(x)[1]>0){
         t <- x %>%
           group_by_(group) %>%
           summarize(m=mean(abund, na.rm=T),
                     v=var(abund, na.rm=T),
                     count=length(abund))

         t <- cbind(t, set=rep(names(parse.list[k]),
                               length(parse.list[[k]][[j]])),
                      block=rep(names(parse.list[[k]][j]),
                                 length(parse.list[[k]][[j]])))

         test <- rbind(test, t)
       }
    }
    test <- mutate(test, se=sqrt(v)/sqrt(count)) %>%
      mutate(lowCI95=m-1.96*se, uppCI95=m+1.96*se) %>%
      mutate(errorbar=abs(uppCI95)-abs(lowCI95))

    test <- data.frame(test, errorbar0=NA, errorbar1=NA)
    test$errorbar0[which(test$errorbar>0)] <- test$m[which(test$errorbar>0)]
    test$errorbar1[which(test$errorbar>0)] <- test$uppCI95[which(test$errorbar>0)]
    test$errorbar1[which(test$errorbar<0)] <- test$m[which(test$errorbar<0)]
    test$errorbar0[which(test$errorbar<0)] <- test$lowCI95[which(test$errorbar<0)]

    names(test) <- c("group", names(test)[-1])

    if (subsample==1){subsample <- F}

    if (subsample!=F){
      if (subsample==T){
        ss <- summary(test$block)
        choose <- which(ss==min(ss))
        test2 <- data.frame()
        for (i in unique(test$block)){
          test3 <- sample_n(test[test$block==i,], ss[choose])
          test2 <- rbind(test2, test3)
        }
      }
      if (is.numeric(subsample) & subsample < 1){
        test2 <- data.frame()
        for (i in as.character(unique(test$block))){
          test3 <- sample_frac(test[test$block==i,], subsample)
          test2 <- rbind(test2, test3)
        }
      }
      if (is.numeric(subsample) & subsample > 1) {stop("subsample is not a ratio between 0 and 1")}
      test <- test2
      rm(test3, test2)
    }


    if(tpl==TRUE){

      pops <- unique(test$block)

      shp <- shps[names(shps)%in%pops]

      path <- paste("coreplots/TPLtest",set,block,group,variable,filename,sep="/")

      if(!dir.exists(paste(path,"boxplots",sep="/"))){dir.create(paste(path,"boxplots",sep="/"), recursive=T)}


      write.csv(test, paste(path,"/tpl_by-",un.set[k], ".csv", sep=""))

      if (length(which(test$m!=0))==0) {next}

      pull <- which(test$v==0 | is.na(test$v))
      if (length(pull)  > 0){test <- test[-pull,]}

      pops <- sort(as.character(pops))

      jpeg(file=paste0(path,"/tpl_by-", un.set[k], ".jpeg"), width=1300, height=1500, res=450)
      print(ggplot(test) + geom_point(aes(m,v, color=block, shape=block), size=2) +
              scale_shape_manual(values=shp) + scale_color_manual(values=clr) +
              stat_smooth(aes(m, v, group=block, color=block), method="lm", formula=y ~ x, se=F, size=0.5, inherit.aes = F, show.legend = F) +
              geom_abline(aes(slope=1, intercept=0), linetype="dashed", color="grey40") +
              ggthemes::theme_few() + theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank()) +
              guides(color=guide_legend(title=NULL, ncol=2), shape=guide_legend(title=NULL, ncol=2)) +
              scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x))) +
              scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x))) +
              ggpmisc::stat_poly_eq(inherit.aes = F, formula=y ~ x, parse=T, label.x = "left", label.y="top", size=3,
                                    aes(x=m, y=v, group=block, color=block, label=paste(stat(eq.label), stat(rr.label), sep="~~~~"))) +
              ggpmisc::stat_fit_glance(inherit.aes=T, method = 'lm', method.args = list(formula = y ~ x), label.x = 'right', label.y = "bottom", size = 3,
                                       aes(x=m, v, group=block, color=block, label = paste0("p = ", signif(..p.value.., digits = 2)))))
      dev.off()

      tpl.model <- data.frame()
      for (i in 1:length(pops)){
        test1 <- subset(test, subset=block==pops[i])

        if (dim(test1)[1] %in% c(0,1,2)){next}

        tpl.seg <- lm(log10(v) ~ log10(m), data=test1)
        print(summary(tpl.seg))
        sum <- summary(tpl.seg)
        mod <- data.frame(set = un.set[k], block = pops[i],
                          interc = tpl.seg$coefficients[[1]], slope = tpl.seg$coefficients[[2]],
                          se = sum[[4]][4], df = sum$df[2], Rse = sum$sigma, mR2 = sum$r.squared,
                          aR2 = sum$adj.r.squared, Fstat = sum$fstatistic[1], Fdf1 = sum$fstatistic[2],
                          Fdf2 = sum$fstatistic[3], t.val = sum[[4]][6], p.val = sum[[4]][8])
        tpl.model <- rbind(tpl.model, mod)
      }

      out <- rbind(out, tpl.model)

      jpeg(paste0(path,"/boxplots/", un.set[k], ".jpg"), 4000, 2300, res=600)
      par(mfrow=c(1,2), cex.axis=1.1, font.axis=2, las=3, xaxt="n")
      boxplot(log10(m) ~ block, data=test, boxwex=0.8,
              boxfill=clr, ylab="log(m)", main=un.set[k])
      axis(1, at=c(1:length(levels(test$block))), labels=FALSE)
      text(c(1:length(levels(test$block))), par("usr")[3], font=2,
           labels=c(levels(test$block)), srt=45, adj=c(1, 1.5), xpd=T, offset=2.5)
      boxplot(log10(v) ~ block, data=test, boxwex=0.8,
              boxfill=clr, ylab="log(v)", main=un.set[k])
      axis(1, at=c(1:length(levels(test$block))), labels=FALSE)
      text(c(1:length(levels(test$block))), par("usr")[3], font=2,
           labels=c(levels(test$block)), srt=45, adj=c(1, 1.5), xpd=T, offset=2.5)
      dev.off()


    }
      out2 <- rbind(out2, test)
  }
  if (tpl==TRUE){
    write.csv(out, paste(path,"/tpl_trendlinedata.csv", sep=""))
  }
  return(out2)
}
