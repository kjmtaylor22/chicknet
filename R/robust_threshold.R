robust.threshold <- function(comm, meta, group, subgroup, select, shapes, colors, filename="ASV_occurrence"){
  
  library(dplyr)
  library(ggplot2)
  
  if (is.null(select)){select <- names(comm)}

  meta <- mutate(meta, !!"group.subgroup" := paste(eval(parse(text=group)), eval(parse(text=subgroup)), sep="."))
  
  grouped <- data.frame(SampleID=row.names(comm), comm) %>% left_join(meta[,c("SampleID", "group.subgroup")]) %>%
    .[,-1] %>% group_by(group.subgroup) %>% summarize_all(funs(mean)) %>% as.data.frame() %>% `row.names<-` (.$group.subgroup) %>%
    .[,-1] %>% apply(MARGIN=1, FUN=ra) %>% t() %>% as.data.frame()
  
  ls.length <- function(x){return(length(unlist(x,F)))}
  
  grp.matches <- function(x){return(length(which(summary(as.factor(unlist(x)),maxsum=1000)==length(x))))}
  
  new <- data.frame()
  matches <- data.frame()
  for (k in c(0.01,seq(0.05,1,0.05))){
    
    common <- core.id(comm[,select], meta, group=group, subgroup=subgroup, margin = k)

    y <- c()
    for (j in names(unlist(common, F))){ 
      if (is.na(unlist(common, F)[[j]])){
        y <- c(y, 0)
      } else {
        y <- c(y, sum(grouped[which(row.names(grouped)==j), unlist(common, F)[[j]]]))
      }
    }
    
    matches <- rbind(matches, data.frame(temp.group=names(unlist(lapply(common, FUN=grp.matches))),
                                         temp.matches=unlist(lapply(common, FUN=grp.matches)), margin=k))
    
    o <- data.frame(group.subgroup=names(unlist(common,F)), margin=k, rel.ab=y,
                    count=unlist(lapply(unlist(common,F), FUN=ls.length)))
    
    new <- rbind(new, o)
    
  }
  
  names(matches)[1] <- group
  
  new <- left_join(new, unique(meta[,c(group, subgroup, "group.subgroup")])) %>%
    left_join(., matches)
  
  y1axis <- subset(new, subset=margin==0.01) %>% mutate(., newcount=count/max(count)*100)
  y2axis <- subset(new, subset=margin==1) %>% mutate(., newperc=paste(round(rel.ab*100,1),"%"))
  y3axis <- subset(new, subset=margin==0.01) %>% mutate(., newmatch=temp.matches/max(count)*100)
  
  #if (file_test("-f", "coreplots/ASV_occurrence.jpeg")==T){
  #  filename <- paste0("ASV_occurrence_",length(list.files("coreplots", "ASV_occurrence"))+1,".jpeg")
  #}
  
  dims <- c(ceiling(sqrt(length(common))),floor(sqrt(length(common))))
  
  
  jpeg(paste0("coreplots/",filename,".jpeg"), width=1100*dims[2]+900, height=900*dims[1]+150, res=300)
  print(ggplot(new) + facet_wrap(vars(eval(parse(text=group))), dir = "v") + theme_bw() + 
    geom_point(aes(margin, scales::rescale(count, c(0,100)), color=eval(parse(text=subgroup)), shape=eval(parse(text=subgroup))), size=3)+
    geom_line(aes(margin, scales::rescale(count, c(0,100)), color=eval(parse(text=subgroup)), linetype="Number of ASVs"))+
    geom_point(aes(margin, rel.ab*100, color=eval(parse(text=subgroup)), shape=eval(parse(text=subgroup))), size=3)+
      geom_line(aes(margin, rel.ab*100, color=eval(parse(text=subgroup)), linetype="Relative Abundance"))+
      geom_line(aes(margin, scales::rescale(temp.matches, c(0,max(y3axis$newmatch))), linetype="Shared ASVs"), color="black")+
      geom_label(data=y1axis, inherit.aes = F, aes(x=margin, y=newcount, label=count, color=eval(parse(text=subgroup))), hjust=1, show.legend = F) +
      geom_label(data=y2axis, inherit.aes = F, aes(x=margin, y=rel.ab*100, label=newperc, color=eval(parse(text=subgroup))), hjust=0, show.legend = F) +
      geom_label(data=y3axis, inherit.aes = F, aes(x=margin, y=newmatch, label=temp.matches), color="black", hjust=1, show.legend = F) +
      theme(strip.text = element_text(face="bold", size=14), legend.text = element_text(size=12),
            axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank()) +
    scale_shape_manual(values=shapes) + scale_color_manual(values=colors) +
    scale_linetype_manual(values=c("Relative Abundance"="dashed","Number of ASVs"="solid", "Shared ASVs"="dotted")) +
    scale_y_continuous(breaks=seq(0, 100, 25), labels=round(seq(0,max(new$count),max(new$count)/4)),
      sec.axis = sec_axis(~ ., breaks=seq(0,100,25),labels = paste0(seq(0,100, 25), "%"), name = "Total relative abundance (%) of detected ASVs")) +
    labs(y="Number of ASVs detected", x="% Hosts in which an ASV is detected", linetype="", color="", shape="") +
    scale_x_continuous(breaks = seq(0,1, 0.1), labels = paste0(seq(0,100, 10), "%"), limits=c(0,1)) +
      coord_cartesian(xlim=c(-0.1,1.2)))
  dev.off()
  
  return(new)

}