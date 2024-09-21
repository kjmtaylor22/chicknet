core.id <- function(comm, meta, tax, group, subgroup, margin){ 
  # where x is a sample-by-OTU table only, rows as samples, columns as otus
  # Row.names should be the sample IDs.
  # Any metadata should be input as 'factors', 
  # and make sure that the ID column in the metadata is labeled "SampleID"
  if (!"dplyr"%in%row.names(installed.packages())) {install.packages("dplyr")}
  library(dplyr)
  
  den <- colnames(comm) 
  
  new <- data.frame(sample=row.names(comm), comm)
  
  id <- grep("sample", names(meta), ignore.case = T)
  if (length(id)==1){
    names(meta)[id] <- "sample"
  } else {
    warning(paste0("Sample ID not found. Using factor column 1 (", names(meta)[1], ") as sample ID. "))
    names(meta)[1] <- "sample"
  }
  
  if (any(!row.names(comm)%in%meta[,which(names(meta)=="sample")])){
    warning("Check that all community samples are accounted for in metadata.")
    warning("May also need to rename sample ID column in metadata.")
    return(NULL)
  }
  
  data <- left_join(new, meta)
  
  rm(new)
  colnames(data) <- c("ID", names(data)[-1])
  
  core.mb <- function(x){
    z1 <- length(x[which(x>0)])
    if (length(margin)==1){
      if (z1/length(x) >= margin){
        return(paste("+"))
      }
    } else {
      if (z1/length(x) >= margin[i]){
        return(paste("+"))
      }
    }
  }
  
  
  out <- as.list(unique(as.character(data[,group])))
  names(out) <- unique(as.character(data[,group]))

  
  for (i in 1:length(out)){
    sub1 <- data[which(data[,group]==out[[i]]),]
    grp <- unique(as.character(sub1[,subgroup]))
    out1 <- as.list(unique(as.character(sub1[,subgroup])))
    
    for (j in 1:length(grp)){
      if (as.character(grp[j])%in%unique(as.character(sub1[,subgroup]))){
        cat(paste(out[[i]]," : ", grp[j], "\n", sep=""))
        sub2 <- sub1[which(sub1[,subgroup]==grp[j]),
                     which(names(sub1)%in%den)]
        z <- apply(sub2, MARGIN=2, FUN=core.mb) %>% unlist()
        out1[[j]] <- append(out1[[j]], names(z))
        out1[[j]] <- out1[[j]][-1]
        names(out1)[[j]] <- as.character(grp[j])
        if (length(out1[[j]])==0){
          out1[[j]] <- NA
        }
      }
    }
    
    out[[i]] <- append(out[[i]], out1)
    out[[i]] <- out[[i]][-1]
  }
  return(out)
}