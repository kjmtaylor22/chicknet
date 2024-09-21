#' @export

biom.as.csv <- function(path){
  if (!"dplyr"%in%row.names(installed.packages())) {
    install.packages("dplyr")
  }
  if (!"biomformat"%in%row.names(installed.packages())){
    #source("https://bioconductor.org/biocLite.R")
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")}
    #biocLite(i, suppressUpdates=T)
    BiocManager::install("biomformat")
  }
  library(dplyr)
  library(biomformat)

  biom <- read_biom(path)

  comm <- biom$data

  den <- function(x){
    y <- x$id
    return(y)
  }

  taxon <- function(x){
    y <- paste0(unlist(x$metadata), collapse="; ")
    return(y)
  }

  tax1 <- lapply(biom$rows, FUN=den)
  tax2 <- lapply(biom$rows, FUN=taxon)

  tax <- data.frame(den.otu=unlist(tax1), taxonomy=unlist(tax2))

  comm <- as.data.frame(comm)

  attach(comm)
  comm <- comm[order(row.names(comm)),]
  detach(comm)

  comm <- t(comm) %>% as.data.frame()

  row.names(comm) <- tax1

  write.csv(comm, "./feature_table.csv")
  write.csv(tax, "taxonomy.csv", row.names=F)
}





bact.tax <- function(taxonomy, database=NULL){
  if (!"dplyr"%in%row.names(installed.packages())) {
    install.packages("dplyr")
  }
  library(dplyr)

  splt <- function(x){
    y <- strsplit(x, split=";", fixed=TRUE)
  }

  trim <- function(x){ #where x is a character string
    if ( is.null(database)==TRUE ){
      stop("Please specify database: SILVA (`silva`) or GreenGenes (`green`)")
    }
    if ( !(database%in%c("silva","green")) ){
      stop("Please specify database: SILVA (`silva`) or GreenGenes (`green`)")
    }
    if (database=="silva"){
      y <- substring(x, first=6)
    }
    if (database=="green"){
      y <- substring(x, first=4)
    }
    return(y)
  }

  silv <- function(x){
    library(stringr)

    for (w in c("^un", "^gut", "^meta", " bacterium", "taxa", "taxon", "[-_()]", "group", "sensu stricto")){
      w1 <- grep(w, x)
      if (w %in% c("^un", "^gut", "^meta", " bacterium", "taxa", "taxon")){
        x[w1] <- ""
      }
      if (w=="[-_()]" & length(w1)>0){
        w2 <- grep("^[A-Z]", x[7])
        if (length(w2)==1){
          if (w1==7){
            x[7] <- ""
          } else {
            x[w1] <- strsplit(x[7], split=" ")[[1]][1]
          }
        }
        if (length(w1)>=1 & length(w2)==0){
          a1 <- grep("Escherichia", x[6])
          a2 <- grep("Shigella", x[6])
          if (length(a1)==0 | length(a2)==0){
            x[w1] <- ""
          }
        }
      }

      if (w=="group" & length(w1)>0){
        for (i in w1){
          w2 <- strsplit(x[i], split=" ")
          w3 <- grep("[A-Z]", w2[[1]])
          if (length(w3) > 1){
            x[i] <- ""
          }
          if (length(w3)==1){
            x[i] <- w2[[1]][1]
          }
        }
      }
      if (w=="sensu stricto" & length(w1)>0){
        for (i in w1){
          x[i] <- strsplit(x[i], split=" ")[[1]][1]
        }
      }
    }

    if (length(x)==7 & x[7]!=""){
      w2 <- strsplit(x[7], split=" ")
      if (length(w2[[1]])==2){
        x[7] <- w2[[1]][2]
      }
      if (length(w2[[1]])>2){
        x[7] <- ""
      }
      if (x[7]=="sp."){
        x[7] <- ""
      }
    }

    if (length(x)==7 & x[6]==""){
      x[7] <- ""
    }


    s <- str_count(x, "[0-9]")
    if (length(s) > 0){
      for (i in 1:length(s)){
        if (s[i]==1){
          x[i] <- strsplit(x[i], split=" ")[[1]][1]
        } else {next}
      }
    }

    if (length(x) > 5){
      w1 <- str_count(x[6], "[A-Z]")
      if (w1==2 & str_count(x[6], "[ ]")==1){
        w3 <- strsplit(x[6], split=" ")[[1]]
        x[6] <- w3[1]
        x[7] <- w3[2] #tolower(w3[2])
      }
    }
    return(x)
  }

  add <- function(x){
    if (length(x)==7){
      y <- x
    } else {
      y <- append(x, rep("", 7-length(x)))
    }
  }

  org <- function(x){
    w <- grep("[1-9]", x)
    if (length(w)==0){
      y <- max(which(x!=""))
      if (y==7){
        z <- paste(x[6], x[7], sep=" ")
      } else {
        z <- paste(x[y])
      }
    } else {
      #z <- paste(x[min(w)-1], paste(x[w], collapse=" "), sep=" ")
      z <- x[min(w)-1]
    }
    pull <- grep("\\[", z)
    if (length(pull) > 0){
      ext <- unlist(strsplit(z, split=""))
      p1 <- grep("\\[", ext)
      p2 <- grep("\\]", ext)
      ext <- ext[-c(p1, p2)]
      z <- paste(ext,collapse="")
    }
    return(z)
  }

  tidy <- function(x){
    pull <- grep("\\[", x, value=TRUE)
    get <- grep("\\[", x)
    if (length(pull)==0){
      return(x)
    } else {
      for (i in length(pull)){
        ext <- unlist(strsplit(pull[i], split=""))
        p1 <- grep("\\[", ext)
        p2 <- grep("\\]", ext)
        ext <- ext[-c(p1, p2)]
        z <- paste(ext, collapse="")
        x[get] <- z
      }
      return(x)
    }
  }

  single <- function(x){
    y <- strsplit(as.character(x), split=" ")
    y1 <- lapply(y, FUN=length)
    y2 <- which(y1==1)
  }

  separated <- sapply(as.character(taxonomy[,2]), FUN=splt)
  trimmed <- lapply(separated, FUN=trim)

  if (database=="silva"){
    trimmed <- lapply(trimmed, FUN=silv)
  }

  evened <- sapply(trimmed, FUN=add)
  #return(evened)

  named <- apply(evened, MARGIN=2, FUN=org)
  names(named) <- NULL
  named <- unlist(named)

  out <- data.frame(den.otu=taxonomy[,1], taxonomy=named)
  out <- mutate(out, tag=paste("asv", row.names(out), sep=""))

  y <- rep(NA, dim(out)[1])
  y[single(out[,2])] <- paste(out[single(out[,2]),2], " (",
                              out[single(out[,2]),3], ")", sep="")
  y[-single(out[,2])] <- paste(out[-single(out[,2]),2])


  out <- data.frame(out, otu.name=y)

  z <- transmute(out, tag.name=paste0(otu.name, " (", tag, ")"))
  z$tag.name[grep("[(]", y)] <- y[grep("[(]", y)]

  out <- data.frame(out, tag.name=z$tag.name)


  evened <- apply(evened, MARGIN=2, FUN=tidy)
  classified <- as.data.frame(t(evened), row.names=1:dim(evened)[2])
  names(classified) <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")


  out <- data.frame(out, classified)

  out$phylum <- as.character(out$phylum)
  out$phylum[out$phylum=="Deinococcus"] <- "Deinococcus-Thermus"
  out$phylum <- as.factor(out$phylum)

  #if genus-species not available, use smallest level possible (but not strain codes;
  #use regular expressions to pull out any names given as alphanumeric codes)
  #and then follow that name with either "sp" or the strain code
  tax <- out
  save(tax, file="./taxonomy.RD")

  return(out)
}
