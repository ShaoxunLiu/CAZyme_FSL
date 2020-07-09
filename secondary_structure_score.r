install.packages("rvest")
install.packages("xml2")
library(xml2)
library(rvest)
library(dplyr)


input <- "~\\ITresults"
ssscore <- function(input = file){
  #get listed IT result 
  ITlist <- list.files(input)
  ITlen <- length(ITlist)
  #construct aas list
  aasList <- matrix(ncol = 1, nrow = ITlen)
  #create raw aa fequency matrix, frequency = type * confidence
  ssf <- matrix(ncol = 4, nrow = 20, data = 0)
  colnames(ssf) <- c("aa", "Helix", "Strand", "Coil")
  ssf[,1] <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  #crawling aa frequency data
  for(a in 1:ITlen){
    #crawlling indentification string
    url <- paste("~\\ITresults\\", ITlist[a],"\\index.html", sep = "")
    web<-read_html(url)
    ssp <- web %>% html_nodes("tr:nth-child(3) tr:nth-child(2) td+ td") %>% html_text()
    aas <- web %>% html_nodes("td~ th+ td") %>% html_text()
    aas <- substr(aas, start = nchar(ssp) + 1, stop = nchar(aas))
    aasList[a,1] <- aas
    ssi <- web %>% html_nodes("div~ div:nth-child(4) > table tr:nth-child(3) td+ td") %>% html_text()
    ssi <- as.character(ssi[1])
    for(b in 1: nchar(aas)){
      sam_aa <- substr(aas, start = b, stop = b)
      sam_sp <- substr(ssp, start = b, stop = b)
      sam_si <- substr(ssi, start = b, stop = b)
      ssf_row <- which(grepl(sam_aa, ssf[,1]))
      if(sam_sp == "H"){
        ssf[ssf_row, 2] <- as.numeric(ssf[ssf_row, 2]) + as.numeric(sam_si)
      }else if(sam_sp == "S"){
        ssf[ssf_row, 3] <- as.numeric(ssf[ssf_row, 3]) + as.numeric(sam_si)
      }else if(sam_sp == "C"){
        ssf[ssf_row, 4] <- as.numeric(ssf[ssf_row, 4]) + as.numeric(sam_si)
      }
    }
  }
  
  ssf_score <- matrix(ncol = 4, nrow = 20, data = 0)
  colnames(ssf_score) <- c("aa", "Helix", "Strand", "Coil")
  ssf_score[,1] <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  for(a in 1:20){
    ave <- mean(c(as.numeric(ssf[a,2:4])))
    sd <- sd(c(as.numeric(ssf[a,2:4])))
    ssf_score[a,2] <- (as.numeric(ssf[a,2]) - ave) / sd
    ssf_score[a,3] <- (as.numeric(ssf[a,3]) - ave) / sd
    ssf_score[a,4] <- (as.numeric(ssf[a,4]) - ave) / sd
  }
  write.csv2(ssf_score, file = paste(input, "\\ssf_score.txt", sep =""))
  write.csv2(aasList, file = paste(input, "\\aasList", sep = ""))
  return(ssf_score)
}

##############################################################################################################################
##HSC score###
##input###
#ssf_score
score <- as.data.frame(ssf_score)
#threshold
thresh <- matrix(nrow = 2, ncol = 3, data = 0)
colnames(thresh) <- c("H", "S", "C")
rownames(thresh) <- c("HighT", "LowT")
thresh[,1] <- c(-0.1,-1.4)
thresh[,2] <- c(-0.1,-1.4)
thresh[,3] <- c(1.10, 1.0)
#consequtive strnd
cons <- 5
#tolerence/missmatch
tol <- 1
#aas
seq <- aasList[1,1]
genSS(seq,ssf_score,thresh,cons,tol)
#####################################3
##execute###
genSS <- function(seq, scoreMatrix = ssf_score, thresholdMatrix = thresh, conseque = cons, tolerence = tol){
  score <- as.data.frame(ssf_score)
  result <- matrix(ncol = 3, nrow = 2, data = 0)
  colnames(result) <- c("H","S", "C")
  rownames(result) <- c("longStrand", "average")
  for(a in 2:4){
    longCount <- 0
    lifeCount <- 0
    seqLen <- 0
    highSum <- 0
    for(b in 1:nchar(seq)){
      achar <- substr(seq, start = b, stop = b)
      arow <- score %>% filter(aa == achar)
      anum <- as.numeric(as.character(arow[1,a]))
      if(anum > thresh[1,a-1]){
        seqLen <- seqLen + 1
        lifeCount <- lifeCount + 1
        highSum <- highSum + 1
      }else if(anum > thresh[2,a-1]){
        if(lifeCount >= tol){
          lifeCount <- 0
          seqLen <- seqLen + 1
        }else if(seqLen >= cons){
          longCount <- longCount + 1
          seqLen <- 0
          lifeCount <- 0
        }else{
          lifeCount <- 0
          seqLen <- 0
        }
      }else{
        if(seqLen >= cons){
          longCount <- longCount + 1
          seqLen <- 0
          lifeCount <- 0        
        }else{
          lifeCount <- 0
          seqLen <- 0        
        }
      }
    }
    result[1,a-1] <- longCount
    result[2,a-1] <- highSum
  }
  return(result)
}



#########cross check to adjust parameter################
#####input####
#####aas.no.
a <- 1
pro <- ITlist[a]
check(input,pro)
file <- input
###########################################
check <- function(file, pro){
  url <- paste(file, "\\", pro,"\\index.html", sep = "")
  web<-read_html(url)
  ssp <- web %>% html_nodes("tr:nth-child(3) tr:nth-child(2) td+ td") %>% html_text()
  H <- 0
  S <- 0
  C <- 0
  for(a in 1:nchar(ssp)){
    if(substr(ssp, start = a, stop = a) == "H"){
      H <- H + 1
    }else if(substr(ssp, start = a, stop = a) == "S"){
      S <- S + 1
    }else{
      C <- C + 1
    }
  }
  result <- matrix(nrow = 1, ncol = 3, data = c(H,S,C))
  return(result)
}

#############gownloading list of SS##############
ssMat <- data.frame(matrix(ncol = 8, nrow = 59))
colnames(ssMat) <- c("aas","sss","Ha","Hl","Sa","Sl","Ca","Cl")
for(a in 1:59){
  url <- paste(file, "\\", ITlist[a],"\\index.html", sep = "")
  web<-read_html(url)
  ssp <- web %>% html_nodes("tr:nth-child(3) tr:nth-child(2) td+ td") %>% html_text()
  aas <- web %>% html_nodes("td~ th+ td") %>% html_text()
  ssMat[a,1] <- aas
  ssMat[a,2] <- ssp
  HH <- 0
  SS <- 0
  CC <- 0
  HL <- 0
  SL <- 0
  CL <- 0
  heapS <- "A"
  heapN <- 0
  for(b in 1:nchar(ssp)){
    aass <- substr(ssp,b,b)
    if(aass == "H"){
      HH <- HH + 1
    }else if(aass=="S"){
      SS <- SS + 1
    }else{
      CC <- CC + 1
    }
    if(aass != heapS && heapN >= 5){
      if(aass == "H"){
        HL <- HL + 1
      }else if(aass=="S"){
        SL <- SL + 1
      }else{
        CL <- CL + 1
      }
      heapS <- aass
      heapN <- 0
    }else if(aass == heapS){
      heapN <- heapN + 1
    }else{
      heapN <- 0
      heapS <- aass
    }
  }
  ssMat[a,3] <- HH
  ssMat[a,4] <- HL
  ssMat[a,5] <- SS
  ssMat[a,6] <- SL
  ssMat[a,7] <- CC
  ssMat[a,8] <- CL
  print(paste("aas", a,"complete"))
}

write.csv2(ssMat, file = "~\\ssMat.txt")





