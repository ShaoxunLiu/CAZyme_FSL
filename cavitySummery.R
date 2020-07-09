library(dplyr)
File <- "~\\MDVset"
DO <- list.files(File)

####construct matrix #####
cavFrags <- aasManue
fifteen <- matrix(ncol = 15, nrow = 59)
cavFrags <- cbind(cavFrags, fifteen)

####initialize loop####
posr <- 5
for(a in 1:59){
  if(a == 5|| a == 6){
    posr <- 6
  }else{
    posr <- 5
  }
  fname <- DO[a]
  cav <- read.table(paste("~/MDVset/", fname, "/cav.pdb", sep = ""), quote="\"", comment.char="", skip = 3, fill = TRUE)
  cav <- cav %>% filter(cav[,1] == "ATOM", cav[,3] == "N")
  ###generate one strand of cavFrag####
  fragList <- vector(length = 100, mode = "character")
  lInd <- 1
  pos <- cav[1,posr]
  frag <- aaTrans(AA321,as.character(cav[1,4]))
  for(b in 1:nrow(cav)){
    new <- cav[b,posr]
    if(new == pos){
      next
    }else if(new < (pos + 2)){
      pos <- new
      frag <- paste(frag, aaTrans(AA321,as.character(cav[b,4])), sep = "")
    }else{
      fragList[lInd] <- frag
      lInd <- lInd + 1
      pos <- cav[b, posr]
      frag <- aaTrans(AA321,as.character(cav[b,4]))
    }
  }
  fragList[lInd] <- frag
  fragLen <- vector(length = 100, mode = "numeric")
  for ( n in 1:100){
    fragLen[n] <- nchar(fragList[n]) / 100
  }
  order  <- as.data.frame(matrix(ncol = 3, nrow = 100, data = c(fragList, 1:100, fragLen)))
  longCav <- order%>%arrange(desc(order$V3))
  cavList <- as.character(longCav[c(1:15), 1])
  cavFrags[a,c(3:17)] <- cavList
}

write.csv2(cavFrags, file = "~/Frags.txt")


#######sugar binding matrix###############
sugarBinding <- matrix(ncol = 5, nrow = 22, data = 0)
sugarBinding[,1] <- as.character(aa321[,3])
colnames(sugarBinding) <- c("aa","1-kestose","Nystose","Raffinose","Stachyose")
for(a in 1 : 59){
  fname <- DO[a]
  seq <- aasManue%>%filter(aasManue$V2 == fname)
  scores <- vector(length = 4, mode = "numeric")
  seq2 <- aasSampleSum%>%filter(aasSampleSum$aas_seq == seq$V1)
  scores <- as.numeric(seq2[1:4])
  for(b in 1:4){
    cav <- read.table(paste("~/MDVset/", fname, "/out_", b, ".pdb", sep = ""), quote="\"", comment.char="", skip = 3, fill = TRUE)
    cav <- cav %>% filter(cav[,1] == "ATOM", cav[,3] == "N")
    for(c in 1:nrow(cav)){
      L <- aaTrans(AA321, as.character(cav[c,4]))
      R <- grep(L, sugarBinding[,1])
      sugarBinding[R,c(2:5)] <- as.numeric(sugarBinding[R,c(2:5)]) + as.numeric(scores)
    }
  }
}

sugarBinding <- as.data.frame(sugarBinding)
sugarBinding <- sugarBinding%>%filter(sugarBinding$Nystose != "0")
for(a in 2:5){
  ave <- mean(as.numeric(sugarBinding[,a]))
  std <- sd(as.numeric(sugarBinding[,a]))
  sugarBinding[,a] <- -(as.numeric(sugarBinding[,a]) - ave) / std
}
write.csv2(sugarBinding, file = "~/IcoSugarBinding.txt")


##########AA3213############
AA321 <- read.csv("~/aa321.txt", row.names=1, sep=";")
aaTrans(AA321, "ALA")
aaTrans <- function(AA321, aa){
  if(nchar(aa) == 1){
    return(as.character(AA321[grepl(aa,AA321$one_letter), 2]))
  } else if(nchar(aa) == 3){
    return(as.character(AA321[grepl(aa,AA321$three_letter), 3]))
  } else{
    stop("invalid aa input")
  }
}
