###preparation of sample set
  ##input 60 AAS sample
  aasList <- read.csv("~/aasList.txt", row.names=1, sep=";")
  ##input label of AAS sample, generates AAS with normalized sugarbinding score
  aasList_labeled <- as.data.frame(matrix(ncol = 5, nrow = 60))
  aasList_labeled[,1] <- aasList[,1]  
  for(a in 1:60){
    aasList_labeled[a,c(2:5)] <- aasSampleSum[grep(substr(as.character(aasList_labeled[a,1]),1,100), as.character(aasSampleSum[,5])),c(1:4)]
  }  

###importing value matrixes
  ##AA exchange matrix
  aaEX <- read.csv("~/aaEX.txt", row.names=1, sep="")
  ##anchor AAS cavity fragments (cavity fragment sequence obtained from 6A search of molegro output)
  Frags <- read.csv("~/Frags.txt", row.names=1, sep=";")
  AAAS_cavity <- data.frame(matrix(nrow = 10, ncol = 16))
  AAAS_cavity[,1] <- aasList[c(1:10), 1]
  for(a in 1 : 10){
    for(b in 1:15){
      AA <- as.character(Frags[grep(as.character(substr(AAAS_cavity[a,1],1,100)), as.character(Frags[,1])), b+2])
      if(nchar(AA) > 1){
        AAAS_cavity[a,b+1] <- AA
      }else{
        AAAS_cavity[a,b+1] <- NA
      }
    }
  }
  ##AA sugar affinity 
    #see append_2 wsSBs and bsSbs was genreated
  
###feature generation part_1 
  ##secondary score generating functions are in append_5
  ##threshold matrix optimization is in append_6
  ##upper bound optimized at -0.02
  
###feature generation part_2
  ##cavity score code and sunction are in append_3
  ##threshold optimization code is in append_4
  ##t = 616 is the optimun threshold
  
###feature generation part_3
  ##whole sequence amino acid substitude search algorithm applied in append_7
  ##highest matching result obtained

###feature generation part_4
  ##appling sugar binding matrix of 4 oligosaccharides' filter on aaEX, optimizing filter ratio in append_8 

###feature generation part_5
  ##summerizing whole seuqence sugar binding average in append_9

###feature generation part_6
  nchar()
  
###feature generation pipeline for any number of substrates in append_10
  
###poisson noise generation for sample aas in append_11
  
###labeling of aas generation and colum binding with feature gerenated in append_12

###search of best Neuralnetwork in append_12
  
  
  
###append_1: read PDB
  ##########generate AA 321###########
  aa321 <- read.delim("~/AA321.txt", header=FALSE, comment.char="#")
  AA321 <- matrix(ncol = 3, nrow = 22)
  colnames(AA321) <- c("name", "three_letter", "one_letter")
  for(a in 1 : 22){
    AA321[a,1] <- as.character(aa321[(a-1)*3+1, 1])
    AA321[a,2] <- toupper(aa321[(a-1)*3+2, 1])
    AA321[a,3] <- as.character(aa321[(a-1)*3+3, 1])
  }
  AA321 <- as.data.frame(AA321)
  write.csv2(AA321, file = "~\AA_321.txt")
  
  ############aa321################
  aaTrans <- function(AA321, aa){
    if(nchar(aa) == 1){
      return(as.character(AA321[grepl(aa,AA321$one_letter), 2]))
    } else if(nchar(aa) == 3){
      return(as.character(AA321[grepl(aa,AA321$three_letter), 3]))
    } else{
      stop("invalid aa input")
    }
  }
  
  
  ###########import pdb file##############
  readPDB <- function(pdbFile, frag = FALSE){
    pdb <- read.csv(pdbFile, sep="", skip = 15, header = FALSE)
    if(!frag){
      point <- 1
      AAS <- character(0)
      for(a in 1:(nrow(pdb)-2)){
        if(pdb[a,6] == point){
          AAS <- paste(AAS, aaTrans(AA321,as.character(pdb[a,4])), sep = "")
          point <- point + 1
        }
      }
      return(AAS)
    }
  }
  
  
  #########get pdb from docking result############
  getPDB <- function(dockFile){
    file <- paste(dockFile, "\\Unnamed_complex.mvdml", sep ="")
    fstr <- readChar(file, 500)
    pos1 <- regexpr("filename=\\\"",fstr) + nchar("filename=\\\"") -1
    pos2 <- regexpr("\\\" >\\n<AtomList count",fstr) -1 
    pdbf <- substr(fstr, pos1, pos2)
    return(pdbf)
  }
  
  
  ##############main############################
  totalFile <- "~/MVD Data"
  totalList <- list.files(totalFile,pattern = "Docking")
  dockRes <- matrix(ncol = 12, nrow = length(MDVfiles))
  code <- c("[00]440_1", "[00]166_1", "[00]439_1", "[00]439_2")
  
  totalList <- list.files(totalFile,pattern = "Docking")
  aasSampleSum <- matrix(ncol = 5, nrow = 60)
  colnames(aasSampleSum) <- c("1-Kestose","Nystose", "Raffinose", "Stachyose", "aas_seq")
  for(a in 1:length(totalList)){
    usef <- paste(totalFile, "\\", totalList[a], sep = "")
    pf <- getPDB(usef)
    samaas <- readPDB(pf)
    if(length(grep(samaas,aasSampleSum[,5])) != 0){
      next
    }
    aasSampleSum[a,5] <- samaas
    for(b in 1:4){
      aasSampleSum[a,b] <- DockingResults[DockingResults$Name == code[b] ,"Energy"]
    }
  }

  
###append_2: AA_sugar binding matrix
  ##import binding site binding score
  bsSBs <- read.csv2("~/IcoSugarBinding.txt", row.names=1)
  ##whole sequence binding score
  wsSBs <- data.frame( matrix(ncol = 5, nrow = 20))
  colnames(wsSBs) <- colnames(bsSBs)
  wsSBs[,1] <- bsSBs[,1]
  wsSBs[c(1:20),c(2:5)] <- 0
    for(a in 1:60){
      aas <- as.character(aasList_labeled[a,1])
      for(b in 1:nchar(aas)){
        aa <- substr(aas,b,b)
        wsSBs[grep(aa,wsSBs[,1]),c(2:5)] <- as.numeric(wsSBs[grep(aa,wsSBs[,1]),c(2:5)]) + as.numeric(aasList_labeled[a,c(2:5)])
      }
    }
  for(a in 2:5){
    ave <- mean(wsSBs[,a])
    std <- sd(wsSBs[,a])
    wsSBs[,a] <- (wsSBs[,a] - ave) / std
  }
  
  
  
###append_3: cavity matching algorism and functions
  ####includes###
  setup <- function(){
    library(dplyr)
    library(compiler)
  }
  ########
  longest_subseq.R <- cmpfun(function(x) {
    P = integer(length(x))
    M = integer(length(x) + 1)
    L = newL = 0
    for (i in seq_along(x) - 1) {
      lo = 1
      hi = L
      while (lo <= hi) {
        mid = (lo + hi)%/%2
        if (x[M[mid + 1] + 1] < x[i + 1]) {
          lo = mid + 1
        } else {
          hi = mid - 1
        }
      }
      newL = lo
      P[i + 1] = M[newL]
      if (newL > L) {
        M[newL + 1] = i
        L = newL
      } else if (x[i + 1] < x[M[newL + 1] + 1]) {
        M[newL + 1] = i
      }
    }
    k = M[L + 1]
    re = integer(L)
    for (i in L:1) {
      re[i] = k + 1
      k = P[k + 1]
    }
    re
  })
  ####function 1###
  aaPatternScore <- function(aaEX, aast, patt){ ####aaEX: matrix of amino acid similarity; aast:string of aas; patt:string of pattern
    n.a <- nchar(aast)
    n.p <- nchar(patt)
    npos <- n.a-n.p+1
    score <- matrix(nrow = 1, ncol = npos, data = 0)
    for(a in 1:npos){
      reps <- 0
      for(b in 1:n.p){
        sub <- substr(aast, a+b-1, a+b-1)
        dir <- substr(patt, b,b)
        reps <- reps + as.numeric(aaEX[dir,sub])
      }
      score[1,a] <- reps / n.p
    }
    return(score)
  }
  
  ###function 2###
  aaCutScore <- function(score, thresh){###score:matrix generated from aaPatternScore; thresh: number from 0 to 1000;
    result <- c()
    for(a in 1:ncol(score)){
      if(score[1,a] > thresh){
        result <- append(result, a, after = length(result))
      }
    }
    return(result)###return: vector of pattern position in aas
  }
  
  ###function 3###
  cavFrag <- function(cavityFrag,aast,max = 15){###cavityFrag: table of integers of start,stop,start,stop.....position of cavity fragments
    ###aast: vector of aas sequences, subject number should be same as nrow of cavityFrag table
    ###max: maximun number of fragments of one aas
    frag <- matrix(nrow = 10, ncol = max, data = NA)
    for(a in 1: 10){
      for(b in 1:max){
        if(!is.na(cavityFrag[a,2*b - 1])){
          start <- cavityFrag[a,2*b - 1]
          stop <- cavityFrag[a,2*b]
          frag[a,b] <- substr(aast[a], start, stop)
        }else{
          break
        }
      }
    }
    return(frag)###return: matrix of string of each fragments of each aas
  }
  
  
  
  ###function 4###
  Cav.Score <- function(aastest, cavPatt, max = 15, thresh = 600){###aastest:string of aas to be tested;
    ###cavPatt: vector of cavity patterns, one row of frag
    ###max: maximun number of fragments of one aas
    ###thresh: threshold of cutting aa distence
    cavPos <- matrix(nrow = max, ncol = 10000, data = NA)
    for(a in 1:max){
      if(!is.na(cavPatt[a])){
        arr <- aaPatternScore(aaEX,aastest,cavPatt[a])
        att <- aaCutScore(arr, thresh)
        if(!is.null(att[b]) ){
          for(b in 1:length(att)){
            cavPos[a,b] <- att[b]
          }
        }
      }else{
        break
      }
    }
    
    ###########assign length weight to frag###########
    lenWeight <- matrix(nrow = 1, ncol = max)
    for(a in 1: max){
      lenWeight[1,a] <- nchar(cavPatt[a])
    }
    
    
    
    #############generate vector of matching cavity sequence###########
    aaslen <- nchar(aastest)
    cazMat <- vector(mode = "numeric", length = aaslen)
    for(a in 1:max){
      for(b in 1:ncol(cavPos)){
        if(!is.na(cavPos[a,b])){
          c <- cavPos[a,b]
          if(cazMat[c] == 0){
            cazMat[c] <- a
          }else if(cazMat[c] > a){
            c <- c + 1
            while(cazMat[c] != 0){
              c <- c + 1
            }
            cazMat[c] <- a
          }else{
            c <- c - 1
            if(c > 0){
              while(cazMat[c] != 0){
                c <- c - 1
                if(c == 0){
                  break
                }
              }
              cazMat[c] <- a
            }
          }
        }
      }
    }
    
    lisPos <- longest_subseq.R(cazMat)
    cavityScore <- 0
    for(a in 1:length(lisPos)){
      if(cazMat[lisPos[a]] != 0){
        cavityScore <- cavityScore + as.numeric(lenWeight[1,cazMat[lisPos[a]]])
      }
    }
    return(cavityScore)###return: raw cavity alignment score, can be normalized with self alignment score
  }
  
  ###function 5###
  seqMatch <- function(seq1, seq2, aaEX){###seq1:string aas to be compared with; seq2: string aas comparing to seq1; aaEX: amino acid distance matrix
    l1 <- nchar(seq1)
    l2 <- nchar(seq2)
    score <- vector("numeric", length = l1 + l2 - 1)
    comp <- matrix(nrow = 2, ncol = (l1 + 2*l2 - 2))
    ######fill in seq1#############
    for(a in l2:(l2+l1-1)){
      comp[1,a] <- substr(seq1, a+1-l2, a+1-l2)
    }
    ######train#####
    max <- 0
    for(a in 1:(l2+l1-1)){
      comp[2,] <- NA
      for(b in a:(a+l2-1)){
        comp[2,b] <- substr(seq2, b+1-a, b+1-a)
      }
      pair <- 0
      ex <- 0
      for(t in 1:(l1 + 2*l2 - 2)){
        if(!is.na(comp[1,t]) && !is.na(comp[2,t])){
          pair <- pair +1
          ex <- ex + aaEX[comp[1,t],comp[2,t]]
        }
      }
      iscore <- ex/l2
      if(iscore > max){
        max <- iscore
      }
    }
    return(max)###return number of maximun average similarity
  }
  
  
###append_4: find optimun threshold  
  cavS <- AAAS_cavity[,c(2:16)]
  threshOpt <- data.frame(matrix(ncol = 4, nrow = 300))
  colnames(threshOpt) <- c("no.", "Scol_std", "rowS_std", "score")
  threshOpt[,1] <- c(501:800)  
  for(t in 501:800){
    normal2 <- vector(length = 10)#generate normalizer
    for(a in 1: 10){
      normal2[a] <- Cav.Score(as.character(aasList[a,1]), cavS[a,], 15, t)
    }
    CSmat <- data.frame(matrix(ncol = 10, nrow = 60, data = 0))#generate cavscore matrix
    for(c in 1 :60){
        useaas <- as.character(aasList_labeled[c,1])
        for(b in 1:10){
          CSmat[c,b] <- (Cav.Score(useaas, cavS[b,], 15, t)) / normal2[b] 
        }
    }
    col_std <- vector(length = 10)
    for(d in 1:10){
      col_std[d] <- sd(CSmat[,d])
    }
    threshOpt[t-500,2] <- sum(col_std)
    threshOpt[t-500,3] <- sd(col_std)
    threshOpt[t-500,4] <- sum(col_std) / (sd(col_std)+1)
    print(paste("Threshold ", t, " has finished", sep = ""))
  }  

  
###append_5: functions for secondary structure score generation
  ##function 1: generate amino acid secondary structure supporting score
  input <- "~\ITresults" #I-Tasser result file
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
      url <- paste(input, ITlist[a],"\\index.html", sep = "")
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
      tt <- as.numeric(ssf[a,2]) + as.numeric(ssf[a,3]) + as.numeric(ssf[a,4])
      ssf[a,2] <- ssf[a,2] / tt
      ssf[a,3] <- ssf[a,3] / tt
      ssf[a,4] <- ssf[a,4] / tt
    }
    for(a in 2:4){
      ave <- mean(c(as.numeric(ssf[1:20,a])))
      sd <- sd(c(as.numeric(ssf[1:20,a])))
      ssf_score[,a] <- (as.numeric(ssf[,a]) - ave) / sd
    }
    write.csv2(ssf_score, file = paste(input, "\\ssf_score.txt", sep =""))
    write.csv2(aasList, file = paste(input, "\\aasList", sep = ""))
    return(ssf_score)
  }
  
  ##function 2: generating secondary structure score matrix 
  genSS <- cmpfun(function(seq, ssf_score, thresh, cons, tol){
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
  })


###append_6: optimization of threshold matrix. threshold matrix consists of higher bound and lower bound. 
  ##importing I-Tasser result
  ssMat <- data.frame(matrix(ncol = 8, nrow = 59))
  colnames(ssMat) <- c("aas","sss","Ha","Hl","Sa","Sl","Ca","Cl")
  for(a in 1:59){
    url <- paste(input, "\\", ITlist[a],"\\index.html", sep = "")
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
  
  ##fixing lower bond and adjust higher bound
  thresh <- matrix(nrow = 2, ncol = 3, data = 0)
  colnames(thresh) <- c("H", "S", "C")
  rownames(thresh) <- c("HighT", "LowT")
  thresh[,1] <- c(-1.0,-1.5)
  thresh[,2] <- c(-1.0,-1.5)
  thresh[,3] <- c(-1.0,-1.5)
  ##initializing scoring matrix
  sshOpt <- data.frame(matrix(ncol = 3, nrow = 251))
  colnames(sshOpt) <- c("H","S","C")
  
  ##Matrix generation
  for(a in 1:251){
    ssa <- data.frame(matrix(ncol = 7, nrow = 59))
    colnames(ssa) <- c("aas","Ha","Hl","Sa","Sl","Ca","Cl")
    ssa[,1] <- ssMat[,1]
    for(b in 1 : 59){
      seq <- as.character(ssa[b,1])
      res <- genSS(seq, ssf_score, thresh, 5 , 1)
      ssa[b,2] <- res[2,1]
      ssa[b,3] <- res[1,1]
      ssa[b,4] <- res[2,2]
      ssa[b,5] <- res[1,2]
      ssa[b,6] <- res[2,3]
      ssa[b,7] <- res[1,3]
    }
    sshOpt[a,1] <- sum( abs(ssa[,2] - ssMat[,3]) / (ssMat[,3] + 1) + sum(ssMat[,3]) / sum(ssMat[,4]) * (abs(ssa[,3] - ssMat[,4]) / (ssMat[,4] + 1)))
    sshOpt[a,2] <- sum( abs(ssa[,4] - ssMat[,5]) / (ssMat[,5] + 1) + sum(ssMat[,5]) / sum(ssMat[,6]) * (abs(ssa[,5] - ssMat[,6]) / (ssMat[,6] + 1)))
    sshOpt[a,3] <- sum( abs(ssa[,6] - ssMat[,7]) / (ssMat[,7] + 1) + sum(ssMat[,7]) / sum(ssMat[,8]) * (abs(ssa[,7] - ssMat[,8]) / (ssMat[,8] + 1)))
    thresh[1,] <- thresh[1,] + 0.01
    print(paste("round",a,"complete"))
  }
  ##high bound optimized at -0.02
  
  ##fixing higher bond and adjust lower bound
  thresh <- matrix(nrow = 2, ncol = 3, data = 0)
  colnames(thresh) <- c("H", "S", "C")
  rownames(thresh) <- c("HighT", "LowT")
  thresh[,1] <- c(-0.02,-1.5)
  thresh[,2] <- c(-0.02,-1.5)
  thresh[,3] <- c(-0.02,-1.5)
  ##initializing scoring matrix
  sshOpt <- data.frame(matrix(ncol = 3, nrow = 251))
  colnames(sshOpt) <- c("H","S","C")
  
  ##Matrix generation
  for(a in 1:251){
    ssa <- data.frame(matrix(ncol = 7, nrow = 59))
    colnames(ssa) <- c("aas","Ha","Hl","Sa","Sl","Ca","Cl")
    ssa[,1] <- ssMat[,1]
    for(b in 1 : 59){
      seq <- as.character(ssa[b,1])
      res <- genSS(seq, ssf_score, thresh, 5 , 1)
      ssa[b,2] <- res[2,1]
      ssa[b,3] <- res[1,1]
      ssa[b,4] <- res[2,2]
      ssa[b,5] <- res[1,2]
      ssa[b,6] <- res[2,3]
      ssa[b,7] <- res[1,3]
    }
    sshOpt[a,1] <- sum( abs(ssa[,2] - ssMat[,3]) / (ssMat[,3] + 1) + sum(ssMat[,3]) / sum(ssMat[,4]) * (abs(ssa[,3] - ssMat[,4]) / (ssMat[,4] + 1)))
    sshOpt[a,2] <- sum( abs(ssa[,4] - ssMat[,5]) / (ssMat[,5] + 1) + sum(ssMat[,5]) / sum(ssMat[,6]) * (abs(ssa[,5] - ssMat[,6]) / (ssMat[,6] + 1)))
    sshOpt[a,3] <- sum( abs(ssa[,6] - ssMat[,7]) / (ssMat[,7] + 1) + sum(ssMat[,7]) / sum(ssMat[,8]) * (abs(ssa[,7] - ssMat[,8]) / (ssMat[,8] + 1)))
    thresh[2,] <- thresh[2,] + 0.01
    print(paste("round",a,"complete"))
  }
  ##lower bound optimized at -1.5
  
  
###append_7: whole sequence amino acid substitution search
  seqMatch <- function(seq1, seq2, aaEX){###seq1:string aas to be compared with; seq2: string aas comparing to seq1; aaEX: amino acid distance matrix
    l1 <- nchar(seq1)
    l2 <- nchar(seq2)
    score <- vector("numeric", length = l1 + l2 - 1)
    comp <- matrix(nrow = 2, ncol = (l1 + 2*l2 - 2))
    ######fill in seq1#############
    for(a in l2:(l2+l1-1)){
      comp[1,a] <- substr(seq1, a+1-l2, a+1-l2)
    }
    ######train#####
    max <- 0
    for(a in 1:(l2+l1-1)){
      comp[2,] <- NA
      for(b in a:(a+l2-1)){
        comp[2,b] <- substr(seq2, b+1-a, b+1-a)
      }
      pair <- 0
      ex <- 0
      for(t in 1:(l1 + 2*l2 - 2)){
        if(!is.na(comp[1,t]) && !is.na(comp[2,t])){
          pair <- pair +1
          ex <- ex + aaEX[comp[1,t],comp[2,t]]
        }
      }
      iscore <- ex/l2
      if(iscore > max){
        max <- iscore
      }
    }
    return(max)###return number of maximun average similarity
  }
  
  
###append_8: optimizing sugarbinding filter ratio
  ##function_1: generate filtered aaEX
  genSSEX <- function(aaEX, sugarBindingCavity, type, rat){
    aaEX_sub <- aaEX
    ave <- mean(sugarBindingCavity[,type])
    for(source in 1:20){
      for(dest in 1:20){
        from <- colnames(aaEX)[source]
        to <- rownames(aaEX)[dest]
        from.val <- sugarBindingCavity[grep(from,sugarBindingCavity[,1]),type]
        to.val <- sugarBindingCavity[grep(to,sugarBindingCavity[,1]),type]
        diff <- (to.val - from.val)
        ratio <- (rat + diff) / rat
        aaEX_sub[dest,source] <- aaEX[dest,source] * ratio
      }
    }
    return(aaEX_sub)
  }

  
  ##operational code
  aaEX_ori <- aaEX
  sbrOpt <- data.frame(matrix(ncol = 5, nrow = 150))
  colnames(sbrOpt) <- c("thresh","1-kestose","nystose","raffinose","stachyose")
  sbrOpt[,1] <- c(1:50)*10
  for(a in 1:50){
    rat <- (a *10)
    for(r in 1:4){
      aaEX <- genSSEX(aaEX,IcoSugarBinding, (1+r), rat)
      normal2 <- vector(length = 10)#generate normalizer
      for(d in 1: 10){
        normal2[d] <- Cav.Score(as.character(aasList[d,1]), cavS[d,], 15, 616)
      }
      CSmat <- data.frame(matrix(ncol = 10, nrow = 60, data = 0))#generate cavscore matrix
      for(c in 1 :60){
        useaas <- as.character(aasList_labeled[c,1])
        for(b in 1:10){
          CSmat[c,b] <- (Cav.Score(useaas, cavS[b,], 15, a)) / normal2[b] 
        }
      }
      col_std <- vector(length = 10)
      for(d in 1:10){
        col_std[d] <- sd(CSmat[,d])
      }
      sbrOpt[a,r+1] <- sum(col_std) / (sd(col_std)+1)
    }
    aaEX <- aaEX_ori
    print(paste("ratio ", rat, " has finished", sep = ""))
  }
  
##ratio optimized at 100
  
###append_9:generating whole seuqence sugar binding score
  wholeSB <- function(useaas, sugarBindingWhole, n){
    SBM <- matrix(nrow = 1, ncol = n, data = 0)
    for(a in 1:nchar(useaas)){
      taa <- substr(useaas,a,a)
      aa <- sugarBindingWhole[grep(taa,sugarBindingWhole[,1]),2:(n+1)]
      SBM <- SBM + aa
    }
    return(SBM/nchar(useaas))
  }
  
  

###append_10 feature generation pipeline
  ##param
  aug_aasList <- aasList
  SBM <- IcoSugarBinding
  SBW <- IcoSugarWhole
  thresh <- matrix(nrow = 2, ncol = 3, data = 0)
  colnames(thresh) <- c("H", "S", "C")
  rownames(thresh) <- c("HighT", "LowT")
  thresh[,1] <- c(-0.02,-1.5)
  thresh[,2] <- c(-0.02,-1.5)
  thresh[,3] <- c(-0.02,-1.5)
  p1_thresh <- thresh
  cavS <- cavS
  anchorList <- data.frame(aasList[1:10,])
  nclass <- 4
  aaEX_ori <- aaEX_ori
  ##execute
  train_data <- featurePipe(SampleSet, p1_thresh, cavS, anchorList, nclass, aaEX_ori, SBM, SBW)
  write.csv2(train_data, file = "~\control_data.txt")
  
  ##function
  featurePipe <- function(aug_aasList, p1_thresh, cavS, anchorList, nclass, aaEX_ori, SBM, SBW){
    naas <- nrow(aug_aasList)
    nach <- nrow(cavS)
    training <- data.frame(matrix(nrow = naas, ncol = (7+nclass+nach*(nclass+2)), data = 0))
    normal2 <- vector(length = nach)
    for(a in 1: nach){
      normal2[a] <- Cav.Score(as.character(anchorList[a,1]), cavS[a,], 15, 616)
    }
    for(a in 1:naas){
      tryCatch({
        useaas <- as.character(aug_aasList[a,1])
        ss_M <- genSS(useaas, ss_s, p1_thresh, 5, 1)
        training[a,1] <- ss_M[1,1]
        training[a,2] <- ss_M[2,1]
        training[a,3] <- ss_M[1,2]
        training[a,4] <- ss_M[2,2]
        training[a,5] <- ss_M[1,3]
        training[a,6] <- ss_M[2,3]
        for(b in 1:nach){
          aaEX <- aaEX_ori
          training[a,6+b] <- (Cav.Score(useaas, cavS[b,], 15, 616)) / normal2[b] * 1000
          training[a,nach+6+b] <- seqMatch(as.character(anchorList[b,1]), useaas, aaEX)
          for(c in 1:nclass){
            aaEX <- genSSEX(aaEX_ori,SBM, (1+r), 100)
            training[a,6+(c+1)*nach+b] <- (Cav.Score(useaas, cavS[b,], 15, 616)) / normal2[b] * 1000
          }
        }
        training[a,c((7+nach*(nclass+2)):(6+nach*(nclass+2)+nclass))] <- wholeSB(useaas,SBW,nclass)
        training[a,(7+nach*(nclass+2)+nclass)] <- nchar(useaas)
      }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
      print(paste("AAS ", a, " completed", sep = ""))
    }
    return(training)
  }
  

###append_11 poisson noise generation
  ##########generating AA distribution using cazyme samples############
  aaDist <- matrix(ncol = 2, nrow = 20)
  aaDist[,1] <- as.character(sugarBindingCavity$name) 
  aaDist[,2] <- 0
  
  for(a in 1 : nrow(as.data.frame(aasList))){
    samaa <- as.character(aasList[a,1])
    for(b in 1:nchar(samaa)){
      aaDist[which(grepl(substr(samaa, start = b, stop = b), aaDist[,1])), 2] <- 
        as.numeric(aaDist[which(grepl(substr(samaa, start = b, stop = b), aaDist[,1])), 2]) + 1
    }
  }
  
  total <- sum(as.numeric(aaDist[,2]))
  for(a in 1 : 20){
    aaDist[a,2] <- as.numeric(aaDist[a,2]) / total
  }
  
  ############genreating poisson noise###############
  pNoise <- 0.1
  aaP <- matrix(ncol = 2, nrow = 20)
  aaP[,1] <- as.character(sugarBindingCavity$name) 
  aaP[,2] <- 0
  
  aaP[,2] <- pNoise * dpois(1,as.numeric(aaDist[,2]))
  for(a in 2:20){
    aaP[a,2] <- sum(as.numeric(aaP[(a-1) : a, 2]))
  }
  
  ##########generating augmented data############
  aug_aasList <- matrix(nrow = 190, ncol = 1, data = character(0))
  for(a in 1 : nrow(aasList)){
    for(b in 1 : 10){
      aug_aasList[(a-1) * 10 + b, 1] <- as.character(aasList[a,1])
      for(c in 1 : nchar(aug_aasList[(a-1) * 10 + b, 1])){
        subl <- which(runif(1) < as.numeric(aaP[,2]))
        if(length(subl) != 0){
          sub <- min(subl)
          substr(aug_aasList[(a-1) * 10 + b, 1], c, c) <- as.character(aaP[sub,1])
        }
      }
    }
  }
  colnames(aug_aasList) <- "V1"
  SampleSet <- rbind(data.frame(aasList), data.frame(aug_aasList))
  str(SampleSet)
  ##SampleSet with 1:60 being original sample, 61:660 being poisson noised sample, 10 in a group
  
  
###append_12: combining training data and label, proceed to best nn training
  library(fields)
  library(lsa)
  library(neuralnet)
  library(dplyr)
  library(pracma)
  library(car)
  ##normalizing training data
    train_data <- read.csv2("~/train_data.txt", row.names=1)
    for(a in 1:71){
      ave <- mean(train_data[,a])
      std <- sd(train_data[,a])
      train_data[,a] <- (train_data[,a] - ave) / std
      train_data[,a] <- sigmoid(train_data[,a])
    }
  ##labeling train_data
    train_class <- read.csv("~/train_class.txt", row.names=1, sep=";")
    label <- data.frame(matrix(ncol = 4, nrow = 660))
    label[1:60,] <- train_class[,2:5]
    colnames(label) <- c("L1","L2","L3","L4")
    for(a in 1 : 10){
      label[(c(1:60)*10 + a + 50),] <- train_class[,2:5]
    }  
    total_data <- cbind(train_data, label)
  
  ##generate formular
    func <- "L1"
    for(a in 2:4){
      func <- paste(func, "+L", a, sep = "")
    }
    func <- paste(func, "~X1", sep = "")
    for(a in 2:71){
      #if(imp[which(as.character(imp$name) == paste("X",a,sep = "")),1] > 0){
      func <- paste(func, "+X", a, sep = "")
      #}
    }
    f <- as.formula(func)
  
  ##generate sample set
    rand <- sample(c(1:60), 10, replace = FALSE)
    mod <- vector(length = 100)
    for(a in 1:10){
      mod[((a-1)*10+1):((a-1)*10+10)] <- c(((rand[a]-1)*10+1):((rand[a]-1)*10+10)) + 60
    }
    mod <- c(mod,(rand))
    test <- total_data[mod,]
    train <- total_data[-mod,]
    
  ##unaugmented control
    rand <- sample(c(1:60), 10, replace = FALSE)
    con_data <- total_data[1:60,]
    test <- con_data[rand,]
    train <- con_data[-rand,]

  
  ##train bnn    
    bnn <- findBest(train, reps = 20, 55, f, hidden = c(100,100,100), 4)
    ##test bnn
    cl <- c(0.5,0.5,0.5,0.5)
    cl <- cutLine(bnn, train, 4)
    predict <- predict(bnn, train)
    test_nn(predict, train,cl)
    predict <- predict(bnn,test)
    test_nn(predict, test, cl)
    getHit(predict, label, 0.6, 10)
    
 
    library(fields)
    library(lsa)
    library(neuralnet)
    library(dplyr)
    library(pracma)
    library(car)
  ##function_1: J-loss calculation
  gen_J <- function(S, Q, fsl_nn, nclass){
    Sample <- neuralnet::compute(fsl_nn,S)
    Predict <- neuralnet::compute(fsl_nn,Q)
    ck <- matrix(nrow = 10, ncol = nclass, data = 0)
    fSk <- matrix(nrow = 1, ncol = nclass, data = 0)
    for(k in 1:nclass){
      Sk <- which(grepl(1, S[,60+k]))
      for(v in 1: nclass){
        fSk[1,v] <- sum(Sample$net.result[Sk,v])/length(Sk)
      }
      ck[k,] <- fSk[1,]
    }
    
    J <- 0
    for(k in 1:nclass){
      Qk <- which(grepl(1, Q[,60+k]))
      if(length(Qk) != 0){
        for(b in 1:length(Qk)){
          Jp <- dist(rbind(ck[k,],Predict$net.result[Qk[b],]))
          ln <- 0
          for(a in 1:nclass){
            if(a != k){
              ln <- ln + exp(-dist(rbind(ck[a,],Predict$net.result[Qk[b],])))
            }
          }
          Jp <- (Jp + log(ln)) / nclass / nrow(Predict$net.result)
          J <- J + Jp
        }
      }
    }
    return(as.numeric(J))
  }

  ##function_2: find best nn according to J-loss
  #str(train)
  #dataSet = train
  #reps = 1
  #leaveOut = 110
  #f = f1
  #hidden =  c(100,100,100)
  #nclass = 4
  findBest <- function(dataSet, reps = 1, leaveOut = 10, f, hidden = c(50,50,50,50), nclass, label){
    std.J <- 0 ##10 in fsl
    best_nn <- NULL
    for(r in 1:reps){
      #a <- c(11:60)
      #s <- sample(a, leaveOut, replace = FALSE)
      ##poisson augment
      lv <- leaveOut/11
      s <- sample(c(11:50), lv, replace = FALSE)
      mod <- vector(length = lv*10)
      for(a in 1:lv){
        mod[((a-1)*10+1):((a-1)*10+10)] <- c(((s[a]+4)*10+1):((s[a]+4)*10+10))
      }
      s <- c(mod,(s))
      test <- dataSet[s,]
      train <- dataSet[-s,]
      nn <- neuralnet(f, data = train,  hidden, act.fct = "logistic", linear.output = FALSE, rep = 1)
      #predict <- predict(nn,test)
      #Js <- test_nn(predict, test)
      #Js <- abs(gen_J(train,test,nn,nclass))
      #cl <- cutLine(nn, train, nclass
      predict <- predict(nn, test)
      Js <- getHitsub(predict, label, 0, 10)
      print(paste("J = ", Js, sep = ""))
      if(Js > std.J ){#<in fsl#
        best_nn <- nn
        std.J <- Js
      }
      print(paste("epoch ", r, " finished"))
    }
    print(paste("J = ", std.J, sep = ""))
    return(best_nn)
  }
  ##function_3: test bnn
  #testset = test ##changing parameter
  #cutline = cl
  #nclass = 4
  test_nn <- function(predict, testset, cutline){
    sc <- 0
    for(a in 1:nrow(testset)){
      for(b in 1:4){
        pre <- 0
        if(predict[a,b] > cutline[b]){
          pre <- 1
        }
        if(pre == testset[a,ncol(testset)-4-1+b]){
          sc <- sc+1
        }
      }
    }
    return(sc / nrow(testset) / 4)
  }
  
  ##function_4: summerize augmented data
  cutLine <- function(nn, traindata, nclass){
    predict <- data.frame(predict(nn, traindata))
    cn <- ncol(traindata)
    op <- cn - nclass
    out <- vector(length = nclass)
    for(a in 1:nclass){
      sum <- sum(traindata[,(op+a)])
      re <- predict%>%arrange(desc(predict[,a]))
      if(sum !=0 ){
        out[a] <- re[sum,a]
      }else{
        out[a] <- 1.0
      }
    }
    return(out)
  }
  
  library(MLmetrics)
  ##function_5: final Hit test
 # thresh <- -1
  #pois <- 10
  getHit <- function(predict, label, thresh, pois){
    for(n in 1 :4){
      predict[,n] <- (predict[,n] - mean(predict[,n]))/sd(predict[,n])  
    }
    ns <- nrow(predict)/(pois+1)
    nc <- ncol(predict)
    result <- data.frame(matrix(ncol = nc, nrow = ns))
    checkm <- data.frame(matrix(ncol = nc, nrow = ns))
    rownames(result) <- rownames(predict[(((ns)*pois+1):(ns*pois+ns)),])
    hit <- 0
    for(a in 1:ns){
      for(b in 1:nc){
        result[a,b] <- mean(predict[c(c((((a-1)*pois)+1):(a*pois)), (a+ns*pois)),b])
        check <- label[as.numeric(rownames(result[a,])),b]
        checkm[a,b] <- check
        if(result[a,b] > thresh){
          result[a,b] <- 1
        }else{
          result[a,b] <- 0
        }
        if(result[a,b] == check){
          hit <- hit + 1
        }
      }
    }
    checkf1 <- c(as.matrix(checkm))
    resultf1 <- c(as.matrix(result))
    precision <- Precision(resultf1, checkf1, positive=NULL)
    recall <- Recall(resultf1, checkf1, positive=NULL)
    F1 <- F1_Score(as.vector(checkf1), as.vector(resultf1), positive = NULL)
    print(paste("F1 score is ", F1, sep = ""))
    acc <- hit / ns / nc
    print(paste("accuracy score is ", acc, sep = ""))
    print(c(precision, recall, F1, acc))
    return(rbind(checkf1, resultf1))
  }
  
getHitsub <- function(predict, label, thresh, pois){
  for(n in 1 :4){
    predict[,n] <- (predict[,n] - mean(predict[,n]))/sd(predict[,n])  
  }
  ns <- nrow(predict)/(pois+1)
  nc <- ncol(predict)
  result <- data.frame(matrix(ncol = nc, nrow = ns))
  checkm <- data.frame(matrix(ncol = nc, nrow = ns))
  rownames(result) <- rownames(predict[(((ns)*pois+1):(ns*pois+ns)),])
  hit <- 0
  for(a in 1:ns){
    for(b in 1:nc){
      result[a,b] <- mean(predict[c(c((((a-1)*pois)+1):(a*pois)), (a+ns*pois)),b])
      check <- label[as.numeric(rownames(result[a,])),b]
      checkm[a,b] <- check
      if(result[a,b] > thresh){
        result[a,b] <- 1
      }else{
        result[a,b] <- 0
      }
      if(result[a,b] == check){
        hit <- hit + 1
      }
    }
  }
  checkf1 <- c(as.matrix(checkm))
  resultf1 <- c(as.matrix(result))
  precision <- Precision(resultf1, checkf1, positive=NULL)
  recall <- Recall(resultf1, checkf1, positive=NULL)
  F1 <- F1_Score(as.vector(checkf1), as.vector(resultf1), positive = NULL)
  acc <- hit / ns / nc
  return(acc)
}
  

getHitsub2 <- function(predict, label, thresh, pois){
  for(n in 1 :4){
    predict[,n] <- (predict[,n] - mean(predict[,n]))/sd(predict[,n])  
  }
  ns <- nrow(predict)/(pois+1)
  nc <- ncol(predict)
  result <- data.frame(matrix(ncol = nc, nrow = ns))
  checkm <- data.frame(matrix(ncol = nc, nrow = ns))
  rownames(result) <- rownames(predict[(((ns)*pois+1):(ns*pois+ns)),])
  hit <- 0
  for(a in 1:ns){
    for(b in 1:nc){
      result[a,b] <- mean(predict[c(c((((a-1)*pois)+1):(a*pois)), (a+ns*pois)),b])
      check <- label[as.numeric(rownames(result[a,])),b]
      checkm[a,b] <- check
      if(result[a,b] > thresh){
        result[a,b] <- 1
      }else{
        result[a,b] <- 0
      }
      if(result[a,b] == check){
        hit <- hit + 1
      }
    }
  }
  checkf1 <- c(as.matrix(checkm))
  resultf1 <- c(as.matrix(result))
  precision <- Precision(resultf1, checkf1, positive=NULL)
  recall <- Recall(resultf1, checkf1, positive=NULL)
  F1 <- F1_Score(as.vector(checkf1), as.vector(resultf1), positive = NULL)
  acc <- hit / ns / nc
  rrr <- c(precision, recall, F1, acc)
  return(rrr)
}
  