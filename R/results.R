install.packages("fields")
install.packages("e1071")
install.packages("lsa")
install.packages("dplyr")
install.packages("pracma")
install.packages("RSNNS")

library(fields)
library(lsa)
library(neuralnet)
library(dplyr)
library(pracma)
library(car)
library(e1071)
library(RSNNS)

file <- "~/Aug_10_feature_matrix.txt"
out <- c()
for(nn in 1:2){
  train_data <-  read.csv(file[1], row.names=1)
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
}  

  for(a in 2:4){
    func <- paste(func, "+L", a, sep = "")
  }
  func <- paste(func, "~X1", sep = "")
  for(a in 2:20){
    #if(imp[which(as.character(imp$name) == paste("X",a,sep = "")),1] > 0){
    func <- paste(func, "+X", a, sep = "")
    #}
  }
  f <- as.formula(func)
  
  
  total_res <- c(0,0,0,0)
  total_mat <- c()
  for(i in 1:50){
    rand <- sample(c(11:60), 10, replace = FALSE)
    mod <- vector(length = 100)
    for(a in 1:10){
      mod[((a-1)*10+1):((a-1)*10+10)] <- c(((rand[a]-1)*10+1):((rand[a]-1)*10+10)) + 60
    }
    mod <- c(mod,(rand))
    test <- total_data[mod,]
    train <- total_data[-mod,]
    
    bnn <- findBest(train, reps = 50, 110, f, hidden = c(100,100,100), 4, label)
    predict <- predict(bnn,test)
    b_res2 <- c(0,0,0,0)
    r_res <- c()
    for(a in -19:19){
        tryCatch({
        t <- a/20
        print(t)
        res <- getHit(predict, label,t, 10)
        res2 <- getHitsub2(predict, label,t, 10)
        if(res2[4] > b_res2[4]){
          b_res2 <- res2
          r_res <- res
        }else if(res2[4] == b_res2[4]){
          if(res2[3] > b_res2[3])
          b_res2 <- res2
          r_res <- res
        }
      },error = function(e) e)
    }
    total_res <- total_res + b_res2
    total_mat <- cbind(total_mat, r_res)
    print(paste("round", i, "complete"))
  }
  total_res <- total_res/50
  out <- c(out, total_res)


library(NeuralNetTools)
library(nnet)

list <- olden(bnn, out_var = "L1", bar_plot = TRUE)
olden1 <- list$data
list <- olden(bnn, out_var = "L2", bar_plot = TRUE)
list <- olden(bnn, out_var = "L3", bar_plot = TRUE)
list <- olden(bnn, out_var = "L4", bar_plot = TRUE)
olden2 <- list$data
olden3 <- list$data
olden4 <- list$data

###nnet
model <- neuralnet(f, data = total_data,  hidden = 20, act.fct = "logistic", linear.output = FALSE, rep = 1)
lek <- lekprofile(model)
str(lek)
lek$data[10000,]


install.packages("pROC")
library(pROC)
plot.roc(roc(res[2,], res[1,]))

######unaug#####
train_data <-  read.csv("~/Aug_10_feature_matrix.txt", row.names=1)
for(a in 1:60){
  train_data[c(((a+5)*10+1):((a+5)*10+10)),] <- train_data[a,]
}


####feature analysis####
core <- c()
corel <- c(ln1234, ln123,ln124, ln134, ln234)
for(i in 1:length(corel)){
  core <- c(core, which(grepl(corel[i], feature_names)))
}
core <- unique(core)



###start###
train_data <-  read.csv("~/Aug_10_feature_matrix.txt", row.names=1)
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
train_data[,-core] <- 0
total_data <- cbind(train_data, label)



func <- "L1"
for(a in 2:4){
  func <- paste(func, "+L", a, sep = "")
}
func <- paste(func, "~X1", sep = "")
for(a in 2:71){
  #if(imp[which(as.character(imp$name) == paste("X",a,sep = "")),1] > 0){
 #if(is.na(match(a,core)) == FALSE){
    func <- paste(func, "+X", a, sep = "")
 # }
  #}
}
f <- as.formula(func)


total_res <- c(0,0,0,0)
total_mat <- c()
for(i in 1:50){
  tryCatch({
    rand <- sample(c(11:60), 10, replace = FALSE)
    mod <- vector(length = 100)
    for(a in 1:10){
      mod[((a-1)*10+1):((a-1)*10+10)] <- c(((rand[a]-1)*10+1):((rand[a]-1)*10+10)) + 60
    }
    mod <- c(mod,(rand))
    test <- total_data[mod,]
    train <- total_data[-mod,]
    
    bnn <- findBest(train, reps = 50, 110, f, hidden = c(100,100,100), 4, label)
    predict <- predict(bnn,test)
    b_res2 <- c(0,0,0,0)
    r_res <- c()
    for(a in -19:19){
      tryCatch({
        t <- a/20
        print(t)
        res <- getHit(predict, label,t, 10)
        res2 <- getHitsub2(predict, label,t, 10)
        if(res2[4] > b_res2[4]){
          b_res2 <- res2
          r_res <- res
        }else if(res2[4] == b_res2[4]){
          if(res2[3] > b_res2[3])
            b_res2 <- res2
          r_res <- res
        }
      },error = function(e) e)
    }
    total_res <- total_res + b_res2
    total_mat <- cbind(total_mat, r_res)
    print(paste("round", i, "complete"))
  },error = function(e) e)
}
total_res <- total_res/50
total_res
out <- c(out, total_res)




same_class <- c("WSAL_02", "SBAR_04", "CFAL_04", "WSAL_04", "CFAL_07", "CFAL_08", "SSSS_01")
all_class <- c("SSSH_01", "SBAK_04", "CFAL_02", "WSAL_01", "SBAR_07", "SBAR_04", "SBAK_07", "SBAN_01", "SBAK_08")

write.csv(im, file = "~\\im.txt")
write.csv(feature_names, file = "~\\feature_names.txt")


im <-  read.csv("~\\im.txt", row.names=1, stringsAsFactors = FALSE)
feature_names <-  read.csv("~\\feature_names", row.names=1, stringsAsFactors = FALSE)
temp <- cbind.data.frame(feature_names, im)
rank <- arrange(temp, desc(mrmr))
rank1 <- rbind(rank[(1:10),],rank[(51:60),], rank[11:20,], rank[41:50,], rank[61:70,], rank[21:30,], rank[31:40,])
core <- c()
for(i in 1:20){
  core <- c(core, grep(as.character(rank1[i+50 ,1]), feature_names[,1]))
}


mvmv <- read.table("~/mvmv.txt",  quote="\"", comment.char="")
mrmr <- c(mvmv[(1:7),3], 0.75)
varimp <- c(mvmv[(8:14),3],0.75) 

plot(mrmr,type = "o", col = "red")

lines(varimp, type = "o", col = "blue")

rank2 <- arrange(mvmv, desc(V4))





