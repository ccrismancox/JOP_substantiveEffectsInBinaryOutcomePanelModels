rm(list=ls())
library(doRNG)
library(doParallel)

library(survival)
library(brglm)
library(lme4)
library(plm)

library(Matrix)
library(spam)

library(ggplot2)
library(reshape2)
library(scales)
library(matrixStats)
library(data.table)
source("sparseFirth.R")


beck.logit <- function(alpha, cml, X, y){
  theta <- c(cml, alpha)
  XB <- drop(X %*% theta)
  return(
    -sum(ifelse(y==1,
                plogis(XB, log.p=T),
                plogis(XB, lower.tail = F, log.p=T))))
}

cl <- makeCluster(30)
registerDoParallel(cl)

set.seed(123456)
n <- c(20, 40, 60, 80, 100)
t <- c(3, 5, 10,25)
g <- c(3.25,2.25)
NT <- expand.grid(n,t,g)
B <- 1000
M <- 50 #cutpoint for numeric issue 

results <- list()
for(i in 1:nrow(NT)){
  N <- NT[i, 1]
  T <- NT[i, 2]
  gamma <- NT[i, 3]
  
  alpha <- rnorm(N, mean=-4)

  
  out <- foreach(b= 1:B, .combine = 'rbind',
                 .packages = c("brglm",  "survival", "lme4", "data.table", "Matrix", "spam", "plm"),
                 .errorhandling = "remove") %dorng%{
                   
                   indexes <- expand.grid(1:T,1:N)
                   
                   c <- rep(alpha, each=T) #this is z in the paper
                   r <- 0
                   data1 <- list(y=0)
                   all.zero <- T
                   while(abs(r)<.3 | sum(data1$y)==0 || all.zero){
                     X <- rnorm(N*T)+ c 
                     e <- rlogis(N*T) + c*gamma
                     r <- cor(X,c)

                     y.linear <- -2*X + e
                     data1 <- data.table(state=factor(indexes[,2]),
                                         year=indexes[,1],
                                         y.linear=y.linear,
                                         X=X)
                     
                     data1$y <- ifelse(data1$y.linear>0,1,0)
                     
                     data1[, sum.y:=sum(y), by=state]
                     data1$cgamma <- c*gamma
                     within.y <- mean(data1[, mean(y), by=state]$V1)
                     data2 <- data1[(sum.y>0 & sum.y<T)]
                     data1[,state2 := ifelse((sum.y==0), N+1, state)]
                     data1[,state2 := ifelse((sum.y==T), N+2, state2)]
                     data1[,state2:=factor(state2)]
                     data2$state <- droplevels(data2$state)
                     ca <- c[data1$sum.y>0 & data1$sum.y<T]
                     prop.all.zero <- 1-length(unique(data2$state))/length(unique(data1$state))
                     if(gamma == g[1]){
                       all.zero <- !((prop.all.zero <1) & (prop.all.zero > 0))
                     }else{
                       all.zero <- FALSE
                     }
                   }
                   #note: state is already a factor variable (see above)
                   
                   t2 <- system.time({g2 <- try(glm(y~X+(state)-1, family=binomial, data=data2, x=T,y=T))})[3]
                   t3 <- system.time({g3 <- try(sp.ffe(y~X+(state)-1,  data=data1))})[3]
                   t4 <- system.time({g4 <- try(clogit(y~X+strata(state), data=data1))})[3]
                   
                   Z <- tapply(data1$X, data1$state, mean) %x% rep(1, T)
                   t5 <- system.time({g5 <- try(glmer(y~X+(1|state)+Z, data=data1, family =binomial))})[3]
                   t6 <- system.time({g6 <- try(glmer(y~X+(1|state), data=data1, family =binomial))})[3]

                   t7 <- system.time({g7 <- try(lm(y~X+state+0,data=data1))})[3]
                   
                   
                   if(class(g4)[1]=="clogit"){
                     t.beck <- system.time({g.beck <- try(optim(par=rep(0, length(unique(data2$state))), 
                                                                fn=beck.logit,method="BFGS",
                                                                cml=g4$coefficients,
                                                                X=sparse.model.matrix(~X+state+0, data=data2),
                                                                y=data2$y))})[3]
                     if(class(g.beck)!="try-error" && g.beck$convergence==0){
                       g.beck <- list(coefficients=c(g4$coefficients, g.beck$par), conv=g.beck$convergence)
                     }else{
                       g.beck <- list(coefficients=rep(NA, 1+length(unique(data2$state))))
                       t.beck <- NA
                     }
                   }else{
                     g.beck <-  list(coefficients=rep(NA, 1+length(unique(data2$state))))
                     t.beck <- NA
                   }
                   t.beck <- t.beck+t4
                   
                   newData <- data1[,mean(X), by=state]; setnames(newData, "V1", "X"); newData$Z <- newData$X
                   newData$state2 <- unique(data1[,list(state, state2)])$state2
                   newData2 <- data2[,mean(X), by=state]; setnames(newData2, "V1", "X")
                   X.mem <- sparse.model.matrix(~X+(state)-1, data=newData)
                   X.mem2 <- sparse.model.matrix(~X+(state)-1, data=newData2)
                   
                   
                   
                   
                   if(class(g2)=="try-error"||!g2$conv||any(is.na(g2$coef))|| abs(g2$coef[1]) > M){
                     g2 <- list(coef=rep(NA,N+1))
                     mar2 <- NA; bmar2 <- NA; ame2 <- NA; t2 <- NA
                   }else{
                     prob.g2 <- g2$fitted
                     prob.g2.new <- g2$fitted #if MEM use predict(g2, newdata=newData2, type="response")
                     mar2 <- mean((prob.g2 -plogis(-2*data2$X + data2$cgamma))^2)
                     bmar2 <-(mean(prob.g2 -plogis(-2*data2$X + data2$cgamma)))
                     ame2 <-(prob.g2.new*(1-prob.g2.new)*(coef(g2)[1]))
                   }
                   if(class(g3)=="try-error"||!(g3$conv==0) || any(is.na(g3$coefficients))||any(abs(g3$coefficients) > M)){
                     g3 <- list(coef=rep(NA,N+1))
                     mar3 <- NA; bmar3 <- NA; ame3 <- NA; t3 <- NA
                   }else{
                     g3$coef.table <- NULL
                     prob.g3 <- g3$fitted
                     prob.g3.new <- g3$fitted #if MEM use plogis(drop(X.mem %*% g3$coef))
                     mar3 <- mean((prob.g3 -plogis(-2*data1$X + data1$cgamma))^2)
                     bmar3 <-(mean(prob.g3 -plogis(-2*data1$X + data1$cgamma)))
                     ame3 <- (prob.g3.new*(1-prob.g3.new)*(coef(g3)[1]))
                   }
                   if(class(g4)=="try-error"||!g4$iter<10|| any(is.na(g4$coef))||any(abs(g4$coef) > M)){
                     g4 <- list(coef=rep(NA,1)); t4 <- NA; t4 <- NA
                   }
                   if(class(g5)=="try-error"||g5@optinfo$conv$opt!=0||anyNA(coef(g5)$state)||abs(unique(coef(g5)$state$X)) > M){
                     beta.CRE <- NA;alpha.CRE <- rep(NA,N);mar5 <- NA;bmar5 <- NA; ame5 <- NA; t5 <- NA
                   }else{
                     beta.CRE <- g5@beta[2]
                     alpha.CRE <- rowSums(unique(cbind(1, Z)) * as.matrix(coef(g5)$state[,-2]))

                     if(any(abs(c(beta.CRE, alpha.CRE))> M)){
                       beta.CRE <- alpha.CRE <- mar5 <- NA;bmar5 <- NA; ame5 <- NA
                     }else{
                       mar5 <- mean((predict(g5, type="response") -plogis(-2*data1$X + c*gamma))^2)
                       bmar5 <-(mean(predict(g5, type="response") -plogis(-2*data1$X + c*gamma)))
                       ame5 <- ((predict(g5, type="response", newdata=data1)*(1-predict(g5, type="response", newdata=data1))*(coef(g5)$state$X[1])))
                     }
                   }
                   if(class(g6)=="try-error"||g6@optinfo$conv$opt!=0||anyNA(coef(g6)$state)||abs(unique(coef(g6)$state$X)) > M){
                     beta.RE <- NA;alpha.RE <- rep(NA,N);mar6 <- NA;bmar6 <- NA; ame6 <- NA; t6 <- NA
                   }else{
                     beta.RE <- g6@beta[2]
                     alpha.RE <- as.matrix(coef(g6)$state[,-2])
                     if(any(abs(c(beta.RE, alpha.RE))> M)){
                       beta.RE <- alpha.RE <- mar6 <- NA;bmar6 <- NA; ame6 <- NA
                     }else{
                       mar6 <- mean((predict(g6, type="response") -plogis(-2*data1$X + data1$cgamma))^2)
                       bmar6 <-(mean(predict(g6, type="response") -plogis(-2*data1$X + data1$cgamma)))
                       ame6 <- ((predict(g6, type="response", newdata=data1)*(1-predict(g6, type="response", newdata=data1))*(coef(g6)$state$X[1])))
                     }
                   }
                   
                   
                   if(class(g.beck)=="try-error"|| any(is.na(g.beck$coef))||any(abs(g.beck$coef) > M)){
                     g.beck <-  list(coef=rep(NA, length(unique(data2$state))+1))
                     mar.beck <- NA; bmar.beck <- NA; ame.beck <- NA; t.beck <- NA
                   }else{
                     X.beck <-  model.matrix(~X+state-1, data=data2)
                     p.beck <- plogis(X.beck %*% g.beck$coef)
                     mar.beck <- mean((p.beck -plogis(-2*data2$X + data2$cgamma))^2)
                     bmar.beck <-(mean(p.beck -plogis(-2*data2$X + data2$cgamma)))
                     p.beck.new <-  plogis(model.matrix(~X+state-1, data=data2) %*% g.beck$coefficients)
                     ame.beck <- (p.beck.new*(1-p.beck.new)*(g.beck$coef[1]))
                   }
                   
                   if(class(g7)=="try-error" || any(is.na(g7$coef))||any(abs(g7$coef) > M)){
                     mar7 <- NA; bmar7 <- NA; ame7 <- NA; t7 <- NA
                   }else{
                     mar7 <- mean((predict(g7) -plogis(-2*data1$X + c*gamma))^2)
                     bmar7 <-(mean(predict(g7) -plogis(-2*data1$X + c*gamma)))
                     ame7 <- coef(g7)[1]
                   }
                   
                   
                   
                   
                   keep <- as.numeric(as.character(unique(data1[!(sum.y==0 | sum.y==T)]$state)))
                   m3 <- mean((gamma*alpha-g3$coef[-1])^2)
                   m2 <- mean((gamma*alpha[keep]-g2$coef[-1])^2)
                   m5 <- mean((gamma*alpha-alpha.CRE)^2)
                   m6 <- mean((gamma*alpha-alpha.RE)^2)
                   mbeck <- mean((gamma*alpha[keep]-g.beck$coef[-1])^2)

                   
                   b3 <-(mean(gamma*alpha-g3$coef[-1]))
                   b2 <-(mean(gamma*alpha[keep]-g2$coef[-1]))
                   b5 <-(mean(gamma*alpha-alpha.CRE))
                   b6 <-(mean(gamma*alpha-alpha.RE))
                   bbeck <-(mean(gamma*alpha[keep]-g.beck$coef[-1]))
                   
                   ame.truth <- (plogis(-2*data1$X + data1$cgamma)*(1-plogis(-2*data1$X + data1$cgamma))*-2)
                   truth.ame <-  mean(ame.truth)
                   
                   ame2a <- sum(ame2)/nrow(data1)
                   ame2 <- mean(ame2)
                   ame3 <- mean(ame3)
                   ame5 <- mean(ame5)
                   ame6 <- mean(ame6)
                   ame.beck <- mean(ame.beck)
                   
                   
                   
                   output <- c(g4$coef[1],
                               g2$coef[1],
                               g3$coef[1],
                               beta.CRE,
                               beta.RE,
                               g.beck$coef[1],
                               m2, m3,  m5, m6, mbeck,
                               b2,  b3,  b5, b6, bbeck,
                               mar2,  mar3,  mar5, mar6, mar7, mar.beck,
                               bmar2,  bmar3,  bmar5, bmar6, bmar7, bmar.beck,
                               ame2, ame2a, ame3, ame5, ame6, ame7, ame.beck, truth.ame,
                               t4, t2, t3,t5, t6, t7, t.beck,
                               within.y, prop.all.zero,
                               r)
                   output
                 }
  colnames(out) <- c("CML", "MLDV",  "FFE", "CRE", "RE", "Beck",
                     "mse.MLDV","mse.FFE", "mse.CRE", "mse.RE", "mse.beck",
                     "bias.MLDV", "bias.FFE", "bias.CRE",  "bias.RE", "bias.beck",
                     "pred.MLDV", "pred.FFE","pred.CRE", "pred.RE", "pred.LPM", "pred.beck",
                     "bias.pred.MLDV", "bias.pred.FFE","bias.pred.CRE", "bias.pred.RE", "bias.pred.LPM", "bias.pred.beck",
                     "ame2", "ame2a", "ame3", "ame5", "ame6", "ame7", "ame.beck", "ame.truth",
                     "time.CML", "time.MLDV",  "time.FFE", "time.CRE", "time.RE", "time.LPM", "time.beck",
                     "Prop.1", "Prop.dropped",
                     "corXc")
  
  out <- as.data.frame(out)                   
  out$bias.mem.MLDV <- out$ame2 - out$ame.truth
  out$bias.mem.MLDV2 <- out$ame2a - out$ame.truth
  out$bias.mem.FFE <- out$ame3 - out$ame.truth
  out$bias.mem.CRE <- out$ame5 - out$ame.truth
  out$bias.mem.RE <- out$ame6 - out$ame.truth
  out$bias.mem.LPM <- out$ame7 - out$ame.truth
  out$bias.mem.beck <- out$ame.beck - out$ame.truth
  
  out$mse.mem.MLDV <- (out$ame2 - out$ame.truth)^2
  out$mse.mem.MLDV2 <- (out$ame2a - out$ame.truth)^2
  out$mse.mem.FFE <-  (out$ame3 - out$ame.truth)^2 
  out$mse.mem.CRE <- (out$ame5 - out$ame.truth)^2
  out$mse.mem.RE <- (out$ame6 - out$ame.truth)^2
  out$mse.mem.LPM <- (out$ame7 - out$ame.truth)^2
  out$mse.mem.beck <- (out$ame.beck - out$ame.truth)^2
  out <- as.matrix(out)
  
  results[[i]] <- out

}
stopCluster(cl)


save.image("MC_main.rdata")



rm(list=ls())
load("MC_main.rdata")
library(ggplot2)
library(extrafont)
library(matrixStats)
library(reshape2) 
library(scales)
library(data.table)

adInfo <- cbind(NT, do.call(rbind,lapply(results, function(x){return(cbind(mean(x[,"Prop.1"]), mean(x[,"Prop.dropped"])))})))
adInfo <- data.table(adInfo)
adInfo[,Var1:=NULL]
colnames(adInfo) <- c("T", "Rare", "ybar", "dropped")

adInfo[,ybar:= round(mean(ybar),2), by=list(T, Rare)]
adInfo[,dropped:=round(mean(dropped),2)*100, by=list(T, Rare)]
adInfo$Rare <- ifelse(adInfo$Rare==g[1], "Rare event", "Non-rare Event")
adInfo$Rare <- factor(adInfo$Rare, levels=c("Rare event", "Non-rare Event"))
adInfo$T <- factor(paste("T =", adInfo$T), levels=paste("T =", t))
adInfo$label <- with(adInfo, paste("textstyle('%'~dropped) == ", dropped,  sep=""))

results <- lapply(results, as.matrix)

p1data <- cbind(NT, 
                do.call(rbind, 
                        lapply(results,
                               function(x){
                                 sqrt(colVars(x[,c(
                                   "CML", "FFE", "CRE", "MLDV", "RE", "Beck")], na.rm=T) +
                                     (colMeans(x[,c(
                                       "CML", "FFE", "CRE", "MLDV", "RE", "Beck")], na.rm=T)+ 2)^2)})))
colnames(p1data)[colnames(p1data)=="FFE"] <- c("PML")
p1data <- melt(p1data, id.vars = c("Var1", "Var2", "Var3"))
colnames(p1data) <- c("N", "T", "Rare", "Estimator", "RMSE")
p1data$Rare <- ifelse(p1data$Rare==g[1], "Rare event", "Non-rare Event")
p1data$Rare <- factor(p1data$Rare, levels=c("Rare event", "Non-rare Event"))
p1data$T <- factor(paste("T =", p1data$T), levels=paste("T =", t))


p2data <- cbind(NT, do.call(rbind, lapply(results,
                                          function(x){ 
                                            abs(colMeans(x[,c( "CML", "FFE", "CRE", "MLDV", "RE", "Beck")]+2, na.rm=T))})))
colnames(p2data)[colnames(p2data)=="FFE"] <- c("PML")
p2data <- melt(p2data, id.vars = c("Var1", "Var2", "Var3"))
colnames(p2data) <- c("N", "T", "Rare", "Estimator", "Bias")
p2data$Rare <- ifelse(p2data$Rare==g[1], "Rare event", "Non-rare Event")
p2data$Rare <- factor(p2data$Rare, levels=c("Rare event", "Non-rare Event"))
p2data$T <- factor(paste("T =", p2data$T),  levels=paste("T =", t))


p3data <- cbind(NT, do.call(rbind, 
                            lapply(results,
                                   function(x){ sqrt(colMeans(x[,c(
                                     "mse.FFE", "mse.CRE", "mse.MLDV", "mse.RE", "mse.beck"
                                   )], na.rm=T))})))
colnames(p3data)[4:8] <- c("PML", "CRE", "MLDV", "RE", "Beck")
p3data <- melt(p3data, id.vars = c("Var1", "Var2", "Var3"))
colnames(p3data) <- c("N", "T", "Rare", "Estimator", "RMSE")
p3data$RMSE <- sqrt(p3data$RMSE)
p3data$Rare <- ifelse(p3data$Rare==g[1], "Rare event", "Non-rare Event")
p3data$Rare <- factor(p3data$Rare, levels=c("Rare event", "Non-rare Event"))
p3data$T <- factor(paste("T =", p3data$T),  levels=paste("T =", t))


p4data <- cbind(NT, do.call(rbind,
                            lapply(results,
                                   function(x){ abs(colMeans(x[,c( #"bias.firth", "bias.cre",
                                     "bias.FFE", "bias.CRE", "bias.MLDV", "bias.RE", "bias.beck"
                                   )], na.rm=T))})))
colnames(p4data)[4:8] <- c("PML", "CRE", "MLDV", "RE", "Beck")
p4data <- melt(p4data, id.vars = c("Var1", "Var2", "Var3"))
colnames(p4data) <- c("N", "T", "Rare", "Estimator", "Bias")
p4data$Rare <- ifelse(p4data$Rare==g[1], "Rare event", "Non-rare Event")
p4data$Rare <- factor(p4data$Rare, levels=c("Rare event", "Non-rare Event"))
p4data$T <- factor(paste("T =", p4data$T), levels=paste("T =", t))


p5data <- cbind(NT, do.call(rbind,
                            lapply(results,
                                   function(x){ colMeans(x[,c( 
                                     "pred.FFE", "pred.CRE", "pred.MLDV", "pred.RE",  "pred.beck", "pred.LPM"
                                   )], na.rm=T)})))
colnames(p5data)[4:9] <- c("PML", "CRE", "MLDV", "RE",  "Beck", "LPM")
p5data <- melt(p5data, id.vars = c("Var1", "Var2", "Var3"))
colnames(p5data) <- c("N", "T", "Rare", "Estimator", "MSE")
p5data$Rare <- ifelse(p5data$Rare==g[1], "Rare event", "Non-rare Event")
p5data$Rare <- factor(p5data$Rare, levels=c("Rare event", "Non-rare Event"))
p5data$T <- factor(paste("T =", p5data$T), levels=paste("T =", t))



p6data <- cbind(NT, do.call(rbind,
                            lapply(results,
                                   function(x){ abs(colMeans(x[,c(
                                     "bias.pred.FFE", "bias.pred.CRE", "bias.pred.MLDV", 
                                     "bias.pred.RE", "bias.pred.beck","bias.pred.LPM"
                                   )], na.rm=T))})))
colnames(p6data)[4:9] <- c("PML", "CRE", "MLDV", "RE",  "Beck","LPM")
p6data <- melt(p6data, id.vars = c("Var1", "Var2", "Var3"))
colnames(p6data) <- c("N", "T", "Rare", "Estimator", "Bias")
p6data$Rare <- ifelse(p6data$Rare==g[1], "Rare event", "Non-rare Event")
p6data$Rare <- factor(p6data$Rare, levels=c("Rare event", "Non-rare Event"))
p6data$T <- factor(paste("T =", p6data$T), levels=paste("T =", t))

p7data <- cbind(NT, do.call(rbind,
                            lapply(results,
                                   function(x){ abs(colMeans(x[,c(
                                     "bias.mem.FFE", "bias.mem.CRE",
                                     "bias.mem.MLDV", "bias.mem.RE",
                                     "bias.mem.beck","bias.mem.LPM"
                                     )], na.rm=T))})))
colnames(p7data)[4:9] <- c("PML", "CRE", "MLDV", "RE", "Beck","LPM")
p7data <- melt(p7data, id.vars = c("Var1", "Var2", "Var3"))
colnames(p7data) <- c("N", "T", "Rare", "Estimator", "Bias")
p7data$Rare <- ifelse(p7data$Rare==g[1], "Rare event", "Non-rare Event")
p7data$Rare <- factor(p7data$Rare, levels=c("Rare event", "Non-rare Event"))
p7data$T <- factor(paste("T =", p7data$T), levels=paste("T =", t))




p8data <- cbind(NT, do.call(rbind,
                            lapply(results,
                                   function(x){colMeans( (x[,c(
                                     "mse.mem.FFE", "mse.mem.CRE", 
                                     "mse.mem.MLDV", "mse.mem.RE",
                                     "mse.mem.beck",  "mse.mem.LPM"
                                     )]), na.rm=T)})))
colnames(p8data)[4:9] <- c("PML", "CRE", "MLDV", "RE", "Beck", "LPM")
p8data <- melt(p8data, id.vars = c("Var1", "Var2", "Var3"))
colnames(p8data) <- c("N", "T", "Rare", "Estimator", "MSE")
p8data$Rare <- ifelse(p8data$Rare==g[1], "Rare event", "Non-rare Event")
p8data$Rare <- factor(p8data$Rare, levels=c("Rare event", "Non-rare Event"))
p8data$T <- factor(paste("T =", p8data$T), levels=paste("T =", t))




p9data <- cbind(NT, do.call(rbind,
                            lapply(results,
                                   function(x){colMeans( (x[,c( 
                                     "time.FFE", "time.CRE",
                                     "time.MLDV", "time.RE",
                                     "time.beck",  "time.LPM")]), na.rm=T)})))
colnames(p9data)[4:9] <- c("PML", "CRE", "MLDV", "RE", "Beck", "LPM")
p9data <- melt(p9data, id.vars = c("Var1", "Var2", "Var3"))
colnames(p9data) <- c("N", "T", "Rare", "Estimator", "Time")
p9data$Rare <- ifelse(p9data$Rare==g[1], "Rare event", "Non-rare Event")
p9data$Rare <- factor(p9data$Rare, levels=c("Rare event", "Non-rare Event"))
p9data$T <- factor(paste("T =", p9data$T), levels=paste("T =", t))



g1 <- ggplot(p1data[! p1data$Estimator %in% c("RE","Beck", "LPM", "PML"),])+ #RMSE in beta
  geom_line(aes(x=N, y=(RMSE), color=Estimator, linetype=Estimator), size=1)+
  facet_grid(Rare~T )+
  geom_text(aes(x=95,y=Inf, label=label),
            data=adInfo, parse=TRUE, vjust=1.1, size=5, hjust=1, check_overlap = TRUE) +
  theme_bw(16)+
  ylab('RMSE')+
  xlab(expression(N))+
  scale_linetype_manual(values=c("dotted", "solid", "dashed", "dotdash"))+
  scale_color_manual(values=hue_pal()(4)[c(2,1,3,4)])+
  theme(legend.position="bottom",
        legend.title = element_text(size=18,
                                    face="bold"),
        legend.text = element_text(size = 18),
        strip.text.y = element_text(size=14),
        legend.key.size = unit(.5,"in"))
print(g1)

g2 <- ggplot(p2data[! p2data$Estimator %in% c("RE","Beck", "LPM", "PML"),])+ #bias in beta
  geom_line(aes(x=N, y=(Bias), color=Estimator, linetype=Estimator), size=1)+
  facet_grid(Rare~T )+
  geom_text(aes(x=95,y=Inf, label=label),
            data=adInfo, parse=TRUE, vjust=1.1, size=5, hjust=.9, check_overlap = TRUE) +
  theme_bw(16)+
  ylab(' Absolute Bias')+
  xlab(expression(N))+
  scale_linetype_manual(values=c("dotted", "solid", "dashed", "dotdash"))+
  scale_color_manual(values=hue_pal()(4)[c(2,1,3,4)])+
  theme(legend.position="bottom",
        legend.title = element_text(size=18,
                                    face="bold"),
        legend.text = element_text(size = 18),
        strip.text.y = element_text(size=14),
        legend.key.size = unit(.5,"in"))
print(g2)

g3 <- ggplot(p3data[! p3data$Estimator %in% c("RE","Beck", "LPM", "PML"),])+ #rmse in c
  geom_line(aes(x=N, y=(RMSE), color=Estimator, linetype=Estimator), size=1)+
  facet_grid(Rare~T)+
  geom_text(aes(x=95,y=Inf, label=label),
            data=adInfo, parse=TRUE, vjust=1.1, size=5, hjust=.9, check_overlap = TRUE) +
  theme_bw(16)+
  ylab('RMSE')+
  xlab(expression(N))+
  scale_linetype_manual(values=c( "solid", "dashed", "dotdash"))+
  scale_color_manual(values=hue_pal()(4)[c(1,3,4)])+
  theme(legend.position="bottom",
        legend.title = element_text(size=18,
                                    face="bold"),
        legend.text = element_text(size = 18),
        strip.text.y = element_text(size=14),
        legend.key.size = unit(.5,"in"))
print(g3)

g5 <- ggplot(p5data[!p5data$Estimator %in% c("RE","Beck", "LPM", "PML"),])+ #rmse in pr
  geom_line(aes(x=N, y=sqrt(MSE), color=Estimator, linetype=Estimator), size=1)+
  facet_grid(Rare~T )+
  geom_text(aes(x=95,y=Inf, label=label),
            data=adInfo, parse=TRUE, vjust=1.1, size=5, hjust=.9, check_overlap = TRUE) +
  theme_bw(16)+
  ylab('RMSE')+ylim(0,.45)+
  xlab(expression(N))+
  scale_linetype_manual(values=c( "solid", "dashed", "dotdash"))+
  scale_color_manual(values=hue_pal()(4)[c(1,3,4)])+
  theme(legend.position="bottom",
        legend.title = element_text(size=18,
                                    face="bold"),
        legend.text = element_text(size = 18),
        strip.text.y = element_text(size=14),
        legend.key.size = unit(.5,"in"))
print(g5)



g6 <- ggplot(p6data[!p6data$Estimator %in% c("RE","Beck", "LPM", "PML"),])+ #bias in pr
  geom_line(aes(x=N, y=Bias, color=Estimator, linetype=Estimator), size=1)+
  facet_grid(Rare~T )+
  geom_text(aes(x=95,y=Inf, label=label),
            data=adInfo, parse=TRUE, vjust=1.1, size=5, hjust=.9, check_overlap = TRUE) +
  theme_bw(16)+
  ylab('Absolute Bias')+ylim(0,.3)+
  xlab(expression(N))+
  scale_linetype_manual(values=c( "solid", "dashed", "dotdash"))+
  scale_color_manual(values=hue_pal()(4)[c(1,3,4)])+
  theme(legend.position="bottom",
        legend.title = element_text(size=18,
                                    face="bold"),
        legend.text = element_text(size = 18),
        strip.text.y = element_text(size=14),
        legend.key.size = unit(.5,"in"))
print(g6)

g7 <- ggplot(p7data[! p7data$Estimator %in% c("RE","Beck", "LPM", "PML"),])+ #bias in marginal
  geom_line(aes(x=N, y=Bias, color=Estimator, linetype=Estimator), size=1)+
  facet_grid(Rare~T )+
  geom_text(aes(x=95,y=Inf, label=label),
            data=adInfo, parse=TRUE, vjust=1.1, size=5, hjust=.9, check_overlap = TRUE) +
  theme_bw(16)+
  ylab('Absolute Bias')+
  xlab(expression(N))+
  scale_linetype_manual(values=c( "solid", "dashed", "dotdash"))+
  scale_color_manual(values=hue_pal()(4)[c(1,3,4)])+
  theme(legend.position="bottom",
        legend.title = element_text(size=18,
                                    face="bold"),
        legend.text = element_text(size = 18),
        strip.text.y = element_text(size=14),
        legend.key.size = unit(.5,"in"))
print(g7)


g8 <- ggplot(p8data[! p8data$Estimator %in% c("RE","Beck", "LPM", "PML"),])+ #rmse in marginal
  geom_line(aes(x=N, y=sqrt(MSE), color=Estimator, linetype=Estimator), size=1)+
  facet_grid(Rare~T )+
  geom_text(aes(x=95,y=Inf, label=label),
            data=adInfo, parse=TRUE, vjust=1.1, size=5, hjust=.9, check_overlap = TRUE) +
  theme_bw(16)+
  ylab('RMSE')+ylim(0,.45)+
  xlab(expression(N))+
  scale_linetype_manual(values=c( "solid", "dashed", "dotdash"))+
  scale_color_manual(values=hue_pal()(4)[c(1,3,4)])+
  theme(legend.position="bottom",
        legend.title = element_text(size=18,
                                    face="bold"),
        legend.text = element_text(size = 18),
        strip.text.y = element_text(size=14),
        legend.key.size = unit(.5,"in"))
print(g8)



h=7; w=9
pdf(file="Figure1.pdf", height=h, width=w); g1; dev.off()
pdf(file="FigureA1.pdf", height=h, width=w); g2; dev.off()
pdf(file="Figure2.pdf", height=h, width=w); g3; dev.off()
pdf(file="Figure3.pdf", height=h, width=w); g5; dev.off()
pdf(file="FigureA2.pdf", height=h, width=w); g6; dev.off()
pdf(file="FigureA3.pdf", height=h, width=w); g7; dev.off()
pdf(file="Figure4.pdf", height=h, width=w); g8; dev.off()




#### APPENDIX B####

final.table <- rbind(c(p1data[p1data$N==100 & p1data$T== "T = 25" & p1data$Rare=="Non-rare Event",]$RMSE, NA),
                     c(NA,sqrt(p5data[p5data$N==100 & p5data$T== "T = 25" & p5data$Rare=="Non-rare Event",]$MSE)),
                     c(NA,sqrt(p8data[p8data$N==100 & p8data$T== "T = 25" & p8data$Rare=="Non-rare Event",]$MSE)))
colnames(final.table) <- c(as.character(p1data[p1data$N==100 & p1data$T== "T = 25" & p1data$Rare=="Non-rare Event",]$Estimator), "LPM")
rownames(final.table) <- c("RMSE in $\\hat\\beta$", 
                           "RMSE in $\\hat{p}$",
                           "RMSE in AME")

print(round((final.table/final.table[,"CRE"])[,-3], 2))






final.table.rare <- rbind(c(p1data[p1data$N==100 & p1data$T== "T = 25" & p1data$Rare=="Rare event",]$RMSE, NA),
                          c(NA,sqrt(p5data[p5data$N==100 & p5data$T== "T = 25" & p5data$Rare=="Rare event",]$MSE)),
                          c(NA,sqrt(p8data[p8data$N==100 & p8data$T== "T = 25" & p8data$Rare=="Rare event",]$MSE)))
colnames(final.table.rare) <-  c(as.character(p1data[p1data$N==100 & p1data$T== "T = 25" & p1data$Rare=="Non-rare Event",]$Estimator), "LPM")
rownames(final.table.rare) <- c("RMSE in $\\hat\\beta$", 
                                "RMSE in $\\hat{p}$",
                                "RMSE in AME")


print(round((final.table.rare/final.table.rare[,"CRE"])[,-3],2))
