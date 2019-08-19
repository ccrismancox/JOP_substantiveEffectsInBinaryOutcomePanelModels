rm(list=ls())
library(doRNG)
library(doParallel)

library(survival)
library(plm)
library(lme4)

library(reshape2)
library(scales)
library(matrixStats)
library(data.table)


cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

set.seed(123456)
n <- 1000
t <- 30
g <- 5.2
NT <- expand.grid(n,t,g)
B <- 1000
M <- 50 #cutpoint for numeric issue 

results <- list()
for(i in 1:nrow(NT)){
  N <- NT[i, 1]
  T <- NT[i, 2]
  gamma <- NT[i, 3] #this is alpha in the paper
  
  alpha <- rnorm(N, mean=-4) #this is z in the paper
  
  
  out <- foreach(b= 1:B, .combine = 'rbind',
                 .packages = c("survival", "lme4", "data.table"),
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
                     data3 <- data1[state2==N+1]
                     data3$state <- droplevels(data3$state)
                     prop.all.zero <- 1-length(unique(data2$state))/length(unique(data1$state))
                     if(gamma == g[1]){
                       all.zero <- !((prop.all.zero <1))
                     }else{
                       all.zero <- FALSE
                     }
                   }
                   #note: state is already a factor variable (see above)
                   
                   data1[,Z := mean(X), by =state]
                   data2[,Z := mean(X), by =state]
                   
                   t2 <- system.time({g2 <- try(glm(y~X+(state)-1, family=binomial, data=data2, x=T,y=T))})[3]

                   Z <- tapply(data1$X, data1$state, mean) %x% rep(1, T)
                   t5 <- system.time({g5 <- try(glmer(y~X+(1|state)+Z, data=data1, family =binomial))})[3]
                   t6 <- system.time({g6 <- try(glmer(y~X+(1|state)+Z, data=data2, family =binomial))})[3]

                   
                   if(class(g2)=="try-error"||!g2$conv||any(is.na(g2$coef))|| abs(g2$coef[1]) > M){
                     g2 <- list(coef=rep(NA,N+1))
                     mar2 <- NA; bmar2 <- NA; ame2 <- NA; t2 <- NA
                   }else{
                     prob.g2 <- g2$fitted
                     prob.g2.new <- g2$fitted #if MEM use predict(g2, newdata=newData2)
                     mar2 <- mean((prob.g2 -plogis(-2*data2$X + data2$cgamma))^2)
                     bmar2 <-(mean(prob.g2 -plogis(-2*data2$X + data2$cgamma)))
                     ame2 <-(prob.g2.new*(1-prob.g2.new)*(coef(g2)[1]))
                   }

                   data2[,Z := mean(X), by =state]
                   data3[,Z := mean(X), by =state]
                   
                   if(class(g5)=="try-error"||g5@optinfo$conv$opt!=0||anyNA(coef(g5)$state)||abs(unique(coef(g5)$state$X)) > M){
                     beta.CRE <- NA;alpha.CRE <- rep(NA,N);mar5 <- NA;bmar5 <- NA; ame5 <- NA; t5 <- NA
                     alpha.CRE.cml <- NA;mar5.cml <- NA;bmar5.cml <- NA; ame5.cml <- NA; ALPHA.cre <- NA
                   }else{
                     beta.CRE <- g5@beta[2]
                     alpha.CRE <- rowSums(unique(cbind(1, data1$Z)) * as.matrix(coef(g5)$state[,-2]))
                     alpha.CRE.cml <- alpha.CRE[unique(as.numeric(as.character(data2$state)))]
                     ALPHA.cre <- alpha.CRE[unique(as.numeric(as.character(data3$state)))]
                     
                     if(any(abs(c(beta.CRE, alpha.CRE))> M)){
                       beta.CRE <- NA;alpha.CRE <- rep(NA,N);mar5 <- NA;bmar5 <- NA; ame5 <- NA; t5 <- NA
                       alpha.CRE.cml <- NA;mar5.cml <- NA;bmar5.cml <- NA; ame5.cml <- NA;ALPHA.cre <- NA
                     }else{
                       mar5.cml <- mean((predict(g5, type="response", newdata=data2) -plogis(-2*data2$X + data2$cgamma))^2)
                       bmar5.cml <-(mean(predict(g5, type="response", newdata=data2) -plogis(-2*data2$X + data2$cgamma)))
                       ame5.cml <- ((predict(g5, type="response", newdata=data2)*(1-predict(g5, type="response", newdata=data2))*(coef(g5)$state$X[1])))
                       mar5 <- mean((predict(g5, type="response", newdata=data3) -plogis(-2*data3$X + data3$cgamma))^2)
                       bmar5 <-(mean(predict(g5, type="response", newdata=data3) -plogis(-2*data3$X + data3$cgamma)))
                       ame5 <- ((predict(g5, type="response", newdata=data3)*(1-predict(g5, type="response", newdata=data3))*(coef(g5)$state$X[1])))
                     }
                   }
                   
                   if(class(g6)=="try-error"||g6@optinfo$conv$opt!=0||anyNA(coef(g6)$state)||abs(unique(coef(g6)$state$X)) > M){
                     beta.CRE6 <- NA;mar6 <- NA;bmar6 <- NA; ame6 <- NA; t6 <- NA;
                     alpha.CRE.cml6 <- NA;mar6.cml <- NA;bmar6.cml <- NA; ame6.cml <- NA;
                   }else{
                     beta.CRE6 <- g6@beta[2]
                     # alpha.CRE6 <- rowSums(unique(data2$Z) * as.matrix(coef(g6)$state[,-2]))
                     alpha.CRE.cml6 <- rowSums(cbind(1,unique(data2$Z)) * as.matrix(coef(g6)$state[,-2]))
                     if(any(abs(c(beta.CRE6, alpha.CRE.cml6))> M)){
                       beta.CRE6 <- NA;mar6 <- NA;bmar6 <- NA; ame6 <- NA; t6 <- NA;
                       alpha.CRE.cml6 <- NA;mar6.cml <- NA;bmar6.cml <- NA; ame6.cml <- NA;
                     }else{
                       mar6.cml <- mean((predict(g6, type="response", newdata=data2) -plogis(-2*data2$X + data2$cgamma))^2)
                       bmar6.cml <-(mean(predict(g6, type="response", newdata=data2) -plogis(-2*data2$X + data2$cgamma)))
                       ame6.cml <- ((predict(g6, type="response", newdata=data2)*(1-predict(g6, type="response", newdata=data2))*(coef(g6)$state$X[1])))
                     }
                   }
                   
                  
                   m2 <- mean((g2$coef[-1]-unique(data2$cgamma))^2)
                   m5.cml <- mean((alpha.CRE.cml-unique(data2$cgamma))^2)
                   m6.cml <- mean((alpha.CRE.cml6-unique(data2$cgamma))^2)
                   m5 <- mean((ALPHA.cre-unique(data3$cgamma))^2)
                   
                   b2 <- mean((g2$coef[-1]-unique(data2$cgamma)))
                   b5.cml <- mean((alpha.CRE.cml-unique(data2$cgamma)))
                   b6.cml <- mean((alpha.CRE.cml6-unique(data2$cgamma)))
                   b5 <- mean((ALPHA.cre-unique(data3$cgamma)))
                   
                   ame.truth.cml <- mean((plogis(-2*data2$X + data2$cgamma)*(1-plogis(-2*data2$X + data2$cgamma))*-2))
                   ame.truth <- mean((plogis(-2*data3$X + data3$cgamma)*(1-plogis(-2*data3$X + data3$cgamma))*-2))
                   
                   ame5.cml <- mean(ame5.cml)
                   ame6.cml <- mean(ame6.cml)
                   
                   ame2 <- mean(ame2)
                   ame5 <- mean(ame5)
                   
                   
                  
                   within.y.90 <- quantile(plogis(-2*data1$X+data1$cgamma), prob=.9)
                   output <- c(m2, m6.cml,  m5.cml,  m5,
                               b2,  b6.cml,  b5.cml,  b5,
                               mar2,  mar6.cml,  mar5.cml,   mar5, 
                               bmar2,  bmar6.cml,  bmar5.cml,   bmar5, 
                               ame2, ame6.cml, ame5.cml,  ame5, ame.truth.cml, ame.truth,
                               within.y, prop.all.zero, within.y.90,
                               r)
                   output
                 }
  colnames(out) <- c("mse.MLDV","mse.rCRE.cml", "mse.CRE.cml", "mse.CRE",
                     "bias.MLDV", "bias.rCRE.cml", "bias.CRE.cml", "bias.CRE",
                     "pred.MLDV", "pred.rCRE.cml","pred.CRE.cml","pred.CRE", 
                     "bias.pred.MLDV", "bias.pred.rCRE.cml","bias.pred.CRE.cml","bias.pred.CRE",
                     "ame2", "ame6.cml", "ame5.cml", "ame5", "ame.truth.cml", "ame.truth",
                     "Prop.1", "Prop.dropped", "within.y.90",
                     "corXc")
  
  out <- as.data.frame(out)                   
  out$bias.mem.MLDV <- out$ame2 - out$ame.truth
  out$bias.mem.CRE <- out$ame5 - out$ame.truth
  
  
  out$bias.mem.rCRE.cml <- out$ame6.cml - out$ame.truth.cml
  out$bias.mem.CRE.cml <- out$ame5.cml - out$ame.truth.cml

  out$mse.mem.MLDV <- (out$ame2 - out$ame.truth)^2
  out$mse.mem.CRE <- (out$ame5 - out$ame.truth)^2
  
  out$mse.mem.rCRE.cml <- (out$ame6.cml - out$ame.truth.cml)^2 
  out$mse.mem.CRE.cml <- (out$ame5.cml - out$ame.truth.cml)^2

  out <- as.matrix(out)
  
  results[[i]] <- out

}
stopCluster(cl)


save.image("MC_highlyCensored.rdata")



rm(list=ls())
load("MC_highlyCensored.rdata")
library(matrixStats)
library(reshape2) 


results <- lapply(results, as.matrix)
loadfonts(quiet=T)

p3data <- cbind(NT, do.call(rbind,
                            lapply(results,
                                   function(x){ sqrt(colMeans(x[,c(
                                     "mse.rCRE.cml", "mse.CRE.cml", "mse.MLDV"
                                   )], na.rm=T))})))
colnames(p3data)[4:6] <- c("rCRE", "CRE", "MLDV")
p3data <- melt(p3data, id.vars = c("Var1", "Var2", "Var3"))
colnames(p3data) <- c("N", "T", "Rare", "Estimator", "RMSE")
p3data$RMSE <- sqrt(p3data$RMSE)
p3data$Rare <- ifelse(p3data$Rare==g[1], "Rare event", "Non-rare Event")
p3data$Rare <- factor(p3data$Rare, levels=c("Rare event", "Non-rare Event"))
p3data$T <- factor(paste("T =", p3data$T),  levels=paste("T =", t))


p5data <- cbind(NT, do.call(rbind,
                            lapply(results,
                                   function(x){ colMeans(x[,c( 
                                     "pred.rCRE.cml", "pred.CRE.cml", "pred.MLDV"
                                   )], na.rm=T)})))
colnames(p5data)[4:6] <- c("rCRE", "CRE", "MLDV")
p5data <- melt(p5data, id.vars = c("Var1", "Var2", "Var3"))
colnames(p5data) <- c("N", "T", "Rare", "Estimator", "MSE")
p5data$Rare <- ifelse(p5data$Rare==g[1], "Rare event", "Non-rare Event")
p5data$Rare <- factor(p5data$Rare, levels=c("Rare event", "Non-rare Event"))
p5data$T <- factor(paste("T =", p5data$T), levels=paste("T =", t))

p8data <- cbind(NT, do.call(rbind,
                            lapply(results,
                                   function(x){colMeans( (x[,c(
                                     "mse.mem.rCRE.cml", "mse.mem.CRE.cml",
                                     "mse.mem.MLDV")]), na.rm=T)})))
colnames(p8data)[4:6] <- c("rCRE", "CRE", "MLDV")
p8data <- melt(p8data, id.vars = c("Var1", "Var2", "Var3"))
colnames(p8data) <- c("N", "T", "Rare", "Estimator", "MSE")
p8data$Rare <- ifelse(p8data$Rare==g[1], "Rare event", "Non-rare Event")
p8data$Rare <- factor(p8data$Rare, levels=c("Rare event", "Non-rare Event"))
p8data$T <- factor(paste("T =", p8data$T), levels=paste("T =", t))




final.table.cml <- rbind(c(p3data$RMSE, NA),
                         sqrt(p5data$MSE),
                         sqrt(p8data$MSE))
colnames(final.table.cml) <-c("rCRE", "CRE", "MLDV")
rownames(final.table.cml) <- c("RMSE in $\\hat{c}$", 
                               "RMSE in $\\hat{p}$",
                               "RMSE in AME")
print(round(final.table.cml,2))