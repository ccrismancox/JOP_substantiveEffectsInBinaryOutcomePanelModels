rm(list=ls())
library(doRNG)
library(doParallel)

library(survival)
library(brglm)
library(lme4)

library(ggplot2)
library(reshape2)
library(scales)
library(matrixStats)
library(data.table)

cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

set.seed(123456)
n <- c(20, 40, 60, 80, 100)
t <- c(3, 5, 10,25)
g <- seq(3.25,2.25, length=5)
b <- c(-3, -2, -1)
NT <- expand.grid(n,t,g,b)
B <- 1000
M <- 50 #cutpoint for numeric issue 

results <- list()
for(i in 1:nrow(NT)){
  N <- NT[i, 1]
  T <- NT[i, 2]
  gamma <- NT[i, 3]
  beta <- NT[i, 4]
  
  alpha <- rnorm(N, mean=-4)

  
  out <- foreach(b= 1:B, .combine = 'rbind',
                 .packages = c( "survival", "lme4", "data.table"),
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

                     y.linear <- beta*X + e
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
                     all.zero <- !((prop.all.zero <1) )
                  
                     
                     print(prop.all.zero)
                   }
                   #note: state is already a factor variable (see above)
                   
                   t2 <- system.time({g2 <- try(glm(y~X+(state)-1, family=binomial, data=data2, x=T,y=T))})[3]
                   t4 <- system.time({g4 <- try(clogit(y~X+strata(state), data=data1))})[3]
                   
                   Z <- tapply(data1$X, data1$state, mean) %x% rep(1, T)
                   t5 <- system.time({g5 <- try(glmer(y~X+(1|state)+Z, data=data1, family =binomial))})[3]

                   
                   t7 <- system.time({g7 <- try(lm(y~X+state+0,data=data1))})[3]
                   
                   
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
                     mar2 <- mean((prob.g2 -plogis(beta*data2$X + data2$cgamma))^2)
                     bmar2 <-(mean(prob.g2 -plogis(beta*data2$X + data2$cgamma)))
                     ame2 <-(prob.g2.new*(1-prob.g2.new)*(coef(g2)[1]))
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
                       mar5 <- mean((predict(g5, type="response") -plogis(beta*data1$X + c*gamma))^2)
                       bmar5 <-(mean(predict(g5, type="response") -plogis(beta*data1$X + c*gamma)))
                       ame5 <- ((predict(g5, type="response", newdata=data1)*(1-predict(g5, type="response", newdata=data1))*(coef(g5)$state$X[1])))
                     }
                   }

                   keep <- as.numeric(as.character(unique(data1[!(sum.y==0 | sum.y==T)]$state)))
                   m2 <- mean((gamma*alpha[keep]-g2$coef[-1])^2)
                   m5 <- mean((gamma*alpha-alpha.CRE)^2)

                   

                   b2 <-(mean(gamma*alpha[keep]-g2$coef[-1]))
                   b5 <-(mean(gamma*alpha-alpha.CRE))

                   ame.truth <- (plogis(beta*data1$X + data1$cgamma)*(1-plogis(beta*data1$X + data1$cgamma))*beta)
                   truth.ame <-  mean(ame.truth)
                   
                   ame2 <- mean(ame2)
                   ame5 <- mean(ame5)

                   
                   
                   output <- c(g4$coef[1],
                               g2$coef[1],
                               beta.CRE,
                               m2,  m5, 
                               b2,  b5,
                               mar2,  mar5,
                               bmar2,   bmar5, 
                               ame2, ame5, truth.ame,
                               t4, t2, t5,
                               within.y, prop.all.zero,
                               r)
                   output
                 }
  colnames(out) <- c("CML", "MLDV",  "CRE", 
                     "mse.MLDV", "mse.CRE",
                     "bias.MLDV",  "bias.CRE",  
                     "pred.MLDV", "pred.CRE", 
                     "bias.pred.MLDV","bias.pred.CRE", 
                     "ame2",  "ame5",  "ame.truth",
                     "time.CML", "time.MLDV",  "time.CRE", 
                     "Prop.1", "Prop.dropped",
                     "corXc")
  
  out <- as.data.frame(out)                   
  out$bias.mem.MLDV <- out$ame2 - out$ame.truth
  out$bias.mem.CRE <- out$ame5 - out$ame.truth

  
  out$mse.mem.MLDV <- (out$ame2 - out$ame.truth)^2
  out$mse.mem.CRE <- (out$ame5 - out$ame.truth)^2
  out <- as.matrix(out)
  
  results[[i]] <- out

}
stopCluster(cl)


save.image("MC_expanded.rdata")



rm(list=ls())
load("MC_expanded.rdata")
library(ggplot2)
library(matrixStats)
library(reshape2) 
library(scales)
library(data.table)

results <- lapply(results, as.matrix)
adInfo <- cbind(NT, do.call(rbind,lapply(results, function(x){return(cbind(mean(x[,"Prop.1"]), mean(x[,"Prop.dropped"])))})))
adInfo <- data.table(adInfo)
adInfo[,Var1:=NULL]
colnames(adInfo) <- c("T", "alpha", "beta", "ybar", "dropped")

adInfo[,ybar:= round(mean(ybar),3), by=list(T, alpha, beta)]
adInfo[,dropped:=round(mean(dropped),2)*100, by=list(T, alpha, beta)]


p1data <- cbind(NT, 
                do.call(rbind, 
                        lapply(results,
                               function(x){
                                 sqrt(colVars(x[,c(
                                   "CML", "CRE", "MLDV")], na.rm=T) +
                                     (colMeans(x[,c(
                                       "CML", "CRE", "MLDV")], na.rm=T)+ 2)^2)})))
p1data <- melt(p1data, id.vars =  paste("Var", 1:ncol(NT), sep=""))
colnames(p1data) <- c("N", "T", "alpha", "beta",  "Estimator", "RMSE")



p5data <- cbind(NT, do.call(rbind,
                            lapply(results,
                                   function(x){ colMeans(x[,c( "pred.CRE", "pred.MLDV"
                                   )], na.rm=T)})))
colnames(p5data)[5:6] <- c( "CRE", "MLDV")
p5data <- melt(p5data, id.vars =paste("Var", 1:ncol(NT), sep=""))
colnames(p5data) <- c("N", "T", "alpha", "beta", "Estimator", "MSE")



p8data <- cbind(NT, do.call(rbind,
                            lapply(results,
                                   function(x){colMeans( (x[,c(
                                     "mse.mem.CRE",  "mse.mem.MLDV"
                                     )]), na.rm=T)})))
colnames(p8data)[5:6] <- c("CRE", "MLDV")
p8data <- melt(p8data, id.vars = paste("Var", 1:ncol(NT), sep=""))
colnames(p8data) <- c("N", "T", "alpha", "beta", "Estimator", "MSE")




adInfo$label <- with(adInfo, paste("textstyle('%'~dropped) == ", dropped,  sep=""))

h=9; w=6
pdf(file="FigureC1.pdf", height=h, width=w)
g1 <- ggplot(p1data)+ #RMSE in beta
  geom_line(aes(x=N, y=(RMSE), color=Estimator, linetype=Estimator), size=1)+
  facet_grid(beta+alpha~T, labeller = label_bquote(rows= atop(beta==.(beta),alpha==.(alpha)), cols=T==.(T)) )+
  geom_text(aes(x=95,y=Inf, label=label),
            data=adInfo, parse=TRUE, vjust=1.1, size=2.5, hjust=1, check_overlap = TRUE) +
  theme_bw(9)+
  ylab('RMSE')+
  xlab(expression(N))+
  scale_linetype_manual(values=c("dotted", "solid", "dashed", "dotdash"))+
  scale_color_manual(values=hue_pal()(4)[c(2,1,3,4)])+
  theme(legend.position="bottom",
        legend.title = element_text(size=10,
                                    face="bold"),
        legend.text = element_text(size = 10),
        strip.text.y = element_text(size=8),
        legend.key.size = unit(.25,"in"))
g1
dev.off()
pdf(file="FigureC2.pdf", height=h, width=w)
g5 <- ggplot(p5data)+ #rmse in pr
  geom_line(aes(x=N, y=sqrt(MSE), color=Estimator, linetype=Estimator), size=1)+
  facet_grid(beta+alpha~T, labeller = label_bquote(rows= atop(beta==.(beta),alpha==.(alpha)), cols=T==.(T)) )+
  geom_text(aes(x=95,y=Inf, label=label),
            data=adInfo, parse=TRUE, vjust=1.1, size=2.5, hjust=1, check_overlap = TRUE) +
  theme_bw(9)+
  ylab('RMSE')+
  xlab(expression(N))+
  scale_linetype_manual(values=c( "solid", "dashed", "dotdash"))+
  scale_color_manual(values=hue_pal()(4)[c(1,3,4)])+
  theme(legend.position="bottom",
        legend.title = element_text(size=10,
                                    face="bold"),
        legend.text = element_text(size = 10),
        strip.text.y = element_text(size=8),
        legend.key.size = unit(.25,"in"))
g5
dev.off()


pdf(file="FigureC3.pdf", height=h, width=w)
g8 <- ggplot(p8data)+ #rmse in marginal
  geom_line(aes(x=N, y=sqrt(MSE), color=Estimator, linetype=Estimator), size=1)+
  facet_grid(beta+alpha~T, labeller = label_bquote(rows= atop(beta==.(beta),alpha==.(alpha)), cols=T==.(T)) )+
  geom_text(aes(x=95,y=Inf, label=label),
            data=adInfo, parse=TRUE, vjust=1.1, size=2.5, hjust=1, check_overlap = TRUE) +
  theme_bw(9)+
  ylab('RMSE')+
  xlab(expression(N))+
  scale_linetype_manual(values=c( "solid", "dashed", "dotdash"))+
  scale_color_manual(values=hue_pal()(4)[c(1,3,4)])+
  theme(legend.position="bottom",
        legend.title = element_text(size=10,
                                    face="bold"),
        legend.text = element_text(size = 10),
        strip.text.y = element_text(size=8),
        legend.key.size = unit(.25,"in"))
g8
dev.off()




