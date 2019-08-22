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
t <- c(3, 5, 10, 25)
g <- c(3.25,2.25)
NT <- expand.grid(n,t,g)
B <- 1000
M <- 50

results <- list()
for(i in 1:nrow(NT)){
  N <- NT[i, 1]
  T <- NT[i, 2]
  gamma <- NT[i, 3]
  
  alpha <- rnorm(N, mean=-4)
  
  out <- foreach(b= 1:B, .combine = 'rbind',
                 .packages = c("survival", "lme4", "data.table"),
                 .errorhandling = "remove") %dorng%{
                   
                   indexes <- expand.grid(1:T,1:N)
                   
                   c <- rep(alpha, each=T) 
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
                     if(gamma == 3.25){
                       all.zero <- !((prop.all.zero <=1) & (prop.all.zero > 0))
                     }else{
                       all.zero <- FALSE
                     }
                   }
                   data1[,state2:=NULL]
                   z2 <- rep(rnorm(N*9, mean=-4), each=T)
                   #note: state is already a factor variable (see above)
                   data1 <- rbind(data1,  #add irrelevant units
                                  data.table(state=factor(rep((N+1):(N+9*N), each=T)),
                                             year=rep(1:T, N*9),
                                             y.linear=-Inf,
                                             X=rnorm((9*N)*T)+z2,
                                             y=0,
                                             sum.y=0,
                                             cgamma=z2*Inf))
                   
                   data2 <- data1[sum.y>0]
                   data2[,state:=droplevels(state)]
                   
                   t2 <- system.time({g2 <- try(glm(y~X+(state)-1, family=binomial, data=data2, x=T,y=T))})[3]
                   t4 <- system.time({g4 <- try(clogit(y~X+strata(state), data=data1))})[3]
                   
                   Z1 <- tapply(data1$X, data1$state, mean) %x% rep(1, T)
                   Z2 <- tapply(data2$X, data2$state, mean) %x% rep(1, T)
                   t5 <- system.time({g5 <- try(glmer(y~X+(1|state)+Z1, data=data1, family =binomial))})[3] #with irrelevant units
                   t6 <- system.time({g6 <- try(glmer(y~X+(1|state)+Z2, data=data2, family =binomial))})[3] #w/o irrelevant units

                   data2[,Z1 := mean(X), by =state]
                   data2[,Z2 := mean(X), by =state]
                   
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
             
                   if(class(g4)=="try-error"||!g4$iter<10|| any(is.na(g4$coef))||any(abs(g4$coef) > M)){
                     g4 <- list(coef=rep(NA,1)); t4 <- NA; t4 <- NA
                   }
                   if(class(g5)=="try-error"||g5@optinfo$conv$opt!=0||anyNA(coef(g5)$state)||abs(unique(coef(g5)$state$X)) > M){
                     beta.CRE <- NA;alpha.CRE <- rep(NA,N);mar5 <- NA;bmar5 <- NA; ame5 <- NA; t5 <- NA
                   }else{
                     beta.CRE <- g5@beta[2]
                     alpha.CRE <- rowSums(unique(cbind(1, Z1)) * as.matrix(coef(g5)$state[,-2]))
                     if(any(abs(c(beta.CRE, alpha.CRE))> M)){
                       beta.CRE <- alpha.CRE <- mar5 <- NA;bmar5 <- NA; ame5 <- NA
                     }else{
                       mar5 <- mean((predict(g5, type="response", newdata=data2) -plogis(-2*data2$X + data2$cgamma))^2)
                       bmar5 <-(mean(predict(g5, type="response", newdata=data2) -plogis(-2*data2$X + data2$cgamma)))
                       ame5 <- ((predict(g5, type="response", newdata=data2)*(1-predict(g5, type="response", newdata=data2))*(coef(g5)$state$X[1])))
                     }
                   }
                  
                   ###NOTE: RE here is the rCRE
                   if(class(g6)=="try-error"||g6@optinfo$conv$opt!=0||anyNA(coef(g6)$state)||abs(unique(coef(g6)$state$X)) > M){
                     beta.RE <- NA;alpha.RE <- rep(NA,N);mar6 <- NA;bmar6 <- NA; ame6 <- NA; t6 <- NA
                   }else{
                     beta.RE <- g6@beta[2]
                     alpha.RE <-rowSums(unique(cbind(1, Z2)) * as.matrix(coef(g6)$state[,-2]))
                     if(any(abs(c(beta.RE, alpha.RE))> M)){
                       beta.RE <- alpha.RE <- mar6 <- NA;bmar6 <- NA; ame6 <- NA
                     }else{
                       mar6 <- mean((predict(g6, type="response", newdata=data2) -plogis(-2*data2$X + data2$cgamma))^2)
                       bmar6 <-(mean(predict(g6, type="response", newdata=data2) -plogis(-2*data2$X + data2$cgamma)))
                       ame6 <- ((predict(g6, type="response", newdata=data2)*(1-predict(g6, type="response", newdata=data2))*(coef(g6)$state$X[1])))
                     }
                   }
                   
                   
                 
                   
                   keep <- as.numeric(as.character(unique(data1[!(sum.y==0 | sum.y==T)]$state)))
                   m2 <- mean((gamma*alpha[keep]-g2$coef[-1])^2)
                   m5 <- mean((gamma*alpha[keep]-alpha.CRE[keep])^2)
                   m6 <- mean((gamma*alpha[keep]-alpha.RE)^2)

                   b2 <-(mean(gamma*alpha[keep]-g2$coef[-1]))
                   b5 <-(mean(gamma*alpha[keep]-alpha.CRE[keep]))
                   b6 <-(mean(gamma*alpha[keep]-alpha.RE))
                   ame.truth <- (plogis(-2*data2$X + data2$cgamma)*(1-plogis(-2*data2$X + data2$cgamma))*-2)
                   
                   
                   truth.ame <-  mean(ame.truth)
                   
                   ame2 <- mean(ame2)
                   ame5 <- mean(ame5)
                   ame6 <- mean(ame6)

                   
                   
                   output <- c(g4$coef[1],
                               g2$coef[1],
                               beta.CRE,
                               beta.RE,
                               m2,   m5, m6,
                               b2,    b5, b6,
                               mar2,    mar5, mar6, 
                               bmar2,   bmar5, bmar6,
                               ame2,  ame5, ame6, truth.ame, 
                               t4, t2, t5, t6,
                               within.y, prop.all.zero, 
                               r)
                   output
                 }
  colnames(out) <- c("CML", "MLDV",  "CRE", "RE",
                     "mse.MLDV", "mse.CRE", "mse.RE",
                     "bias.MLDV", "bias.CRE",  "bias.RE",
                     "pred.MLDV", "pred.CRE", "pred.RE",
                     "bias.pred.MLDV","bias.pred.CRE", "bias.pred.RE", 
                     "ame2", "ame5", "ame6", "ame.truth",
                     "time.CML", "time.MLDV",   "time.CRE", "time.RE", 
                     "Prop.1", "Prop.dropped",
                     "corXc")
  out <- as.data.frame(out)
  out$bias.mem.MLDV <- out$ame2 - out$ame.truth
  out$bias.mem.CRE <- out$ame5 - out$ame.truth
  out$bias.mem.RE <- out$ame6 - out$ame.truth
  

  out$mse.mem.MLDV <- (out$ame2 - out$ame.truth)^2
  out$mse.mem.CRE <- (out$ame5 - out$ame.truth)^2
  out$mse.mem.RE <- (out$ame6 - out$ame.truth)^2
  results[[i]] <- as.matrix(out)

}
stopCluster(cl)


save.image("MCset_irrelevant0s.rdata")



rm(list=ls())
library(ggplot2)
library(matrixStats)
library(reshape2) 
library(scales)
library(stringr)
library(data.table)
load("MCset_irrelevant0s.rdata")

for(i in 1:length(results)){
  results[[i]][,"Prop.dropped"] <- (results[[i]][,"Prop.dropped"]*NT[i,1] + 9*NT[i,1])/(10*NT[i,1])
}
adInfo <- cbind(NT, do.call(rbind,lapply(results, function(x){return(cbind(mean(x[,"Prop.1"]), mean(x[,"Prop.dropped"])))})))
adInfo <- data.table(adInfo)
colnames(adInfo) <- c("N", "T", "Rare", "ybar", "dropped")
adInfo[,N:=NULL]
NT[,1] <- NT[,1]+9*NT[,1]

adInfo[,ybar:= round(mean(ybar),2), by=list(T, Rare)]
adInfo[,dropped:=round(mean(dropped),2)*100, by=list(T, Rare)]
adInfo$Rare <- ifelse(adInfo$Rare==g[1], "Rare event", "Non-rare Event")
adInfo$Rare <- factor(adInfo$Rare, levels=c("Rare event", "Non-rare Event"))
adInfo$T <- factor(paste("T =", adInfo$T), levels=paste("T =", t))
adInfo$label <- with(adInfo, paste("textstyle('%'~dropped) == ", dropped,  sep=""))



results <- lapply(results, as.matrix)
loadfonts(quiet=T)
p1data <- cbind(NT, 
                do.call(rbind, 
                        lapply(results,
                               function(x){
                                 sqrt(colVars(x[,c(
                                   "CML", "CRE", "MLDV", "RE")], na.rm=T) +
                                     (colMeans(x[,c(
                                       "CML", "CRE", "MLDV", "RE")], na.rm=T)+ 2)^2)})))
colnames(p1data)[colnames(p1data)=="RE"] <- c("restricted-CRE")
p1data <- melt(p1data, id.vars = c("Var1", "Var2", "Var3"))
colnames(p1data) <- c("N", "T", "Rare", "Estimator", "RMSE")
p1data$Rare <- ifelse(p1data$Rare==g[1], "Rare event", "Non-rare Event")
p1data$Rare <- factor(p1data$Rare, levels=c("Rare event", "Non-rare Event"))
p1data$T <- factor(paste("T =", p1data$T), levels=paste("T =", t))


p5data <- cbind(NT, do.call(rbind,
                            lapply(results,
                                   function(x){ colMeans(x[,c(
                                     "pred.CRE", "pred.MLDV", "pred.RE"
                                   )], na.rm=T)})))
colnames(p5data)[4:6] <- c("CRE", "MLDV", "restricted-CRE")
p5data <- melt(p5data, id.vars = c("Var1", "Var2", "Var3"))
colnames(p5data) <- c("N", "T", "Rare", "Estimator", "MSE")
p5data$Rare <- ifelse(p5data$Rare==g[1], "Rare event", "Non-rare Event")
p5data$Rare <- factor(p5data$Rare, levels=c("Rare event", "Non-rare Event"))
p5data$T <- factor(paste("T =", p5data$T), levels=paste("T =", t))



p8data <- cbind(NT, do.call(rbind,
                            lapply(results,
                                   function(x){colMeans( (x[,c( 
                                     "mse.mem.CRE", 
                                     "mse.mem.MLDV", "mse.mem.RE")]), na.rm=T)})))
colnames(p8data)[4:6] <- c("CRE",  "MLDV", "restricted-CRE")
p8data <- melt(p8data, id.vars = c("Var1", "Var2", "Var3"))
colnames(p8data) <- c("N", "T", "Rare", "Estimator", "MSE")
p8data$Rare <- ifelse(p8data$Rare==g[1], "Rare event", "Non-rare Event")
p8data$Rare <- factor(p8data$Rare, levels=c("Rare event", "Non-rare Event"))
p8data$T <- factor(paste("T =", p8data$T), levels=paste("T =", t))



g1 <- ggplot(p1data)+ #RMSE in beta
  geom_line(aes(x=N, y=(RMSE), color=Estimator, linetype=Estimator), size=1)+
  facet_grid(Rare~T )+
  geom_text(aes(x=950,y=Inf, label=label),
            data=adInfo, parse=TRUE, vjust=1.1, size=5, hjust=.9, check_overlap = TRUE) +
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



g5 <- ggplot(p5data)+ #rmse in pr
  geom_line(aes(x=N, y=sqrt(MSE), color=Estimator, linetype=Estimator), size=1)+
  facet_grid(Rare~T )+
  geom_text(aes(x=950,y=Inf, label=label),
            data=adInfo, parse=TRUE, vjust=1.1, size=5, hjust=.9, check_overlap = TRUE) +
  theme_bw(16)+
  ylab('RMSE')+ ylim(0, .45)+
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


g8 <- ggplot(p8data)+ #rmse in marginal
  geom_line(aes(x=N, y=sqrt(MSE), color=Estimator, linetype=Estimator), size=1)+
  facet_grid(Rare~T )+
  geom_text(aes(x=950,y=Inf, label=label),
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
print(g8)


h=7; w=9.5
pdf(file="FigureD5.pdf", height=h, width=w); g1; dev.off()
pdf(file="FigureD6.pdf", height=h, width=w); g5; dev.off()
pdf(file="FigureD7.pdf", height=h, width=w); g8; dev.off()
