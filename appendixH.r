###Code taken from:
# Cook, Scott; Hays, Jude; Franzese, Rob, 2018, "Replication Data for: Fixed Effects in Rare Events Data:
# A Penalized Maximum Likelihood Solution", https://doi.org/10.7910/DVN/H7QR3Q, Harvard Dataverse, V1, 
# UNF:6:sFrVI667SmdP3yNpG9WFZA== [fileUNF] 
###Updated by Crisman-Cox for Appendix H and added parallelization



library(lme4)
library(survival)
library(MASS)
library(brglm)
library(doParallel)
library(doRNG)

rm(list = ls())

trials = 1000
n <-cbind(50,100)
t <-cbind(20,50)

cor <- cbind(0.25,0.5)
sigma_b <- 1
sigma_w <- cbind(1,2)
beta <- 1

sim_par <- data.frame(expand.grid(n, t, cor, sigma_w))
pc <- length(sim_par[,1])

rmse.beta <- matrix(data=NA,nrow=pc,ncol=5)
rmse.mfx <-  matrix(data=NA,nrow=pc,ncol=4)
av_cens <- matrix(data=NA,nrow=pc,ncol=1)


cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)
for(j in 1:pc){
  set.seed(123456)
  nt <- sim_par[j,1]*sim_par[j,2]
  tscs_df <- expand.grid(time = 1:sim_par[j,2], unit = 1:sim_par[j,1])
  
  omega = matrix(c(1,sim_par[j,3],sim_par[j,3],1),nrow=2)
  tau = diag(c(sqrt(1),sqrt(sigma_b)))
  covar = t(tau)%*%omega%*%tau  
  u = mvrnorm(sim_par[j,1], c(-4,0), covar)
  
  alpha <- u[,1] 
  x.bar <- u[,2]
  
  for (i in 1:nt) {
    tscs_df$alpha[tscs_df$unit == i] <- alpha[i]
    tscs_df$x.bar[tscs_df$unit == i] <- x.bar[i]
  }
  
  
  tscs_df$x <- rnorm(nt,tscs_df$x.bar,sim_par[j,4]) 
  
  z <- tscs_df$alpha + beta*tscs_df$x      
  pr <- 1/(1+exp(-z))         
  
  beta1.est<- matrix(data=NA,nrow=trials,ncol=5)
  mfx<- matrix(data=NA,nrow=trials,ncol=4)
  cens<- matrix(data=NA,nrow=trials,ncol=1)
  
  ests <- foreach(k = 1:trials, .combine = "rbind",
                  .packages = c("survival", "lme4", "brglm")) %dorng% {
                    
                    tscs_df$y <- rbinom(nt,1,pr)      
                    
                    y_panel <- aggregate(y ~ unit, tscs_df, mean)

                    y_panel$censor <- ifelse(y_panel$y == 0, 1, 0)
                    y_panel$uncensor <- ifelse(y_panel$y > 0, 1, 0)
                    y_panel <-cbind(y_panel, y_panel$censor*u)
                    y_panel <-cbind(y_panel, y_panel$uncensor*u)
                    
                    for (i in 1:nrow(y_panel)) {
                      tscs_df$censor[tscs_df$unit == y_panel$unit[i]] <- y_panel$censor[i]
                    }
                    
                    tscs_df$obs.id <- tscs_df$unit 
                    tscs_df$uncensor <- 1- tscs_df$censor
                    tscs_df$obs.id_un <- tscs_df$obs.id*tscs_df$uncensor
                    tscs_df$z <- unlist(by(tscs_df, tscs_df$unit, function(x){rep(mean(x$x), each=sim_par[j,2])}))
                    tscs_df2 <- tscs_df[ which(tscs_df$uncensor=='1'), ]
                    
                    
                    naive <- glm(y ~ x, data = tscs_df, family = binomial(link = "logit"))
                    uncfe <- glm(y ~ x + factor(unit), data = tscs_df2, family = binomial(link = "logit"))
                    confe <- clogit(y ~ x + strata(unit), data = tscs_df)
                    re    <- glmer(y ~ x + z + (1|unit), data = tscs_df, family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa"),nAGQ = 10)
                    pmlfe <- brglm(y ~ x + factor(obs.id_un), family=binomial(link="logit"), method = "brglm.fit", p1 = T, data = tscs_df)
                    
                    
                    mfx.naive <- mean(predict(naive,type="response")*(1-predict(naive,type="response"))*summary(naive)$coeff[2])
                    mfx.uncfe <- mean(predict(uncfe,type="response")*(1-predict(uncfe,type="response"))*summary(uncfe)$coeff[2])
                    mfx.re <- mean(predict(re,type="response")*(1-predict(re,type="response"))*summary(re)$coeff[2])
                    mfx.pmlfe <- mean(predict(pmlfe,type="response")*(1-predict(pmlfe,type="response"))*summary(pmlfe)$coeff[2])
                    
                    c(cbind(summary(naive)$coeff[2], summary(uncfe)$coeff[2], summary(confe)$coeff[1], summary(re)$coeff[2], summary(pmlfe)$coeff[2]),
                      cbind(mfx.naive,mfx.uncfe,mfx.re,mfx.pmlfe),
                      mean(y_panel$censor))
                  }
  beta1.est <- ests[,1:5]
  mfx <- ests[,6:9]
  cens <- ests[,10]
  
  
  rmse.beta[j,] <- mapply(function(m) mean(sqrt((1-beta1.est[,m])^2), na.rm=T),seq(1:5))
  true.mfx <- mean(pr*(1-pr)*beta)
  rmse.mfx[j,] <- mapply(function(m) mean(sqrt((true.mfx-mfx[,m])^2), na.rm=T),seq(1:4))
  av_cens[j,] <- mean(cens)
  
  rm(beta1.est,mfx,cens)
  
}
stopCluster(cl)
save.image("MC_pml.rdata")

rm(list=ls())
load("MC_pml.rdata")
library(xtable)
#Table 1
table1_n50_t20 <- mapply(function(m) rmse.beta[5,5]/rmse.beta[5,m],cbind(1,4,2,3))                       
table1_n100_t20 <- mapply(function(m) rmse.beta[6,5]/rmse.beta[6,m],cbind(1,4,2,3))                       
table1_n50_t50 <- mapply(function(m) rmse.beta[7,5]/rmse.beta[7,m],cbind(1,4,2,3))                       
table1_n100_t50 <- mapply(function(m) rmse.beta[8,5]/rmse.beta[8,m],cbind(1,4,2,3))                       
table1_results <- round(rbind(table1_n50_t20,table1_n100_t20, table1_n50_t50, table1_n100_t50),2)
table1_results <- cbind(c("T=20, T=50", "T=50, N=100","T=20, T=100","T=50, N=100"), table1_results)
colnames(table1_results) <- c(" ", "Pooled", "CRE", "MLDV", "CML")
print(xtable(table1_results, align='llcccc'), sanitize.text.function=function(x){x}, include.rownames=FALSE)

#Table 2 
table2_n50_t20 <- mapply(function(m) rmse.mfx[5,4]/rmse.mfx[5,m],cbind(1,3,2))                       
table2_n100_t20 <- mapply(function(m) rmse.mfx[6,4]/rmse.mfx[6,m],cbind(1,3,2))                       
table2_n50_t50 <- mapply(function(m) rmse.mfx[7,4]/rmse.mfx[7,m],cbind(1,3,2))                       
table2_n100_t50 <- mapply(function(m) rmse.mfx[8,4]/rmse.mfx[8,m],cbind(1,3,2))                       
table2_results <- round(rbind(table2_n50_t20,table2_n100_t20, table2_n50_t50, table2_n100_t50),2)
table2_results <- cbind(table2_results,round(av_cens[5:8,1]*100))
table2_results <- cbind(c("T=20, T=50", "T=50, N=100","T=20, T=100","T=50, N=100"), table2_results)
colnames(table2_results) <- c(" ", "Pooled", "CRE", "MLDV", "Censoring (\\%)")
print(xtable(table2_results, align='llcccc'), sanitize.text.function=function(x){x}, include.rownames=FALSE)


#Table 3 
table3_n50_t20 <- mapply(function(m) rmse.mfx[1,4]/rmse.mfx[1,m],cbind(1,3,2))                       
table3_n100_t20 <- mapply(function(m) rmse.mfx[2,4]/rmse.mfx[2,m],cbind(1,3,2))                       
table3_n50_t50 <- mapply(function(m) rmse.mfx[3,4]/rmse.mfx[3,m],cbind(1,3,2))                       
table3_n100_t50 <- mapply(function(m) rmse.mfx[4,4]/rmse.mfx[4,m],cbind(1,3,2))                       
table3_results <- round(rbind(table3_n50_t20,table3_n100_t20, table3_n50_t50, table3_n100_t50),2)
table3_results <- cbind(table3_results,round(av_cens[1:4,1]*100))
table3_results <- cbind(c("T=20, T=50", "T=50, N=100","T=20, T=100","T=50, N=100"), table3_results)
colnames(table3_results) <- c(" ", "Pooled", "CRE", "MLDV", "Censoring (\\%)")
print(xtable(table3_results, align='llcccc'), sanitize.text.function=function(x){x}, include.rownames=FALSE)


#Table 4 
table4_n50_t20 <- mapply(function(m) rmse.mfx[13,4]/rmse.mfx[13,m],cbind(1,3,2))                       
table4_n100_t20 <- mapply(function(m) rmse.mfx[14,4]/rmse.mfx[14,m],cbind(1,3,2))                       
table4_n50_t50 <- mapply(function(m) rmse.mfx[15,4]/rmse.mfx[15,m],cbind(1,3,2))                       
table4_n100_t50 <- mapply(function(m) rmse.mfx[16,4]/rmse.mfx[16,m],cbind(1,3,2))                       
table4_results <- round(rbind(table4_n50_t20,table4_n100_t20, table4_n50_t50, table4_n100_t50),2)
table4_results <- cbind(table4_results,round(av_cens[13:16,1]*100))
table4_results <- cbind(c("T=20, T=50", "T=50, N=100","T=20, T=100","T=50, N=100"), table4_results)
colnames(table4_results) <- c(" ", "Pooled", "CRE", "MLDV", "Censoring (\\%)")
print(xtable(table4_results, align='llcccc'), sanitize.text.function=function(x){x}, include.rownames=FALSE)


#Table 5
table5_n50_t20 <- mapply(function(m) rmse.mfx[9,4]/rmse.mfx[9,m],cbind(1,3,2))                       
table5_n100_t20 <- mapply(function(m) rmse.mfx[10,4]/rmse.mfx[10,m],cbind(1,3,2))                       
table5_n50_t50 <- mapply(function(m) rmse.mfx[11,4]/rmse.mfx[11,m],cbind(1,3,2))                       
table5_n100_t50 <- mapply(function(m) rmse.mfx[12,4]/rmse.mfx[12,m],cbind(1,3,2))                       
table5_results <- round(rbind(table5_n50_t20,table5_n100_t20, table5_n50_t50, table5_n100_t50),2)
table5_results <- cbind(table5_results, round(av_cens[9:12,1]*100))
table5_results <- cbind(c("T=20, T=50", "T=50, N=100","T=20, T=100","T=50, N=100"), table5_results)
colnames(table5_results) <- c(" ", "Pooled", "CRE", "MLDV", "Censoring (\\%)")
print(xtable(table5_results, align='llcccc'), sanitize.text.function=function(x){x}, include.rownames=FALSE)

