library(survival)
library(lme4)
library(readstata13)
library(data.table)
library(MASS)
library(lmtest)
library(stargazer)
library(ggplot2)
library(extrafont)
rm(list=ls())

loadfonts(quiet = TRUE)
# source("sparseFirth.R")
# source("xtlogit_stan.r") 


emw <- data.table(read.dta13("temp-micro.dta"))
emw[,proUse:=mean(progov),by=dcode] #to match their coding
cml.obs <- copy(emw)
cml.obs[,mean.y:=mean(protest01, na.rm=T), by=dcode]
cml.obs[,sum.y:=sum(protest01, na.rm=T), by=dcode]
cml.obs <- cml.obs[mean.y>0 & mean.y < 1]


### COLUMN 4####
t4 <- system.time({m4 <- clogit(protest01~remit + remit:proUse+cellphone+ lage+ education+ wealth+ male+ 
                                  employment+ travel+strata(dcode),  
                                data=cml.obs,x=T)})
coeftest(m4)

# t4.pml <- system.time(m4.pml <- sp.ffe(protest01~remit + cellphone+ lage+ education+ wealth+ male+ 
#                                          employment+ travel+remitXpro+factor(dcode)+0, data=emw))
# head(m4.pml$coef.table,9)


#### CRE #####
emw[,`:=`(remit.bar = mean(remit,na.rm=T),
          cellphone.bar = mean(cellphone,na.rm=T),
          lage.bar = mean(lage,na.rm=T),
          education.bar = mean(education,na.rm=T),
          wealth.bar = mean(wealth,na.rm=T),
          male.bar = mean(male,na.rm=T),
          employment.bar = mean(employment,na.rm=T),
          travel.bar = mean(travel,na.rm=T),
          remitXpro.bar = mean(remitXpro,na.rm=T)),
    by = dcode]
cml.obs[,`:=`(remit.bar = mean(remit,na.rm=T),
          cellphone.bar = mean(cellphone,na.rm=T),
          lage.bar = mean(lage,na.rm=T),
          education.bar = mean(education,na.rm=T),
          wealth.bar = mean(wealth,na.rm=T),
          male.bar = mean(male,na.rm=T),
          employment.bar = mean(employment,na.rm=T),
          travel.bar = mean(travel,na.rm=T),
          remitXpro.bar = mean(remitXpro,na.rm=T)),
    by = dcode]
mle.cre.out <- glmer(protest01 ~ remit + remit:proUse + cellphone + lage +
                       education + wealth + male + employment+ travel +
                       remit.bar + remitXpro.bar +cellphone.bar + lage.bar +
                       education.bar + wealth.bar + male.bar + employment.bar +
                       travel.bar + (1|dcode),
                     data=emw, family=binomial(), nAGQ = 12,
                     control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e5)))

mle.rcre.out <- glmer(protest01 ~ remit + remit:proUse + cellphone + lage +
                       education + wealth + male + employment+ travel +
                       remit.bar + remitXpro.bar +cellphone.bar + lage.bar +
                       education.bar + wealth.bar + male.bar + employment.bar +
                       travel.bar + (1|dcode),
                     data=cml.obs, family=binomial(), nAGQ = 12,
                     control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e5)))

# #time without compile
# y1 <- with(emw, cbind(protest01,dcode))
# X1 <- with(emw, cbind(remit, cellphone, lage, education ,wealth, male, employment, travel, remitXpro))
# Xbar1 <- with(emw, cbind(remit.bar, cellphone.bar, lage.bar, education.bar, 
#                          wealth.bar, male.bar, employment.bar, travel.bar, remitXpro.bar))
# 
# XX <- na.omit(cbind(y1, X1,Xbar1))
# X1 <- XX[, colnames(X1)]
# Xbar1 <- XX[, colnames(Xbar1)]
# y1 <- XX[, colnames(y1)]
# group1 <- y1[,"dcode"]
# y1 <- y1[,"protest01"]
# 
# t4.stan.cre <- system.time({m4.stan.cre <- CRE(y=y1, X=X1, Xbar=Xbar1, group.id = group1, 
#                                                chains=4, cores=4, seed=1234567, iter=5000, thin=4)})
# print(m4.stan.cre$info$Rhat)
# summary(m4.stan.cre$FE$Rhat)
# summary(m4.stan.cre$RE$Rhat)
# summary(m4.stan.cre$info$MCMC)
# beta.cre <- m4.stan.cre$FE$coefficients[1:9]
# se.cre <- sqrt(diag(m4.stan.cre$FE$vcov)[1:9])
# z.cre <- beta.cre/se.cre
# p.cre <- 2*pnorm(abs(z.cre), lower=F)
# round(cbind(beta.cre,se.cre,z.cre,p.cre),3)





t4.glm2 <- system.time(m4.glm2 <- glm(protest01~remit + remit:proUse + cellphone+ lage+ education+ wealth+ male+ 
                                        employment+ travel+factor(dcode)+0, 
                                      data=cml.obs,family=binomial,x=T))

cre.margins <- summary(margins(mle.cre.out, variables="remit", at=list(proUse=c(.18,.8)), data=as.data.frame(emw)))[c("lower", "upper", "AME")]*5
cre.cml.margins <- summary(margins(mle.cre.out, variables="remit", at=list(proUse=c(.18,.8)), data=as.data.frame(cml.obs)))[c("lower", "upper", "AME")]*5
rcre.margins <- summary(margins(mle.rcre.out, variables="remit", at=list(proUse=c(.18,.8)),  data=as.data.frame(cml.obs)))[c("lower", "upper", "AME")]*5
mldv.margins <- summary(margins(m4.glm2, variables="remit",   at=list(proUse=c(.18,.8)), data=as.data.frame(cml.obs)))[c("lower", "upper", "AME")]*5

margins.out <- rbind.data.frame(cre.margins, cre.cml.margins, rcre.margins, mldv.margins)
colnames(margins.out) <- c("lo", "hi", "Estimate")
margins.out$x <-  c("Opposition\n18% Progovernment",
                    "Stronghold\n80% Progovernment")
margins.out$Estimator <- rep(c("CRE-AME", "CRE-cAME", "rCRE", "MLDV"), each=2)
save.image("emwOut.rdata")

rm(list=ls())
load("emwOut.rdata")
library(ggplot2)
library(MASS)
library(lme4)
library(survival)
library(stargazer)
library(extrafont)
loadfonts(quiet=TRUE)
set.seed(1)

##### EMW's "Marginal Effects"#####
m4.boot <- 	mvrnorm(1000, m4$coef, vcov(m4)) 

# 
# lo = quantile(emw$mean_progov, probs=.1) #lo progov
# hi = quantile(emw$mean_progov, probs=.9) #hi progov
lo=.18
hi=.8
r = 5  # linear interpolation (moving up five units)

marg.remitO <- rowMeans(t(plogis(m4$x %*% t(m4.boot)))*(m4.boot[,"remit"] + m4.boot[,"remit:proUse"]*lo))*r
marg.remitS <- rowMeans(t(plogis(m4$x %*% t(m4.boot)))*(m4.boot[,"remit"] + m4.boot[,"remit:proUse"]*hi))*r

plot2.dataORIGINGAL <- rbind.data.frame(c(quantile(marg.remitO, probs=c(0.025,.975)), mean(marg.remitO)),      
                                        c(quantile(marg.remitS, probs=c(0.025,.975)), mean(marg.remitS)))
colnames(plot2.dataORIGINGAL) <- c("lo", "hi", "Estimate")
plot2.dataORIGINGAL$x <- c("Opposition\n18% Progovernment",
                           "Stronghold\n80% Progovernment")


##### What I think EMW intended ####
marg.remitO <- rowMeans(t(dlogis(m4$x %*% t(m4.boot)))*(m4.boot[,"remit"] + m4.boot[,"remit:proUse"]*lo))*r
marg.remitS <- rowMeans(t(dlogis(m4$x %*% t(m4.boot)))*(m4.boot[,"remit"] + m4.boot[,"remit:proUse"]*hi))*r

plot2.dataFIXED1 <- rbind.data.frame(c(quantile(marg.remitO, probs=c(0.025,.975)), mean(marg.remitO)),      
                                     c(quantile(marg.remitS, probs=c(0.025,.975)), mean(marg.remitS)))
colnames(plot2.dataFIXED1) <- c("lo", "hi", "Estimate")
plot2.dataFIXED1$x <- c("Opposition\n18% Progovernment",
                        "Stronghold\n80% Progovernment")
plot2.dataFIXED1 <- rbind(plot2.dataFIXED1, plot2.dataORIGINGAL)
plot2.dataFIXED1$Result <- rep(c("Without typo", "With typo"),each=2)
emw.plot1 <- ggplot(plot2.dataFIXED1)+
  geom_pointrange(aes(x=x, y=Estimate, ymin=lo,ymax=hi,color=Result, shape=Result, linetype=Result),
                  position=position_dodge(width=.5), size=1)+
  geom_abline(aes(intercept=0, slope=0), linetype="dashed")+
  xlab("")+
  ylab("Marginal effect on\nthe probability of protest")+
  theme_bw(16)+
  theme(text=element_text(family="CM Roman"),
        legend.position = "bottom")

print(emw.plot1)
h=5; w=6

ggsave(emw.plot1, file="/home/cox/Dropbox/xtlogit/figures/emw2.pdf", height=h, width=w)
embed_fonts(file="/home/cox/Dropbox/xtlogit/figures/emw2.pdf")

# 
# 
# 
# ## Marginals from PML and CRE
# pml4.boot <- 	mvrnorm(1000, m4.pml$coefficients, m4.pml$vcov)
# 
# 
# lr.model <- sp.ffe(protest01~remit + cellphone+ lage+ education+ wealth+ male+ employment+ travel+remitXpro+0, data=emw)
# #LR test of whether the c_i are all 0?
# pchisq(2*(m4.pml$logLik-lr.model$logLik), lower=F,
#        df=length(m4.pml$coefficients)-length(lr.model$coefficients)) 


# X <- sparse.model.matrix(protest01~remit + cellphone+ lage+ education+ wealth+ male+ employment+ travel+remitXpro+factor(dcode)+0, data=emw)
# 
# marg.remitO <- rowMeans(t(dlogis(as.matrix(X %*% t(pml4.boot))))*(pml4.boot[,"remit"] + pml4.boot[,"remitXpro"]*lo))*r
# marg.remitS <- rowMeans(t(dlogis(as.matrix(X %*% t(pml4.boot))))*(pml4.boot[,"remit"] + pml4.boot[,"remitXpro"]*hi))*r
# 
# plot2.data <- rbind.data.frame(c(quantile(marg.remitO, probs=c(0.025,.975)), mean(marg.remitO)),      
#                                c(quantile(marg.remitS, probs=c(0.025,.975)), mean(marg.remitS)))
# colnames(plot2.data) <- c("lo", "hi", "Estimate")
# plot2.data$x <- c("Opposition\n18% Progovernment",
#                   "Stronghold\n80% Progovernment")



# ## CRE
# m4.stan.boot <- do.call(rbind,m4.stan.cre$FE$MCMC)
# marg.remitO <- rowMeans(t(dlogis(as.matrix(X %*% t(m4.stan.boot))))*(m4.stan.boot[,"remit"] + m4.stan.boot[,"remitXpro"]*lo))*r
# marg.remitS <- rowMeans(t(dlogis(as.matrix(X %*% t(m4.stan.boot))))*(m4.stan.boot[,"remit"] + m4.stan.boot[,"remitXpro"]*hi))*r
# 
# plot2.data.stan <- rbind.data.frame(c(quantile(marg.remitO, probs=c(0.025,.975)), mean(marg.remitO)),      
#                                     c(quantile(marg.remitS, probs=c(0.025,.975)), mean(marg.remitS)))
# colnames(plot2.data.stan) <- c("lo", "hi", "Estimate")
# plot2.data.stan$x <- c("Opposition\n18% Progovernment",
#                        "Stronghold\n80% Progovernment")
# 
# plot2.data.stan <- rbind(plot2.data.stan, plot2.data, plot2.dataFIXED1[1:2,-5])
# plot2.data.stan$Estimator <- factor(rep(c("CRE", "PML", "EMW"), each=2), levels=c("EMW", "PML", "CRE"))


## MLDV for comparison
# ml4.glm2.boot <- 	mvrnorm(1000, m4.glm2$coefficients, vcov(m4.glm2))
# 
# 
# X.cml <- m4.glm2$x
# 
# marg.remitO <- rowMeans(t(dlogis(as.matrix(X.cml %*% t(ml4.glm2.boot))))*(ml4.glm2.boot[,"remit"] + ml4.glm2.boot[,"remitXpro"]*lo))*r
# marg.remitS <- rowMeans(t(dlogis(as.matrix(X.cml %*% t(ml4.glm2.boot))))*(ml4.glm2.boot[,"remit"] + ml4.glm2.boot[,"remitXpro"]*hi))*r
# 
# plot2.dataGLM <- rbind.data.frame(c(quantile(marg.remitO, probs=c(0.025,.975)), mean(marg.remitO)),
#                                   c(quantile(marg.remitS, probs=c(0.025,.975)), mean(marg.remitS)))
# colnames(plot2.dataGLM) <- c("lo", "hi", "Estimate")
# plot2.dataGLM$x <- c("Opposition\n18% Progovernment",
#                      "Stronghold\n80% Progovernment")
# plot2.dataGLM$Estimator <- "MLDV"
# 
# 
# 
# margins.full.est <- as.numeric(read.table("../stataFiles/emw_fullSample_margin.txt", header = TRUE)[1,])
# margins.full.se <-sqrt(diag(as.matrix(read.table("../stataFiles/emw_fullSample_margin_vcov.txt", header = TRUE)[1:2,1:2])))
# margins.full.z <- margins.full.est/margins.full.se
# margins.full.p <- pnorm(abs(margins.full.z), lower=F)*2
# margins.full.low <- margins.full.est-qnorm(.975)*margins.full.se
# margins.full.high <- margins.full.est+qnorm(.975)*margins.full.se
# cre.full.margins <- data.frame(lo=margins.full.low,
#                                hi=margins.full.high,
#                                Estimate=margins.full.est)
# cre.full.margins$x <- c("Opposition\n18% Progovernment",
#                      "Stronghold\n80% Progovernment")
# cre.full.margins$Estimator <- "CRE-AME"
# print(cre.full.margins)
# 
# margins.cml.est <- as.numeric(read.table("../stataFiles/emw_CMLSample_margin.txt", header = TRUE)[1,])
# margins.cml.se <-sqrt(diag(as.matrix(read.table("../stataFiles/emw_CMLSample_margin_vcov.txt", header = TRUE)[1:2,1:2])))
# margins.cml.z <- margins.cml.est/margins.cml.se
# margins.cml.p <- pnorm(abs(margins.cml.z), lower=F)*2
# margins.cml.low <- margins.cml.est-qnorm(.975)*margins.cml.se
# margins.cml.high <- margins.cml.est+qnorm(.975)*margins.cml.se
# cre.cml.margins <-data.frame(lo=margins.cml.low,
#                              hi=margins.cml.high,
#                              Estimate=margins.cml.est)
# cre.cml.margins$x <- c("Opposition\n18% Progovernment",
#                         "Stronghold\n80% Progovernment")
# cre.cml.margins$Estimator <- "rCRE"
# print(cre.cml.margins)
# 
# margins.full2.est <- as.numeric(read.table("../stataFiles/emw_fullSample_margin2.txt", header = TRUE)[1,])
# margins.full2.se <-sqrt(diag(as.matrix(read.table("../stataFiles/emw_fullSample_margin_vcov2.txt", header = TRUE)[1:2,1:2]))) #lost something 
# margins.full2.z <- margins.full2.est/margins.full2.se
# margins.full2.p <- pnorm(abs(margins.full2.z), lower=F)*2
# margins.full2.low <- margins.full2.est-qnorm(.975)*margins.full2.se
# margins.full2.high <- margins.full2.est+qnorm(.975)*margins.full2.se
# cre.full2.margins <-data.frame(lo=margins.full2.low,
#                              hi=margins.full2.high,
#                              Estimate=margins.full2.est)
# cre.full2.margins$x <- c("Opposition\n18% Progovernment",
#                        "Stronghold\n80% Progovernment")
# cre.full2.margins$Estimator <- "CRE-cAME"
# print(cre.full2.margins)



# plot2.dataGLM <- data.table(rbind(plot2.dataGLM,cre.cml.margins, cre.full2.margins, cre.full.margins))
margins.out$Estimator <- factor(margins.out$Estimator,
                         levels=c("CRE-AME", "CRE-cAME", "rCRE", "MLDV"))
full.plot <- ggplot(margins.out)+
  geom_pointrange(aes(x=x, y=Estimate, ymin=lo,ymax=hi,color=Estimator, shape=Estimator, linetype=Estimator), position=position_dodge(width=.5), size=1)+
  geom_abline(aes(intercept=0, slope=0), linetype="dashed")+
  xlab("")+
  ylab("Marginal effect on\nthe probability of protest")+
  theme_bw(16)+
  theme(text=element_text(family="CM Roman"),
        legend.position = "bottom")

print(full.plot)


ggsave(full.plot, file="/home/cox/Dropbox/xtlogit/figures/emw2019.pdf", height=1.25*h, width=1.5*w)
embed_fonts(file="/home/cox/Dropbox/xtlogit/figures/emw2019.pdf")



# beta.cre <- as.numeric(read.table("../stataFiles/emw_fullSample_est.txt", header = TRUE)[1,])[1:9]
# beta.cre.restricted <- as.numeric(read.table("../stataFiles/emw_CMLSample_est.txt", header = TRUE)[1,])[1:9]
# names(beta.cre.restricted) <- c("remit", "remitXpro" ,"cellphone", "lage" ,"education" ,"wealth" ,"male" ,"employment" ,"travel")
# names(beta.cre) <-c("remit", "remitXpro" ,"cellphone", "lage" ,"education" ,"wealth" ,"male" ,"employment" ,"travel")
# 
# 
# vcov.cre <- as.matrix(read.table("../stataFiles/emw_fullSample_vcov.txt", header = TRUE))[1:9, 1:9]
# vcov.cre.restricted <- as.matrix(read.table("../stataFiles/emw_CMLSample_vcov.txt", header = TRUE))[1:9, 1:9]
# rownames(vcov.cre) <- colnames(vcov.cre) <- names(beta.cre)
# rownames(vcov.cre.restricted) <- colnames(vcov.cre.restricted) <- names(beta.cre.restricted)
# se.cre <- sqrt(diag(vcov.cre))
# z.cre <- beta.cre/se.cre
# p.cre <- 2*pnorm(abs(z.cre), lower=F)
# round(cbind(beta.cre,se.cre,z.cre,p.cre),2)
# 
# se.cre.restricted <- sqrt(diag(vcov.cre.restricted))
# z.cre.restricted <- beta.cre.restricted/se.cre.restricted
# p.cre.restricted <- 2*pnorm(abs(z.cre.restricted), lower=F)
# round(cbind(beta.cre.restricted,se.cre.restricted,z.cre.restricted,p.cre.restricted),3)
# library(stargazer)

stargazer(m4,m4.glm2,m4.glm2, m4.glm2,
          # covariate.labels = 
          omit=c("factor*", "*.bar", "Constant"),
          omit.stat=c("all"), 
          title="Replicating ",
          label="tab.dirtypool",
          star.cutoffs = c(NA),
          digits=2,
          no.space=TRUE,
          coef=list(NULL, NULL, 
                    summary(mle.rcre.out)$coefficients[names(m4$coef), 1],
                    summary(mle.cre.out)$coefficients[names(m4$coef), 1]),
          se=list(NULL, NULL, 
                  summary(mle.rcre.out)$coefficients[names(m4$coef), 2],
                  summary(mle.cre.out)$coefficients[names(m4$coef), 2]))

# 
# (m4$coef-beta.cre[names(m4$coef)]) %*% solve(vcov(m4)- vcov.cre[names(m4$coef),names(m4$coef)]) %*% (m4$coef-beta.cre[names(m4$coef)])
# (m4$coef-m4.glm2$coef[names(m4$coef)]) %*% solve(vcov(m4)- vcov(m4.glm2)[names(m4$coef),names(m4$coef)]) %*% (m4$coef-m4.glm2$coef[names(m4$coef)])
