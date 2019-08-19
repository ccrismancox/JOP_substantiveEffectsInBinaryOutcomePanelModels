library(survival)

library(foreign)
library(margins)
library(lmtest)
library(data.table)
library(MASS)
library(stringr)
library(ggplot2)
library(extrafont)

rm(list=ls())

disp <- data.table(read.dta("disputes1.dta"))
disp[, dyadid.full:=str_pad(as.character(dyadid), 6, side="left", pad="0")]
disp[, max.contig := max(contig), by=dyadid]
disp[, ccode1:=str_sub(dyadid.full, 0, 3)]
disp[, ccode2:=str_sub(dyadid.full, 4, -1)]
disp <- disp[nyears >=20]
disp[,sum.y:=sum(dispute), by=dyadid]
disp[,dyadid2:=ifelse(sum.y==0,-100,dyadid)]
disp2 <- disp[sum.y >0]
major <- c("002", "200", "220", "365", "710")
disp.relevant <- subset(disp, ccode1 %in% major | ccode2 %in% major | max.contig ==1  | sum.y>0)
disp.relevant2 <- subset(disp, ccode1 %in% major | ccode2 %in% major | max.contig ==1)




###CML####
#clogit drops 2,877 groups due to no war
t1 <- system.time({m1 <- clogit(dispute ~contig+ lcaprat+ mingrow+ allied +mindem+ mineco1+strata(dyadid), data=disp2)})
print(coeftest(m1))


#### CRE #####
disp[,`:=`(contig.bar = mean(contig),
           lcaprat.bar = mean(lcaprat),
           mingrow.bar = mean(mingrow),
           allied.bar = mean(allied),
           mindem.bar = mean(mindem),
           mineco1.bar = mean(mineco1)),
     by = dyadid]
disp2[,`:=`(contig.bar = mean(contig),
           lcaprat.bar = mean(lcaprat),
           mingrow.bar = mean(mingrow),
           allied.bar = mean(allied),
           mindem.bar = mean(mindem),
           mineco1.bar = mean(mineco1)),
     by = dyadid]
disp.relevant[,`:=`(contig.bar = mean(contig),
            lcaprat.bar = mean(lcaprat),
            mingrow.bar = mean(mingrow),
            allied.bar = mean(allied),
            mindem.bar = mean(mindem),
            mineco1.bar = mean(mineco1)),
      by = dyadid]
disp.relevant2[,`:=`(contig.bar = mean(contig),
                    lcaprat.bar = mean(lcaprat),
                    mingrow.bar = mean(mingrow),
                    allied.bar = mean(allied),
                    mindem.bar = mean(mindem),
                    mineco1.bar = mean(mineco1)),
              by = dyadid]
cre1 <- glmer(dispute ~contig+ lcaprat+ mingrow+ allied +mindem+ mineco1+
              contig.bar + lcaprat.bar + mingrow.bar + allied.bar + mindem.bar + mineco1.bar+
                (1|dyadid), family=binomial, data=disp,  nAGQ = 12,
              control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=500000)))
rcml1 <- glmer(dispute ~contig+ lcaprat+ mingrow+ allied +mindem+ mineco1+
                contig.bar + lcaprat.bar + mingrow.bar + allied.bar + mindem.bar + mineco1.bar+
                (1|dyadid), family=binomial, data=disp2,  nAGQ = 12,
              control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=500000)))
rcml.relevant <- glmer(dispute ~contig+ lcaprat+ mingrow+ allied +mindem+ mineco1+
                 contig.bar + lcaprat.bar + mingrow.bar + allied.bar + mindem.bar + mineco1.bar+
                 (1|dyadid), family=binomial, data=disp.relevant,  nAGQ = 12,
               control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=500000)))
rcml.relevant2 <- glmer(dispute ~contig+ lcaprat+ mingrow+ allied +mindem+ mineco1+
                         contig.bar + lcaprat.bar + mingrow.bar + allied.bar + mindem.bar + mineco1.bar+
                         (1|dyadid), family=binomial, data=disp.relevant2,  nAGQ = 12,
                       control=glmerControl(optimizer = "bobyqa", optCtrl=list(maxfun=500000)))
cre.margins <- summary(margins(cre1, variables=names(m1$coef),  data=as.data.frame(disp)))
cre.cml.margins <- summary(margins(cre1, variables=names(m1$coef), data=as.data.frame(disp2)))
rcre.margins <- summary(margins(rcml1, variables=names(m1$coef),  data=as.data.frame(disp2)))
cre.relevant.margins <- summary(margins(cre1, variables=names(m1$coef),  data=as.data.frame(disp.relevant)))
rcre.relevant.margins <- summary(margins(rcml.relevant, variables=names(m1$coef),  data=as.data.frame(disp.relevant)))
cre.relevant.margins2 <- summary(margins(cre1, variables=names(m1$coef),  data=as.data.frame(disp.relevant2)))
rcre.relevant.margins2 <- summary(margins(rcml.relevant2, variables=names(m1$coef),  data=as.data.frame(disp.relevant2)))



mldv <- glm( dispute~0+contig+ lcaprat+ mingrow+ allied +mindem+ mineco1+factor(dyadid),data=disp2, family=binomial, x=T, y=T)
mldv.margins <- summary(margins(mldv, variables=names(mldv$coef[1:6])))

##constants
Z1 <- unique(cbind(1, with(disp, cbind(dyadid, contig.bar, lcaprat.bar, mingrow.bar, allied.bar, mindem.bar, mineco1.bar))))[,-2]
alpha1 <- data.table(dyadid=unique(disp$dyadid), constant=rowSums(Z1*coef(cre1)$dyadid[,c(1, 8:13)]))
Z1 <- unique(cbind(1, with(disp2, cbind(dyadid, contig.bar, lcaprat.bar, mingrow.bar, allied.bar, mindem.bar, mineco1.bar))))[,-2]
ralpha1 <- data.table(dyadid=unique(disp2$dyadid), constant=rowSums(Z1*coef(rcml1)$dyadid[,c(1, 8:13)]))
mldv.alpha <-data.table(dyadid=unique(disp2$dyadid), constant= mldv$coef[-c(1:6)])
Z1 <- unique(cbind(1, with(disp.relevant, cbind(dyadid, contig.bar, lcaprat.bar, mingrow.bar, allied.bar, mindem.bar, mineco1.bar))))[,-2]
alpha.relevant <- data.table(dyadid=unique(disp.relevant$dyadid), constant=rowSums(Z1*coef(rcml.relevant)$dyadid[,c(1, 8:13)]))
Z1 <- unique(cbind(1, with(disp.relevant2, cbind(dyadid, contig.bar, lcaprat.bar, mingrow.bar, allied.bar, mindem.bar, mineco1.bar))))[,-2]
alpha.relevant2 <- data.table(dyadid=unique(disp.relevant2$dyadid), constant=rowSums(Z1*coef(rcml.relevant2)$dyadid[,c(1, 8:13)]))

save.image("dirtypoolout.rdata")



#Major powers: China, russia, france, usa, united kingdom


#time without compile
# X1 <- with(disp, cbind(contig, lcaprat, mingrow, allied ,mindem, mineco1))
# Xbar1 <- with(disp, cbind(contig.bar, lcaprat.bar, mingrow.bar, allied.bar ,mindem.bar, mineco1.bar))
# set.seed(1)
# t3.stan <- system.time({m3.stan <- CRE(y=disp$dispute, X=X1, Xbar=Xbar1, group.id = disp$dyadid, 
#                                        prior_var = 1, beta_var=2, psi_var=2,
#                                        chains=4,cores=4, seed=1,iter=5000, thin=5)})
# print(m3.stan$info$Rhat)
# summary(m3.stan$FE$Rhat)
# summary(m3.stan$RE$Rhat)
# summary(m3.stan$info$MCMC)
# beta.cre <- m3.stan$FE$coefficients[1:6]
# se.cre <- sqrt(diag(m3.stan$FE$vcov)[1:6])
# z.cre <- beta.cre/se.cre
# p.cre <- 2*pnorm(abs(z.cre), lower=F)
# round(cbind(beta.cre,se.cre,z.cre,p.cre),3)
# 
# alpha.CRE.stan <- m3.stan$FE$coefficients[-c(1:6)]



####MLDV####
## Sparsity matters from the start
# N <- length(unique(disp$dyadid))
# set.seed(1)
# theta <- c(m1$coef, rnorm(N)*.01) #start with consistent beta
# t2 <- system.time({xout <- sp.ffe(f1, data=disp, start=theta, optim.control=list(trace=TRUE, maxit=200))})
# print(xout$coef.table[1:6,])
# 
# f1 <- dispute~0+contig+ lcaprat+ mingrow+ allied +mindem+ mineco1+factor(dyadid)

# set.seed(1)
# ffe.boot <- mvrnorm(50, mu=xout$coefficients, Sigma=xout$vcov)
# m3.stan.boot <- do.call(rbind, m3.stan$FE$MCMC)
# mldv.boot <- mvrnorm(50, mldv$coefficients, Sigma=vcov(mldv))
# 
# X <- sparse.model.matrix( ~ 0 + contig + lcaprat + mingrow + allied + mindem + 
#                             mineco1 + factor(dyadid), data=disp)
# X2 <- sparse.model.matrix( ~ 0 + contig + lcaprat + mingrow + allied + mindem + 
#                              mineco1 + factor(dyadid), data=disp2)
# 
# 
# 
# mldv.ame.ci <- sapply(1:6, function(x){quantile(rowMeans(t(dlogis(as.matrix(X2 %*% t(mldv.boot))))*mldv.boot[,x]), c(0.025, .5,.975))})
# m3.stan.ame.ci <- sapply(1:6, function(x){quantile(rowMeans(t(dlogis(as.matrix(X %*% t(m3.stan.boot))))*m3.stan.boot[,x]), c(0.025, .5,.975))})
# ffe.ame.ci <-sapply(1:6, function(x){quantile(rowMeans(t(dlogis(as.matrix(X %*% t(ffe.boot))))*ffe.boot[,x]), c(0.025, .5,.975))})
# colnames(mldv.ame.ci) <- colnames(m3.stan.ame.ci) <- colnames(ffe.ame.ci) <- colnames(ffe.boot)[1:6]
# 
# 
# ALPHA <- cbind(xout$coeff[-c(1:6)], alpha.CRE.stan)
# colnames(ALPHA ) <- c("PML", "CRE")
# pairs(ALPHA)
# cor(ALPHA)
# 
# ALPHA2 <- ALPHA[(unique(disp$dyadid) %in% unique(disp2$dyadid)),]
# pairs(ALPHA2)
# cor(ALPHA2)
# 
# 
# 
#
# cre.fitted <- plogis(as.matrix(xout$X %*% m3.stan$FE$coefficients))
# beta.cre <- as.numeric(read.table("../stataFiles/dirtypool_fullSample_est.txt", header = TRUE)[1,])[1:6]
# beta.cre.restricted <- as.numeric(read.table("../stataFiles/dirtypool_CMLSample_est.txt", header = TRUE)[1,])[1:6]
# names(beta.cre.restricted) <- c("allied", "contig","mingrow","mindem","lcaprat", "mineco1")
# names(beta.cre) <- c("allied", "contig","mingrow","mindem","lcaprat", "mineco1")
# 
# 
# vcov.cre <- as.matrix(read.table("../stataFiles/dirtypool_fullSample_vcov.txt", header = TRUE))[1:6, 1:6]
# vcov.cre.restricted <- as.matrix(read.table("../stataFiles/dirtypool_CMLSample_vcov.txt", header = TRUE))[1:6, 1:6]
# rownames(vcov.cre) <- colnames(vcov.cre) <- names(beta.cre)
# rownames(vcov.cre.restricted) <- colnames(vcov.cre.restricted) <- names(beta.cre.restricted)
# se.cre <- sqrt(diag(vcov.cre))
# z.cre <- beta.cre/se.cre
# p.cre <- 2*pnorm(abs(z.cre), lower=F)
# round(cbind(beta.cre,se.cre,z.cre,p.cre),3)
# 
# se.cre.restricted <- sqrt(diag(vcov.cre.restricted))
# z.cre.restricted <- beta.cre.restricted/se.cre.restricted
# p.cre.restricted <- 2*pnorm(abs(z.cre.restricted), lower=F)
# round(cbind(beta.cre.restricted,se.cre.restricted,z.cre.restricted,p.cre.restricted),3)
# 
# 
# 
# 
# 
# margins.full.est <- as.numeric(read.table("../stataFiles/dirtypool_fullSample_margin.txt", header = TRUE)[1,])
# margins.full.se <-sqrt(diag(as.matrix(read.table("../stataFiles/dirtypool_fullSample_margin_vcov.txt", header = TRUE))))
# margins.full.z <- margins.full.est/margins.full.se
# margins.full.p <- pnorm(abs(margins.full.z), lower=F)*2
# margins.full.low <- margins.full.est-qnorm(.975)*margins.full.se
# margins.full.high <- margins.full.est+qnorm(.975)*margins.full.se
# cre.full.margins <- data.frame(factor=c("allied", "contig","mingrow","mindem","lcaprat", "mineco1") ,
#                                AME=margins.full.est, 
#                                SE=margins.full.se, 
#                                z=margins.full.z,
#                                p=margins.full.p,
#                                lower=margins.full.low,
#                                upper=margins.full.high)
# print(cre.full.margins)
# 
# margins.cml.est <- as.numeric(read.table("../stataFiles/dirtypool_CMLSample_margin.txt", header = TRUE)[1,])
# margins.cml.se <-sqrt(diag(as.matrix(read.table("../stataFiles/dirtypool_CMLSample_margin_vcov.txt", header = TRUE))))
# margins.cml.z <- margins.cml.est/margins.cml.se
# margins.cml.p <- pnorm(abs(margins.cml.z), lower=F)*2
# margins.cml.low <- margins.cml.est-qnorm(.975)*margins.cml.se
# margins.cml.high <- margins.cml.est+qnorm(.975)*margins.cml.se
# cre.cml.margins <- data.frame(factor=c("allied", "contig","mingrow","mindem","lcaprat", "mineco1") ,
#                               AME=margins.cml.est, 
#                               SE=margins.cml.se, 
#                               z=margins.cml.z,
#                               p=margins.cml.p,
#                               lower=margins.cml.low,
#                               upper=margins.cml.high)
# print(cre.cml.margins)
# 
# margins.full2.est <- as.numeric(read.table("../stataFiles/dirtypool_fullSample_margin2.txt", header = TRUE)[1,])
# margins.full2.se <-sqrt(diag(as.matrix(read.table("../stataFiles/dirtypool_fullSample_margin_vcov2.txt", header = TRUE))))
# margins.full2.z <- margins.full2.est/margins.full2.se
# margins.full2.p <- pnorm(abs(margins.full2.z), lower=F)*2
# margins.full2.low <- margins.full2.est-qnorm(.975)*margins.full2.se
# margins.full2.high <- margins.full2.est+qnorm(.975)*margins.full2.se
# cre.full2.margins <- data.frame(factor=c("allied", "contig","mingrow","mindem","lcaprat", "mineco1") ,
#                                 AME=margins.full2.est, 
#                                 SE=margins.full2.se, 
#                                 z=margins.full2.z,
#                                 p=margins.full2.p,
#                                 lower=margins.full2.low,
#                                 upper=margins.full2.high)
# print(cre.full2.margins)


rm(list=ls())
library(ggplot2)
library(MASS)
library(lme4)
library(survival)
library(stargazer)
load("dirtypoolout.rdata")



# p1 <- as.data.table(rbind.data.frame(summary(mldv.margins), cre.full.margins,cre.cml.margins,cre.full2.margins))
p1 <- as.data.table(rbind.data.frame(mldv.margins, cre.margins,
                                     cre.cml.margins, rcre.margins,
                                     cre.relevant.margins, rcre.relevant.margins))
# colnames(p1) <- c("lo","Estimate", "hi")
p1$variable <- factor(p1$factor,
                      levels=c("mindem", "mineco1","contig","lcaprat","mingrow","allied"),
                      labels=c("Demo.", "Trade", "Contiguity", "Cap. Rat. (log)", "Growth", "Ally"))
p1$Estimator <- factor(rep(c("MLDV", "CRE-AME",  "CRE-cAME", "rCRE",  "CRE-PR", "pAME"), each=6),
                       levels=c("CRE-AME", "CRE-PR","pCRE", "CRE-cAME", "rCRE", "MLDV"))
p1$Estimator <- factor(rep(c("MLDV", "CRE",  "CRE", "rCRE",  "CRE", "pCRE"), each=6),
                       levels=c("CRE","pCRE", "rCRE", "MLDV"))
p1$Effect <- c(rep("conditional AME",  6),
               rep("full-sample AME", 6),
               rep("conditional AME", 12),
               rep("politically relevant AME", 12))
p1$Effect <- factor(p1$Effect, levels=c("full-sample AME", "politically relevant AME", "conditional AME"))
g2 <- ggplot(p1)+
  geom_pointrange(aes(x=variable, y=AME,ymin=lower,ymax=upper, color=Estimator, shape=Estimator, linetype=Estimator), 
                  position=position_dodge(width=.75), size=.8)+
  ylab("Average marginal effects")+
  xlab("")+
  facet_wrap(~Effect,ncol=1)+
  geom_abline(aes(intercept=0, slope=0), linetype="dashed")+
  theme_bw(14)+
  theme(text=element_text(family="CM Roman"),
        legend.position = "bottom")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, 5.5, 1),alpha=.3)

print(g2)

ggsave(g2, file="/home/cox/Dropbox/xtlogit/figures/dirtypoolmem.pdf", height=10, width=8)
embed_fonts(file="/home/cox/Dropbox/xtlogit/figures/dirtypoolmem.pdf")

library(stargazer)

stargazer(m1,mldv, mldv, mldv, mldv,
          # covariate.labels = 
          omit=c("factor*", "*.bar", "Constant"),
          omit.stat=c("all"), 
          title="Replicating ",
          label="tab.dirtypool",
          star.cutoffs = c(NA),
          digits=2,
          no.space=TRUE,
          coef=list(NULL, NULL, 
                    summary(rcml1)$coef[names(m1$coef),1],
                    summary(rcml.relevant)$coef[names(m1$coef),1],
                    summary(cre1)$coef[names(m1$coef),1]),
          se=list(NULL, NULL, 
                  summary(rcml1)$coef[names(m1$coef),2],
                  summary(rcml.relevant)$coef[names(m1$coef),2],
                  summary(cre1)$coef[names(m1$coef),2]))


#FOR GKY2
p1 <- as.data.table(rbind.data.frame( cre.relevant.margins2, rcre.relevant.margins2))
# colnames(p1) <- c("lo","Estimate", "hi")
p1$variable <- factor(p1$factor,
                      levels=c("mindem", "mineco1","contig","lcaprat","mingrow","allied"),
                      labels=c("Demo.", "Trade", "Contiguity", "Cap. Rat. (log)", "Growth", "Ally"))
p1$Estimator <- factor(rep(c("CRE", "pCRE"), each=6),
                       levels=c("CRE","pCRE"))
# p1$Effect <- c(rep("conditional AME",  6),
#                rep("full-sample AME", 6),
#                rep("conditional AME", 12),
#                rep("politically relevant AME", 12))
# p1$Effect <- factor(p1$Effect, levels=c("full-sample AME", "politically relevant AME", "conditional AME"))
g2 <- ggplot(p1)+
  geom_pointrange(aes(x=variable, y=AME,ymin=lower,ymax=upper, color=Estimator, shape=Estimator, linetype=Estimator), 
                  position=position_dodge(width=.75), size=.8)+
  ylab("Average marginal effects")+
  xlab("")+
  # facet_wrap(~Effect,ncol=1)+
  geom_abline(aes(intercept=0, slope=0), linetype="dashed")+
  theme_bw(14)+
  theme(text=element_text(family="CM Roman"),
        legend.position = "bottom")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, 5.5, 1),alpha=.3)

print(g2)

ggsave(g2, file="/home/cox/Dropbox/xtlogit/figures/dirtypoolmem2.pdf", height=6, width=8)
embed_fonts(file="/home/cox/Dropbox/xtlogit/figures/dirtypoolmem2.pdf")

library(stargazer)

stargazer(m1, m1,
          # covariate.labels = 
          omit=c("factor*", "*.bar", "Constant"),
          omit.stat=c("all"), 
          title="Replicating ",
          label="tab.dirtypool",
          star.cutoffs = c(NA),
          digits=2,
          no.space=TRUE,
          coef=list(summary(rcml.relevant2)$coef[names(m1$coef),1],
                    summary(cre1)$coef[names(m1$coef),1]),
          se=list(summary(rcml.relevant2)$coef[names(m1$coef),2],
                  summary(cre1)$coef[names(m1$coef),2]))


