library(readstata13)
library(data.table)

library(survival)
library(lme4)
library(margins)
library(MASS)

library(ggplot2)
library(xtable)
library(stargazer)
rm(list=ls())



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

h=5; w=6


margins.out$Estimator <- factor(margins.out$Estimator,
                         levels=c("CRE-AME", "CRE-cAME", "rCRE", "MLDV"))
full.plot <- ggplot(margins.out)+
  geom_pointrange(aes(x=x, y=Estimate, ymin=lo,ymax=hi,color=Estimator, shape=Estimator, linetype=Estimator), position=position_dodge(width=.5), size=1)+
  geom_abline(aes(intercept=0, slope=0), linetype="dashed")+
  xlab("")+
  ylab("Marginal effect on\nthe probability of protest")+
  theme_bw(16)+
  theme(legend.position = "bottom")

print(full.plot)


ggsave(full.plot, file="Figure5.pdf", height=1.25*h, width=1.5*w)

print(stargazer(m4,m4.glm2,m4.glm2, m4.glm2,
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
                  summary(mle.cre.out)$coefficients[names(m4$coef), 2])))

##### From EMW's code for marginal effects#####
set.seed(1)
m4.boot <- 	mvrnorm(1000, m4$coef, vcov(m4)) 
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
  theme(legend.position = "bottom")

print(emw.plot1)


ggsave(emw.plot1, file="FigureJ1.pdf", height=h, width=w)
