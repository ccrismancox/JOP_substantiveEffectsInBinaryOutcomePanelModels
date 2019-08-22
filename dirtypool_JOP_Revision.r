library(foreign)
library(data.table)
library(stringr)

library(survival)
library(lme4)
library(margins)

library(ggplot2)
library(xtable)
library(stargazer)
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





###CML####
#clogit drops 2,877 groups due to no war
t1 <- system.time({m1 <- clogit(dispute ~contig+ lcaprat+ mingrow+ allied +mindem+ mineco1+strata(dyadid), data=disp2)})

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

cre.margins <- summary(margins(cre1, variables=names(m1$coef),  data=as.data.frame(disp)))
cre.cml.margins <- summary(margins(cre1, variables=names(m1$coef), data=as.data.frame(disp2)))
rcre.margins <- summary(margins(rcml1, variables=names(m1$coef),  data=as.data.frame(disp2)))
cre.relevant.margins <- summary(margins(cre1, variables=names(m1$coef),  data=as.data.frame(disp.relevant)))
rcre.relevant.margins <- summary(margins(rcml.relevant, variables=names(m1$coef),  data=as.data.frame(disp.relevant)))



mldv <- glm( dispute~0+contig+ lcaprat+ mingrow+ allied +mindem+ mineco1+factor(dyadid),data=disp2, family=binomial, x=T, y=T)
mldv.margins <- summary(margins(mldv, variables=names(mldv$coef[1:6])))





####Present Results####
p1 <- as.data.table(rbind.data.frame(mldv.margins, cre.margins,
                                     cre.cml.margins, rcre.margins,
                                     cre.relevant.margins, rcre.relevant.margins))

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
  theme(legend.position = "bottom")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, 5.5, 1),alpha=.3)

print(g2)

ggsave(g2, file="Figure7.pdf", height=10, width=8)


stargazer(m1,mldv, mldv, mldv, mldv,
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
