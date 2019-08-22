library(data.table)

library(survival)
library(lme4)
library(margins)

library(ggplot2)
library(xtable)
library(stargazer)
rm(list=ls())

#This file is not on dataverse per the terms and conditions
#of the NAES
#email the author for this file: ccrismancox@gmail.com
load("goldmanData.rdata")


#Recode the data for "first differences" using his Stata file has a guide
waves <- merge(wave1, wave3, by="RKEY")

waves <- waves[firstdemHC_1 != firstdemHC_3]
waves[, FFWperceptions := FFWperceptionsZ1_3-FFWperceptionsZ1_1]
waves[, perceivedFFB := perceivedFFBZ1_3-perceivedFFBZ1_1]
waves[, pidstrength := pidstrengthZ1_3-pidstrengthZ1_1]

waves[, campcontact := campcontactZ1_3-campcontactZ1_1]
waves[, polinterest := polinterest_3-polinterest_1]
waves[, viableDem := viableDem_3-viableDem_1]
waves[, electableDem := electableDem_3-electableDem_1]
waves[, HCBOJEissue6 := HCBOJEissue6_3-HCBOJEissue6_1]
waves[, ideoproxHCBOJE := ideoproxHCBOJE_3-ideoproxHCBOJE_1]

setnames(waves, 
         c("RKEY",
           "firstdemHC_3", "FFWperceptions", "perceivedFFB", "pidstrength",
           "black_1", "educZ1_1", "incomeZ1_1", "ageZ1_1", "oldsouth_1", "campcontact",
           "polinterest", "viableDem", "HCBOJEissue6","ideoproxHCBOJE", "male_1"),
         c("RKEY" , "firstdemHC", "FFWpercep", "FFBpercep", "pidstrength",
           "black", "educ", "income", "age", "south", "contact", "interest", "viable", "issues", "ideo", "male"))

#Replicate Goldman's Table 5 (sanity check)
model5.male <- glm(firstdemHC~FFWpercep + FFWperceptionsZ1_1 + FFBpercep + perceivedFFBZ1_1 + 
                     pidstrength + pidstrengthZ1_1 + black + educ + income + 
                     age + south + contact + campcontactZ1_1 + interest + 
                     polinterest_1 + viable + viableDem_1 + issues + HCBOJEissue6_1 +
                     ideo + ideoproxHCBOJE_1, x=T,
                   subset =(male==1 ), data=waves, family=binomial)

# rebuild build the data for two period fixed effects
wave1.subset <- wave1[,list(RKEY, firstdemHC_1, FFWperceptionsZ1_1,perceivedFFBZ1_1,pidstrengthZ1_1, black_1,
                            educZ1_1, incomeZ1_1, ageZ1_1,oldsouth_1, campcontactZ1_1,polinterest_1,
                            viableDem_1,HCBOJEissue6_1,ideoproxHCBOJE_1,male_1)]
wave3.subset <- wave3[,list(RKEY, firstdemHC_3, FFWperceptionsZ1_3,perceivedFFBZ1_3,pidstrengthZ1_3, black_3,
                            educZ1_3, incomeZ1_3, ageZ1_3,oldsouth_3, campcontactZ1_3,polinterest_3,
                            viableDem_3,HCBOJEissue6_3,ideoproxHCBOJE_3, male_3)]

setnames(wave1.subset, c("RKEY", "firstdemHC", "FFWpercep", "FFBpercep", "pidstrength", "black", "educ", "income", "age", "south",
                         "contact", "interest", "viable", "issues", "ideo", "male"))
setnames(wave3.subset, c("RKEY", "firstdemHC", "FFWpercep", "FFBpercep", "pidstrength", "black", "educ", "income", "age", "south",
                         "contact", "interest", "viable", "issues", "ideo", "male"))
wave1.subset[,wave:=1]
wave3.subset[,wave:=3]
wave3.subset[,age:=ifelse(is.na(age), NA, wave1.subset$age)]
wave13 <- rbind(wave1.subset, wave3.subset)
wave13 <- na.omit(wave13)
wave13[,nr:= length(firstdemHC), by=RKEY]
wave13 <- wave13[nr==2]
wave13[, mean.y := mean(firstdemHC, na.rm = T), by="RKEY"]
wave13.cml <- wave13[mean.y < 1 & mean.y >0]

#CML
c1 <- clogit(firstdemHC ~ FFWpercep + FFBpercep + pidstrength +
               contact + interest + viable +  issues + ideo +
               strata(RKEY),  subset =(male==1),
             data=wave13.cml)

#MLDV
glm1 <- glm(firstdemHC ~ FFWpercep + FFBpercep + pidstrength +
              contact + interest + viable +  issues + ideo +
              factor(RKEY)+0,  subset =(male==1), family=binomial, x=T,
            data=wave13.cml)

###CRE###
wave13[,`:=`( FFWpercep.bar = mean(FFWpercep, na.rm=T),
              FFBpercep.bar = mean(FFBpercep, na.rm=T),
              pidstrength.bar = mean(pidstrength, na.rm=T),
              contact.bar = mean(contact, na.rm=T),
              interest.bar = mean(interest, na.rm=T),
              viable.bar = mean(viable, na.rm=T),
              issues.bar = mean(issues, na.rm=T),
              ideo.bar = mean(ideo, na.rm=T)),
       by=RKEY]
###rCRE###
wave13.cml[,`:=`( FFWpercep.bar = mean(FFWpercep, na.rm=T),
              FFBpercep.bar = mean(FFBpercep, na.rm=T),
              pidstrength.bar = mean(pidstrength, na.rm=T),
              contact.bar = mean(contact, na.rm=T),
              interest.bar = mean(interest, na.rm=T),
              viable.bar = mean(viable, na.rm=T),
              issues.bar = mean(issues, na.rm=T),
              ideo.bar = mean(ideo, na.rm=T)),
       by=RKEY]


cre1 <- glmer(firstdemHC ~ FFWpercep + FFBpercep + pidstrength +
                contact + interest + viable +  issues + ideo+
                FFWpercep.bar + FFBpercep.bar + pidstrength.bar +
                contact.bar + interest.bar + viable.bar +  issues.bar + ideo.bar+
                (1|RKEY), family=binomial, data=wave13, subset=(male==1), nAGQ = 12,
              control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e5)))
rcre1 <- glmer(firstdemHC ~ FFWpercep + FFBpercep + pidstrength +
                contact + interest + viable +  issues + ideo+
                FFWpercep.bar + FFBpercep.bar + pidstrength.bar +
                contact.bar + interest.bar + viable.bar +  issues.bar + ideo.bar+
                (1|RKEY), family=binomial, data=wave13.cml, subset=(male==1), nAGQ = 12,
              control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=5e5)))

#####Marginals####
mldv.margins <- summary(margins(glm1, variables=names(c1$coefficients)))
cre.margins <- summary(margins(cre1, variables=names(c1$coefficients),  data=as.data.frame(wave13[male==1])))
cre.cml.margins <- summary(margins(cre1, variables=names(c1$coefficients), data=as.data.frame(wave13.cml[male==1])))
rcre.margins <- summary(margins(rcre1, variables=names(c1$coefficients),  data=as.data.frame(wave13.cml[male==1])))

male.out <- rbind.data.frame(mldv.margins, rcre.margins, cre.cml.margins, cre.margins)[, c("factor", "AME", "lower", "upper")] 
male.out <- male.out[male.out$factor=="FFWpercep",][,-1]
colnames(male.out) <- c("est", "lo", "hi")
male.out$Estimator <- c( "MLDV","rCRE", "CRE-cAME", "CRE-AME")



### Report the results ### 
stargazer(c1,glm1,glm1, glm1,
          omit=c("factor*", "*.bar", "Constant"),
          omit.stat=c("all"), 
          title="Replicating Goldsmith",
          label="tab.gold",
          star.cutoffs = c(NA),
          digits=2,
          no.space=TRUE,
          coef=list(NULL, NULL,
                    summary(rcre1)$coef[names(c1$coef),1],
                    summary(cre1)$coef[names(c1$coef),1]),
          se=list(NULL, NULL,
                  summary(rcre1)$coef[names(c1$coef),2],
                  summary(cre1)$coef[names(c1$coef),2]))



p1 <- as.data.table( rbind.data.frame(mldv.margins, rcre.margins, cre.cml.margins, cre.margins))
p1$variable <- factor(p1$factor,
                      levels=names(c1$coefficients)[8:1],
                      labels=c("Gender Fav.", "Racial Fav.", "Party ID", "Contact",
                               "Interest", "Viability", "Issue Agreement", "Ideo. Agreement")[8:1])
p1$Estimator <- factor(rep(c("MLDV",  "rCRE", "CRE-cAME", "CRE-AME"), each=8),
                       levels=c("CRE-AME", "CRE-PR","pCRE", "CRE-cAME", "rCRE", "MLDV"))

g2 <- ggplot(p1)+
  geom_pointrange(aes(x=variable, y=AME,ymin=lower,ymax=upper, color=Estimator, shape=Estimator, linetype=Estimator), 
                  position=position_dodge(width=.9), size=.8)+
  ylab("Average marginal effects")+
  xlab("")+
  geom_abline(aes(intercept=0, slope=0), linetype="dashed")+
  theme_bw(16)+
  theme(legend.position = "bottom")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, 7.5, 1),alpha=.3)
print(g2)

ggsave(g2, file="Figure6.pdf", height=9, width=8)
