library(survival)
library(readstata13)
library(data.table)
library(MASS)
library(lmtest)
library(xtable)
library(lme4)
library(margins)
rm(list=ls())


source("sparseFirth.R")
source("xtlogit_stan.r")

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

#Replicate Goldman's Table 5
model5.male <- glm(firstdemHC~FFWpercep + FFWperceptionsZ1_1 + FFBpercep + perceivedFFBZ1_1 + 
                     pidstrength + pidstrengthZ1_1 + black + educ + income + 
                     age + south + contact + campcontactZ1_1 + interest + 
                     polinterest_1 + viable + viableDem_1 + issues + HCBOJEissue6_1 +
                     ideo + ideoproxHCBOJE_1, x=T,
                   subset =(male==1 ), data=waves, family=binomial)
#Replicate Goldman's idea correctly
model5.correct <- glm(firstdemHC~FFWpercep +  FFBpercep +  
                     pidstrength +   contact +  interest + 
                      viable +  issues +ideo+0, x=T,
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

c1 <- clogit(firstdemHC ~ FFWpercep + FFBpercep + pidstrength +
               contact + interest + viable +  issues + ideo +
               strata(RKEY),  subset =(male==1),
             data=wave13.cml)


glm1 <- glm(firstdemHC ~ FFWpercep + FFBpercep + pidstrength +
              contact + interest + viable +  issues + ideo +
              factor(RKEY)+0,  subset =(male==1), family=binomial, x=T,
            data=wave13.cml)

# sp1 <- sp.ffe(firstdemHC ~ FFWpercep + FFBpercep + pidstrength +
#                 contact + interest + viable +  issues + ideo +
#                 factor(RKEY)+0,  subset =(male==1),
#               data=wave13, optim.control=list(trace=T))

wave13[,`:=`( FFWpercep.bar = mean(FFWpercep, na.rm=T),
              FFBpercep.bar = mean(FFBpercep, na.rm=T),
              pidstrength.bar = mean(pidstrength, na.rm=T),
              contact.bar = mean(contact, na.rm=T),
              interest.bar = mean(interest, na.rm=T),
              viable.bar = mean(viable, na.rm=T),
              issues.bar = mean(issues, na.rm=T),
              ideo.bar = mean(ideo, na.rm=T)),
       by=RKEY]
wave13.cml[,`:=`( FFWpercep.bar = mean(FFWpercep, na.rm=T),
              FFBpercep.bar = mean(FFBpercep, na.rm=T),
              pidstrength.bar = mean(pidstrength, na.rm=T),
              contact.bar = mean(contact, na.rm=T),
              interest.bar = mean(interest, na.rm=T),
              viable.bar = mean(viable, na.rm=T),
              issues.bar = mean(issues, na.rm=T),
              ideo.bar = mean(ideo, na.rm=T)),
       by=RKEY]


#CRE
#Gets close but not close enough. Gradient still > 1e-5
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
goldman <- summary(margins(model5.male, variables=names(c1$coefficients)))
cml.margins <- summary(margins(model5.correct, variables=names(c1$coefficients)))
mldv.margins <- summary(margins(glm1, variables=names(c1$coefficients)))
cre.margins <- summary(margins(cre1, variables=names(c1$coefficients),  data=as.data.frame(wave13[male==1])))
cre.cml.margins <- summary(margins(cre1, variables=names(c1$coefficients), data=as.data.frame(wave13.cml[male==1])))
rcre.margins <- summary(margins(rcre1, variables=names(c1$coefficients),  data=as.data.frame(wave13.cml[male==1])))

save.image("genderOut.rdata")

rm(list=ls())
library(ggplot2)
library(MASS)
library(lme4)
library(survival)
library(xtable)
library(stargazer)
load("genderOut.rdata")
# summary(margins(cre1, variables="FFWpercep", data=as.data.frame(wave13[male==1])))
# y1 <- with(wave13[male==1], cbind(firstdemHC,RKEY))
# X1 <- with(wave13[male==1], cbind(FFWpercep, FFBpercep, pidstrength, contact, interest, viable,
#                                   issues, ideo))
# Xbar1 <- with(wave13[male==1], cbind(FFWpercep.bar, FFBpercep.bar, pidstrength.bar, contact.bar,
#                                      interest.bar, viable.bar, issues.bar, ideo.bar))
# 
# 
# XX <- na.omit(cbind(y1, X1,Xbar1))
# X1 <- XX[, colnames(X1)]
# Xbar1 <- XX[, colnames(Xbar1)]
# y1 <- XX[, colnames(y1)]
# group1 <- y1[,"RKEY"]
# y1 <- y1[,"firstdemHC"]
# 
# #Sanity checks
# length(unique(group1)) 
# dim(unique(cbind(group1, Xbar1)))
# 
# 
# model.male <- CRE(y=y1, X=X1, Xbar=Xbar1, group.id = group1, 
#                   chains=4, cores=4, seed=1234567, iter=1500, thin=5)
# print(model.male$info$Rhat)
# summary(model.male$FE$Rhat)
# summary(model.male$RE$Rhat)
# summary(model.male$info$MCMC)
# beta.cre <- model.male$FE$coefficients[1:8]
# se.cre <- sqrt(diag(model.male$FE$vcov)[1:8])
# z.cre <- beta.cre/se.cre
# p.cre <- 2*pnorm(abs(z.cre), lower=F)
# round(cbind(beta.cre,se.cre,z.cre,p.cre),3)

# 
# beta.cre <- as.numeric(read.table("../stataFiles/gender_fullSample_est.txt", header = TRUE)[1,])[1:8]
# names(beta.cre) <- c(names(glm1$coef)[1:8])
# vcov.cre <- as.matrix(read.table("../stataFiles/gender_fullSample_vcov.txt", header = TRUE))[1:8, 1:8]
# rownames(vcov.cre) <- colnames(vcov.cre) <- names(beta.cre)
# 
# 
# se.cre <- sqrt(diag(vcov.cre))
# z.cre <- beta.cre/se.cre
# p.cre <- 2*pnorm(abs(z.cre), lower=F)
# round(cbind(beta.cre,se.cre,z.cre,p.cre),3)
# 
# 
# 
# 
# 
# 
# beta.cre.restricted <- as.numeric(read.table("../stataFiles/gender_cmlSample_est.txt", header = TRUE)[1,])[1:8]
# names(beta.cre.restricted) <- c(names(glm1$coef)[1:8])
# vcov.cre.restricted <- as.matrix(read.table("../stataFiles/gender_cmlSample_vcov.txt", header = TRUE))[1:8, 1:8]
# rownames(vcov.cre.restricted) <- colnames(vcov.cre.restricted) <- names(beta.cre.restricted)
# 
# 
# se.cre.restricted <- sqrt(diag(vcov.cre.restricted))
# z.cre.restricted <- beta.cre.restricted/se.cre.restricted
# p.cre.restricted <- 2*pnorm(abs(z.cre.restricted), lower=F)
# round(cbind(beta.cre.restricted,se.cre.restricted,z.cre.restricted,p.cre.restricted),3)
# 
# set.seed(1)
# 
# 
# 
# margins.full.est <- as.numeric(read.table("../stataFiles/gender_fullSample_margin.txt", header = TRUE)[1,])
# margins.full.se <-sqrt(as.numeric(read.table("../stataFiles/gender_fullSample_margin_vcov.txt", header = TRUE)[1,]))
# margins.full.z <- margins.full.est/margins.full.se
# margins.full.p <- pnorm(abs(margins.full.z), lower=F)*2
# margins.full.low <- margins.full.est-qnorm(.975)*margins.full.se
# margins.full.high <- margins.full.est+qnorm(.975)*margins.full.se
# cre.full.margins <- data.frame(factor="FFWpercep", 
#                                AME=margins.full.est, 
#                                SE=margins.full.se, 
#                                z=margins.full.z,
#                                p=margins.full.p,
#                                lower=margins.full.low,
#                                upper=margins.full.high)
# print(cre.full.margins)
# 
# margins.cml.est <- as.numeric(read.table("../stataFiles/gender_cmlSample_margin.txt", header = TRUE)[1,])
# margins.cml.se <-sqrt(as.numeric(read.table("../stataFiles/gender_cmlSample_margin_vcov.txt", header = TRUE)[1,]))
# margins.cml.z <- margins.cml.est/margins.cml.se
# margins.cml.p <- pnorm(abs(margins.cml.z), lower=F)*2
# margins.cml.low <- margins.cml.est-qnorm(.975)*margins.cml.se
# margins.cml.high <- margins.cml.est+qnorm(.975)*margins.cml.se
# cre.cml.margins <- data.frame(factor="FFWpercep", 
#                               AME=margins.cml.est, 
#                               SE=margins.cml.se, 
#                               z=margins.cml.z,
#                               p=margins.cml.p,
#                               lower=margins.cml.low,
#                               upper=margins.cml.high)
# print(cre.cml.margins)
# 
# margins.full2.est <- as.numeric(read.table("../stataFiles/gender_fullSample_margin2.txt", header = TRUE)[1,])
# margins.full2.se <-sqrt(as.numeric(read.table("../stataFiles/gender_fullSample_margin_vcov2.txt", header = TRUE)[1,]))
# margins.full2.z <- margins.full2.est/margins.full2.se
# margins.full2.p <- pnorm(abs(margins.full2.z), lower=F)*2
# margins.full2.low <- margins.full2.est-qnorm(.975)*margins.full2.se
# margins.full2.high <- margins.full2.est+qnorm(.975)*margins.full2.se
# cre.full2.margins <- data.frame(factor="FFWpercep", 
#                                 AME=margins.full2.est, 
#                                 SE=margins.full2.se, 
#                                 z=margins.full2.z,
#                                 p=margins.full2.p,
#                                 lower=margins.full2.low,
#                                 upper=margins.full2.high)
# print(cre.full2.margins)
# 


male.out <- rbind.data.frame(mldv.margins, rcre.margins, cre.cml.margins, cre.margins)[, c("factor", "AME", "lower", "upper")] #what to do here?
male.out <- male.out[male.out$factor=="FFWpercep",][,-1]
colnames(male.out) <- c("est", "lo", "hi")
# male.out$sex <- "Male"
male.out$Estimator <- c( "MLDV","rCRE", "CRE-cAME", "CRE-AME")



### Report the results ### 


male.tab <- rbind(as.character(round(male.out$est, 2)),
                  paste("(", round(male.out$lo, 2),
                        ", ", round(male.out$hi, 2), ")", sep=""))
rownames(male.tab) <- c("Est.","95\\% CI")
colnames(male.tab) <- male.out$Estimator
print(male.tab) #table 3
print(xtable(male.tab,
             align="rcccc",
             label='tab.ame.gender',
             caption="The AME of perceived gender favoritism on the probability of supporting Clinton among male democrats."), 
      include.rownames=T, caption.placement="top",
      sanitize.text.function=function(x){x})

stargazer(c1,glm1,glm1, glm1,
          # covariate.labels = 
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
# colnames(p1) <- c("lo","Estimate", "hi")
p1$variable <- factor(p1$factor,
                      levels=names(c1$coefficients)[8:1],
                      labels=c("Gender Fav.", "Racial Fav.", "Party ID", "Contact",
                               "Interest", "Viability", "Issue Agreement", "Ideo. Agreement")[8:1])
p1$Estimator <- factor(rep(c("MLDV",  "rCRE", "CRE-cAME", "CRE-AME"), each=8),
                       levels=c("CRE-AME", "CRE-PR","pCRE", "CRE-cAME", "rCRE", "MLDV"))
# p1$Estimator <- factor(rep(c("MLDV", "CRE",  "CRE", "rCRE",  "CRE", "pCRE"), each=8),
#                        levels=c("CRE","pCRE", "rCRE", "MLDV"))
# p1$Effect <- c(rep("cAME",  6),
#                rep("AME", 6),
#                rep("cAME", 12),
#                rep("pAME", 12))
# p1$Effect <- factor(p1$Effect, levels=c("AME", "pAME", "cAME"))
g2 <- ggplot(p1)+
  geom_pointrange(aes(x=variable, y=AME,ymin=lower,ymax=upper, color=Estimator, shape=Estimator, linetype=Estimator), 
                  position=position_dodge(width=.9), size=.8)+
  ylab("Average marginal effects")+
  xlab("")+
  # facet_wrap(~Effect,ncol=1)+
  geom_abline(aes(intercept=0, slope=0), linetype="dashed")+
  theme_bw(16)+
  theme(text=element_text(family="CM Roman"),
        legend.position = "bottom")+
  coord_flip()+
  geom_vline(xintercept=seq(1.5, 7.5, 1),alpha=.3)
print(g2)

ggsave(g2, file="/home/cox/Dropbox/xtlogit/figures/genderAME.pdf", height=9, width=8)
embed_fonts(file="/home/cox/Dropbox/xtlogit/figures/genderAME.pdf")
