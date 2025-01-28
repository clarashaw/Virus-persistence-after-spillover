setwd("~/PostDoc Computer/SpilloverCharacteristics")

require("boot")
require("lme4")
require("zoo")
library(dplyr)
library(tidyr)
library(MASS)
library(AICcmodavg)
library(ggplot2)
########################################################################################################
data<-read.csv("spill.char.dataset1.csv")

#is prevalence different across strains and blocks?
prev.model<-glm(best.p~Strain+block, data=data, family="quasibinomial")
drop1(prev.model, test="F")

#is y different across strains and blocks?
y.model<-glm(best.y~Strain+block, data=data, family="quasibinomial")
drop1(y.model, test="F")

#is intensity different across strains and blocks?
int.model<-lm(median.ct~Strain+block, data=data)
drop1(int.model, test="F")
#Strains are different, but blocks are not different

###############################################
#Transmission models
library(performance)

data.small<-filter(data, !is.na(median.ct))
m1 <- glmer.nb(trans.ability ~ best.p + median.ct + best.y + rel.susc + (1|block)+(1|Strain), data = data.small)
#singularity warning
summary(m1)
drop1(m1, test="Chisq", na.rm=TRUE)
#best.p and best.y are significant
r2_nakagawa(m1)

#Single factor models: (use full dataset when possible)
prev.model <- glmer.nb(trans.ability ~ best.p + (1|Strain)+(1|block), data = data)
summary(prev.model)
drop1(prev.model, test="Chisq")

best.p=seq(0.001,0.999,.001)    #Set range for plotting model predictions
trans.ability<-exp((best.p)*fixef(prev.model)[2]+fixef(prev.model)[1])
trans.ability.l<-exp((best.p)*(fixef(prev.model)[2]-coef(summary(prev.model))[4])+fixef(prev.model)[1]-coef(summary(prev.model))[3])
trans.ability.h<-exp((best.p)*(fixef(prev.model)[2]+coef(summary(prev.model))[4])+fixef(prev.model)[1]+coef(summary(prev.model))[3])
df.prev.predict<-as.data.frame(cbind(best.p, trans.ability, trans.ability.l, trans.ability.h))
r2_nakagawa(prev.model)

#Intensity
intensity.model <- glmer.nb(trans.ability ~ median.ct + (1|Strain)+(1|block), data = data.small)
#failed to converge
drop1(intensity.model, test="Chisq")
summary(intensity.model)
r2_nakagawa(intensity.model)

median.ct=seq(18,37,.01)    #Set range for plotting model predictions
trans.ability<-exp((median.ct)*fixef(intensity.model)[2]+fixef(intensity.model)[1])
trans.ability.l<-exp((median.ct)*(fixef(intensity.model)[2]-coef(summary(intensity.model))[4])+fixef(intensity.model)[1]-coef(summary(intensity.model))[3])
trans.ability.h<-exp((median.ct)*(fixef(intensity.model)[2]+coef(summary(intensity.model))[4])+fixef(intensity.model)[1]+coef(summary(intensity.model))[3])
df.intensity.predict<-as.data.frame(cbind(median.ct, trans.ability, trans.ability.l, trans.ability.h))

#Shedding
shedding.model <- glmer.nb(trans.ability ~ best.y + (1|block)+(1|Strain), data = data)
#singular
drop1(shedding.model, test="Chisq")
r2_nakagawa(shedding.model)

#singular, P<0.001
best.y=seq(0.001,0.999,.001)    #Set range for plotting model predictions
trans.ability<-exp((best.y)*fixef(shedding.model)[2]+fixef(shedding.model)[1])
trans.ability.l<-exp((best.y)*(fixef(shedding.model)[2]-coef(summary(shedding.model))[4])+fixef(shedding.model)[1]-coef(summary(shedding.model))[3])
trans.ability.h<-exp((best.y)*(fixef(shedding.model)[2]+coef(summary(shedding.model))[4])+fixef(shedding.model)[1]+coef(summary(shedding.model))[3])
df.shedding.predict<-as.data.frame(cbind(best.y, trans.ability, trans.ability.l, trans.ability.h))

#rel.susc
rel.susc.model <- glmer.nb(trans.ability ~ rel.susc + (1|block)+(1|Strain), data = data)
drop1(rel.susc.model, test="Chisq")
#P=0.038
r2_nakagawa(rel.susc.model)

#failed to converge
rel.susc=seq(0.001,0.2,.002)    #Set range for plotting model predictions
trans.ability<-exp((rel.susc)*fixef(rel.susc.model)[2]+fixef(rel.susc.model)[1])
trans.ability.l<-exp((rel.susc)*(fixef(rel.susc.model)[2]-coef(summary(rel.susc.model))[4])+fixef(rel.susc.model)[1]-coef(summary(rel.susc.model))[3])
trans.ability.h<-exp((rel.susc)*(fixef(rel.susc.model)[2]+coef(summary(rel.susc.model))[4])+fixef(rel.susc.model)[1]+coef(summary(rel.susc.model))[3])
df.rel.susc.predict<-as.data.frame(cbind(rel.susc, trans.ability, trans.ability.l, trans.ability.h))

