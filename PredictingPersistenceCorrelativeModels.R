setwd("~/PostDoc Computer/SpilloverCharacteristics/PLOS Biology/Code")

require("boot")
require("lme4")
require("zoo")
library(dplyr)
library(tidyr)
library(MASS)
library(AICcmodavg)
library(ggplot2)
library(performance)

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

###############################################
#Transmission models

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


library(ggeffects)
library(scales)
hex <- hue_pal()(9)
data$strain.overall <- factor(data$strain.overall, levels = c("ZF1092","JU1428","VX80","JU1857","JU1873","JU724","JU4093","SB454"))

F3a<-ggplot(data=data, aes(best.p, trans.ability))+
  geom_jitter(data=data, aes(best.p, trans.ability, color=strain.overall, shape=block), size=2, alpha=0.8, height = 0.15)+
  ylab("Passages until loss")+xlab("Maximum likelihood \nprevalence (p)")+theme_bw()+
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits=c(0,10))+
  theme(axis.text = element_text(color="black"))+
  scale_color_manual(name="Strain", values=hex)+
  geom_line(data=df.prev.predict, aes(best.p, trans.ability))+
  geom_line(data=df.prev.predict, aes(best.p, trans.ability.h), lty="dashed")+
  geom_line(data=df.prev.predict, aes(best.p, trans.ability.l), lty="dashed")+
  scale_shape_manual(values=c(16,15,18,17))
  

#Intensity
intensity.model <- glmer.nb(trans.ability ~ median.ct + (1|Strain)+(1|block), data = data.small)
summary(intensity.model)
#failed to converge
drop1(intensity.model, test="Chisq")
summary(intensity.model)
r2_nakagawa(intensity.model)

median.ct=seq(18,37,.01)    #Set range for plotting model predictions
trans.ability<-exp((median.ct)*fixef(intensity.model)[2]+fixef(intensity.model)[1])
trans.ability.l<-exp((median.ct)*(fixef(intensity.model)[2]-coef(summary(intensity.model))[4])+fixef(intensity.model)[1]-coef(summary(intensity.model))[3])
trans.ability.h<-exp((median.ct)*(fixef(intensity.model)[2]+coef(summary(intensity.model))[4])+fixef(intensity.model)[1]+coef(summary(intensity.model))[3])
df.intensity.predict<-as.data.frame(cbind(median.ct, trans.ability, trans.ability.l, trans.ability.h))

data.small$strain.overall <- factor(data.small$strain.overall, levels = c("ZF1092","JU1428","VX80","JU1857","JU1873","JU724","JU4093","SB454"))

F3c<-ggplot(data=data.small, aes(median.ct, trans.ability))+
  ylab("Passages until loss")+
  xlab(expression(atop("Infection intensity",paste("(corrected median Ct)"~(mu)))))+
  theme_bw()+
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits=c(0,10))+
  theme(axis.text = element_text(color="black"))+scale_color_manual(name="Strain", values=hex)+
  geom_jitter(data=data.small, aes(median.ct, trans.ability, color=strain.overall, shape=block), size=2, alpha=0.8, height=0.15)+
  scale_x_reverse()+
  geom_line(data=df.intensity.predict, aes(median.ct, trans.ability), color="gray")+
  geom_line(data=df.intensity.predict, aes(median.ct, trans.ability.h), lty="dashed", color="gray")+
  geom_line(data=df.intensity.predict, aes(median.ct, trans.ability.l), lty="dashed", color="gray")+
  scale_shape_manual(values=c(16,15,18,17))


#Shedding
shedding.model <- glmer.nb(trans.ability ~ best.y + (1|block)+(1|Strain), data = data)
summary(shedding.model)
#singular
drop1(shedding.model, test="Chisq")
r2_nakagawa(shedding.model)

#singular, P<0.001
best.y=seq(0.001,0.999,.001)    #Set range for plotting model predictions
trans.ability<-exp((best.y)*fixef(shedding.model)[2]+fixef(shedding.model)[1])
trans.ability.l<-exp((best.y)*(fixef(shedding.model)[2]-coef(summary(shedding.model))[4])+fixef(shedding.model)[1]-coef(summary(shedding.model))[3])
trans.ability.h<-exp((best.y)*(fixef(shedding.model)[2]+coef(summary(shedding.model))[4])+fixef(shedding.model)[1]+coef(summary(shedding.model))[3])
df.shedding.predict<-as.data.frame(cbind(best.y, trans.ability, trans.ability.l, trans.ability.h))

F3b<-ggplot(data=data, aes(best.y, trans.ability))+
  geom_jitter(data=data, aes(best.y, trans.ability, color=strain.overall, shape=block),size=2, alpha=0.8, height=0.15)+
  ylab("Passages until loss")+xlab("Maximum likelihood \nshedding ability (y)")+theme_bw()+
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits=c(0,10))+
  theme(axis.text = element_text(color="black"))+
  geom_line(data=df.shedding.predict, aes(best.y, trans.ability))+
  geom_line(data=df.shedding.predict, aes(best.y, trans.ability.h), lty="dashed")+
  geom_line(data=df.shedding.predict, aes(best.y, trans.ability.l), lty="dashed")+
  scale_shape_manual(values=c(16,15,18,17))

#rel.susc
rel.susc.model <- glmer.nb(trans.ability ~ rel.susc + (1|block)+(1|Strain), data = data)
summary(rel.susc.model)
drop1(rel.susc.model, test="Chisq")
#P=0.038
r2_nakagawa(rel.susc.model)

#failed to converge
rel.susc=seq(0.001,0.2,.002)    #Set range for plotting model predictions
trans.ability<-exp((rel.susc)*fixef(rel.susc.model)[2]+fixef(rel.susc.model)[1])
trans.ability.l<-exp((rel.susc)*(fixef(rel.susc.model)[2]-coef(summary(rel.susc.model))[4])+fixef(rel.susc.model)[1]-coef(summary(rel.susc.model))[3])
trans.ability.h<-exp((rel.susc)*(fixef(rel.susc.model)[2]+coef(summary(rel.susc.model))[4])+fixef(rel.susc.model)[1]+coef(summary(rel.susc.model))[3])
df.rel.susc.predict<-as.data.frame(cbind(rel.susc, trans.ability, trans.ability.l, trans.ability.h))

F3d<-ggplot(data=data, aes(rel.susc, trans.ability))+
  geom_jitter(data=data, aes(rel.susc, trans.ability,color=strain.overall, shape=block),size=2, alpha=0.8, height=0.15)+
  ylab("Passages until loss")+xlab("Relative susceptibility (s)")+theme_bw()+
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits = c(0,10))+
  theme(axis.text = element_text(color="black"))+scale_color_manual(name="Strain", values=hex)+
  geom_line(data=df.rel.susc.predict, aes(rel.susc, trans.ability))+
  geom_line(data=df.rel.susc.predict, aes(rel.susc, trans.ability.h), lty="dashed")+
  geom_line(data=df.rel.susc.predict, aes(rel.susc, trans.ability.l), lty="dashed")+
  scale_shape_manual(name="Block",values=c(16,15,18,17))

Correlations<-plot_grid(F3a+theme(legend.position = "none"), 
                        F3b+theme(legend.position = "none", axis.title.y = element_blank()),
                        F3c+theme(legend.position="none", axis.title.y = element_blank()),
                        F3d+theme(axis.title.y = element_blank()), ncol=4, rel_widths = c(1.1,1,1,1.5),
                        labels=c("A","B","C","D"), align="h", axis="b")
save_plot("Figure3.tiff", Correlations, base_height = 5, base_width = 10)

library(car)
vif(m1)
library("grid")
library("ggplotify")
library("corrplot")

#Make a correlations plot
intensity<-data.small$median.ct
prevalence<-data.small$best.p
y<-data.small$best.y
rel.susc<-data.small$rel.susc

Corplot<-cbind(intensity, prevalence, y, rel.susc)
Corplot<-as.data.frame(Corplot)
plot(Corplot)

M<-cor(Corplot, use="complete.obs")
corrplot(M, method="number")

cor.grob<-as.ggplot(~plot(Corplot))+geom_text(x=0.46, y=.87, label="-0.36", color="red")+
  geom_text(x=0.67, y=.87, label="-0.37", color="red")+
  geom_text(x=0.87, y=.87, label="-0.23", color="red")+
  geom_text(x=0.67, y=.67, label="0.63", color="red")+
  geom_text(x=0.87, y=.67, label="0.44", color="red")+
  geom_text(x=0.87, y=.47, label="0.27", color="red")
save_plot("SpillChar.CorrelationPlot.jpg", cor.grob, base_height = 6.1, base_width = 8)

