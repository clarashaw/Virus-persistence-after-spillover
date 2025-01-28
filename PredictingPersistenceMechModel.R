setwd("~/PostDoc Computer/SpilloverCharacteristics")

library("boot")
library("lme4")
library(dplyr)
library(tidyr)
library(AICcmodavg)
library(ggplot2)

############################################################################################################
#mechanistic model
#Function to calculate the probability of successful passage from the proportion infected, intensity of infection,
#susceptibility, and shedding

ProbOfPassage=function(p, y, Cj, avgi, sig)        
{
  TotalWorms=20
  ProbOfNoTrans=0
  for (i in 0:TotalWorms)
  {
    ProbOfNoTrans=ProbOfNoTrans + dbinom(i,TotalWorms,p)*dbinom(0,i,y)*((1-pnorm(Cj,avgi,sig))^y)
  }
  return(1-ProbOfNoTrans)
}


Data=read.csv("spill.char.dataset2.csv")    

##############################################################################################
#ct cutoff=25
Data25=read.csv("spill.char.dataset2.csv")    
Data25$infected=as.numeric(Data25$Ct<25)    #Recode the data so that successful passage requires a CT < 25

#but then we have to cut some data with the following loop (the updata dataset will be DATA.1):
DATA.1<-NULL
LINE<-unique(Data25$Line)
for (i in LINE){
  thisline<-filter(Data25, Line==i) #Look at 1 passage line at a time
  thislineinfected<-filter(thisline, infected==1) #These are the ones that are infected with the ct<25 criterion
  thislineuninfected<-filter(thisline, infected==0) #These are the ones that are uninfected
  thisline.updated<-rbind(thislineinfected, thislineuninfected[1,]) #Keep the infected ones and the first uninfected one
  DATA.1<-rbind(DATA.1, thisline.updated)
}

DATA.1<-as.data.frame(DATA.1)
DATA.1<-filter(DATA.1, !is.na(Line)) #clean things up - these come from the lines that don't have uninfected passages
Data25<-DATA.1 #rename

Passage.best = ProbOfPassage(Data25$best.p, Data25$best.y, a+log2(Data25$rel.susc), Data25$pre.ct, sig)      #Maximum likleihood estimate for the probability of successful passage
Weight=.8                                            #Factor to adjust passage estimates by to avoid 0% and 100% probability estimates.  Larger values weight the data more.
MechEstimate=Passage.best*Weight+ 0.5*(1-Weight)     #Slightly adjust the probability of passage to avoid 0% and 100%
LogitMechPassage= logit(MechEstimate)
Data25$LogitMechPassage<-LogitMechPassage
Data25$MechEstimate<-MechEstimate 

#Data25p1<-filter(Data25, Passage==1)

Model_0 =glm(infected ~ 0 + LogitMechPassage, family="binomial", data=Data25)    #Run null GLM model without random effects
#Model_0 =glm(infected ~ 0 + LogitMechPassage, family="binomial", data=Data25p1)    #Run null GLM model without random effects

Model_Null =glmer(infected ~ 0 + (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset = LogitMechPassage, family="binomial", data=Data25)    #Run null GLM model with random effects
Model_Null0 =glmer(infected ~ 0 + (1|Line) + (1|Strain) + (1|Passage) + (1|block), family="binomial", data=Data25)    #Run null GLM model with random effects

# Calculate R-squared
deviance <- summary(Model_0)$deviance

null_deviance <- summary(Model_0)$null.deviance
rsquared <- 1 - (deviance / null_deviance) #0.38

r<-1-(deviance1/devianceNull)
null_deviance <- summary(Model_Null)$null.deviance
rsquared <- 1 - (deviance / null_deviance) #0.38

Model_Null =glmer(infected ~ 0 + LogitMechPassage + (1|Line) + (1|Strain) + (1|Passage) + (1|block), family="binomial", data=Data25)    #Run null GLM model with random effects
#Model_Null =glmer(infected ~ 0 + LogitMechPassage + (1|Line) + (1|Strain) + (1|block), family="binomial", data=Data25p1)    #Just passage 1 Run null GLM model with random effects

summary(Model_Null) 

R2_LogitMechPassage<- var(Data25$LogitMechPassage*fixef(Model_Null)[1])/(var(Data25$LogitMechPassage*fixef(Model_Null)[1])+summary(Model_Null)$varcor$Line[1]+summary(Model_Null)$varcor$Strain[1]+summary(Model_Null)$varcor$Passage[1]+summary(Model_Null)$varcor$block[1]+(pi^2)/3)
R2_Line<- summary(Model_Null)$varcor$Line[1]/(var(Data25$LogitMechPassage*fixef(Model_Null)[1])+summary(Model_Null)$varcor$Line[1]+summary(Model_Null)$varcor$Strain[1]+summary(Model_Null)$varcor$Passage[1]+summary(Model_Null)$varcor$block[1]+(pi^2)/3)
R2_Strain<- summary(Model_Null)$varcor$Strain[1]/(var(Data25$LogitMechPassage*fixef(Model_Null)[1])+summary(Model_Null)$varcor$Line[1]+summary(Model_Null)$varcor$Strain[1]+summary(Model_Null)$varcor$Passage[1]+summary(Model_Null)$varcor$block[1]+(pi^2)/3)
R2_Passage<- summary(Model_Null)$varcor$Passage[1]/(var(Data25$LogitMechPassage*fixef(Model_Null)[1])+summary(Model_Null)$varcor$Line[1]+summary(Model_Null)$varcor$Strain[1]+summary(Model_Null)$varcor$Passage[1]+summary(Model_Null)$varcor$block[1]+(pi^2)/3)
R2_block<- summary(Model_Null)$varcor$block[1]/(var(Data25$LogitMechPassage*fixef(Model_Null)[1])+summary(Model_Null)$varcor$Line[1]+summary(Model_Null)$varcor$Strain[1]+summary(Model_Null)$varcor$Passage[1]+summary(Model_Null)$varcor$block[1]+(pi^2)/3)
R2_cond<- (var(Data25$LogitMechPassage*fixef(Model_Null)[1])+summary(Model_Null)$varcor$Line[1]+summary(Model_Null)$varcor$Strain[1]+summary(Model_Null)$varcor$Passage[1]+summary(Model_Null)$varcor$block[1])/(var(Data25$LogitMechPassage*fixef(Model_Null)[1])+summary(Model_Null)$varcor$Line[1]+summary(Model_Null)$varcor$Strain[1]+summary(Model_Null)$varcor$Passage[1]+summary(Model_Null)$varcor$block[1]+(pi^2)/3)

library(MuMIn)
r.squaredGLMM(Model_Null)
r.squaredGLMM(Model_0)
#This is for plotting the predicted lines on the graph
XRange=seq(0.001,0.999,.001)    #Set range for plotting model predictions
pred<-inv.logit((logit(XRange))*fixef(Model_Null)[1])
pred.l<-inv.logit((logit(XRange))*(fixef(Model_Null)[1]-coef(summary(Model_Null))[2]))
pred.h<-inv.logit((logit(XRange))*(fixef(Model_Null)[1]+coef(summary(Model_Null))[2]))
df.pred<-as.data.frame(cbind(XRange, pred, pred.l, pred.h))

#The points are the mech model estimate on the x and the actual result on the y 
p3<-ggplot(Data25, aes(MechEstimate, infected, color=factor(Passage)))+geom_jitter(height = 0.08, width=0.02)+
  theme_bw()+
  theme(axis.text = element_text(color="black"))+
  geom_line(data=df.pred, aes(XRange, pred), color="blue")+
  geom_line(data=df.pred, aes(XRange, pred.l), color="blue", lty="dashed")+
  geom_line(data=df.pred, aes(XRange, pred.h), color="blue", lty="dashed")+
  geom_abline(slope=1, intercept = 0, color="red")+
  scale_color_discrete(name="Passage")+
  xlab("Mechanistic model passage estimate")+
  ylab("Passage success")
save_plot("Figure3.tiff", p3, base_height = 4, base_width = 4.5)

############ Various versions of the model including and excluding spillover characteristics############
#Cutoff<25
#LogitMechPassage (from mechanistic model) is used as the offset... actually, here it's not because I took out the "offset=", but you could put that back in if you want.
#Data25$median.ct<-(40-Data25$median.ct)/40

OFFSET=TRUE
if (OFFSET==TRUE)
{
  Model1111 = glmer(infected ~ 0 + best.y + best.p + median.ct + rel.susc + (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset=LogitMechPassage, family="binomial", data=Data25)
  Model1110 = glmer(infected ~ 0 + best.y + best.p + median.ct +            (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset=LogitMechPassage, family="binomial", data=Data25)
  
  Model1101 = glmer(infected ~ 0 + best.y + best.p +             rel.susc + (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset=LogitMechPassage, family="binomial", data=Data25)
  Model1100 = glmer(infected ~ 0 + best.y + best.p +                        (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset=LogitMechPassage, family="binomial", data=Data25)
  
  Model1011 = glmer(infected ~ 0 + best.y +          median.ct + rel.susc + (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset=LogitMechPassage , family="binomial", data=Data25)
  Model1010 = glmer(infected ~ 0 + best.y +          median.ct +            (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset=LogitMechPassage, family="binomial", data=Data25)
  
  Model0111 = glmer(infected ~ 0 +          best.p + median.ct + rel.susc + (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset=LogitMechPassage , family="binomial", data=Data25)
  Model0110 = glmer(infected ~ 0 +          best.p + median.ct +            (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset=LogitMechPassage , family="binomial", data=Data25)
  
  Model1001 = glmer(infected ~ 0 + best.y +                      rel.susc + (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset=LogitMechPassage , family="binomial", data=Data25)
  Model1000 = glmer(infected ~ 0 + best.y +                                 (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset=LogitMechPassage , family="binomial", data=Data25)
  
  Model0101 = glmer(infected ~ 0 +          best.p +             rel.susc + (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset=LogitMechPassage , family="binomial", data=Data25)
  Model0100 = glmer(infected ~ 0 +          best.p +                        (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset=LogitMechPassage , family="binomial", data=Data25)
  
  Model0011 = glmer(infected ~ 0 +                   median.ct + rel.susc + (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset=LogitMechPassage , family="binomial", data=Data25)
  Model0010 = glmer(infected ~ 0 +                   median.ct +            (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset=LogitMechPassage , family="binomial", data=Data25)
  
  Model0001 = glmer(infected ~ 0 +                               rel.susc + (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset=LogitMechPassage , family="binomial", data=Data25)
  Model0000 = glmer(infected ~ 0 +                                          (1|Line) + (1|Strain) + (1|Passage) + (1|block), offset=LogitMechPassage , family="binomial", data=Data25)
  
}else {
  Model1111 = glmer(infected ~ 0 + best.y + best.p + median.ct + rel.susc + (1|Line) + (1|Strain) + (1|Passage) + (1|block) + LogitMechPassage, family="binomial", data=Data25)
  Model1110 = glmer(infected ~ 0 + best.y + best.p + median.ct +            (1|Line) + (1|Strain) + (1|Passage) + (1|block) + LogitMechPassage, family="binomial", data=Data25)
  
  Model1101 = glmer(infected ~ 0 + best.y + best.p +             rel.susc + (1|Line) + (1|Strain) + (1|Passage) + (1|block) + LogitMechPassage, family="binomial", data=Data25)
  Model1100 = glmer(infected ~ 0 + best.y + best.p +                        (1|Line) + (1|Strain) + (1|Passage) + (1|block) + LogitMechPassage, family="binomial", data=Data25)
  
  Model1011 = glmer(infected ~ 0 + best.y +          median.ct + rel.susc + (1|Line) + (1|Strain) + (1|Passage) + (1|block) + LogitMechPassage , family="binomial", data=Data25)
  Model1010 = glmer(infected ~ 0 + best.y +          median.ct +            (1|Line) + (1|Strain) + (1|Passage) + (1|block) + LogitMechPassage, family="binomial", data=Data25)
  
  Model0111 = glmer(infected ~ 0 +          best.p + median.ct + rel.susc + (1|Line) + (1|Strain) + (1|Passage) + (1|block) + LogitMechPassage , family="binomial", data=Data25)
  Model0110 = glmer(infected ~ 0 +          best.p + median.ct +            (1|Line) + (1|Strain) + (1|Passage) + (1|block) + LogitMechPassage , family="binomial", data=Data25)
  
  Model1001 = glmer(infected ~ 0 + best.y +                      rel.susc + (1|Line) + (1|Strain) + (1|Passage) + (1|block) + LogitMechPassage , family="binomial", data=Data25)
  Model1000 = glmer(infected ~ 0 + best.y +                                 (1|Line) + (1|Strain) + (1|Passage) + (1|block) + LogitMechPassage , family="binomial", data=Data25)
  
  Model0101 = glmer(infected ~ 0 +          best.p +             rel.susc + (1|Line) + (1|Strain) + (1|Passage) + (1|block) + LogitMechPassage , family="binomial", data=Data25)
  Model0100 = glmer(infected ~ 0 +          best.p +                        (1|Line) + (1|Strain) + (1|Passage) + (1|block) + LogitMechPassage , family="binomial", data=Data25)
  
  Model0011 = glmer(infected ~ 0 +                   median.ct + rel.susc + (1|Line) + (1|Strain) + (1|Passage) + (1|block) + LogitMechPassage , family="binomial", data=Data25)
  Model0010 = glmer(infected ~ 0 +                   median.ct +            (1|Line) + (1|Strain) + (1|Passage) + (1|block) + LogitMechPassage , family="binomial", data=Data25)
  
  Model0001 = glmer(infected ~ 0 +                               rel.susc + (1|Line) + (1|Strain) + (1|Passage) + (1|block) + LogitMechPassage , family="binomial", data=Data25)
  Model0000 = glmer(infected ~ 0 +                                          (1|Line) + (1|Strain) + (1|Passage) + (1|block) + LogitMechPassage , family="binomial", data=Data25)
}

Models=list(Model1111, Model1110, Model1101, Model1100, Model1011, Model1010, 
            Model0111, Model0110, Model1001, Model1000, Model0101, Model0100, 
            Model0011, Model0010, Model0001, Model0000)
print(lapply(Models, AIC))      #print out the AIC values for all versions of the model

AICvalues<-unlist(lapply(Models, AIC))      #print out the AIC values for all versions of the model
min(AICvalues)

#model weights
M1111<-exp(-(AICvalues[1]-min(AICvalues))/2)
M1110<-exp(-(AICvalues[2]-min(AICvalues))/2)
M1101<-exp(-(AICvalues[3]-min(AICvalues))/2)
M1100<-exp(-(AICvalues[4]-min(AICvalues))/2)
M1011<-exp(-(AICvalues[5]-min(AICvalues))/2)
M1010<-exp(-(AICvalues[6]-min(AICvalues))/2)
M0111<-exp(-(AICvalues[7]-min(AICvalues))/2)
M0110<-exp(-(AICvalues[8]-min(AICvalues))/2)
M1001<-exp(-(AICvalues[9]-min(AICvalues))/2)
M1000<-exp(-(AICvalues[10]-min(AICvalues))/2)
M0101<-exp(-(AICvalues[11]-min(AICvalues))/2)
M0100<-exp(-(AICvalues[12]-min(AICvalues))/2)
M0011<-exp(-(AICvalues[13]-min(AICvalues))/2)
M0010<-exp(-(AICvalues[14]-min(AICvalues))/2)
M0001<-exp(-(AICvalues[15]-min(AICvalues))/2)
M0000<-exp(-(AICvalues[16]-min(AICvalues))/2)

sumall<-M1111+M1110+M1101+M1100+M1011+M1010+
  M0111+M0110+M1001+M1000+M0101+M0100+
  M0011+M0010+M0001+M0000

ct.weight<-(exp(-(AICvalues[1]-min(AICvalues))/2)+exp(-(AICvalues[2]-min(AICvalues))/2)+exp(-(AICvalues[5]-min(AICvalues))/2)+exp(-(AICvalues[6]-min(AICvalues))/2)+
  exp(-(AICvalues[7]-min(AICvalues))/2)+exp(-(AICvalues[8]-min(AICvalues))/2)+exp(-(AICvalues[13]-min(AICvalues))/2)+exp(-(AICvalues[14]-min(AICvalues))/2))/
  (exp(-(AICvalues[1]-min(AICvalues))/2)+exp(-(AICvalues[2]-min(AICvalues))/2)+exp(-(AICvalues[3]-min(AICvalues))/2)+exp(-(AICvalues[4]-min(AICvalues))/2)+
  exp(-(AICvalues[5]-min(AICvalues))/2)+exp(-(AICvalues[6]-min(AICvalues))/2)+exp(-(AICvalues[7]-min(AICvalues))/2)+exp(-(AICvalues[8]-min(AICvalues))/2)+
  exp(-(AICvalues[9]-min(AICvalues))/2)+exp(-(AICvalues[10]-min(AICvalues))/2)+exp(-(AICvalues[11]-min(AICvalues))/2)+exp(-(AICvalues[12]-min(AICvalues))/2)+
  exp(-(AICvalues[13]-min(AICvalues))/2)+exp(-(AICvalues[14]-min(AICvalues))/2)+exp(-(AICvalues[15]-min(AICvalues))/2)+exp(-(AICvalues[16]-min(AICvalues))/2))

p.weight<-(exp(-(AICvalues[1]-min(AICvalues))/2)+exp(-(AICvalues[2]-min(AICvalues))/2)+exp(-(AICvalues[3]-min(AICvalues))/2)+exp(-(AICvalues[4]-min(AICvalues))/2)+
  exp(-(AICvalues[7]-min(AICvalues))/2)+exp(-(AICvalues[8]-min(AICvalues))/2)+exp(-(AICvalues[11]-min(AICvalues))/2)+exp(-(AICvalues[12]-min(AICvalues))/2))/
  (exp(-(AICvalues[1]-min(AICvalues))/2)+exp(-(AICvalues[2]-min(AICvalues))/2)+exp(-(AICvalues[3]-min(AICvalues))/2)+exp(-(AICvalues[4]-min(AICvalues))/2)+
     exp(-(AICvalues[5]-min(AICvalues))/2)+exp(-(AICvalues[6]-min(AICvalues))/2)+exp(-(AICvalues[7]-min(AICvalues))/2)+exp(-(AICvalues[8]-min(AICvalues))/2)+
     exp(-(AICvalues[9]-min(AICvalues))/2)+exp(-(AICvalues[10]-min(AICvalues))/2)+exp(-(AICvalues[11]-min(AICvalues))/2)+exp(-(AICvalues[12]-min(AICvalues))/2)+
     exp(-(AICvalues[13]-min(AICvalues))/2)+exp(-(AICvalues[14]-min(AICvalues))/2)+exp(-(AICvalues[15]-min(AICvalues))/2)+exp(-(AICvalues[16]-min(AICvalues))/2))
  
y.weight<-(M1111+M1110+M1101+M1011+M1100+M1001+M1010+M1000)/sumall
p.weight<- (M1111+M1110+M1101+M1100+M0111+M0110+M0101+M0100)/sumall
ct.weight<- (M1111+M1110+M1011+M0111+M0110+M0011+M1010+M0010)/sumall
TCID50.weight<-(M1111+M1101+M1011+M0111+M1001+M0011+M0101+M0001)/sumall

#model averaged parameter estimates
Models.P=list(Model1111, Model1110, Model1101, Model1100, Model1011, Model1010,
            Model0111, Model0110, Model1001, Model1000, Model0101, Model0100,
            Model0011, Model0010, Model0001)
modavg(Models.P,parm="best.p")
modavg(Models.P,parm="best.y")
modavg(Models.P,parm="rel.susc")
modavg(Models.P,parm="median.ct")
modavg(Models.P,parm="LogitMechPassage")

#Analysis of best model - not sure how this works with the offset. Also, had to do na.rm for ct. 
#Probably a better way to do this would be to filter the data for all.
R2_p<- var(Data25$best.p*fixef(Model0110)[1])/(var(Data25$best.p*fixef(Model0110)[1])+var(Data25$median.ct*fixef(Model0110)[2], na.rm = TRUE)+summary(Model0110)$varcor$Line[1]+summary(Model0110)$varcor$Strain[1]+summary(Model0110)$varcor$Passage[1]+summary(Model0110)$varcor$block[1]+(pi^2)/3)
R2_ct<- var(Data25$median.ct*fixef(Model0110)[2], na.rm = TRUE)/(var(Data25$best.p*fixef(Model0110)[1])+var(Data25$median.ct*fixef(Model0110)[2], na.rm = TRUE)+summary(Model0110)$varcor$Line[1]+summary(Model0110)$varcor$Strain[1]+summary(Model0110)$varcor$Passage[1]+summary(Model0110)$varcor$block[1]+(pi^2)/3)
R2_Line<- summary(Model0110)$varcor$Line[1]/(var(Data25$best.p*fixef(Model0110)[1])+var(Data25$median.ct*fixef(Model0110)[2], na.rm = TRUE)+summary(Model0110)$varcor$Line[1]+summary(Model0110)$varcor$Strain[1]+summary(Model0110)$varcor$Passage[1]+summary(Model0110)$varcor$block[1]+(pi^2)/3)
R2_Strain<- summary(Model0110)$varcor$Strain[1]/(var(Data25$best.p*fixef(Model0110)[1])+var(Data25$median.ct*fixef(Model0110)[2], na.rm = TRUE)+summary(Model0110)$varcor$Line[1]+summary(Model0110)$varcor$Strain[1]+summary(Model0110)$varcor$Passage[1]+summary(Model0110)$varcor$block[1]+(pi^2)/3)
R2_Passage<- summary(Model0110)$varcor$Passage[1]/(var(Data25$best.p*fixef(Model0110)[1])+var(Data25$median.ct*fixef(Model0110)[2], na.rm = TRUE)+summary(Model0110)$varcor$Line[1]+summary(Model0110)$varcor$Strain[1]+summary(Model0110)$varcor$Passage[1]+summary(Model0110)$varcor$block[1]+(pi^2)/3)
R2_block<- summary(Model0110)$varcor$block[1]/(var(Data25$best.p*fixef(Model0110)[1])+var(Data25$median.ct*fixef(Model0110)[2], na.rm = TRUE)+summary(Model0110)$varcor$Line[1]+summary(Model0110)$varcor$Strain[1]+summary(Model0110)$varcor$Passage[1]+summary(Model0110)$varcor$block[1]+(pi^2)/3)


r.squaredGLMM(Model0110)


