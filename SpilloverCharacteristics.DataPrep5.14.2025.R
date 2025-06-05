#combine FCX, FCY, FCZ, FDI datasets for spillover characteristics experiment
setwd("~/PostDoc Computer/SpilloverCharacteristics/PLOS Biology/Code")

library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(cowplot)
library(rstan)

#upload and combine transmission data
#FCX, FXY, FCZ, FDI are codes I used to refer to the experimental blocks
FCX<-read.csv("FCXexperiment.population.csv", header=TRUE)
FCX$block<-rep("FCX", length(FCX$Strain))
FCY<-read.csv("FCYexperiment.population.csv", header=TRUE)
FCY$block<-rep("FCY", length(FCY$Strain))
FCZ<-read.csv("FCZexperiment.population.csv", header=TRUE)
FCZ$block<-rep("FCZ", length(FCZ$Strain))
FDI<-read.csv("FDIexperiment.population.csv", header=TRUE)
FDI$block<-rep("FDI", length(FDI$Strain))

all.data<-rbind(FCX, FCY, FCZ, FDI)

#Correct some errors in the datasheets
all.data$Species[all.data$Species=="sp.66"]<-"sp.65"
all.data$Species[all.data$Species=="sp.66 "]<-"sp.65"
all.data$Species[all.data$Species=="sp.65 "]<-"sp.65"
all.data$Species[all.data$Species=="sp.25 "]<-"sp.25"

all.data<-filter(all.data, Strain!="JU1331") #included this strain in FCX, but decided not to pursue further
all.data<-filter(all.data, !is.na(Passage))

#Remove negative control C. elegans data
all.data<-filter(all.data, Line!="FCX1")
all.data<-filter(all.data, Line!="FCY1")
all.data<-filter(all.data, Line!="FCZ1")
all.data<-filter(all.data, Line!="FDI1")

#Remove data from lines that became contaminated and had to be dropped
all.data<-filter(all.data, Line!="FDI20") #contaminated
all.data<-filter(all.data, Line!="FDI19") #contaminated
all.data<-filter(all.data, Line!="FDI22") #contaminated
all.data<-filter(all.data, Line!="FDI16") #contaminated
all.data<-filter(all.data, Line!="FDI21") #contaminated
all.data<-filter(all.data, Line!="FDI34") #just stopped growing
all.data<-filter(all.data, Sample!="FCY17.5") #continued infinite cts
all.data<-filter(all.data, Sample!="FCY17.6") #continued infinite cts
all.data<-filter(all.data, Sample!="FDI25.8") #false positive
all.data<-filter(all.data, Sample!="FDI24.4") #false positive
all.data<-filter(all.data, Sample!="FCX28.5") #continued infinite
all.data<-filter(all.data, Sample!="FCX28.4") #false positive
all.data<-filter(all.data, Sample!="FCX27.3") #continued infinite
all.data<-filter(all.data, Sample!="FCY11.2") #continued infinite
all.data<-filter(all.data, Sample!="FCY30.2") #continued infinite

#Graphing transmission - add species name to strains
all.data$NAME[all.data$Strain=="ZF1092"]<-"C. sp. 25\n ZF1092"
all.data$NAME[all.data$Strain=="JU1873"]<-"C. wallacei\n JU1873"
all.data$NAME[all.data$Strain=="JU1428"]<-"C. tropicalis\n JU1428"
all.data$NAME[all.data$Strain=="VX80"]<-"C. latens\n VX80"
all.data$NAME[all.data$Strain=="JU1857"]<-"C. macrosperma\n JU1857"
all.data$NAME[all.data$Strain=="JU4093"]<-"C. sp. 65\n JU4093"
all.data$NAME[all.data$Strain=="JU724"]<-"C. latens\n JU724"
all.data$NAME[all.data$Strain=="SB454"]<-"C. sulstoni\n SB454"
all.data$NAME[all.data$Strain=="JU1580"]<-"C. elegans\n JU1580"

#Figure out order to put strains in from least to most passages
orders<-filter(all.data, Ct!=Inf)
totals<-all.data%>%group_by(Strain)%>%summarise(number=n())
orders<-orders%>%group_by(Strain)%>%summarise(number=n())

#reorder the factors
all.data$Strain<-factor(all.data$Strain, levels = c("ZF1092","JU1428","VX80","JU1857","JU1873","JU724","JU4093","SB454","JU1580"))
all.data$NAME<-factor(all.data$NAME, levels = c("C. sp. 25\n ZF1092","C. tropicalis\n JU1428","C. latens\n VX80","C. macrosperma\n JU1857","C. wallacei\n JU1873","C. latens\n JU724","C. sp. 65\n JU4093","C. sulstoni\n SB454","C. elegans\n JU1580"))

Figure1<-ggplot(all.data, aes(Passage, Ct, group=Line, color=Strain, shape=block))+
  geom_point(fill="white", size=3)+facet_wrap(~NAME, nrow=1)+geom_line()+theme_bw()+
  scale_shape_manual(values=c(21,22,23,24))+
  scale_y_reverse()+scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10))+
  theme(legend.position = "none")+theme(strip.text = element_text(face="italic"))+
  theme(strip.background = element_rect(fill="white"))+
  theme(axis.text = element_text(color="black"))+
  theme(panel.grid.minor = element_blank())+
  ylab("Orsay virus RNA1 cycle threshold (Ct)")

###########################################################################
#Organize data into two dataframes:
  #"farthest" has the passages until loss.
  #"each.passage" has whether a passage was successful or not 

#assign time until loss (dataframe 1)
all.data.inf<-filter(all.data, Ct!=Inf)
Passage0<-filter(all.data, Passage==0)
P0<-filter(Passage0, Ct==Inf)
all.data.inf<-rbind(P0, all.data.inf) #combine all P0 plates with infected plates
ordered.data<-all.data.inf[order(all.data.inf$Line),]
farthest<-all.data.inf%>%group_by(Line, Strain, Species, block)%>%summarise(trans.ability=max(Passage))
hist(farthest$trans.ability)

#designate passages as infected or not (0 or 1) (dataframe 2)
each.passage<-select(all.data, Sample,Passage, Line, Ct) 
each.passage$infected[each.passage$Ct==Inf]<-0
each.passage$infected[each.passage$Ct!=Inf]<-1
#"each.passage has 0 for a line not passaging and 1 for it passaging at each passage"

###############################################################################
###upload prevalence and intensity data
FCX.s<-read.csv("FCX.strips.csv", header=TRUE)
FCX.s$block<-rep("FCX", length(FCX.s$Strain))
FCX.s<-filter(FCX.s, X!=247) #these were mis-assigned
FCX.s<-filter(FCX.s, X!=254)
FCX.s<-filter(FCX.s, X!=261)
FCX.s<-filter(FCX.s, X!=268)
FCX.s<-filter(FCX.s, X!=276)
FCX.s<-filter(FCX.s, X!=284)
FCX.s<-filter(FCX.s, X!=289)
FCX.s<-filter(FCX.s, X!=292)
FCX.s<-filter(FCX.s, X!=297)
FCX.s<-filter(FCX.s, X!=300)
FCY.s<-read.csv("FCY.strips.csv", header=TRUE)
FCY.s$block<-rep("FCY", length(FCY.s$Strain)) #found mistake line 382, A2 FCY28 - changed from 0 to 2 worms
FCY.s<-filter(FCY.s, X!=394)
FCY.s<-filter(FCY.s, X!=405)
FCZ.s<-read.csv("FCZ.strips.csv", header=TRUE)
FCZ.s$block<-rep("FCZ", length(FCZ.s$Strain))
FDI.s<-read.csv("FDI.strips.csv", header=TRUE)
FDI.s$block<-rep("FDI", length(FDI.s$Strain))

Passage.Line<-unique(all.data$Line)
strips<-bind_rows(FCX.s, FCY.s, FCZ.s, FDI.s)
check<-anti_join(strips, all.data, by="Sample") #This all checks out- appropriate rows will be deleted by the next line of code
strips$Species[strips$Species=="sp.66"]<-"sp.65" #correct mistake
strips$Species[strips$Species=="sp.66 "]<-"sp.65" #correct mistake
strips$Species[strips$Species=="sp.65 "]<-"sp.65" #correct mistake
strips$Species[strips$Species=="sp.25 "]<-"sp.25" #correct mistake
strips<-filter(strips, Sample%in%Passage.Line)

#Designate samples as infected (1) or not (0)
strips$inf[strips$Ct==Inf]<-0
strips$inf[is.na(strips$inf)]<-1
write.csv(strips, "strip.ct.data.csv") #make dataset of combined strip data

#First group strip data by number infected of the different categories (0,1,2,or 3 worms)
strip.groups<-strips%>%group_by(Sample, WormNumber, Species, Strain, block)%>%summarise(inf.number=sum(inf), number.wells=length(inf))
write.csv(strip.groups, "strip.extractions.csv")

samplesize<-strip.groups%>%group_by(Sample)%>%summarise(totalwells=sum(number.wells))
#Look at control wells with 0 worms
zero.strips<-filter(strips, WormNumber==0)

###################################################################################
#Load shedding data
shedding.X<-read.csv("FCX.shedding.csv",header=TRUE)
shedding.Y<-read.csv("FCY.shedding.csv",header=TRUE)
shedding.Z<-read.csv("FCZ.shedding.csv",header=TRUE)
shedding.I<-read.csv("FDI.shedding.csv",header=TRUE)

shedding<-rbind(shedding.X, shedding.Y, shedding.Z, shedding.I)
check<-anti_join(shedding, all.data, by="Sample") #This all checks out- appropriate rows will be deleted by the next line of code
shedding<-filter(shedding, Sample%in%Passage.Line) #dropping some- make sure they're not just mis-la
write.csv(shedding, "sheddingdata.csv") #combined shedding datasheet

#Loop to estimate maximum likelihood prevalence and shedding together
#Estimate each at the plate level
#This gives maximum likelihood prevalence estimate and maximum likelihood chance that an infected worm makes a plate glow
LINE<-unique(all.data$Line)
ESTIMATE<-NULL
for (i in LINE){
  shed.data<-filter(shedding,Sample==i) #look at shedding from 1 plate at a time
  prev.data<-filter(strip.groups, Sample==i) #look at strip tube results for 1 plate at a time
  single<-prev.data$inf.number[2] #this is how many times a single came up infected
  double<-prev.data$inf.number[3] #This is how many times a double extraction came up infected
  triple<-prev.data$inf.number[4] #This is how many times a triple extraction came up infected
  single[is.na(single)]<-0
  double[is.na(double)]<-0
  triple[is.na(triple)]<-0 #Sometimes I didn't have any triples, so put in 0 so we don't have NAs.
  single.tries<-prev.data$number.wells[2] #this is how many wells were devoted to singles
  double.tries<-prev.data$number.wells[3] #this is how many wells were devoted to doubles
  triple.tries<-prev.data$number.wells[4] #this is how many wells were devoted to triples
  single.tries[is.na(single.tries)]<-0 
  double.tries[is.na(double.tries)]<-0 
  triple.tries[is.na(triple.tries)]<-0 #Sometimes I didn't have any triples, so put in 0 so we don't have NAs

  CALC<-NULL #Make an empty matrix for each plate
  for(p in seq(0,1,0.01)){ #infection prevalence ranges from 0-1
    for (y in seq(0,1,0.01)){ #This is the probability that 1 worm causes glowing (ranges from 0-1)
      NOGLOW<-NULL #For each p and y, figure out the probability of 0-15 worms being infected 
      for (N in seq(0,15,1)){
        outcome<-(dbinom(N,15,p))*((1-y)^N)
        NOGLOW<-rbind(outcome,NOGLOW)
      } 
      chance.NoGlow<-sum(NOGLOW)
      chance.NoGlow<-round(chance.NoGlow, digits = 8)
      loglike<-dbinom(shed.data$Plates.Glowing,shed.data$Plates.Total, 1-chance.NoGlow, log=TRUE)+
        dbinom(single,single.tries,1-dbinom(0,1,p), log=TRUE)+
        dbinom(double, double.tries, 1-dbinom(0,2,p), log=TRUE)+
        dbinom(triple, triple.tries, 1-dbinom(0,3,p), log=TRUE)
      record<-c(p, y, loglike)
      CALC<-rbind(CALC, record)
    }
  }
  
  CALC<-as.data.frame(CALC)
  colnames(CALC)<-c("p","y","loglike")
  max.like<-CALC%>%filter(loglike>(max(loglike)-1.92))
  best.y<-max(max.like$y[max.like$loglike==max(max.like$loglike)]) #should this be a median?
  l.y<-min(max.like$y)
  h.y<-max(max.like$y)
  best.p<-max(max.like$p[max.like$loglike==max(max.like$loglike)]) #should this be a median?
  l.p<-min(max.like$p)
  h.p<-max(max.like$p)
  output2<-c(i, best.y,l.y,h.y,best.p,l.p,h.p)
  ESTIMATE<-rbind(output2, ESTIMATE)
}

ESTIMATE<-as.data.frame(ESTIMATE)
colnames(ESTIMATE)<-c("Sample","best.y","l.y","h.y", "best.p","l.p","h.p")

ESTIMATE$best.y<-as.numeric(ESTIMATE$best.y)
ESTIMATE$best.p<-as.numeric(ESTIMATE$best.p)
ESTIMATE$l.y<-as.numeric(ESTIMATE$l.y)
ESTIMATE$h.y<-as.numeric(ESTIMATE$h.y)
ESTIMATE$l.p<-as.numeric(ESTIMATE$l.p)
ESTIMATE$h.p<-as.numeric(ESTIMATE$h.p)
#118 sample lines

#Consolidate the metadata about each line
INFO<-unique(all.data[c("Line","Strain","Species","block")])
est.info<-full_join(ESTIMATE, INFO, by=c("Sample"="Line"))

#Make graphs of prevalence (p) and shedding (y) across strains (Figure 2 b,c)
est.info.graph = group_by(est.info, Strain, block) %>% mutate(n = n())
est.info.graph$Strain<-factor(est.info.graph$Strain, levels = c("ZF1092","JU1428","VX80","JU1857","JU1873","JU724","JU4093","SB454", "JU1580"))

F2b<-ggplot(est.info.graph, aes(Strain, best.p, shape=block, color=Strain))+
  geom_point(position=position_jitterdodge(jitter.width = 0.15))+
  geom_boxplot(data=subset(est.info.graph, n>1), alpha=0)+
  theme_bw()+
  scale_shape_manual(values=c(21,22,23,24))+
  ylab("Maximum likelihood prevalence \nestimate (p)")+
  theme(legend.position = "none")+
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5))
F2c<-ggplot(est.info.graph, aes(Strain, best.y, shape=block, color=Strain))+
  geom_point(position=position_jitterdodge(jitter.width = 0.15))+
  scale_shape_manual(values=c(21,22,23,24))+
  geom_boxplot(data=subset(est.info.graph, n>1), alpha=0)+
  theme_bw()+ylab("Maximum likelihood shedding \nability (y)")+theme(legend.position = "none")+
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5))

###############################################################################
#Estimate intensity using prevalence to correct.
#Make a dataframe of the chances that only 1 worm will be infected in 2 worm sets and 3 worm sets
CUTOFFS<-NULL
for (i in seq(0,1,by=0.01)){
  est3<-dbinom(1,3,i)/(1-dbinom(0,3,i))
  est2<-dbinom(1,2,i)/(1-dbinom(0,2,i))
  cutoffs<-c(i,est3, est2)
  CUTOFFS<-rbind(cutoffs, CUTOFFS)
}

CUTOFFS<-as.data.frame(CUTOFFS)
colnames(CUTOFFS)<-c("prevalence", "est3", "est2")
CUTOFFS$est3[101]<-1 #this is just so code doesn't break later
CUTOFFS$est2[101]<-1 #this is just so code doesn't break later

#Use infected dataframe and cut out those with a high chance (>30%) of more than one worm
infected<-filter(strips, inf==1)
infected.small<-filter(infected, WormNumber!=0)
estlines<-unique(est.info$Sample)
int.lines<-unique(infected.small$Line)
check<-anti_join(est.info, infected.small, by=c("Sample"="Line")) #This all checks out- deleted lines don't have any infections for intensity to be measured
check<-anti_join(infected.small,est.info, by=c("Line"="Sample")) #checks out
LINE<-intersect(estlines, int.lines)

NEW<-NULL
for (i in LINE){
  this.line<-filter(infected.small,Line==i)
  check.line<-filter(est.info, Sample==i)
  prev<-check.line$best.p[1]
  check.cutoff<-filter(CUTOFFS,prevalence==prev)
  if(check.cutoff$est3[1]<0.70){
    new<-filter(this.line, WormNumber!=3)
  } else (new<-this.line)
  if(check.cutoff$est2[1]<0.70){
    new<-filter(new, WormNumber!=2)
  }else(new<-this.line)
  NEW<-rbind(new, NEW)
}

NEW.intensities<-as.data.frame(NEW) #Length =442

#compare uncorrected and corrected intensities
intensity.data<-NEW.intensities%>%group_by(Strain, Sample, Species, block)%>%summarise(median.ct=median(Ct))
intensity.data.uncor<-infected.small%>%group_by(Strain, Sample, Species, block)%>%summarise(median.ct.uncor=median(Ct))
comb.data<-full_join(intensity.data, intensity.data.uncor, by=c("Strain","Species","block","Sample"))
n_o<-full_join(comb.data, est.info, by=c("Strain","Species","block","Sample"))

#Make supplemental figure to compare corrected and uncorrected intensities
s6.1<-ggplot(n_o, aes(median.ct.uncor, median.ct))+geom_point()+scale_y_reverse()+scale_x_reverse()+theme_bw()+
  xlab("Uncorrected infection intensity (median Ct)")+ylab("Corrected infection intensity (median Ct)")
s6.2<-ggplot(n_o, aes(best.p, median.ct))+geom_point(shape=21, fill="Black", color="Red")+
  geom_point(data=n_o, aes(best.p, median.ct.uncor), shape=1)+theme_bw()+
  xlab("Estimated infection prevalence")+ylab("Infection intensity (median Ct)")
s6<-plot_grid(s6.1, s6.2, labels = c("A","B"))
save_plot("SupplementalIntensities.jpg", s6, base_height = 5)

#make figures for intensity
intensity.data$Strain<-factor(intensity.data$Strain, levels = c("ZF1092","JU1428","VX80","JU1857","JU1873","JU724","JU4093","SB454", "JU1580"))
intensity.graph = group_by(intensity.data, Strain, block) %>% mutate(n = n())
F2d<-ggplot(intensity.data, aes(Strain, median.ct, shape=block, color=Strain))+
  geom_point(position=position_jitterdodge(jitter.width = 0.15))+
  scale_shape_manual(values=c(21,22,23,24))+
  geom_boxplot(data=subset(intensity.graph, n>1), alpha=0)+
  theme_bw()+ylab(expression("    Infection intensity \n(corrected median Ct)"~(mu)))+theme(legend.position = "none")+scale_y_reverse()+
  theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5))

#Predict Ct values from prevalence data - use all correct intensities datapoints (not just means)
pre.ct<-full_join(NEW.intensities, est.info, by=c("Sample", "Strain","Species", "block")) #length = 454
pre.ct$Strain<-factor(pre.ct$Strain, levels = c("ZF1092","JU1428","VX80","JU1857","JU1873","JU724","JU4093","SB454", "JU1580"))
pre.ct<-filter(pre.ct, !is.na(Ct)) #these all check out

ggplot(pre.ct, aes(best.p, Ct, shape=factor(WormNumber), color=Strain))+
  geom_point()+
  theme_bw()+
  xlab("Estimated infection prevalence")+ylab("Infection intensity (Ct)")

model<-lmer(data=pre.ct, Ct~best.p+(1|Strain))
summary(model)
sd(resid(model)) #5.24

pre.ct$pre.ct<-predict(model) #add in predicted effects of NA (with strain random effect prevalence)
fixef(model)
ranef(model)

#changed 5/16 because it wasn't working
intensity.data<-pre.ct%>%group_by(Strain, Sample, Species, block)%>%summarise(median.ct=median(Ct), pre.ct=mean(pre.ct))

#Put intensity, prevalence, and shedding data together
prev.int<-full_join(intensity.data,est.info, by=c("Sample","Species","Strain","block"))
#fill in nas for pre.ct
prev.int$pre.ct[107]<-fixef(model)[1]+fixef(model)[2]*prev.int$best.p[107]+ranef(model)$Strain[1,1]
prev.int$pre.ct[108]<-fixef(model)[1]+fixef(model)[2]*prev.int$best.p[108]+ranef(model)$Strain[1,1]
prev.int$pre.ct[109]<-fixef(model)[1]+fixef(model)[2]*prev.int$best.p[109]+ranef(model)$Strain[1,1]
prev.int$pre.ct[110]<-fixef(model)[1]+fixef(model)[2]*prev.int$best.p[110]+ranef(model)$Strain[1,1]
prev.int$pre.ct[115]<-fixef(model)[1]+fixef(model)[2]*prev.int$best.p[115]+ranef(model)$Strain[1,1]
prev.int$pre.ct[111]<-fixef(model)[1]+fixef(model)[2]*prev.int$best.p[111]+ranef(model)$Strain[7,1]
prev.int$pre.ct[114]<-fixef(model)[1]+fixef(model)[2]*prev.int$best.p[114]+ranef(model)$Strain[7,1]
prev.int$pre.ct[118]<-fixef(model)[1]+fixef(model)[2]*prev.int$best.p[118]+ranef(model)$Strain[7,1]
prev.int$pre.ct[112]<-fixef(model)[1]+fixef(model)[2]*prev.int$best.p[112]+ranef(model)$Strain[5,1]
prev.int$pre.ct[113]<-fixef(model)[1]+fixef(model)[2]*prev.int$best.p[113]+ranef(model)$Strain[2,1]
prev.int$pre.ct[116]<-fixef(model)[1]+fixef(model)[2]*prev.int$best.p[116]+ranef(model)$Strain[6,1]
prev.int$pre.ct[117]<-fixef(model)[1]+fixef(model)[2]*prev.int$best.p[117]+ranef(model)$Strain[4,1]

s7<-ggplot(prev.int, aes(median.ct, pre.ct, color=Strain, shape=block))+geom_point()+
  theme_bw()+
  xlab("Corrected median Ct")+ylab("Predicted Ct")
save_plot("SupplementalFig7.jpg", s7, base_height = 4)


#combine intensity, prevalence, and shedding data with "farthest" transmission datasheet
dataset1<-full_join(farthest, prev.int, by=c("Line"="Sample","Strain","Species","block"))

#Single correlation plots for supplement
dataset1.exp<-filter(dataset1, Strain!="JU1580")
dataset1.exp$Strain <- factor(dataset1.exp$Strain, levels = c("ZF1092","JU1428","VX80","JU1857","JU1873","JU724","JU4093","SB454","JU1580"))

s.1<-ggplot(dataset1, aes(best.p, trans.ability, color=Strain))+
  geom_point()+facet_wrap(~Strain, nrow=3)+geom_smooth(method="lm")+
  theme_bw()+ylab("Passages until loss")+
  xlab("Maximum likelihood prevalence")+
  theme(axis.text = element_text(color="black"))
save_plot("Supplement1.jpg", s.1, base_height = 5)

s.3<-ggplot(dataset1, aes(best.y, trans.ability, color=Strain))+
  geom_point()+facet_wrap(~Strain, nrow=3)+geom_smooth(method="lm")+
  theme_bw()+ylab("Passages until loss")+
  xlab("Maximum likelihood shedding ability")+
  theme(axis.text = element_text(color="black"))
save_plot("Supplement3.jpg", s.3, base_height = 5)

s.2<-ggplot(dataset1, aes(median.ct, trans.ability, color=Strain))+
  geom_point()+facet_wrap(~Strain, nrow=3)+geom_smooth(method="lm")+
  theme_bw()+ylab("Passages until loss")+
  scale_x_reverse()+xlab("Average infection intensity (median Ct)")+
  theme(axis.text = element_text(color="black"))
save_plot("Supplement2.jpg", s.2, base_height = 5)


####################################################################################
#Add TCID50 data; Run the following for supplement, otherwise skip to line 326
#don't have the TCID50 data for JU1428.old, VX80.new

blocks.strains<-read.csv("Strainsandblocks.csv")
TC<-read.csv("TCID50.csv")
TC<-select(TC, strain, TCID50.est, conf.low, conf.high, block, strain.1)

#rename
TC$strain[TC$strain=="ZF1092.old"]<-"ZF1092.1"
TC$strain[TC$strain=="ZF1092.new"]<-"ZF1092.2"
TC$strain[TC$strain=="JU1873.old"]<-"JU1873.1"
TC$strain[TC$strain=="JU1873.new"]<-"JU1873.2"
TC$strain[TC$strain=="JU1428.new"]<-"JU1428.2"
TC$strain[TC$strain=="VX80.old"]<-"VX80.1"
TC$strain[TC$strain=="JU1857.old"]<-"JU1857.1"
TC$strain[TC$strain=="JU1857.new"]<-"JU1857.2"
TC$strain[TC$strain=="JU4093.old"]<-"JU4093.1"
TC$strain[TC$strain=="JU4093.new"]<-"JU4093.2"
TC$strain[TC$strain=="JU724.old"]<-"JU724.1"
TC$strain[TC$strain=="JU724.new"]<-"JU724.2"
TC$strain[TC$strain=="SB454.old"]<-"SB454.1"
TC$strain[TC$strain=="SB454.new"]<-"SB454.2"

TC$strain<-factor(TC$strain, levels = c("ZF1092.1","ZF1092.2","JU1428.2","VX80.1","JU1857.1","JU1857.2","JU1873.1","JU1873.2","JU724.1","JU724.2","JU4093.1","JU4093.2","SB454.1","SB454.2", "JU1580"))
TC$strain.1<-factor(TC$strain.1, levels = c("ZF1092","JU1428","VX80","JU1857","JU1873","JU724","JU4093","SB454", "JU1580"))
TC$Group<-c(1,1,1,1,1,2,
            1,1,1,1,1,3,
            1,1,1,4,1,1)
pd <- position_dodge(width=0.4)

TCID50graph<-ggplot(TC, aes(strain, TCID50.est, color=strain.1, group=Group))+
  geom_point(position=pd)+
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), position=pd)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))+
  scale_y_log10()+ylab("TCID50 estimate")+xlab("Strain/thaw line")+
  theme(legend.position = "none")+theme(axis.text.x = element_text(color="black"))+
  theme(axis.title = element_text(size=14))+theme(axis.text = element_text(size=13))
save_plot("SpillChar.Supp.S5.jpg", TCID50graph, base_height = 7)

##############################################################################################################################################
#Add TCID50 data to analysis dataframe (combined)
blocks.strains<-read.csv("Strainsandblocks.csv")
TC.com<-read.csv("TCID50.combined.csv")

TC.com$Strain[TC.com$Strain=="JU724.old"]<-"JU724.1"
TC.com$Strain[TC.com$Strain=="JU724.new"]<-"JU724.2"

TC.com$strain <- factor(TC.com$Strain, levels = c("ZF1092","JU1428","VX80","JU1857","JU1873","JU724.1","JU724.2","JU4093","SB454","JU1580"))
TC.com$Strain <- factor(TC.com$strain.group, levels = c("ZF1092","JU1428","VX80","JU1857","JU1873","JU724","JU4093","SB454","JU1580"))

#divide TCID50 by TCID50 of JU1580 for a measure of relative susceptibility
TC.com$rel.susc<-TC.com$max.like/TC.com$max.like[2] #This is the value for JU1580
TC.com$rel.susc.low<-TC.com$conf.low/TC.com$max.like[2] #This is the value for JU1580
TC.com$rel.susc.high<-TC.com$conf.high/TC.com$max.like[2] #This is the value for JU1580

F2e<-ggplot(TC.com, aes(strain, rel.susc, color=Strain))+geom_point()+
  geom_errorbar(aes(ymin=rel.susc.low, ymax=rel.susc.high))+theme_bw()+
  scale_shape_manual(values=c(21,22,23,24))+
  scale_y_log10()+ylab("Relative suscetibility (s)")+theme(legend.position = "none")+
  xlab("Strain")+theme(axis.text = element_text(color="black"))+
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5))

#Designate different JU724 thaws in dataset1
dataset1$Strain<-as.character(dataset1$Strain)
for (i in 1:length(dataset1$Strain)){
  if (dataset1$Strain[i]=="JU724") {
    if(dataset1$block[i]=="FCX"){
    dataset1$Strain[i]<-"JU724.1"
  } else if (dataset1$block[i]=="FCY"){
    dataset1$Strain[i]<-"JU724.1"
  } else if (dataset1$block[i]=="FCZ"){
    dataset1$Strain[i]<-"JU724.2"
  } else if (dataset1$block[i]=="FDI"){
    dataset1$Strain[i]<-"JU724.2"}
  else {(dataset1$Strain[i]<-dataset1$Strain[i])}
  } else {(dataset1$Strain[i]<-dataset1$Strain[i])}
}

colnames(TC.com)<-c("col1","strain.overall", "number", "TCID50","max.like.value", "sub.val", "conf.low", "conf.high", "strain.group", "Strain", "rel.susc", "rel.susc.low", "rel.susc.high")  
TC.com.small<-select(TC.com, strain.overall, TCID50, conf.low, conf.high, strain.group, Strain, rel.susc, rel.susc.low, rel.susc.high)

#Add TCID50 data to dataset1
dataset1.all<-full_join(dataset1, TC.com.small, by=c("Strain"))

#calculate p for the paper
p_paper =function(p, y, Cj, avgi, sig)        
{
  return(p*y * pnorm(Cj,avgi,sig))
}

a<-40 #Ct value of JU1580 for transmission (if it's detectable, then it transmits).
sig<-5.25 #sd(residuals from fitted lmer) line 299 SpilloverCharacteristics.DataPrep9.22.2023
dataset1.all$p.paper=p_paper(dataset1.all$best.p, dataset1.all$best.y, a+log2(dataset1.all$rel.susc), dataset1.all$pre.ct, sig)
dataset1.all$Strain <- factor(dataset1.all$Strain, levels = c("ZF1092","JU1428","VX80","JU1857","JU1873","JU724.1","JU724.2","JU4093","SB454","JU1580"))
dataset1.all$strain.group <- factor(dataset1.all$strain.group, levels = c("ZF1092","JU1428","VX80","JU1857","JU1873","JU724","JU4093","SB454","JU1580"))

#rename for graphs
dataset1.all$NAME[dataset1.all$strain.group=="ZF1092"]<-"C. sp. 25\n ZF1092"
dataset1.all$NAME[dataset1.all$strain.group=="JU1873"]<-"C. wallacei\n JU1873"
dataset1.all$NAME[dataset1.all$strain.group=="JU1428"]<-"C. tropicalis\n JU1428"
dataset1.all$NAME[dataset1.all$strain.group=="VX80"]<-"C. latens\n VX80"
dataset1.all$NAME[dataset1.all$strain.group=="JU1857"]<-"C. macrosperma\n JU1857"
dataset1.all$NAME[dataset1.all$strain.group=="JU4093"]<-"C. sp. 65\n JU4093"
dataset1.all$NAME[dataset1.all$strain.group=="JU724"]<-"C. latens\n JU724"
dataset1.all$NAME[dataset1.all$strain.group=="SB454"]<-"C. sulstoni\n SB454"
dataset1.all$NAME[dataset1.all$strain.group=="JU1580"]<-"C. elegans\n JU1580"

dataset1.all$NAME <- factor(dataset1.all$NAME, levels = c("C. sp. 25\n ZF1092","C. tropicalis\n JU1428",
                                                          "C. latens\n VX80","C. macrosperma\n JU1857",
                                                          "C. wallacei\n JU1873","C. latens\n JU724",
                                                          "C. sp. 65\n JU4093","C. sulstoni\n SB454","C. elegans\n JU1580"))

F2f<-ggplot(dataset1.all, aes(NAME, p.paper, shape=block, color=strain.group))+
  geom_point(position=position_jitterdodge(jitter.width = 0.15))+
  theme_bw()+
  scale_shape_manual(values=c(21,22,23,24))+
  geom_boxplot(data=filter(dataset1.all, Strain!="JU1580"), alpha=0)+
  ylab("Probability that a worm can \ninitiate infection on new plate (q)")+
  theme(legend.position = "none")+
  theme(axis.text = element_text(color="black"))+
  xlab("")


summaryplots<-plot_grid(F2b,F2c,F2d,F2e, labels=c("B","C","D","E"),ncol=4)
step2<-plot_grid(Figure1, summaryplots, F2f, labels=c("A", "", "F"), nrow=3)
save_plot("SpillChar.Big.jpg", step2, base_height = 10.5)
save_plot("Figure1.tiff", step2, base_height = 10.5)

dataset1.exp<-filter(dataset1.all, Strain!="JU1580")
write.csv(dataset1.exp, "spill.char.dataset1.csv")

#Prepare dataset for mechanistic model
dataset2<-full_join(dataset1.exp, each.passage, by="Line")
dataset2<-filter(dataset2, Passage!=0)
dataset2<-filter(dataset2, !is.na(Strain)) #This takes out the positive controls
dataset2$Passage<-as.character(dataset2$Passage)

write.csv(dataset2, "spill.char.dataset2.csv")
