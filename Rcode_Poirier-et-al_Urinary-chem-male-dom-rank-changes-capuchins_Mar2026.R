

# R code associated to the manuscript submitted to Chemical Senses:
# Longitudinal variation in urinary chemicals associated with 
# dominance rank changes in wild male capuchin monkeys
#
# Authors: 
# Alice C. Poirier, Nelle K. Kulick, Megan Petersdorf, Marlen Kücklich, Suheidy Romero Morales, 
# Brigitte M. Weiß, Claudia Wiesner, Nicholas Chapoy, Wendy Téllez Arias, Anja Widdig, 
# Jacinta C. Beehner, Amanda D. Melin, Katharine M. Jack




# Libraries  ----
library(devtools)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(car)
library(lme4)
library(glmmTMB)
library(ggpubr)
library(randomForest)
library(DHARMa)
library(performance)
library(predictmeans)
library(sjPlot)
library(emmeans)
library(grid)



# 1) Datasets ----

## Chem dataset
# Dataset_CompList_Poirier-et-al_Urinary-chem-and-male-dom-rank-changes-capuchins_Feb2026.csv
AMR.chem<- read.csv(file.choose(), header=T, na.strings="NA", sep=";")
AMR.chem[,c(1:2)]<- lapply(AMR.chem[,c(1:2)], factor)
str(AMR.chem) # 276 filenames, 103 comp

## Sample dataset
# Dataset_SampleList_Poirier-et-al_Urinary-chem-and-male-dom-rank-changes-capuchins_Feb2026.csv
AMR.spl<- read.csv(file.choose(), header=T, na.strings="NA", sep=";")
AMR.spl[,c(1:10)]<- lapply(AMR.spl[,c(1:10)], factor)
str(AMR.spl) # 276 filenames

### z-transform seasonality variable for full dataset (hereafter: all) ----
AMR.spl.all<- mutate(AMR.spl, z.Seasonality=as.vector(scale(Seasonality)))
AMR.spl.all[,c(1:10)]<- lapply(AMR.spl.all[,c(1:10)], factor)

### Merge and add number of comp per sample (Ncomp) and sum of all peak areas (SumArea) to chem and spl datasets ----
AMR.chem.all<- merge(AMR.spl.all, AMR.chem, by="Filename")
AMR.chem.all[,c(1:10,18)]<- lapply(AMR.chem.all[,c(1:10,18)], factor)
AMR.chem.allN<- data.frame(AMR.chem.all %>% group_by(Filename) %>% reframe(Ncomp=n()))
AMR.chem.allA<- data.frame(AMR.chem.all %>% group_by(Filename) %>% reframe(SArea=sum(Area)))

AMR.chem.all<- merge(AMR.chem.all, AMR.chem.allN, by="Filename")
AMR.chem.all<- merge(AMR.chem.all, AMR.chem.allA, by="Filename")
AMR.chem.all[,c(1:10,18)]<- lapply(AMR.chem.all[,c(1:10,18)], factor)

AMR.spl.all<- merge(AMR.spl.all, AMR.chem.allN, by="Filename")
AMR.spl.all<- merge(AMR.spl.all, AMR.chem.allA, by="Filename")
AMR.spl.all[,c(1:10)]<- lapply(AMR.spl.all[,c(1:10)], factor)

rm(list=c("AMR.chem.allN", "AMR.chem.allA"))



### Same steps for Within-AMR dataset (hereafter: C)  ----
AMR.spl.C<- subset(AMR.spl, !is.na(AMR_C3))
AMR.spl.C[,c(1:10)]<- lapply(AMR.spl.C[,c(1:10)], factor)

AMR.spl.C<- mutate(AMR.spl.C, z.Seasonality=as.vector(scale(Seasonality)))
AMR.spl.C[,c(1:10)]<- lapply(AMR.spl.C[,c(1:10)], factor)

AMR.chem.C<- merge(AMR.spl.C, AMR.chem, by="Filename")
AMR.chem.C[,c(1:10,18)]<- lapply(AMR.chem.C[,c(1:10,18)], factor)
AMR.chem.CN<- data.frame(AMR.chem.C %>% group_by(Filename) %>% reframe(Ncomp=n()))
AMR.chem.CA<- data.frame(AMR.chem.C %>% group_by(Filename) %>% reframe(SArea=sum(Area)))

AMR.chem.C<- merge(AMR.chem.C, AMR.chem.CN, by="Filename")
AMR.chem.C<- merge(AMR.chem.C, AMR.chem.CA, by="Filename")
AMR.chem.C[,c(1:10,18)]<- lapply(AMR.chem.C[,c(1:10,18)], factor)

AMR.spl.C<- merge(AMR.spl.C, AMR.chem.CN, by="Filename")
AMR.spl.C<- merge(AMR.spl.C, AMR.chem.CA, by="Filename")
AMR.spl.C[,c(1:10)]<- lapply(AMR.spl.C[,c(1:10)], factor)

rm(list=c("AMR.chem.CN", "AMR.chem.CA"))





## Descriptive stats (used for Suppl.Table S1) ----

# number of compounds per sample
range(AMR.spl.all$Ncomp)
mean(AMR.spl.all$Ncomp)
sd(AMR.spl.all$Ncomp)

range(AMR.spl.C$Ncomp)
mean(AMR.spl.C$Ncomp)
sd(AMR.spl.C$Ncomp)

# number of samples per comp
AMR.chem.all %>% count(Comp, name="Nspl")

# number of samples per MRC
AMR.spl.all %>% count(MRC, name="Nmrc")

AMR.spl.all.mrc.AMR.male<- data.frame(AMR.spl.all %>% group_by(MRC, AMR, Male) %>% tally())
arrange(AMR.spl.all.mrc.AMR.male, AMR, Male)

AMR.spl.all.mrc.male<- data.frame(AMR.spl.all %>% group_by(MRC, Male) %>% tally())
arrange(AMR.spl.all.mrc.male, Male)

# number of samples per male
AMR.spl.all.male<-AMR.spl.all %>% count(Male, name="Nspl")
head(AMR.spl.all.male)
range(AMR.spl.all.male$Nspl)
mean(AMR.spl.all.male$Nspl)
sd(AMR.spl.all.male$Nspl)

AMR.spl.male<- data.frame(AMR.spl.all %>% group_by(Male) %>% tally())
mean(AMR.spl.male$n)
sd(AMR.spl.male$n)

# Number of samples per male rank cat and period
data.frame(AMR.spl.all %>% group_by(MRC, AMR_A3) %>% tally())

data.frame(AMR.spl.C %>% group_by(MRC, AMR_C3) %>% tally())

# male age during AMR
range(AMR.spl.C[which(AMR.spl.C$MRC=="Ascending_alpha"),]$Age)
range(AMR.spl.C[which(AMR.spl.C$MRC=="Descending_alpha"),]$Age)
range(AMR.spl.C[which(AMR.spl.C$MRC=="Stable_subordinate"),]$Age)

range(AMR.spl.C[which(AMR.spl.C$MRC=="Ascending_alpha"|AMR.spl.C$MRC=="Descending_alpha"),]$Age)
mean(AMR.spl.C[which(AMR.spl.C$MRC=="Ascending_alpha"|AMR.spl.C$MRC=="Descending_alpha"),]$Age)
sd(AMR.spl.C[which(AMR.spl.C$MRC=="Ascending_alpha"|AMR.spl.C$MRC=="Descending_alpha"),]$Age)


### Hormonal data (hereafter: T) per male ----
AMR.spl.all.T<- subset(AMR.spl.all, !is.na(T.conc))
AMR.spl.all.T[,c(1:10)]<- lapply(AMR.spl.all.T[,c(1:10)], factor)

AMR.spl.all.T.male<-data.frame(AMR.spl.all.T %>% group_by(Male) %>% tally())
range(AMR.spl.all.T.male$n)
mean(AMR.spl.all.T.male$n)
sd(AMR.spl.all.T.male$n)

### Behavior data (hereafter: B) per male ----
AMR.spl.all.B<- subset(AMR.spl.all, !is.na(UW_rate))
AMR.spl.all.B[,c(1:10)]<- lapply(AMR.spl.all.B[,c(1:10)], factor)

AMR.spl.all.B.male<-data.frame(AMR.spl.all.B %>% group_by(Male) %>% tally())
range(AMR.spl.all.B.male$n)
mean(AMR.spl.all.B.male$n)
sd(AMR.spl.all.B.male$n)




## VEGAN package: RPA standardization, vegdist and MDS per category, Shannon index ----
# full dataset (all)
AMR.chem.all<- AMR.chem.all %>% arrange(Filename, Area)
AMR.chem.RPA.all<- data.frame(AMR.chem.all %>% group_by(Filename) %>% reframe(Area=Area, RPA=round((Area/sum(Area))*100, digits=4), logRPA=round(log(((Area/sum(Area))*100)+1), digits=4), asinRPA=round((log(asin(sqrt(((Area/sum(Area))*100)/100))+0.01)), digits=4)))
AMR.chem.RPA.all<- AMR.chem.RPA.all %>% arrange(Filename, Area)
AMR.chem.std.all<- cbind(AMR.chem.all, AMR.chem.RPA.all[,-(1:3)])
AMR.chem.mat.all<- xtabs(logRPA~Filename+Comp, AMR.chem.std.all)
AMR.chem.vdim.all<-vegdist(AMR.chem.mat.all, method="bray")
AMR.chem.nmds.all<- metaMDS(AMR.chem.vdim.all, k=3, autotransform=F, noshare=F, wascores=F, expand=T, trymax=200, trace=F)
AMR.chem.nmds.sco.all<- scores(AMR.chem.nmds.all)
AMR.chem.nmds.all$stress
Shannon<- as.vector(diversity(AMR.chem.mat.all))
AMR.spl.all<- AMR.spl.all %>% arrange(Filename)
AMR.chem.res.all<- cbind(AMR.spl.all, AMR.chem.nmds.sco.all, Shannon)


# Within-AMR dataset (C)
AMR.chem.C<- AMR.chem.C %>% arrange(Filename, Area)
AMR.chem.RPA.C<- data.frame(AMR.chem.C %>% group_by(Filename) %>% reframe(Area=Area, RPA=round((Area/sum(Area))*100, digits=4), logRPA=round(log(((Area/sum(Area))*100)+1), digits=4), asinRPA=round((log(asin(sqrt(((Area/sum(Area))*100)/100))+0.01)), digits=4)))
AMR.chem.RPA.C<- AMR.chem.RPA.C %>% arrange(Filename, Area)
AMR.chem.std.C<- cbind(AMR.chem.C, AMR.chem.RPA.C[,-(1:3)])
AMR.chem.mat.C<- xtabs(logRPA~Filename+Comp, AMR.chem.std.C)
AMR.chem.vdim.C<-vegdist(AMR.chem.mat.C, method="bray")
AMR.chem.nmds.C<- metaMDS(AMR.chem.vdim.C, k=3, autotransform=F, noshare=F, wascores=F, expand=T, trymax=200, trace=F)
AMR.chem.nmds.sco.C<- vegan::scores(AMR.chem.nmds.C)
AMR.chem.nmds.C$stress
Shannon<- as.vector(diversity(AMR.chem.mat.C))
AMR.spl.C<- AMR.spl.C %>% arrange(Filename)
AMR.chem.res.C<- cbind(AMR.spl.C, AMR.chem.nmds.sco.C, Shannon)






# 2) Start STATS here ----

# 2a) Random Slope Model approach  ----
## Check structure ----
vif(lm(asinRPA ~ MRC+AMR_A3+z.Seasonality+Group, data=AMR.chem.std.all))
# OK, group = 2.56
vif(lm(asinRPA ~ MRC+AMR_C3+z.Seasonality+Group, data=AMR.chem.std.C))
# OK, group = 1.59

xtabs(~AMR+MRC, AMR.chem.std.all)                # ok
xtabs(~AMR+AMR_A3, AMR.chem.std.all)             # ok
xtabs(~AMR+Group, AMR.chem.std.all)              # many 0
xtabs(~Male+AMR_A3, AMR.chem.std.all)            # many 0
xtabs(~Male+Group, AMR.chem.std.all)             # many 0
xtabs(~Batch+MRC, AMR.chem.std.all)              # ok
xtabs(~Batch+AMR_A3, AMR.chem.std.all)           # ok
xtabs(~Batch+Group, AMR.chem.std.all)            # ok
xtabs(~Handler+MRC, AMR.chem.std.all)            # ok
xtabs(~Handler+AMR_A3, AMR.chem.std.all)         # ok
xtabs(~Handler+Group, AMR.chem.std.all)          # many 0

xtabs(~AMR+MRC, AMR.chem.std.C)                # ok
xtabs(~AMR+AMR_C3, AMR.chem.std.C)             # ok
xtabs(~AMR+Group, AMR.chem.std.C)              # many 0
xtabs(~Male+AMR_C3, AMR.chem.std.C)            # ok
xtabs(~Male+Group, AMR.chem.std.C)             # many 0
xtabs(~Batch+MRC, AMR.chem.std.C)              # ok
xtabs(~Batch+AMR_C3, AMR.chem.std.C)           # ok
xtabs(~Batch+Group, AMR.chem.std.C)            # ok
xtabs(~Handler+MRC, AMR.chem.std.C)            # ok
xtabs(~Handler+AMR_C3, AMR.chem.std.C)         # ok
xtabs(~Handler+Group, AMR.chem.std.C)          # many 0



### RSM Full dataset ----
# null model
Null.mod.all=lmer(asinRPA ~ MRC*AMR_A3 +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Handler) +(1|Batch) +(1|Date) +(1|Filename)
                  +(1+z.Seasonality||Comp) +(1|Comp:Group) 
                  +(1|AMR:MRC*AMR_A3) +(1|Handler:MRC*AMR_A3) +(1|Batch:MRC*AMR_A3) +(1|Batch:Group),
                  data=AMR.chem.std.all, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Null.mod.all)

# full model
Full.mod.all=lmer(asinRPA ~ MRC*AMR_A3 +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Handler) +(1|Batch) +(1|Date) +(1|Filename)
                  +(1|Comp:MRC*AMR_A3) +(1+z.Seasonality||Comp) +(1|Comp:Group)
                  +(1|AMR:MRC*AMR_A3) +(1|Handler:MRC*AMR_A3) +(1|Batch:MRC*AMR_A3) +(1|Batch:Group),
                  data=AMR.chem.std.all, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Full.mod.all)
anova(Null.mod.all, Full.mod.all, test="Chisq") #***

# reduced no MRC
Red.mod.all.mrc=lmer(asinRPA ~ MRC*AMR_A3 +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Handler) +(1|Batch) +(1|Date) +(1|Filename)
                     +(1|Comp:AMR_A3) +(1+z.Seasonality||Comp) +(1|Comp:Group)
                     +(1|AMR:MRC*AMR_A3) +(1|Handler:MRC*AMR_A3) +(1|Batch:MRC*AMR_A3) +(1|Batch:Group),
                     data=AMR.chem.std.all, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Red.mod.all.mrc)
anova(Red.mod.all.mrc, Full.mod.all, test="Chisq") #NS

# reduced no AMR timing
Red.mod.all.tim=lmer(asinRPA ~ MRC*AMR_A3 +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Handler) +(1|Batch) +(1|Date) +(1|Filename)
                     +(1|Comp:MRC) +(1+z.Seasonality||Comp) +(1|Comp:Group)
                     +(1|AMR:MRC*AMR_A3) +(1|Handler:MRC*AMR_A3) +(1|Batch:MRC*AMR_A3) +(1|Batch:Group),
                     data=AMR.chem.std.all, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Red.mod.all.tim)
anova(Red.mod.all.tim, Full.mod.all, test="Chisq") #***



### RSM Within-AMR dataset ----
# null model
Null.mod.C=lmer(asinRPA ~ MRC*AMR_C3 +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Handler) +(1|Batch) +(1|Date) +(1|Filename)
                +(1+z.Seasonality||Comp) +(1|Comp:Group) 
                +(1|AMR:MRC*AMR_C3) +(1|Handler:MRC*AMR_C3) +(1|Batch:MRC*AMR_C3) +(1|Batch:Group),
                data=AMR.chem.std.C, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Null.mod.C)

# full model
Full.mod.C=lmer(asinRPA ~ MRC*AMR_C3 +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Handler) +(1|Batch) +(1|Date) +(1|Filename)
                +(1|Comp:MRC*AMR_C3) +(1+z.Seasonality||Comp) +(1|Comp:Group)
                +(1|AMR:MRC*AMR_C3) +(1|Handler:MRC*AMR_C3) +(1|Batch:MRC*AMR_C3) +(1|Batch:Group),
                data=AMR.chem.std.C, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Full.mod.C)
anova(Null.mod.C, Full.mod.C, test="Chisq") #***

# reduced no MRC
Red.mod.C.mrc=lmer(asinRPA ~ MRC*AMR_C3 +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Handler) +(1|Batch) +(1|Date) +(1|Filename)
                   +(1|Comp:AMR_C3) +(1+z.Seasonality||Comp) +(1|Comp:Group)
                   +(1|AMR:MRC*AMR_C3) +(1|Handler:MRC*AMR_C3) +(1|Batch:MRC*AMR_C3) +(1|Batch:Group),
                   data=AMR.chem.std.C, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Red.mod.C.mrc)
anova(Red.mod.C.mrc, Full.mod.C, test="Chisq") #NS

# reduced no AMR timing
Red.mod.C.tim=lmer(asinRPA ~ MRC*AMR_C3 +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Handler) +(1|Batch) +(1|Date) +(1|Filename)
                   +(1|Comp:MRC) +(1+z.Seasonality||Comp) +(1|Comp:Group)
                   +(1|AMR:MRC*AMR_C3) +(1|Handler:MRC*AMR_C3) +(1|Batch:MRC*AMR_C3) +(1|Batch:Group),
                   data=AMR.chem.std.C, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Red.mod.C.tim)
anova(Red.mod.C.tim, Full.mod.C, test="Chisq") #***





# 2b) Diversity models ----
## Plot tools ----
MRC.order.all<-c("Ascending_alpha","Descending_alpha","Stable_subordinate")

AMR.chem.res.all<- transform(AMR.chem.res.all, MRC=factor(MRC, levels=MRC.order.all))
AMR.chem.res.C<- transform(AMR.chem.res.C, MRC=factor(MRC, levels=MRC.order.all))

levels(AMR.chem.res.all$MRC) <- c("Ascending alphas", "Descending alphas", "Stable subordinates")
levels(AMR.chem.res.C$MRC) <- c("Ascending alphas", "Descending alphas", "Stable subordinates")



## Richness  ----
vif(lm(Ncomp ~ MRC+AMR_A3+Group, data=AMR.chem.res.all))
# ok, vif group ID = 2.69
vif(lm(Ncomp ~ MRC+AMR_C3+Group, data=AMR.chem.res.C))
# ok, vif group ID = 1.62


### Full dataset ----
# Null model
ModRic.allN=glmmTMB(Ncomp ~ z.Seasonality +Group +(1|AMR) +(1|Date) +(1|Male) +(1|Handler) +(1|Batch),
                   data=AMR.chem.res.all, family=poisson(link="log"), control=glmmTMBControl())
summary(ModRic.allN)

# Full model
ModRic.all=glmmTMB(Ncomp ~ MRC*AMR_A3 +z.Seasonality +Group +(1|AMR) +(1|Date) +(1|Male) +(1|Handler) +(1|Batch),
                    data=AMR.chem.res.all, family=poisson(link="log"), control=glmmTMBControl())
summary(ModRic.all)
anova(ModRic.allN, ModRic.all, test="Chisq") #NS

emm.ModRic.all<- data.frame(emmeans(ModRic.all, list(pairwise ~ MRC*AMR_A3), adjust="tukey"))
emmeans(ModRic.all, list(pairwise ~ MRC*AMR_A3), adjust="tukey")


### Within-AMR dataset ----
# Null model
ModRic.C3.N=glmmTMB(Ncomp ~ z.Seasonality +Group +(1|AMR) +(1|Date) +(1|Male) +(1|Handler) +(1|Batch),
                    data=AMR.chem.res.C, family=poisson(link="log"), control=glmmTMBControl())
summary(ModRic.C3.N)

# Full model
ModRic.C3=glmmTMB(Ncomp ~ MRC*AMR_C3 +z.Seasonality +Group +(1|AMR) +(1|Date) +(1|Male) +(1|Handler) +(1|Batch),
                   data=AMR.chem.res.C, family=poisson(link="log"), control=glmmTMBControl())
summary(ModRic.C3)
anova(ModRic.C3.N, ModRic.C3, test="Chisq") #NS

emm.ModRic.C3<- data.frame(emmeans(ModRic.C3, list(pairwise ~ MRC*AMR_C3), adjust="tukey"))
emmeans(ModRic.C3, list(pairwise ~ MRC*AMR_C3), adjust="tukey")




## Total intensity  ----
### Full dataset ----
# Null model
ModSar.allN=lmer(log(SArea) ~ z.Seasonality +Group +(1|AMR) +(1|Date) +(1|Male) +(1|Handler) +(1|Batch),
                 data=AMR.chem.res.all, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))
summary(ModSar.allN) 

# Full model
ModSar.all=lmer(log(SArea) ~ MRC*AMR_A3 +z.Seasonality +Group +(1|AMR) +(1|Date) +(1|Male) +(1|Handler) +(1|Batch),
                   data=AMR.chem.res.all, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))
summary(ModSar.all) 
anova(ModSar.allN, ModSar.all, test="Chisq") #**

emm.ModSar.all<- data.frame(emmeans(ModSar.all, list(pairwise ~ MRC*AMR_A3), adjust="tukey"))
emmeans(ModSar.all, list(pairwise ~ MRC*AMR_A3), adjust="tukey")


### Within-AMR dataset ----
# Null model
ModSar.C3.N=lmer(log(SArea) ~ z.Seasonality +Group +(1|AMR) +(1|Date) +(1|Male) +(1|Handler) +(1|Batch),
                 data=AMR.chem.res.C, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))
summary(ModSar.C3.N) # seasonality

# Full model
ModSar.C3=lmer(log(SArea) ~ MRC*AMR_C3 +z.Seasonality +Group +(1|AMR) +(1|Date) +(1|Male) +(1|Handler) +(1|Batch),
                   data=AMR.chem.res.C, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))
summary(ModSar.C3) 
anova(ModSar.C3.N, ModSar.C3, test="Chisq") #NS

emm.ModSar.C3<- data.frame(emmeans(ModSar.C3, list(pairwise ~ MRC*AMR_C3), adjust="tukey"))
emmeans(ModSar.C3, list(pairwise ~ MRC*AMR_C3), adjust="tukey")



## Shannon index ----
### Full dataset ----
# Null model
ModSha.allN=glmer(Shannon ~ z.Seasonality +Group +(1|AMR) +(1|Date) +(1|Male) +(1|Handler) +(1|Batch),
                  data=AMR.chem.res.all, family=Gamma(link='identity'),
                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))
summary(ModSha.allN) #

# Full model
ModSha.all=glmer(Shannon ~ MRC*AMR_A3 +z.Seasonality +Group +(1|AMR) +(1|Date) +(1|Male) +(1|Handler) +(1|Batch),
                data=AMR.chem.res.all, family=Gamma(link='identity'),
                control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))
summary(ModSha.all) # 
anova(ModSha.allN, ModSha.all, test="Chisq") #NS

emm.ModSha.all<- data.frame(emmeans(ModSha.all, list(pairwise ~ MRC*AMR_A3), adjust="tukey"))
emmeans(ModSha.all, list(pairwise ~ MRC*AMR_A3), adjust="tukey")


### Within-AMR dataset ----
# Null model
ModSha.C3.N=glmer(Shannon ~ z.Seasonality +Group +(1|AMR) +(1|Date) +(1|Male) +(1|Handler) +(1|Batch),
                  data=AMR.chem.res.C, family=Gamma(link='identity'),
                  control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))
summary(ModSha.C3.N) # ns

# Full model
ModSha.C3=glmer(Shannon ~ MRC*AMR_C3 +z.Seasonality +Group +(1|AMR) +(1|Date) +(1|Male) +(1|Handler) +(1|Batch),
                 data=AMR.chem.res.C, family=Gamma(link='identity'),
                 control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F))
summary(ModSha.C3) 
anova(ModSha.C3.N, ModSha.C3, test="Chisq") #*

emm.ModSha.C3<- data.frame(emmeans(ModSha.C3, list(pairwise ~ MRC*AMR_C3), adjust="tukey"))
emmeans(ModSha.C3, list(pairwise ~ MRC*AMR_C3), adjust="tukey")




## Figure 2: Plots diversity indices ----
pRic.all<-
  plot_model(ModRic.all, type="pred", terms=c("AMR_A3","MRC"), bias_correction=T, 
             colors=c("#F8766D","blue","#00BA38"), title="a. Richness",
             axis.title=c("AMR period (full dataset)","Number of compounds"), legend.title="Male rank category",
             show.legend=F, line.size=1.5, dot.size=5, dodge=.7)

pRic.C3<-
  plot_model(ModRic.C3, type="pred", terms=c("AMR_C3","MRC"), bias_correction=T, 
             colors=c("#F8766D","blue","#00BA38"), title="d. Richness",
             axis.title=c("AMR period (within-AMR dataset)","Number of compounds"), legend.title="Male rank category",
             show.legend=F, line.size=1.5, dot.size=5, dodge=.7)

pSar.all<-
  plot_model(ModSar.all, type="pred", terms=c("AMR_A3","MRC"), bias_correction=T, 
             colors=c("#F8766D","blue","#00BA38"), title="b. Total intensity",
             axis.title=c("AMR period (full dataset)","log(Sum of peak areas)"), legend.title="Male rank category",
             show.legend=F, line.size=1.5, dot.size=5, dodge=.7)

pSar.C3<-
  plot_model(ModSar.C3, type="pred", terms=c("AMR_C3","MRC"), bias_correction=T, 
             colors=c("#F8766D","blue","#00BA38"), title="e. Total intensity",
             axis.title=c("AMR period (within-AMR dataset)","log(Sum of peak areas)"), legend.title="Male rank category",
             show.legend=F, line.size=1.5, dot.size=5, dodge=.7)

pSha.all<-
  plot_model(ModSha.all, type="pred", terms=c("AMR_A3","MRC"), bias_correction=T, 
             colors=c("#F8766D","blue","#00BA38"), title="c. Shannon index",
             axis.title=c("AMR period (full dataset)","Shannon's H"), legend.title="Male rank category",
             show.legend=T, line.size=1.5, dot.size=5, dodge=.7)

pSha.C3<-
  plot_model(ModSha.C3, type="pred", terms=c("AMR_C3","MRC"), bias_correction=T, 
             colors=c("#F8766D","blue","#00BA38"), title="f. Shannon index",
             axis.title=c("AMR period (within-AMR dataset)","Shannon's H"), legend.title="Male rank category",
             show.legend=T, line.size=1.5, dot.size=5, dodge=.7)


#dev.new()
set_theme(base = theme_classic(base_size=10))

# Plot Figure 2
grid.draw(rbind(cbind(ggplotGrob(pRic.all), ggplotGrob(pSar.all), ggplotGrob(pSha.all), size="first"),
                cbind(ggplotGrob(pRic.C3), ggplotGrob(pSar.C3), ggplotGrob(pSha.C3), size="first")))





# 3) Compounds most influenced by predictors ----

# 3a) RSM slope estimates ----

## Full dataset ----
### MRC: not run since this predictor is non-signif in the RSM (section 2a)

### AMR timing ----
# 0+ then 1| necessary to have all three slopes, otherwise only 2 are shown 
ranef.mod.all.tim=lmer(asinRPA ~ AMR_A3 +(0+AMR_A3|Comp) +(1|Comp) +(1|AMR) +(1|Date) +(1|Male) +(1|Handler) +(1|Batch) +(1|Filename),
                          data=AMR.chem.std.all, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(ranef.mod.all.tim) #

ranef.table.all.tim=round(ranef(ranef.mod.all.tim)$Comp, 2)
ranef.table.all.tim$Comp<- row.names(ranef.table.all.tim)
head(ranef.table.all.tim)

colnames(ranef.table.all.tim)<- c("Intercept","SlopeBEF","SlopeDUR","SlopeAFT","Comp")
ranef.slopes.all.tim<- mutate(ranef.table.all.tim, Abs.slopeBEF=abs(SlopeBEF),Abs.slopeDUR=abs(SlopeDUR), Abs.slopeAFT=abs(SlopeAFT), 
                                  Sum.abs.slope=Abs.slopeBEF+Abs.slopeDUR+Abs.slopeAFT,
                                  BEF_DUR=abs(SlopeBEF-SlopeDUR), BEF_AFT=abs(SlopeBEF-SlopeAFT), DUR_AFT=abs(SlopeDUR-SlopeAFT), 
                                  Abs.sum.diff=BEF_DUR+BEF_AFT+DUR_AFT)
ranef.slopes.all.tim<- ranef.slopes.all.tim %>% arrange(desc(Abs.sum.diff))

ranef.slopesDiff.all.tim<- c(ranef.slopes.all.tim$BEF_DUR, ranef.slopes.all.tim$BEF_AFT, ranef.slopes.all.tim$DUR_AFT)
ranef.threshold.all.tim<- mean(ranef.slopesDiff.all.tim)+2*sd(ranef.slopesDiff.all.tim)
ranef.threshold.all.tim

ranef.infl.all.tim<- ranef.slopes.all.tim %>% filter((BEF_DUR>=ranef.threshold.all.tim |BEF_AFT>=ranef.threshold.all.tim 
                                                      |DUR_AFT>=ranef.threshold.all.tim) %>% replace_na(T))
nrow(ranef.infl.all.tim) #13
ranef.infl.all.tim


## Within-AMR dataset ----
### MRC: not run since this predictor was NS in the RSM (section 2a)

### AMR timing ----
# 0+ then 1| necessary to have all three slopes, otherwise only 2 are shown 
ranef.mod.C.tim=lmer(asinRPA ~ AMR_C3 +(0+AMR_C3|Comp) +(1|Comp) +(1|AMR) +(1|Date) +(1|Male) +(1|Handler) +(1|Batch) +(1|Filename),
                       data=AMR.chem.std.C, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(ranef.mod.C.tim) #

ranef.table.C.tim=round(ranef(ranef.mod.C.tim)$Comp, 2)
ranef.table.C.tim$Comp<- row.names(ranef.table.C.tim)
head(ranef.table.C.tim)

colnames(ranef.table.C.tim)<- c("Intercept","SlopeSTA","SlopeINS","SlopeEND","Comp")
ranef.slopes.C.tim<- mutate(ranef.table.C.tim, Abs.slopeSTA=abs(SlopeSTA),Abs.slopeINS=abs(SlopeINS), Abs.slopeEND=abs(SlopeEND), 
                              Sum.abs.slope=Abs.slopeSTA+Abs.slopeINS+Abs.slopeEND,
                              STA_INS=abs(SlopeSTA-SlopeINS), STA_END=abs(SlopeSTA-SlopeEND), INS_END=abs(SlopeINS-SlopeEND), 
                              Abs.sum.diff=STA_INS+STA_END+INS_END)
ranef.slopes.C.tim<- ranef.slopes.C.tim %>% arrange(desc(Abs.sum.diff))

ranef.slopesDiff.C.tim<- c(ranef.slopes.C.tim$STA_INS, ranef.slopes.C.tim$STA_END, ranef.slopes.C.tim$INS_END)
ranef.threshold.C.tim<- mean(ranef.slopesDiff.C.tim)+2*sd(ranef.slopesDiff.C.tim)
ranef.threshold.C.tim

ranef.infl.C.tim<- ranef.slopes.C.tim %>% filter((STA_INS>=ranef.threshold.C.tim |STA_END>=ranef.threshold.C.tim 
                                                      |INS_END>=ranef.threshold.C.tim) %>% replace_na(T))
nrow(ranef.infl.C.tim) #12
ranef.infl.C.tim




# 3b) Random Forest ----
### Full dataset - AMR period ----
AMR.chem.res.all<-arrange(AMR.chem.res.all, Filename)

RF.mod.all.tim<-randomForest::tuneRF(x=AMR.chem.mat.all, y=AMR.chem.res.all$AMR_A3, 
                                         strata=AMR.chem.res.all$Male, mtryStart=8, ntreeTry=20000, stepFactor=1, improve=0.0001, trace=FALSE, plot=FALSE, doBest=TRUE)
RF.mod.all.tim

RF.mod.all.tim.imp<- data.frame(rownames(RF.mod.all.tim$importance), RF.mod.all.tim$importance) 
colnames(RF.mod.all.tim.imp)=c("Peak", "MeanDecreaseGini")
Gini.thresh.all.tim<- mean(RF.mod.all.tim.imp$MeanDecreaseGini)+2*sd(RF.mod.all.tim.imp$MeanDecreaseGini)
Gini.infl.all.tim<- RF.mod.all.tim.imp[RF.mod.all.tim.imp$MeanDecreaseGini>=Gini.thresh.all.tim,]
nrow(Gini.infl.all.tim) # 6 comps
arrange(Gini.infl.all.tim, -MeanDecreaseGini)


### Full dataset - MRC ----
RF.mod.all.mal<-randomForest::tuneRF(x=AMR.chem.mat.all, y=AMR.chem.res.all$MRC, 
                                         strata=AMR.chem.res.all$Male, mtryStart=8, ntreeTry=20000, stepFactor=1, improve=0.0001, trace=FALSE, plot=FALSE, doBest=TRUE)
RF.mod.all.mal

RF.mod.all.mal.imp<- data.frame(rownames(RF.mod.all.mal$importance), RF.mod.all.mal$importance) 
colnames(RF.mod.all.mal.imp)=c("Peak", "MeanDecreaseGini")
Gini.thresh.all.mal<- mean(RF.mod.all.mal.imp$MeanDecreaseGini)+2*sd(RF.mod.all.mal.imp$MeanDecreaseGini)
Gini.infl.all.mal<- RF.mod.all.mal.imp[RF.mod.all.mal.imp$MeanDecreaseGini>=Gini.thresh.all.mal,]
nrow(Gini.infl.all.mal) # 6 comps
arrange(Gini.infl.all.mal, -MeanDecreaseGini)


### Within-AMR dataset - AMR period ----
AMR.chem.res.C<-arrange(AMR.chem.res.C, Filename)

RF.mod.C.tim<-randomForest::tuneRF(x=AMR.chem.mat.C, y=AMR.chem.res.C$AMR_C3, 
                                        strata=AMR.chem.res.C$Male, mtryStart=8, ntreeTry=20000, stepFactor=1, improve=0.0001, trace=FALSE, plot=FALSE, doBest=TRUE)
RF.mod.C.tim

RF.mod.C.tim.imp<- data.frame(rownames(RF.mod.C.tim$importance), RF.mod.C.tim$importance) 
colnames(RF.mod.C.tim.imp)=c("Peak", "MeanDecreaseGini")
Gini.thresh.C.tim<- mean(RF.mod.C.tim.imp$MeanDecreaseGini)+2*sd(RF.mod.C.tim.imp$MeanDecreaseGini)
Gini.infl.C.tim<- RF.mod.C.tim.imp[RF.mod.C.tim.imp$MeanDecreaseGini>=Gini.thresh.C.tim,]
nrow(Gini.infl.C.tim) # 3 comps
arrange(Gini.infl.C.tim, -MeanDecreaseGini)


### Within-AMR dataset - MRC ----
RF.mod.C.mal<-randomForest::tuneRF(x=AMR.chem.mat.C, y=AMR.chem.res.C$MRC, 
                                        strata=AMR.chem.res.C$Male, mtryStart=8, ntreeTry=20000, stepFactor=1, improve=0.0001, trace=FALSE, plot=FALSE, doBest=TRUE)
RF.mod.C.mal

RF.mod.C.mal.imp<- data.frame(rownames(RF.mod.C.mal$importance), RF.mod.C.mal$importance) 
colnames(RF.mod.C.mal.imp)=c("Peak", "MeanDecreaseGini")
Gini.thresh.C.mal<- mean(RF.mod.C.mal.imp$MeanDecreaseGini)+2*sd(RF.mod.C.mal.imp$MeanDecreaseGini)
Gini.infl.C.mal<- RF.mod.C.mal.imp[RF.mod.C.mal.imp$MeanDecreaseGini>=Gini.thresh.C.mal,]
nrow(Gini.infl.C.mal) # 5 comps
arrange(Gini.infl.C.mal, -MeanDecreaseGini)




# 3c) Presence/absence ----
AMR.chem.std.all.comp<-data.frame(AMR.chem.std.all %>% group_by(Comp) %>% summarise(MAreaC=mean(Area)))
AMR.chem.std.all.comp<-arrange(AMR.chem.std.all.comp, -MAreaC)
head(AMR.chem.std.all.comp)

# Keep only the 75 most abundant compounds (MArea>1% highest MArea)
AMR.chem.std.all.comp75<- subset(AMR.chem.std.all.comp[c(1:75),]) 
AMR.chem.std.all.comp75$Comp<- as.factor(AMR.chem.std.all.comp75$Comp)
head(AMR.chem.std.all.comp75)

Plist<- as.factor(sort(as.vector(AMR.chem.std.all.comp75$Comp)))





# 4) Single compound differences ----

### Keep only comp selected by RF, slopes and manually (n = 38)  ----
Pselec38.pk<- c("M043.P03","M045.P06","M045.P07","M055.P05","M057.P05",
              "M059.P04.3","M071.P01","M074.P01","N074.P02","N091.P01",
              "N108.P01","Q069.P02","Q070.P04","Q097.P03","Q105.P04",
              "Q125.P01","Q128.P02","R042.P02","R057.P06","R135.P01",
              "R137.P01","S055.P06","S058.P11","S084.P01","S085.P03",
              "S113.P01","S128.P01","S133.P01","S133.P05","S136.P01",
              "S138.P01","S149.P01","S150.P01","S165.P04","V043.P07",
              "V045.P04","W117.P03","W137.P01")

Pselec38.num<- c("C002","C014","C017","C022","C016",
                "C025","C001","C010","C037","C032",
                "C036","C064","C066","C063","C062",
                "C054","C051","C073","C076","C074",
                "C069","C097","C098","C079","C090",
                "C077","C089","C078","C083","C085",
                "C086","C088","C082","C096","C012",
                "C011","C070","C072")

Pselec38.names<- c("Acetic acid","1,2-Propanediol","Hexan-2-ol","Unidentified comp 13","2,3-Butanediol",
                   "N-Methylformamide","2-Methyl-3-buten-2-ol","Unidentified comp 06","Unidentified comp 18","Ethylbenzene",
                   "2,5-Dimethylpyrazine","Unidentified comp 29","Unidentified comp 31","1H-Pyrrole-2,5-dione","Acetophenone",
                   "Unidentified comp 25","Unk. nitrogen-containing comp 07","Unidentified comp 34","Unidentified comp 36","Unidentified comp 35",
                   "3,7-Dimethyl-2,3,3a,4,5,6-hexahydro-1-benzofuran","Unidentified comp 50","Unidentified comp 51","2-Cyclopentylcyclopentanone","Unidentified comp 44",
                   "3-Methyl-4-propylpyrrolidine-2,5-dione","Unidentified comp 43","3-Ethylacetophenone","4-Ethylacetophenone","4-Methoxy-3,5-dimethylaniline",
                   "Unidentified comp 40","Unidentified comp 42","3-(Diethylamino)phenol","Unidentified comp 49","Unidentified comp 07",
                   "4-Methylpentan-2-ol","Unidentified comp 33","3,4-Dihydroxyacetophenone")

Pselec38<- as.data.frame(cbind(Pselec38.pk, Pselec38.num, Pselec38.names))
names(Pselec38)<- c("Comp","Comp.order","Comp.name")


# 38 comp of interest
AMR.chem.std.all.P<- subset(AMR.chem.std.all, Comp %in% Pselec38.pk)
AMR.chem.std.all.P<- merge(AMR.chem.std.all.P, Pselec38, by="Comp")
AMR.chem.std.all.P[,c(1:11,25,26)]<- lapply(AMR.chem.std.all.P[,c(1:11,25,26)], factor)

AMR.chem.std.C.P<- subset(AMR.chem.std.C, Comp %in% Pselec38.pk)
AMR.chem.std.C.P<- merge(AMR.chem.std.C.P, Pselec38, by="Comp")
AMR.chem.std.C.P[,c(1:11,25,26)]<- lapply(AMR.chem.std.C.P[,c(1:11,25,26)], factor)



### Plot tools ----
MRC.order.all<-c("Ascending_alpha","Descending_alpha","Stable_subordinate")

AMR.chem.std.all.P<- transform(AMR.chem.std.all.P, MRC=factor(MRC, levels=MRC.order.all))
AMR.chem.std.C.P<- transform(AMR.chem.std.C.P, MRC=factor(MRC, levels=MRC.order.all))

levels(AMR.chem.std.all.P$MRC) <- c("AA", "DA", "SS")
levels(AMR.chem.std.C.P$MRC) <- c("AA", "DA", "SS")



##Individual compounds models + plots ----
### Full dataset ----
AMR.chem.std.all.P38<- split(AMR.chem.std.all.P, AMR.chem.std.all.P$Comp.order)
head(AMR.chem.std.all.P38[1]) #list


#
AMR.chem.std.all.01.butenol<- data.frame(AMR.chem.std.all.P38[1])
AMR.chem.std.all.01.butenol[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.01.butenol[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.01.butenol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.01=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.01.butenol)
summary(glm.modA.01) 
Anova(glm.modA.01) # 0.544774
Aplot.01<-
  plot_model(glm.modA.01, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="2-Methyl-3-buten-2-ol", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.02.acetic<- data.frame(AMR.chem.std.all.P38[2])
AMR.chem.std.all.02.acetic[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.02.acetic[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.02.acetic)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.02=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.02.acetic)
summary(glm.modA.02) 
Anova(glm.modA.02) # 0.0001987
Aplot.02<-
  plot_model(glm.modA.02, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Acetic acid", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.03.unk<- data.frame(AMR.chem.std.all.P38[3])
AMR.chem.std.all.03.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.03.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.03.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.03=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.03.unk)
summary(glm.modA.03) 
Anova(glm.modA.03) # NA
Aplot.03<-
  plot_model(glm.modA.03, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 06", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.04.pentanol<- data.frame(AMR.chem.std.all.P38[4])
AMR.chem.std.all.04.pentanol[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.04.pentanol[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.04.pentanol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.04=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.04.pentanol)
summary(glm.modA.04) 
Anova(glm.modA.04) # 0.3857
Aplot.04<-
  plot_model(glm.modA.04, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="4-Methylpentan-2-ol", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.05.unk<- data.frame(AMR.chem.std.all.P38[5])
AMR.chem.std.all.05.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.05.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.05.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.05=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.05.unk)
summary(glm.modA.05) 
Anova(glm.modA.05) # 0.1466
Aplot.05<-
  plot_model(glm.modA.05, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 07", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.06.propanediol<- data.frame(AMR.chem.std.all.P38[6])
AMR.chem.std.all.06.propanediol[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.06.propanediol[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.06.propanediol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.06=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.06.propanediol)
summary(glm.modA.06) 
Anova(glm.modA.06) # 0.2304
Aplot.06<-
  plot_model(glm.modA.06, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="1,2-Propanediol", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.07.butanediol<- data.frame(AMR.chem.std.all.P38[7])
AMR.chem.std.all.07.butanediol[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.07.butanediol[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.07.butanediol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.07=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.07.butanediol)
summary(glm.modA.07) 
Anova(glm.modA.07) # 0.8756
Aplot.07<-
  plot_model(glm.modA.07, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="2,3-Butanediol", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.08.hexanol<- data.frame(AMR.chem.std.all.P38[8])
AMR.chem.std.all.08.hexanol[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.08.hexanol[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.08.hexanol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.08=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.08.hexanol)
summary(glm.modA.08) 
Anova(glm.modA.08) # 0.40969
Aplot.08<-
  plot_model(glm.modA.08, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Hexan-2-ol", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.09.unk<- data.frame(AMR.chem.std.all.P38[9])
AMR.chem.std.all.09.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.09.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.09.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.09=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.09.unk)
summary(glm.modA.09) 
Anova(glm.modA.09) # 0.77935
Aplot.09<-
  plot_model(glm.modA.09, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 13", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.10.formamide<- data.frame(AMR.chem.std.all.P38[10])
AMR.chem.std.all.10.formamide[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.10.formamide[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.10.formamide)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.10=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.10.formamide)
summary(glm.modA.10) 
Anova(glm.modA.10) # 0.08383
Aplot.10<-
  plot_model(glm.modA.10, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="N-Methylformamide", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.11.benzene<- data.frame(AMR.chem.std.all.P38[11])
AMR.chem.std.all.11.benzene[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.11.benzene[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.11.benzene)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.11=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.11.benzene)
summary(glm.modA.11) 
Anova(glm.modA.11) #0.26255
Aplot.11<-
  plot_model(glm.modA.11, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Ethylbenzene", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.12.pyrazine<- data.frame(AMR.chem.std.all.P38[12])
AMR.chem.std.all.12.pyrazine[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.12.pyrazine[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.12.pyrazine)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.12=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.12.pyrazine)
summary(glm.modA.12) 
Anova(glm.modA.12) # 0.51962
Aplot.12<-
  plot_model(glm.modA.12, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="2,5-Dimethylpyrazine", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.13.unk<- data.frame(AMR.chem.std.all.P38[13])
AMR.chem.std.all.13.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.13.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.13.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.13=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.13.unk)
summary(glm.modA.13) 
Anova(glm.modA.13) # 0.000002454
Aplot.13<-
  plot_model(glm.modA.13, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 18", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.14.unk<- data.frame(AMR.chem.std.all.P38[14])
AMR.chem.std.all.14.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.14.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.14.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.14=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.14.unk)
summary(glm.modA.14) 
Anova(glm.modA.14) #NA
Aplot.14<-
  plot_model(glm.modA.14, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unk. nitrogen-containing comp 07", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.15.unk<- data.frame(AMR.chem.std.all.P38[15])
AMR.chem.std.all.15.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.15.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.15.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.15=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.15.unk)
summary(glm.modA.15) 
Anova(glm.modA.15) # 0.3038
Aplot.15<-
  plot_model(glm.modA.15, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 25", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.16.acetophenone<- data.frame(AMR.chem.std.all.P38[16])
AMR.chem.std.all.16.acetophenone[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.16.acetophenone[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.16.acetophenone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.16=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.16.acetophenone)
summary(glm.modA.16) 
Anova(glm.modA.16) # NA
Aplot.16<-
  plot_model(glm.modA.16, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Acetophenone", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.17.pyrrole<- data.frame(AMR.chem.std.all.P38[17])
AMR.chem.std.all.17.pyrrole[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.17.pyrrole[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.17.pyrrole)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.17=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.17.pyrrole)
summary(glm.modA.17) 
Anova(glm.modA.17) # 0.6571
Aplot.17<-
  plot_model(glm.modA.17, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="1H-Pyrrole-2,5-dione", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.18.unk<- data.frame(AMR.chem.std.all.P38[18])
AMR.chem.std.all.18.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.18.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.18.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.18=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.18.unk)
summary(glm.modA.18) 
Anova(glm.modA.18) # 0.0267022
Aplot.18<-
  plot_model(glm.modA.18, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 29", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.19.unk<- data.frame(AMR.chem.std.all.P38[19])
AMR.chem.std.all.19.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.19.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.19.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.19=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.19.unk)
summary(glm.modA.19) 
Anova(glm.modA.19) # NA
Aplot.19<-
  plot_model(glm.modA.19, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 31", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.20.benzofuran<- data.frame(AMR.chem.std.all.P38[20])
AMR.chem.std.all.20.benzofuran[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.20.benzofuran[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.20.benzofuran)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.20=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.20.benzofuran)
summary(glm.modA.20) 
Anova(glm.modA.20) # 0.026995
Aplot.20<-
  plot_model(glm.modA.20, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="3,7-Dimethyl-(...)-hexahydro-1-benzofuran", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.21.unk<- data.frame(AMR.chem.std.all.P38[21])
AMR.chem.std.all.21.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.21.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.21.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.21=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.21.unk)
summary(glm.modA.21) 
Anova(glm.modA.21) # NA
Aplot.21<-
  plot_model(glm.modA.21, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 33", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.22.hydroxyacetophenone<- data.frame(AMR.chem.std.all.P38[22])
AMR.chem.std.all.22.hydroxyacetophenone[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.22.hydroxyacetophenone[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.22.hydroxyacetophenone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.22=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.22.hydroxyacetophenone)
summary(glm.modA.22) 
Anova(glm.modA.22) # 0.5552
Aplot.22<-
  plot_model(glm.modA.22, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="3,4-Dihydroxyacetophenone", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.23.unk<- data.frame(AMR.chem.std.all.P38[23])
AMR.chem.std.all.23.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.23.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.23.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.23=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.23.unk)
summary(glm.modA.23) 
Anova(glm.modA.23) # NA
Aplot.23<-
  plot_model(glm.modA.23, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 34", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.24.unk<- data.frame(AMR.chem.std.all.P38[24])
AMR.chem.std.all.24.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.24.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.24.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.24=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.24.unk)
summary(glm.modA.24) 
Anova(glm.modA.24) # 0.01429
Aplot.24<-
  plot_model(glm.modA.24, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 35", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.25.unk<- data.frame(AMR.chem.std.all.P38[25])
AMR.chem.std.all.25.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.25.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.25.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.25=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.25.unk)
summary(glm.modA.25) 
Anova(glm.modA.25) # NA
Aplot.25<-
  plot_model(glm.modA.25, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 36", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.26.pyrrolidine<- data.frame(AMR.chem.std.all.P38[26])
AMR.chem.std.all.26.pyrrolidine[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.26.pyrrolidine[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.26.pyrrolidine)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.26=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.26.pyrrolidine)
summary(glm.modA.26) 
Anova(glm.modA.26) # 0.73894
Aplot.26<-
  plot_model(glm.modA.26, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="3-Methyl-4-propylpyrrolidine-2,5-dione", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.27.3acetophenone<- data.frame(AMR.chem.std.all.P38[27])
AMR.chem.std.all.27.3acetophenone[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.27.3acetophenone[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.27.3acetophenone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.27=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.27.3acetophenone)
summary(glm.modA.27) 
Anova(glm.modA.27) # 0.1050
Aplot.27<-
  plot_model(glm.modA.27, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="3-Ethylacetophenone", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.28.cyclopentanone<- data.frame(AMR.chem.std.all.P38[28])
AMR.chem.std.all.28.cyclopentanone[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.28.cyclopentanone[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.28.cyclopentanone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.28=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.28.cyclopentanone)
summary(glm.modA.28) 
Anova(glm.modA.28) # 0.315385
Aplot.28<-
  plot_model(glm.modA.28, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="2-Cyclopentylcyclopentanone", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.29.aminophenol<- data.frame(AMR.chem.std.all.P38[29])
AMR.chem.std.all.29.aminophenol[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.29.aminophenol[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.29.aminophenol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.29=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.29.aminophenol)
summary(glm.modA.29) 
Anova(glm.modA.29) # 0.60441
Aplot.29<-
  plot_model(glm.modA.29, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="3-(Diethylamino)phenol", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.30.4acetophenone<- data.frame(AMR.chem.std.all.P38[30])
AMR.chem.std.all.30.4acetophenone[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.30.4acetophenone[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.30.4acetophenone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.30=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.30.4acetophenone)
summary(glm.modA.30) 
Anova(glm.modA.30) # 0.54061
Aplot.30<-
  plot_model(glm.modA.30, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="4-Ethylacetophenone", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.31.aniline<- data.frame(AMR.chem.std.all.P38[31])
AMR.chem.std.all.31.aniline[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.31.aniline[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.31.aniline)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.31=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.31.aniline)
summary(glm.modA.31) 
Anova(glm.modA.31) # 0.0382
Aplot.31<-
  plot_model(glm.modA.31, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="4-Methoxy-3,5-dimethylaniline", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.32.unk<- data.frame(AMR.chem.std.all.P38[32])
AMR.chem.std.all.32.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.32.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.32.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.32=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.32.unk)
summary(glm.modA.32) 
Anova(glm.modA.32) # 0.04116
Aplot.32<-
  plot_model(glm.modA.32, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 40", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.33.unk<- data.frame(AMR.chem.std.all.P38[33])
AMR.chem.std.all.33.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.33.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.33.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.33=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.33.unk)
summary(glm.modA.33) 
Anova(glm.modA.33) # 0.25885
Aplot.33<-
  plot_model(glm.modA.33, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 42", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.34.unk<- data.frame(AMR.chem.std.all.P38[34])
AMR.chem.std.all.34.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.34.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.34.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.34=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.34.unk)
summary(glm.modA.34) 
Anova(glm.modA.34) # 0.9961
Aplot.34<-
  plot_model(glm.modA.34, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 43", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.35.unk<- data.frame(AMR.chem.std.all.P38[35])
AMR.chem.std.all.35.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.35.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.35.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.35=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.35.unk)
summary(glm.modA.35) 
Anova(glm.modA.35) # 0.08564
Aplot.35<-
  plot_model(glm.modA.35, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 44", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.36.unk<- data.frame(AMR.chem.std.all.P38[36])
AMR.chem.std.all.36.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.36.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.36.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.36=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.36.unk)
summary(glm.modA.36) 
Anova(glm.modA.36) # 0.4951
Aplot.36<-
  plot_model(glm.modA.36, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 49", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.37.unk<- data.frame(AMR.chem.std.all.P38[37])
AMR.chem.std.all.37.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.37.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.37.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.37=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.37.unk)
summary(glm.modA.37) 
Anova(glm.modA.37) # NA
Aplot.37<-
  plot_model(glm.modA.37, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 50", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.all.38.unk<- data.frame(AMR.chem.std.all.P38[38])
AMR.chem.std.all.38.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.all.38.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.all.38.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modA.38=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.all.38.unk)
summary(glm.modA.38) 
Anova(glm.modA.38) # 0.01716
Aplot.38<-
  plot_model(glm.modA.38, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 51", show.legend=F, line.size=1, dot.size=3, dodge=.7)




### Within-AMR dataset  ----
AMR.chem.std.C.P38<- split(AMR.chem.std.C.P, AMR.chem.std.C.P$Comp.order)
head(AMR.chem.std.C.P38[1]) #list


#
AMR.chem.std.C.01.butenol<- data.frame(AMR.chem.std.C.P38[1])
AMR.chem.std.C.01.butenol[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.01.butenol[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.01.butenol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.01=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.01.butenol)
summary(glm.modC.01) 
Anova(glm.modC.01) # 0.42906
Cplot.01<-
  plot_model(glm.modC.01, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="2-Methyl-3-buten-2-ol", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.02.acetic<- data.frame(AMR.chem.std.C.P38[2])
AMR.chem.std.C.02.acetic[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.02.acetic[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.02.acetic)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.02=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.02.acetic)
summary(glm.modC.02) 
Anova(glm.modC.02) # 0.0394
Cplot.02<-
  plot_model(glm.modC.02, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Acetic acid", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.03.unk<- data.frame(AMR.chem.std.C.P38[3])
AMR.chem.std.C.03.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.03.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.03.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.03=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.03.unk)
summary(glm.modC.03) 
Anova(glm.modC.03) # < 0.00000001
Cplot.03<-
  plot_model(glm.modC.03, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 06", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.04.pentanol<- data.frame(AMR.chem.std.C.P38[4])
AMR.chem.std.C.04.pentanol[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.04.pentanol[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.04.pentanol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.04=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.04.pentanol)
summary(glm.modC.04) 
Anova(glm.modC.04) # 0.6678
Cplot.04<-
  plot_model(glm.modC.04, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="4-Methylpentan-2-ol", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.05.unk<- data.frame(AMR.chem.std.C.P38[5])
AMR.chem.std.C.05.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.05.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.05.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.05=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.05.unk)
summary(glm.modC.05) 
Anova(glm.modC.05) # 0.3808
Cplot.05<-
  plot_model(glm.modC.05, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 07", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.06.propanediol<- data.frame(AMR.chem.std.C.P38[6])
AMR.chem.std.C.06.propanediol[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.06.propanediol[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.06.propanediol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.06=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.06.propanediol)
summary(glm.modC.06) 
Anova(glm.modC.06) # 0.7117
Cplot.06<-
  plot_model(glm.modC.06, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="1,2-Propanediol", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.07.butanediol<- data.frame(AMR.chem.std.C.P38[7])
AMR.chem.std.C.07.butanediol[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.07.butanediol[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.07.butanediol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.07=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.07.butanediol)
summary(glm.modC.07) 
Anova(glm.modC.07) # 0.8112
Cplot.07<-
  plot_model(glm.modC.07, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="2,3-Butanediol", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.08.hexanol<- data.frame(AMR.chem.std.C.P38[8])
AMR.chem.std.C.08.hexanol[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.08.hexanol[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.08.hexanol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.08=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.08.hexanol)
summary(glm.modC.08) 
Anova(glm.modC.08) # 0.8124
Cplot.08<-
  plot_model(glm.modC.08, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Hexan-2-ol", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.09.unk<- data.frame(AMR.chem.std.C.P38[9])
AMR.chem.std.C.09.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.09.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.09.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.09=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.09.unk)
summary(glm.modC.09) 
Anova(glm.modC.09) # 0.8355
Cplot.09<-
  plot_model(glm.modC.09, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 13", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.10.formamide<- data.frame(AMR.chem.std.C.P38[10])
AMR.chem.std.C.10.formamide[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.10.formamide[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.10.formamide)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.10=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.10.formamide)
summary(glm.modC.10) 
Anova(glm.modC.10) # 0.003546
Cplot.10<-
  plot_model(glm.modC.10, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="N-Methylformamide", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.11.benzene<- data.frame(AMR.chem.std.C.P38[11])
AMR.chem.std.C.11.benzene[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.11.benzene[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.11.benzene)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.11=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.11.benzene)
summary(glm.modC.11) 
Anova(glm.modC.11) # 0.8103
Cplot.11<-
  plot_model(glm.modC.11, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Ethylbenzene", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.12.pyrazine<- data.frame(AMR.chem.std.C.P38[12])
AMR.chem.std.C.12.pyrazine[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.12.pyrazine[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.12.pyrazine)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.12=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.12.pyrazine)
summary(glm.modC.12) 
Anova(glm.modC.12) # 0.046934
Cplot.12<-
  plot_model(glm.modC.12, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="2,5-Dimethylpyrazine", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.13.unk<- data.frame(AMR.chem.std.C.P38[13])
AMR.chem.std.C.13.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.13.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.13.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.13=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.13.unk)
summary(glm.modC.13) 
Anova(glm.modC.13) # 0.8251016
Cplot.13<-
  plot_model(glm.modC.13, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 18", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.14.unk<- data.frame(AMR.chem.std.C.P38[14])
AMR.chem.std.C.14.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.14.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.14.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.14=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.14.unk)
summary(glm.modC.14) 
Anova(glm.modC.14) # 0.06214
Cplot.14<-
  plot_model(glm.modC.14, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unk. nitrogen-containing comp 07", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.15.unk<- data.frame(AMR.chem.std.C.P38[15])
AMR.chem.std.C.15.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.15.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.15.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.15=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.15.unk)
summary(glm.modC.15) 
Anova(glm.modC.15) # 0.02333
Cplot.15<-
  plot_model(glm.modC.15, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 25", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.16.acetophenone<- data.frame(AMR.chem.std.C.P38[16])
AMR.chem.std.C.16.acetophenone[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.16.acetophenone[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.16.acetophenone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.16=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.16.acetophenone)
summary(glm.modC.16) 
Anova(glm.modC.16) # NA
Cplot.16<-
  plot_model(glm.modC.16, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Acetophenone", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.17.pyrrole<- data.frame(AMR.chem.std.C.P38[17])
AMR.chem.std.C.17.pyrrole[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.17.pyrrole[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.17.pyrrole)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.17=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.17.pyrrole)
summary(glm.modC.17) 
Anova(glm.modC.17) # 0.2367
Cplot.17<-
  plot_model(glm.modC.17, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="1H-Pyrrole-2,5-dione", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.18.unk<- data.frame(AMR.chem.std.C.P38[18])
AMR.chem.std.C.18.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.18.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.18.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.18=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.18.unk)
summary(glm.modC.18) 
Anova(glm.modC.18) # <0.000001
Cplot.18<-
  plot_model(glm.modC.18, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 29", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.19.unk<- data.frame(AMR.chem.std.C.P38[19])
AMR.chem.std.C.19.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.19.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.19.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.19=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.19.unk)
summary(glm.modC.19) 
Anova(glm.modC.19) # NA
Cplot.19<-
  plot_model(glm.modC.19, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 31", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.20.benzofuran<- data.frame(AMR.chem.std.C.P38[20])
AMR.chem.std.C.20.benzofuran[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.20.benzofuran[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.20.benzofuran)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.20=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.20.benzofuran)
summary(glm.modC.20) 
Anova(glm.modC.20) # 0.21409
Cplot.20<-
  plot_model(glm.modC.20, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="3,7-Dimethyl-(...)-hexahydro-1-benzofuran", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.21.unk<- data.frame(AMR.chem.std.C.P38[21])
AMR.chem.std.C.21.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.21.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.21.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.21=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.21.unk)
summary(glm.modC.21) 
Anova(glm.modC.21) # 0.862414
Cplot.21<-
  plot_model(glm.modC.21, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 33", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.22.hydroxyacetophenone<- data.frame(AMR.chem.std.C.P38[22])
AMR.chem.std.C.22.hydroxyacetophenone[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.22.hydroxyacetophenone[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.22.hydroxyacetophenone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.22=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.22.hydroxyacetophenone)
summary(glm.modC.22) 
Anova(glm.modC.22) # 0.4947
Cplot.22<-
  plot_model(glm.modC.22, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="3,4-Dihydroxyacetophenone", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.23.unk<- data.frame(AMR.chem.std.C.P38[23])
AMR.chem.std.C.23.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.23.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.23.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.23=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.23.unk)
summary(glm.modC.23) 
Anova(glm.modC.23) # 0.1411
Cplot.23<-
  plot_model(glm.modC.23, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 34", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.24.unk<- data.frame(AMR.chem.std.C.P38[24])
AMR.chem.std.C.24.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.24.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.24.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.24=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.24.unk)
summary(glm.modC.24) 
Anova(glm.modC.24) # 0.1200579
Cplot.24<-
  plot_model(glm.modC.24, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 35", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.25.unk<- data.frame(AMR.chem.std.C.P38[25])
AMR.chem.std.C.25.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.25.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.25.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.25=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.25.unk)
summary(glm.modC.25) 
Anova(glm.modC.25) # NA
Cplot.25<-
  plot_model(glm.modC.25, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 36", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.26.pyrrolidine<- data.frame(AMR.chem.std.C.P38[26])
AMR.chem.std.C.26.pyrrolidine[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.26.pyrrolidine[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.26.pyrrolidine)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.26=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.26.pyrrolidine)
summary(glm.modC.26) 
Anova(glm.modC.26) # 0.23879
Cplot.26<-
  plot_model(glm.modC.26, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="3-Methyl-4-propylpyrrolidine-2,5-dione", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.27.3acetophenone<- data.frame(AMR.chem.std.C.P38[27])
AMR.chem.std.C.27.3acetophenone[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.27.3acetophenone[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.27.3acetophenone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.27=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.27.3acetophenone)
summary(glm.modC.27) 
Anova(glm.modC.27) # 0.7409
Cplot.27<-
  plot_model(glm.modC.27, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="3-Ethylacetophenone", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.28.cyclopentanone<- data.frame(AMR.chem.std.C.P38[28])
AMR.chem.std.C.28.cyclopentanone[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.28.cyclopentanone[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.28.cyclopentanone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.28=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.28.cyclopentanone)
summary(glm.modC.28) 
Anova(glm.modC.28) # 0.94055
Cplot.28<-
  plot_model(glm.modC.28, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="2-Cyclopentylcyclopentanone", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.29.aminophenol<- data.frame(AMR.chem.std.C.P38[29])
AMR.chem.std.C.29.aminophenol[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.29.aminophenol[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.29.aminophenol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.29=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.29.aminophenol)
summary(glm.modC.29) 
Anova(glm.modC.29) # 0.95351
Cplot.29<-
  plot_model(glm.modC.29, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="3-(Diethylamino)phenol", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.30.4acetophenone<- data.frame(AMR.chem.std.C.P38[30])
AMR.chem.std.C.30.4acetophenone[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.30.4acetophenone[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.30.4acetophenone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.30=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.30.4acetophenone)
summary(glm.modC.30) 
Anova(glm.modC.30) # 0.84851
Cplot.30<-
  plot_model(glm.modC.30, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="4-Ethylacetophenone", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.31.aniline<- data.frame(AMR.chem.std.C.P38[31])
AMR.chem.std.C.31.aniline[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.31.aniline[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.31.aniline)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.31=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.31.aniline)
summary(glm.modC.31) 
Anova(glm.modC.31) # 0.000000003637
Cplot.31<-
  plot_model(glm.modC.31, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="4-Methoxy-3,5-dimethylaniline", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.32.unk<- data.frame(AMR.chem.std.C.P38[32])
AMR.chem.std.C.32.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.32.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.32.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.32=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.32.unk)
summary(glm.modC.32) 
Anova(glm.modC.32) # 0.7952
Cplot.32<-
  plot_model(glm.modC.32, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 40", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.33.unk<- data.frame(AMR.chem.std.C.P38[33])
AMR.chem.std.C.33.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.33.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.33.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.33=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.33.unk)
summary(glm.modC.33) 
Anova(glm.modC.33) # 0.06893
Cplot.33<-
  plot_model(glm.modC.33, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 42", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.34.unk<- data.frame(AMR.chem.std.C.P38[34])
AMR.chem.std.C.34.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.34.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.34.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.34=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.34.unk)
summary(glm.modC.34) 
Anova(glm.modC.34) # 0.04518
Cplot.34<-
  plot_model(glm.modC.34, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 43", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.35.unk<- data.frame(AMR.chem.std.C.P38[35])
AMR.chem.std.C.35.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.35.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.35.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.35=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.35.unk)
summary(glm.modC.35) 
Anova(glm.modC.35) # 0.02700
Cplot.35<-
  plot_model(glm.modC.35, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 44", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.36.unk<- data.frame(AMR.chem.std.C.P38[36])
AMR.chem.std.C.36.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.36.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.36.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.36=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.36.unk)
summary(glm.modC.36) 
Anova(glm.modC.36) # 0.0008905
Cplot.36<-
  plot_model(glm.modC.36, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 49", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.37.unk<- data.frame(AMR.chem.std.C.P38[37])
AMR.chem.std.C.37.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.37.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.37.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.37=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.37.unk)
summary(glm.modC.37) 
Anova(glm.modC.37) # NA
Cplot.37<-
  plot_model(glm.modC.37, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 50", show.legend=F, line.size=1, dot.size=3, dodge=.7)


#
AMR.chem.std.C.38.unk<- data.frame(AMR.chem.std.C.P38[38])
AMR.chem.std.C.38.unk[,c(1:13,31,32)]<- lapply(AMR.chem.std.C.38.unk[,c(1:13,31,32)], factor)
names(AMR.chem.std.C.38.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_A3","AMR_period","Handler","Batch","Column","DOB","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","n.Date","z.Age","z.Seasonality","RT","Area","Ncomp","SArea","meanArea","RPA","logRPA","asinRPA","Comp.order","Comp.name")

glm.modC.38=lmer(logRPA ~ MRC*AMR_period +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.C.38.unk)
summary(glm.modC.38) 
Anova(glm.modC.38) # NA
Cplot.38<-
  plot_model(glm.modC.38, type="pred", terms=c("AMR_period","MRC"), bias_correction=T,
             title="Unidentified comp 51", show.legend=F, line.size=1, dot.size=3, dodge=.7)


save.image(file="current_CapuAMR.RData")




### Suppl. Figure S1 ----
#dev.new()
set_theme(base = theme_classic(base_size=6))

# Full dataset
grid.draw(rbind(cbind(ggplotGrob(Aplot.01), ggplotGrob(Aplot.02), ggplotGrob(Aplot.03), ggplotGrob(Aplot.04), ggplotGrob(Aplot.05), ggplotGrob(Aplot.06), ggplotGrob(Aplot.07), size="first"),
                cbind(ggplotGrob(Aplot.08), ggplotGrob(Aplot.09), ggplotGrob(Aplot.10), ggplotGrob(Aplot.11), ggplotGrob(Aplot.12), ggplotGrob(Aplot.13), ggplotGrob(Aplot.14), size="first"),
                cbind(ggplotGrob(Aplot.15), ggplotGrob(Aplot.16), ggplotGrob(Aplot.17), ggplotGrob(Aplot.18), ggplotGrob(Aplot.19), ggplotGrob(Aplot.20), ggplotGrob(Aplot.21), size="first"),
                cbind(ggplotGrob(Aplot.22), ggplotGrob(Aplot.23), ggplotGrob(Aplot.24), ggplotGrob(Aplot.25), ggplotGrob(Aplot.26), ggplotGrob(Aplot.27), ggplotGrob(Aplot.28), size="first"),
                cbind(ggplotGrob(Aplot.29), ggplotGrob(Aplot.30), ggplotGrob(Aplot.31), ggplotGrob(Aplot.32), ggplotGrob(Aplot.33), ggplotGrob(Aplot.34), ggplotGrob(Aplot.35), size="first"),
                cbind(ggplotGrob(Aplot.36), ggplotGrob(Aplot.37), ggplotGrob(Aplot.38), ggplotGrob(Aplot.38), ggplotGrob(Aplot.38), ggplotGrob(Aplot.38), ggplotGrob(Aplot.38), size="first")))

# Within-AMR dataset
grid.draw(rbind(cbind(ggplotGrob(Cplot.01), ggplotGrob(Cplot.02), ggplotGrob(Cplot.03), ggplotGrob(Cplot.04), ggplotGrob(Cplot.05), ggplotGrob(Cplot.06), ggplotGrob(Cplot.07), size="first"),
                cbind(ggplotGrob(Cplot.08), ggplotGrob(Cplot.09), ggplotGrob(Cplot.10), ggplotGrob(Cplot.11), ggplotGrob(Cplot.12), ggplotGrob(Cplot.13), ggplotGrob(Cplot.14), size="first"),
                cbind(ggplotGrob(Cplot.15), ggplotGrob(Cplot.16), ggplotGrob(Cplot.17), ggplotGrob(Cplot.18), ggplotGrob(Cplot.19), ggplotGrob(Cplot.20), ggplotGrob(Cplot.21), size="first"),
                cbind(ggplotGrob(Cplot.22), ggplotGrob(Cplot.23), ggplotGrob(Cplot.24), ggplotGrob(Cplot.25), ggplotGrob(Cplot.26), ggplotGrob(Cplot.27), ggplotGrob(Cplot.28), size="first"),
                cbind(ggplotGrob(Cplot.29), ggplotGrob(Cplot.30), ggplotGrob(Cplot.31), ggplotGrob(Cplot.32), ggplotGrob(Cplot.33), ggplotGrob(Cplot.34), ggplotGrob(Cplot.35), size="first"),
                cbind(ggplotGrob(Cplot.36), ggplotGrob(Cplot.37), ggplotGrob(Cplot.38), ggplotGrob(Cplot.38), ggplotGrob(Cplot.38), ggplotGrob(Cplot.38), ggplotGrob(Cplot.38), size="first")))







# 5) Fecal androgen concentration (T) ----
## T datasets ----
# Full and within-AMR complete datasets
AMR.chem.std.T.all<- subset(AMR.chem.std.all, !is.na(T.conc))
AMR.chem.std.T.all<- mutate(AMR.chem.std.T.all, logT.conc=log(T.conc))
AMR.chem.std.T.all[,c(1:10,18)]<- lapply(AMR.chem.std.T.all[,c(1:10,18)], factor)
AMR.chem.std.T.all<- AMR.chem.std.T.all %>% arrange(Comp)

AMR.chem.std.T.C<- subset(AMR.chem.std.C, !is.na(T.conc))
AMR.chem.std.T.C<- mutate(AMR.chem.std.T.C, logT.conc=log(T.conc))
AMR.chem.std.T.C[,c(1:10,18)]<- lapply(AMR.chem.std.T.C[,c(1:10,18)], factor)
AMR.chem.std.T.C<- AMR.chem.std.T.C %>% arrange(Comp)

# 38 selected comp T datasets
AMR.chem.std.T.all.P<- subset(AMR.chem.std.all.P, !is.na(T.conc))
AMR.chem.std.T.all.P<- mutate(AMR.chem.std.T.all.P, logT.conc=log(T.conc))
AMR.chem.std.T.all.P[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.P[,c(1:11,25,26)], factor)
AMR.chem.std.T.all.P<- AMR.chem.std.T.all.P %>% arrange(Comp)



## T RSM   ----
### Full dataset ----
# Full model
Full.modT.all=lmer(asinRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler) +(1|Filename)
                      +(1+logT.conc+z.Seasonality||Comp) +(1|Comp:Group) +(1|Batch:Group),
                      data=AMR.chem.std.T.all, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Full.modT.all) # 

# Reduced model
Red.modT.all=lmer(asinRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler) +(1|Filename)
                     +(1+z.Seasonality||Comp) +(1|Comp:Group) +(1|Batch:Group),
                     data=AMR.chem.std.T.all, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Red.modT.all)
anova(Red.modT.all, Full.modT.all, test="Chisq") #NS



### Within-AMR dataset ----
# Full model
Full.modT.C=lmer(asinRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler) +(1|Filename)
                   +(1+logT.conc+z.Seasonality||Comp) +(1|Comp:Group) +(1|Batch:Group),
                   data=AMR.chem.std.T.C, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Full.modT.C) # 

# Reduced model
Red.modT.C=lmer(asinRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler) +(1|Filename)
                  +(1+z.Seasonality||Comp) +(1|Comp:Group) +(1|Batch:Group),
                  data=AMR.chem.std.T.C, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Red.modT.C)
anova(Red.modT.C, Full.modT.C, test="Chisq") #NS p=1.000


save.image(file="current_CapuAMR.RData")





## T GLMM Comp of interest ----
# 38 comps
AMR.chem.std.T.all.indiv<- split(AMR.chem.std.T.all.P, AMR.chem.std.T.all.P$Comp.order)
head(AMR.chem.std.T.all.indiv[1])


AMR.chem.std.T.all.01.butenol<- data.frame(AMR.chem.std.T.all.indiv[1])
AMR.chem.std.T.all.01.butenol[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.01.butenol[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.01.butenol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.01=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.01.butenol)
summary(glm.modT.01) 
Anova(glm.modT.01) # 0.82754
Tplot.01<-
  plot_model(glm.modT.01, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="2-Methyl-3-buten-2-ol", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.02.acetic<- data.frame(AMR.chem.std.T.all.indiv[2])
AMR.chem.std.T.all.02.acetic[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.02.acetic[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.02.acetic)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.02=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.02.acetic)
summary(glm.modT.02) 
Anova(glm.modT.02) # 0.51372
Tplot.02<-
  plot_model(glm.modT.02, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Acetic acid", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.03.unk<- data.frame(AMR.chem.std.T.all.indiv[3])
AMR.chem.std.T.all.03.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.03.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.03.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.03=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.03.unk)
summary(glm.modT.03) 
Anova(glm.modT.03) # 0.157782
Tplot.03<-
  plot_model(glm.modT.03, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 06", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.04.pentanol<- data.frame(AMR.chem.std.T.all.indiv[4])
AMR.chem.std.T.all.04.pentanol[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.04.pentanol[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.04.pentanol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.04=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.04.pentanol)
summary(glm.modT.04) 
Anova(glm.modT.04) # 0.8966
Tplot.04<-
  plot_model(glm.modT.04, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="4-Methylpentan-2-ol", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.05.unk<- data.frame(AMR.chem.std.T.all.indiv[5])
AMR.chem.std.T.all.05.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.05.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.05.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.05=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.05.unk)
summary(glm.modT.05) 
Anova(glm.modT.05) # 0.90643
Tplot.05<-
  plot_model(glm.modT.05, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 07", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.06.propanediol<- data.frame(AMR.chem.std.T.all.indiv[6])
AMR.chem.std.T.all.06.propanediol[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.06.propanediol[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.06.propanediol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.06=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.06.propanediol)
summary(glm.modT.06) 
Anova(glm.modT.06) # 0.8959
Tplot.06<-
  plot_model(glm.modT.06, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="1,2-Propanediol", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.07.butanediol<- data.frame(AMR.chem.std.T.all.indiv[7])
AMR.chem.std.T.all.07.butanediol[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.07.butanediol[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.07.butanediol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.07=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.07.butanediol)
summary(glm.modT.07) 
Anova(glm.modT.07) # 0.7858
Tplot.07<-
  plot_model(glm.modT.07, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="2,3-Butanediol", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.08.hexanol<- data.frame(AMR.chem.std.T.all.indiv[8])
AMR.chem.std.T.all.08.hexanol[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.08.hexanol[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.08.hexanol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.08=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.08.hexanol)
summary(glm.modT.08) 
Anova(glm.modT.08) # 0.5089
Tplot.08<-
  plot_model(glm.modT.08, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Hexan-2-ol", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.09.unk<- data.frame(AMR.chem.std.T.all.indiv[9])
AMR.chem.std.T.all.09.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.09.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.09.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.09=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.09.unk)
summary(glm.modT.09) 
Anova(glm.modT.09) # 0.16040
Tplot.09<-
  plot_model(glm.modT.09, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 13", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.10.formamide<- data.frame(AMR.chem.std.T.all.indiv[10])
AMR.chem.std.T.all.10.formamide[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.10.formamide[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.10.formamide)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.10=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.10.formamide)
summary(glm.modT.10) 
Anova(glm.modT.10) # 0.1333
Tplot.10<-
  plot_model(glm.modT.10, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="N-Methylformamide", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.11.benzene<- data.frame(AMR.chem.std.T.all.indiv[11])
AMR.chem.std.T.all.11.benzene[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.11.benzene[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.11.benzene)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.11=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.11.benzene)
summary(glm.modT.11) 
Anova(glm.modT.11) #* 0.04019
Tplot.11<-
  plot_model(glm.modT.11, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Ethylbenzene", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.12.pyrazine<- data.frame(AMR.chem.std.T.all.indiv[12])
AMR.chem.std.T.all.12.pyrazine[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.12.pyrazine[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.12.pyrazine)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.12=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.12.pyrazine)
summary(glm.modT.12) 
Anova(glm.modT.12) # 0.0998985
Tplot.12<-
  plot_model(glm.modT.12, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="2,5-Dimethylpyrazine", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.13.unk<- data.frame(AMR.chem.std.T.all.indiv[13])
AMR.chem.std.T.all.13.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.13.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.13.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.13=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.13.unk)
summary(glm.modT.13) 
Anova(glm.modT.13) # 0.5192
Tplot.13<-
  plot_model(glm.modT.13, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 18", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.14.unk<- data.frame(AMR.chem.std.T.all.indiv[14])
AMR.chem.std.T.all.14.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.14.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.14.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.14=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.14.unk)
summary(glm.modT.14) 
Anova(glm.modT.14) # 0.99893
Tplot.14<-
  plot_model(glm.modT.14, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unk. nitrogen-containing comp 07", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.15.unk<- data.frame(AMR.chem.std.T.all.indiv[15])
AMR.chem.std.T.all.15.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.15.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.15.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.15=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.15.unk)
summary(glm.modT.15) 
Anova(glm.modT.15) # 0.1660
Tplot.15<-
  plot_model(glm.modT.15, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 25", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.16.acetophenone<- data.frame(AMR.chem.std.T.all.indiv[16])
AMR.chem.std.T.all.16.acetophenone[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.16.acetophenone[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.16.acetophenone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.16=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.16.acetophenone)
summary(glm.modT.16) 
Anova(glm.modT.16) # 0.9988
Tplot.16<-
  plot_model(glm.modT.16, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Acetophenone", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.17.pyrrole<- data.frame(AMR.chem.std.T.all.indiv[17])
AMR.chem.std.T.all.17.pyrrole[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.17.pyrrole[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.17.pyrrole)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.17=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.17.pyrrole)
summary(glm.modT.17) 
Anova(glm.modT.17) # 0.6332
Tplot.17<-
  plot_model(glm.modT.17, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="1H-Pyrrole-2,5-dione", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.18.unk<- data.frame(AMR.chem.std.T.all.indiv[18])
AMR.chem.std.T.all.18.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.18.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.18.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.18=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.18.unk)
summary(glm.modT.18) 
Anova(glm.modT.18) # 0.27862
Tplot.18<-
  plot_model(glm.modT.18, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 29", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.19.unk<- data.frame(AMR.chem.std.T.all.indiv[19])
AMR.chem.std.T.all.19.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.19.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.19.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.19=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.19.unk)
summary(glm.modT.19) 
Anova(glm.modT.19) #* 0.0000000000000002
Tplot.19<-
  plot_model(glm.modT.19, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 31", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.20.benzofuran<- data.frame(AMR.chem.std.T.all.indiv[20])
AMR.chem.std.T.all.20.benzofuran[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.20.benzofuran[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.20.benzofuran)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.20=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.20.benzofuran)
summary(glm.modT.20) 
Anova(glm.modT.20) # 0.65363
Tplot.20<-
  plot_model(glm.modT.20, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="3,7-Dimethyl-(...)-hexahydro-1-benzofuran", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.21.unk<- data.frame(AMR.chem.std.T.all.indiv[21])
AMR.chem.std.T.all.21.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.21.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.21.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.21=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.21.unk)
summary(glm.modT.21) 
Anova(glm.modT.21) # 0.54811
Tplot.21<-
  plot_model(glm.modT.21, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 33", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.22.hydroxyacetophenone<- data.frame(AMR.chem.std.T.all.indiv[22])
AMR.chem.std.T.all.22.hydroxyacetophenone[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.22.hydroxyacetophenone[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.22.hydroxyacetophenone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.22=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.22.hydroxyacetophenone)
summary(glm.modT.22) 
Anova(glm.modT.22) #  0.4374
Tplot.22<-
  plot_model(glm.modT.22, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="3,4-Dihydroxyacetophenone", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.23.unk<- data.frame(AMR.chem.std.T.all.indiv[23])
AMR.chem.std.T.all.23.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.23.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.23.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.23=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.23.unk)
summary(glm.modT.23) 
Anova(glm.modT.23) # 0.2269
Tplot.23<-
  plot_model(glm.modT.23, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 34", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.24.unk<- data.frame(AMR.chem.std.T.all.indiv[24])
AMR.chem.std.T.all.24.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.24.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.24.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.24=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.24.unk)
summary(glm.modT.24) 
Anova(glm.modT.24) # 0.836766
Tplot.24<-
  plot_model(glm.modT.24, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 35", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.25.unk<- data.frame(AMR.chem.std.T.all.indiv[25])
AMR.chem.std.T.all.25.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.25.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.25.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.25=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.25.unk)
summary(glm.modT.25) 
Anova(glm.modT.25) # 0.5107
Tplot.25<-
  plot_model(glm.modT.25, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 36", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.26.pyrrolidine<- data.frame(AMR.chem.std.T.all.indiv[26])
AMR.chem.std.T.all.26.pyrrolidine[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.26.pyrrolidine[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.26.pyrrolidine)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.26=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.26.pyrrolidine)
summary(glm.modT.26) 
Anova(glm.modT.26) # 0.6197
Tplot.26<-
  plot_model(glm.modT.26, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="3-Methyl-4-propylpyrrolidine-2,5-dione", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.27.3acetophenone<- data.frame(AMR.chem.std.T.all.indiv[27])
AMR.chem.std.T.all.27.3acetophenone[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.27.3acetophenone[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.27.3acetophenone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.27=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.27.3acetophenone)
summary(glm.modT.27) 
Anova(glm.modT.27) # 0.2780
Tplot.27<-
  plot_model(glm.modT.27, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="3-Ethylacetophenone", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.28.cyclopentanone<- data.frame(AMR.chem.std.T.all.indiv[28])
AMR.chem.std.T.all.28.cyclopentanone[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.28.cyclopentanone[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.28.cyclopentanone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.28=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.28.cyclopentanone)
summary(glm.modT.28) 
Anova(glm.modT.28) # 0.4085
Tplot.28<-
  plot_model(glm.modT.28, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="2-Cyclopentylcyclopentanone", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.29.aminophenol<- data.frame(AMR.chem.std.T.all.indiv[29])
AMR.chem.std.T.all.29.aminophenol[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.29.aminophenol[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.29.aminophenol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.29=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.29.aminophenol)
summary(glm.modT.29) 
Anova(glm.modT.29) # 0.6004
Tplot.29<-
  plot_model(glm.modT.29, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="3-(Diethylamino)phenol", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.30.4acetophenone<- data.frame(AMR.chem.std.T.all.indiv[30])
AMR.chem.std.T.all.30.4acetophenone[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.30.4acetophenone[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.30.4acetophenone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.30=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.30.4acetophenone)
summary(glm.modT.30) 
Anova(glm.modT.30) # 0.5301
Tplot.30<-
  plot_model(glm.modT.30, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="4-Ethylacetophenone", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.31.aniline<- data.frame(AMR.chem.std.T.all.indiv[31])
AMR.chem.std.T.all.31.aniline[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.31.aniline[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.31.aniline)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.31=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.31.aniline)
summary(glm.modT.31) 
Anova(glm.modT.31) # 0.1500
Tplot.31<-
  plot_model(glm.modT.31, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="4-Methoxy-3,5-dimethylaniline", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.32.unk<- data.frame(AMR.chem.std.T.all.indiv[32])
AMR.chem.std.T.all.32.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.32.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.32.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.32=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.32.unk)
summary(glm.modT.32) 
Anova(glm.modT.32) # 0.9635
Tplot.32<-
  plot_model(glm.modT.32, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 40", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.33.unk<- data.frame(AMR.chem.std.T.all.indiv[33])
AMR.chem.std.T.all.33.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.33.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.33.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.33=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.33.unk)
summary(glm.modT.33) 
Anova(glm.modT.33) # 0.3249
Tplot.33<-
  plot_model(glm.modT.33, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 42", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.34.unk<- data.frame(AMR.chem.std.T.all.indiv[34])
AMR.chem.std.T.all.34.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.34.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.34.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.34=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.34.unk)
summary(glm.modT.34) 
Anova(glm.modT.34) # 0.46182
Tplot.34<-
  plot_model(glm.modT.34, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 43", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.35.unk<- data.frame(AMR.chem.std.T.all.indiv[35])
AMR.chem.std.T.all.35.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.35.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.35.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.35=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.35.unk)
summary(glm.modT.35) 
Anova(glm.modT.35) # 0.8211
Tplot.35<-
  plot_model(glm.modT.35, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 44", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.36.unk<- data.frame(AMR.chem.std.T.all.indiv[36])
AMR.chem.std.T.all.36.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.36.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.36.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.36=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.36.unk)
summary(glm.modT.36) 
Anova(glm.modT.36) # 0.157036
Tplot.36<-
  plot_model(glm.modT.36, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 49", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.37.unk<- data.frame(AMR.chem.std.T.all.indiv[37])
AMR.chem.std.T.all.37.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.37.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.37.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.37=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.37.unk)
summary(glm.modT.37) 
Anova(glm.modT.37) # 0.08379
Tplot.37<-
  plot_model(glm.modT.37, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 50", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.T.all.38.unk<- data.frame(AMR.chem.std.T.all.indiv[38])
AMR.chem.std.T.all.38.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.T.all.38.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.T.all.38.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logT.conc")

glm.modT.38=lmer(logRPA ~ logT.conc +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.T.all.38.unk)
summary(glm.modT.38) 
Anova(glm.modT.38) # 0.9398658
Tplot.38<-
  plot_model(glm.modT.38, type="pred", terms=c("logT.conc"), bias_correction=T,
             title="Unidentified comp 51", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)




## T FDR adjustments on p-values ----
Tpval<- c(0.82754, 0.51372, 0.157782, 0.8966, 0.90643, 0.8959, 0.7858, 0.5089, 0.16040, 0.1333, 0.04019, 0.0998985, 0.5192, 0.99893, 0.1660, 0.9988, 0.6332, 0.27862, 0.0000000000000002, 0.65363, 0.54811, 0.4374, 0.2269, 0.836766, 0.5107, 0.6197, 0.2780, 0.4085, 0.6004, 0.5301, 0.1500, 0.9635, 0.3249, 0.46182, 0.8211, 0.157036, 0.08379, 0.9398658)

Tpval.adjust<- round(p.adjust(Tpval, method="fdr", n=length(Tpval)), digits=3)
Tpval.adjust





## Suppl. Figure S2: T gridplots ----
empty.plot<-   ggplot() +
  theme_void(base_size=8) +
  geom_text(aes(0,0,label='')) +
  xlab(NULL)

#dev.new()
set_theme(base = theme_classic(base_size=8))

grid.draw(rbind(cbind(ggplotGrob(Tplot.01), ggplotGrob(Tplot.02), ggplotGrob(Tplot.03), ggplotGrob(Tplot.04), ggplotGrob(Tplot.05), ggplotGrob(Tplot.06), ggplotGrob(Tplot.07), size="first"),
                cbind(ggplotGrob(Tplot.08), ggplotGrob(Tplot.09), ggplotGrob(Tplot.10), ggplotGrob(Tplot.11), ggplotGrob(Tplot.12), ggplotGrob(Tplot.13), ggplotGrob(Tplot.14), size="first"),
                cbind(ggplotGrob(Tplot.15), ggplotGrob(Tplot.16), ggplotGrob(Tplot.17), ggplotGrob(Tplot.18), ggplotGrob(Tplot.19), ggplotGrob(Tplot.20), ggplotGrob(Tplot.21), size="first"),
                cbind(ggplotGrob(Tplot.22), ggplotGrob(Tplot.23), ggplotGrob(Tplot.24), ggplotGrob(Tplot.25), ggplotGrob(Tplot.26), ggplotGrob(Tplot.27), ggplotGrob(Tplot.28), size="first"),
                cbind(ggplotGrob(Tplot.29), ggplotGrob(Tplot.30), ggplotGrob(Tplot.31), ggplotGrob(Tplot.32), ggplotGrob(Tplot.33), ggplotGrob(Tplot.34), ggplotGrob(Tplot.35), size="first"),
                cbind(ggplotGrob(Tplot.36), ggplotGrob(Tplot.37), ggplotGrob(Tplot.38), ggplotGrob(empty.plot), ggplotGrob(empty.plot), ggplotGrob(empty.plot), ggplotGrob(empty.plot), size="first")))






# 6) Olfactory behavior (B) ----
## B datasets ----
AMR.chem.std.B.all<- subset(AMR.chem.std.all, !is.na(UW_rate))
AMR.chem.std.B.all<- mutate(AMR.chem.std.B.all, logUW=log(UW_rate+1))
AMR.chem.std.B.all[,c(1:10,18)]<- lapply(AMR.chem.std.B.all[,c(1:10,18)], factor)
AMR.chem.std.B.all<- AMR.chem.std.B.all %>% arrange(Comp)

AMR.chem.std.B.C<- subset(AMR.chem.std.C, !is.na(UW_rate))
AMR.chem.std.B.C<- mutate(AMR.chem.std.B.C, logUW=log(UW_rate+1))
AMR.chem.std.B.C[,c(1:10,18)]<- lapply(AMR.chem.std.B.C[,c(1:10,18)], factor)
AMR.chem.std.B.C<- AMR.chem.std.B.C %>% arrange(Comp)

# 38 selected comp B datasets
AMR.chem.std.B.all.P<- subset(AMR.chem.std.all.P, !is.na(UW_rate))
AMR.chem.std.B.all.P<- mutate(AMR.chem.std.B.all.P, logUW=log(UW_rate+1))
AMR.chem.std.B.all.P[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.P[,c(1:11,25,26)], factor)
AMR.chem.std.B.all.P<- AMR.chem.std.B.all.P %>% arrange(Comp)




## B RSM ----
### Full dataset ----
# Full model
Full.modB.all=lmer(asinRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler) +(1|Filename)
                   +(1+logUW+z.Seasonality||Comp) +(1|Comp:Group) +(1|Batch:Group),
                   data=AMR.chem.std.B.all, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Full.modB.all) # 

# Reduced model
Red.modB.all=lmer(asinRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler) +(1|Filename)
                  +(1+z.Seasonality||Comp) +(1|Comp:Group) +(1|Batch:Group),
                  data=AMR.chem.std.B.all, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Red.modB.all)
anova(Red.modB.all, Full.modB.all, test="Chisq") #NS p=1.000



### Within-AMR dataset ----
# Full model
Full.modB.C=lmer(asinRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler) +(1|Filename)
                 +(1+logUW+z.Seasonality||Comp) +(1|Comp:Group) +(1|Batch:Group),
                 data=AMR.chem.std.B.C, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Full.modB.C) # 

# Reduced model
Red.modB.C=lmer(asinRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler) +(1|Filename)
                +(1+z.Seasonality||Comp) +(1|Comp:Group) +(1|Batch:Group),
                data=AMR.chem.std.B.C, REML=F, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5), calc.derivs=F), verbose=3)
summary(Red.modB.C)
anova(Red.modB.C, Full.modB.C, test="Chisq") #NS







## B GLMM Comp of interest ----
# 38 comps
AMR.chem.std.B.all.indiv<- split(AMR.chem.std.B.all.P, AMR.chem.std.B.all.P$Comp.order)
head(AMR.chem.std.B.all.indiv[1])


AMR.chem.std.B.all.01.butenol<- data.frame(AMR.chem.std.B.all.indiv[1])
AMR.chem.std.B.all.01.butenol[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.01.butenol[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.01.butenol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.01=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.01.butenol)
summary(glm.modB.01) 
Anova(glm.modB.01) # 0.7486
Bplot.01<-
  plot_model(glm.modB.01, type="pred", terms=c("logUW"), bias_correction=T,
             title="2-Methyl-3-buten-2-ol", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.02.acetic<- data.frame(AMR.chem.std.B.all.indiv[2])
AMR.chem.std.B.all.02.acetic[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.02.acetic[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.02.acetic)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.02=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.02.acetic)
summary(glm.modB.02) 
Anova(glm.modB.02) # 0.1684854
Bplot.02<-
  plot_model(glm.modB.02, type="pred", terms=c("logUW"), bias_correction=T,
             title="Acetic acid", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.03.unk<- data.frame(AMR.chem.std.B.all.indiv[3])
AMR.chem.std.B.all.03.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.03.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.03.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.03=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.03.unk)
summary(glm.modB.03) 
Anova(glm.modB.03) # 0.1056
Bplot.03<-
  plot_model(glm.modB.03, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 06", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.04.pentanol<- data.frame(AMR.chem.std.B.all.indiv[4])
AMR.chem.std.B.all.04.pentanol[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.04.pentanol[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.04.pentanol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.04=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.04.pentanol)
summary(glm.modB.04) 
Anova(glm.modB.04) # 0.3833
Bplot.04<-
  plot_model(glm.modB.04, type="pred", terms=c("logUW"), bias_correction=T,
             title="4-Methylpentan-2-ol", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.05.unk<- data.frame(AMR.chem.std.B.all.indiv[5])
AMR.chem.std.B.all.05.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.05.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.05.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.05=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.05.unk)
summary(glm.modB.05) 
Anova(glm.modB.05) # 0.1344
Bplot.05<-
  plot_model(glm.modB.05, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 07", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.06.propanediol<- data.frame(AMR.chem.std.B.all.indiv[6])
AMR.chem.std.B.all.06.propanediol[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.06.propanediol[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.06.propanediol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.06=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.06.propanediol)
summary(glm.modB.06) 
Anova(glm.modB.06) # 0.1109
Bplot.06<-
  plot_model(glm.modB.06, type="pred", terms=c("logUW"), bias_correction=T,
             title="1,2-Propanediol", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.07.butanediol<- data.frame(AMR.chem.std.B.all.indiv[7])
AMR.chem.std.B.all.07.butanediol[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.07.butanediol[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.07.butanediol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.07=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.07.butanediol)
summary(glm.modB.07) 
Anova(glm.modB.07) # 0.9142
Bplot.07<-
  plot_model(glm.modB.07, type="pred", terms=c("logUW"), bias_correction=T,
             title="2,3-Butanediol", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.08.hexanol<- data.frame(AMR.chem.std.B.all.indiv[8])
AMR.chem.std.B.all.08.hexanol[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.08.hexanol[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.08.hexanol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.08=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.08.hexanol)
summary(glm.modB.08) 
Anova(glm.modB.08) # 0.5919
Bplot.08<-
  plot_model(glm.modB.08, type="pred", terms=c("logUW"), bias_correction=T,
             title="Hexan-2-ol", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.09.unk<- data.frame(AMR.chem.std.B.all.indiv[9])
AMR.chem.std.B.all.09.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.09.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.09.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.09=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.09.unk)
summary(glm.modB.09) 
Anova(glm.modB.09) # 0.0000000000000002
Bplot.09<-
  plot_model(glm.modB.09, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 13", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.10.formamide<- data.frame(AMR.chem.std.B.all.indiv[10])
AMR.chem.std.B.all.10.formamide[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.10.formamide[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.10.formamide)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.10=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.10.formamide)
summary(glm.modB.10) 
Anova(glm.modB.10) # 0.2860522
Bplot.10<-
  plot_model(glm.modB.10, type="pred", terms=c("logUW"), bias_correction=T,
             title="N-Methylformamide", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.11.benzene<- data.frame(AMR.chem.std.B.all.indiv[11])
AMR.chem.std.B.all.11.benzene[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.11.benzene[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.11.benzene)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.11=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.11.benzene)
summary(glm.modB.11) 
Anova(glm.modB.11) #* 0.07081
Bplot.11<-
  plot_model(glm.modB.11, type="pred", terms=c("logUW"), bias_correction=T,
             title="Ethylbenzene", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.12.pyrazine<- data.frame(AMR.chem.std.B.all.indiv[12])
AMR.chem.std.B.all.12.pyrazine[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.12.pyrazine[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.12.pyrazine)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.12=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.12.pyrazine)
summary(glm.modB.12) 
Anova(glm.modB.12) # 0.7480
Bplot.12<-
  plot_model(glm.modB.12, type="pred", terms=c("logUW"), bias_correction=T,
             title="2,5-Dimethylpyrazine", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.13.unk<- data.frame(AMR.chem.std.B.all.indiv[13])
AMR.chem.std.B.all.13.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.13.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.13.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.13=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.13.unk)
summary(glm.modB.13) 
Anova(glm.modB.13) # 0.05153
Bplot.13<-
  plot_model(glm.modB.13, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 18", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.14.unk<- data.frame(AMR.chem.std.B.all.indiv[14])
AMR.chem.std.B.all.14.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.14.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.14.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.14=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.14.unk)
summary(glm.modB.14) 
Anova(glm.modB.14) # 0.3646556
Bplot.14<-
  plot_model(glm.modB.14, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unk. nitrogen-containing comp 07", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.15.unk<- data.frame(AMR.chem.std.B.all.indiv[15])
AMR.chem.std.B.all.15.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.15.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.15.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.15=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.15.unk)
summary(glm.modB.15) 
Anova(glm.modB.15) # 0.006876
Bplot.15<-
  plot_model(glm.modB.15, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 25", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.16.acetophenone<- data.frame(AMR.chem.std.B.all.indiv[16])
AMR.chem.std.B.all.16.acetophenone[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.16.acetophenone[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.16.acetophenone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.16=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.16.acetophenone)
summary(glm.modB.16) 
Anova(glm.modB.16) # 0.9534
Bplot.16<-
  plot_model(glm.modB.16, type="pred", terms=c("logUW"), bias_correction=T,
             title="Acetophenone", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.17.pyrrole<- data.frame(AMR.chem.std.B.all.indiv[17])
AMR.chem.std.B.all.17.pyrrole[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.17.pyrrole[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.17.pyrrole)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.17=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.17.pyrrole)
summary(glm.modB.17) 
Anova(glm.modB.17) # 0.1716
Bplot.17<-
  plot_model(glm.modB.17, type="pred", terms=c("logUW"), bias_correction=T,
             title="1H-Pyrrole-2,5-dione", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.18.unk<- data.frame(AMR.chem.std.B.all.indiv[18])
AMR.chem.std.B.all.18.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.18.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.18.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.18=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.18.unk)
summary(glm.modB.18) 
Anova(glm.modB.18) # 0.0003855
Bplot.18<-
  plot_model(glm.modB.18, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 29", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.19.unk<- data.frame(AMR.chem.std.B.all.indiv[19])
AMR.chem.std.B.all.19.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.19.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.19.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.19=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.19.unk)
summary(glm.modB.19) 
Anova(glm.modB.19) # 0.8731
Bplot.19<-
  plot_model(glm.modB.19, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 31", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.20.benzofuran<- data.frame(AMR.chem.std.B.all.indiv[20])
AMR.chem.std.B.all.20.benzofuran[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.20.benzofuran[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.20.benzofuran)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.20=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.20.benzofuran)
summary(glm.modB.20) 
Anova(glm.modB.20) # 0.8547
Bplot.20<-
  plot_model(glm.modB.20, type="pred", terms=c("logUW"), bias_correction=T,
             title="3,7-Dimethyl-(...)-hexahydro-1-benzofuran", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.21.unk<- data.frame(AMR.chem.std.B.all.indiv[21])
AMR.chem.std.B.all.21.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.21.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.21.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.21=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.21.unk)
summary(glm.modB.21) 
Anova(glm.modB.21) # 0.85350
Bplot.21<-
  plot_model(glm.modB.21, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 33", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.22.hydroxyacetophenone<- data.frame(AMR.chem.std.B.all.indiv[22])
AMR.chem.std.B.all.22.hydroxyacetophenone[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.22.hydroxyacetophenone[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.22.hydroxyacetophenone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.22=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.22.hydroxyacetophenone)
summary(glm.modB.22) 
Anova(glm.modB.22) #  0.3511
Bplot.22<-
  plot_model(glm.modB.22, type="pred", terms=c("logUW"), bias_correction=T,
             title="3,4-Dihydroxyacetophenone", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.23.unk<- data.frame(AMR.chem.std.B.all.indiv[23])
AMR.chem.std.B.all.23.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.23.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.23.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.23=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.23.unk)
summary(glm.modB.23) 
Anova(glm.modB.23) # 0.04654
Bplot.23<-
  plot_model(glm.modB.23, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 34", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.24.unk<- data.frame(AMR.chem.std.B.all.indiv[24])
AMR.chem.std.B.all.24.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.24.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.24.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.24=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.24.unk)
summary(glm.modB.24) 
Anova(glm.modB.24) # 0.4971
Bplot.24<-
  plot_model(glm.modB.24, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 35", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.25.unk<- data.frame(AMR.chem.std.B.all.indiv[25])
AMR.chem.std.B.all.25.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.25.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.25.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.25=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.25.unk)
summary(glm.modB.25) 
Anova(glm.modB.25) # 0.0000000000000984
Bplot.25<-
  plot_model(glm.modB.25, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 36", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.26.pyrrolidine<- data.frame(AMR.chem.std.B.all.indiv[26])
AMR.chem.std.B.all.26.pyrrolidine[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.26.pyrrolidine[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.26.pyrrolidine)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.26=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.26.pyrrolidine)
summary(glm.modB.26) 
Anova(glm.modB.26) # 0.8583
Bplot.26<-
  plot_model(glm.modB.26, type="pred", terms=c("logUW"), bias_correction=T,
             title="3-Methyl-4-propylpyrrolidine-2,5-dione", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.27.3acetophenone<- data.frame(AMR.chem.std.B.all.indiv[27])
AMR.chem.std.B.all.27.3acetophenone[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.27.3acetophenone[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.27.3acetophenone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.27=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.27.3acetophenone)
summary(glm.modB.27) 
Anova(glm.modB.27) # 0.8625
Bplot.27<-
  plot_model(glm.modB.27, type="pred", terms=c("logUW"), bias_correction=T,
             title="3-Ethylacetophenone", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.28.cyclopentanone<- data.frame(AMR.chem.std.B.all.indiv[28])
AMR.chem.std.B.all.28.cyclopentanone[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.28.cyclopentanone[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.28.cyclopentanone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.28=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.28.cyclopentanone)
summary(glm.modB.28) 
Anova(glm.modB.28) # 0.1955
Bplot.28<-
  plot_model(glm.modB.28, type="pred", terms=c("logUW"), bias_correction=T,
             title="2-Cyclopentylcyclopentanone", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.29.aminophenol<- data.frame(AMR.chem.std.B.all.indiv[29])
AMR.chem.std.B.all.29.aminophenol[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.29.aminophenol[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.29.aminophenol)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.29=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.29.aminophenol)
summary(glm.modB.29) 
Anova(glm.modB.29) # 0.2066
Bplot.29<-
  plot_model(glm.modB.29, type="pred", terms=c("logUW"), bias_correction=T,
             title="3-(Diethylamino)phenol", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.30.4acetophenone<- data.frame(AMR.chem.std.B.all.indiv[30])
AMR.chem.std.B.all.30.4acetophenone[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.30.4acetophenone[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.30.4acetophenone)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.30=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.30.4acetophenone)
summary(glm.modB.30) 
Anova(glm.modB.30) # 0.2941
Bplot.30<-
  plot_model(glm.modB.30, type="pred", terms=c("logUW"), bias_correction=T,
             title="4-Ethylacetophenone", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.31.aniline<- data.frame(AMR.chem.std.B.all.indiv[31])
AMR.chem.std.B.all.31.aniline[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.31.aniline[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.31.aniline)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.31=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.31.aniline)
summary(glm.modB.31) 
Anova(glm.modB.31) # 0.9927
Bplot.31<-
  plot_model(glm.modB.31, type="pred", terms=c("logUW"), bias_correction=T,
             title="4-Methoxy-3,5-dimethylaniline", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.32.unk<- data.frame(AMR.chem.std.B.all.indiv[32])
AMR.chem.std.B.all.32.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.32.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.32.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.32=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.32.unk)
summary(glm.modB.32) 
Anova(glm.modB.32) # 0.2251
Bplot.32<-
  plot_model(glm.modB.32, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 40", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.33.unk<- data.frame(AMR.chem.std.B.all.indiv[33])
AMR.chem.std.B.all.33.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.33.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.33.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.33=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.33.unk)
summary(glm.modB.33) 
Anova(glm.modB.33) # 0.9552
Bplot.33<-
  plot_model(glm.modB.33, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 42", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.34.unk<- data.frame(AMR.chem.std.B.all.indiv[34])
AMR.chem.std.B.all.34.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.34.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.34.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.34=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.34.unk)
summary(glm.modB.34) 
Anova(glm.modB.34) # 0.000007113
Bplot.34<-
  plot_model(glm.modB.34, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 43", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.35.unk<- data.frame(AMR.chem.std.B.all.indiv[35])
AMR.chem.std.B.all.35.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.35.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.35.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.35=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Date) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.35.unk)
summary(glm.modB.35) 
Anova(glm.modB.35) # 0.7610
Bplot.35<-
  plot_model(glm.modB.35, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 44", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.36.unk<- data.frame(AMR.chem.std.B.all.indiv[36])
AMR.chem.std.B.all.36.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.36.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.36.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.36=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.36.unk)
summary(glm.modB.36) 
Anova(glm.modB.36) # 0.5996
Bplot.36<-
  plot_model(glm.modB.36, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 49", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.37.unk<- data.frame(AMR.chem.std.B.all.indiv[37])
AMR.chem.std.B.all.37.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.37.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.37.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.37=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.37.unk)
summary(glm.modB.37) 
Anova(glm.modB.37) # 0.08415
Bplot.37<-
  plot_model(glm.modB.37, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 50", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)


#
AMR.chem.std.B.all.38.unk<- data.frame(AMR.chem.std.B.all.indiv[38])
AMR.chem.std.B.all.38.unk[,c(1:11,25,26)]<- lapply(AMR.chem.std.B.all.38.unk[,c(1:11,25,26)], factor)
names(AMR.chem.std.B.all.38.unk)<- c("Comp","Filename","Date","Group","Male","MRC","AMR","AMR_period","AMR_C3","Handler","Batch","Age","Seasonality","T.conc","Focal_sec","UW_count","UW_rate","z.Seasonality","RT","Area","Ncomp","SArea","logRPA","asinRPA","Comp.order","Comp.name","logUW")

glm.modB.38=lmer(logRPA ~ logUW +z.Seasonality +Group +(1|AMR) +(1|Male) +(1|Batch) +(1|Handler), 
                 data=AMR.chem.std.B.all.38.unk)
summary(glm.modB.38) 
Anova(glm.modB.38) # 0.22985
Bplot.38<-
  plot_model(glm.modB.38, type="pred", terms=c("logUW"), bias_correction=T,
             title="Unidentified comp 51", show.legend=F, line.size=1, dot.size=3, dodge=.7, show.data=T,  jitter=0.0003)





## B FDR adjustments on p-values ----
Bpval<- c(0.7486, 0.1684854, 0.1056, 0.3833, 0.1344, 0.1109, 0.9142, 0.5919, 0.0000000000000002, 0.2860522, 0.07081, 0.7480, 0.05153, 0.3646556, 0.006876, 0.9534, 0.1716, 0.0003855, 0.8731, 0.8547, 0.85350, 0.3511, 0.04654, 0.4971, 0.0000000000000984, 0.8583, 0.8625, 0.1955, 0.2066, 0.2941, 0.9927, 0.2251, 0.9552, 0.000007113, 0.7610, 0.5996, 0.08415, 0.22985)

Bpval.adjust<- round(p.adjust(Bpval, method="fdr", n=length(Bpval)), digits=3)
Bpval.adjust





## Suppl. Figure S4: B gridplots ----
empty.plot<-   ggplot() +
  theme_void(base_size=8) +
  geom_text(aes(0,0,label='')) +
  xlab(NULL)

#dev.new()
set_theme(base = theme_classic(base_size=8))

grid.draw(rbind(cbind(ggplotGrob(Bplot.01), ggplotGrob(Bplot.02), ggplotGrob(Bplot.03), ggplotGrob(Bplot.04), ggplotGrob(Bplot.05), ggplotGrob(Bplot.06), ggplotGrob(Bplot.07), size="first"),
                cbind(ggplotGrob(Bplot.08), ggplotGrob(Bplot.09), ggplotGrob(Bplot.10), ggplotGrob(Bplot.11), ggplotGrob(Bplot.12), ggplotGrob(Bplot.13), ggplotGrob(Bplot.14), size="first"),
                cbind(ggplotGrob(Bplot.15), ggplotGrob(Bplot.16), ggplotGrob(Bplot.17), ggplotGrob(Bplot.18), ggplotGrob(Bplot.19), ggplotGrob(Bplot.20), ggplotGrob(Bplot.21), size="first"),
                cbind(ggplotGrob(Bplot.22), ggplotGrob(Bplot.23), ggplotGrob(Bplot.24), ggplotGrob(Bplot.25), ggplotGrob(Bplot.26), ggplotGrob(Bplot.27), ggplotGrob(Bplot.28), size="first"),
                cbind(ggplotGrob(Bplot.29), ggplotGrob(Bplot.30), ggplotGrob(Bplot.31), ggplotGrob(Bplot.32), ggplotGrob(Bplot.33), ggplotGrob(Bplot.34), ggplotGrob(Bplot.35), size="first"),
                cbind(ggplotGrob(Bplot.36), ggplotGrob(Bplot.37), ggplotGrob(Bplot.38), ggplotGrob(empty.plot), ggplotGrob(empty.plot), ggplotGrob(empty.plot), ggplotGrob(empty.plot), size="first")))





# 7) Check model assumptions ----
mymod<- Full.mod.all

#library(DHARMa)
simulationOutput<- simulateResiduals(fittedModel=mymod, n=1000, refit=F, use.u=F)
testResiduals(simulationOutput) #tests unif+ouliers+disp

#dev.new()
#plot(simulationOutput)

#library(predictmeans)
residplot(mymod)



# end of script.
