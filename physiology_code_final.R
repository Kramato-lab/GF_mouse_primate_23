
library(nlme)
library(car)
library(mosaic)

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Projects/GF_Mouse_Human_v_Primate/GF_Mouse_Human/First_experiment/Final_analysis")

data<-read.csv("Phys_combined_corrected_adults.csv")

percent_weight<-lm(percent_weight_gain~Species, data=data)
summary(aov(percent_weight))
TukeyHSD(aov(percent_weight))

fat <- lm(Fat ~ Species, data = data)
summary(aov(fat))
TukeyHSD(fat)

lw <- lm(Liver_norm ~ Species, data = data)
summary(aov(lw))
TukeyHSD(aov(lw))

pw <- lm(Pancreas_norm ~ Species, data = data)
summary(aov(pw))
TukeyHSD(pw)

bw <- lm(Brain_norm ~ Species, data = data)
summary(aov(bw))
TukeyHSD(aov(bw))

gly<-lm(glycogen~Species, data=data)
summary(aov(gly))
TukeyHSD(aov(gly))

data2<-read.csv("adult_food_intake.csv")
food<-lm(Food_per_body~Species, data=data2)
summary(aov(food))
TukeyHSD(aov(food))

gluc <- lm(Glucose ~ Species, data = data)
summary(aov(gluc))
TukeyHSD(aov(gluc))

AUC <- lm(AUC ~ Species, data = data)
summary(aov(AUC))
TukeyHSD(aov(AUC))

ALP <- lm(ALP ~ Species, data = data)
summary(aov(ALP))
TukeyHSD(aov(ALP))

logAST <- lm(logAST ~ Species, data = data)
summary(aov(logAST))
TukeyHSD(aov(logAST))

logALT <- lm(logALT ~ Species, data = data)
summary(aov(logALT))
TukeyHSD(aov(logALT))

chol <- lm(Cholesterol ~ Species, data = data)
summary(aov(chol))
TukeyHSD(aov(chol))

tri <- lm(logTriglycerides ~ Species, data = data)
summary(aov(tri))
TukeyHSD(aov(tri))

HDL <- lm(HDL ~ Species, data = data)
summary(aov(HDL))
TukeyHSD(aov(HDL))


AA <- lm(AA ~ Species, data = data)
summary(aov(AA))
TukeyHSD(aov(AA))

PA <- lm(PA ~ Species, data = data)
summary(aov(PA))
TukeyHSD(aov(PA))

BA <- lm(BA ~ Species, data = data)
summary(aov(BA))
TukeyHSD(aov(BA))

VA <- lm(VA ~ Species, data = data)
summary(aov(VA))
TukeyHSD(aov(VA))

Total<- lm(Total_SCFA ~ Species, data = data)
summary(aov(Total))
TukeyHSD(aov(Total))

perAA <- lm(Per_AA ~ Species, data = data)
summary(aov(perAA))
TukeyHSD(aov(perAA))

perPA <- lm(Per_PA ~ Species, data = data)
summary(aov(perPA))
TukeyHSD(aov(perPA))

perBA <- lm(Per_BA ~ Species, data = data)
summary(aov(perBA))
TukeyHSD(aov(perBA))

perVA <- lm(Per_VA ~ Species, data = data)
summary(aov(perVA))
TukeyHSD(aov(perVA))

