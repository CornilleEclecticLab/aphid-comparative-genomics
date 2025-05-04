---
#title: "Chemosensory gene Random effects test"
#original script modify by: "Sergio Olvera"
---
  
library("readxl")
library("reshape2")
library("ggplot2")
library("tidyr")
library("ggpubr")
library("lme4")
library("MASS")



setwd("pathway")

All_sp<-read_excel("file.xlsx")
All_sp <- data.frame(All_sp)
head(All_sp)
data_mod<- melt(All_sp, id.vars= 'Type', 
                measure.vars='dNdS')

p<- ggplot(data_mod) +
  geom_boxplot(aes(x=All_sp$Type, y=All_sp$dNdS, fill=Type))+
  labs(subtitle="dNdS distribution")

print(p)

p + theme(axis.text.x = element_text(face="bold", color="black", 
                                     size=6, angle=45),
          axis.text.y = element_text(face="bold", color="black", 
                                     size=6, angle=45)) + scale_color_brewer(palette="Dark2") + theme(legend.position="bottom") + 
  theme(legend.title = element_text(size = 5), 
        legend.text = element_text(face="bold", color="black", size = 5)) + guides(fill=guide_legend(nrow=4, byrow=TRUE))

print(p)

mod.pois1 <- glm(dNdS~Types,family="quasipoisson", data=All_sp)
mod.pois1
summary(mod.pois1)

anova_pois<- anova(mod.pois1,test = "F")
anova_pois
summary(anova_pois)

glmer.nb <- lmer(dNdS ~ Species + (1 | Type), data = All_sp)
summary(glmer.nb)
anova_glmer<- anova(glmer.nb,test="F")
summary(anova_glmer)