---
#title: "Chemosensory gene dNdS GLM test"
#original script modify by: "Sergio Olvera"
---
library("readxl")
library(reshape2)
library(ggplot2)
library("tidyr")

setwd("")

All_sp<-read_excel("file.xlsx")
All_sp <- data.frame(All_sp)
head(All_sp)
data_mod<- melt(All_sp, id.vars= 'Species', 
                measure.vars='dNdS')

p<- ggplot(data_mod) +
  geom_boxplot(aes(x=All_sp$Species, y=All_sp$dNdS, fill=Species))+
  labs(subtitle="dNdS distribution")

All_sp$Species <- factor(All_sp$Species, levels=c('Species'))

print(p)

p + theme(axis.text.x = element_text(face="bold", color="black", 
                                     size=10, angle=45),
          axis.text.y = element_text(face="bold", color="black", 
                                     size=10, angle=45)) + scale_color_brewer(palette="Dark2")

p<- ggplot(data_mod) +
  geom_boxplot(aes(x=All_sp$Species, y=All_sp$dNdS, color=Species))+
  labs(subtitle="OG and Chemogenes dNdS distribution") + theme(axis.text.x = element_text(face="bold", color="black", 
                                                                                          size=10, angle=45),
                                                               axis.text.y = element_text(face="bold", color="black", 
                                                                                          size=10, angle=45))

p+scale_fill_manual(values=c("#DF0101", "#FF9326"))


mod.pois1 <- glm(dNdS~Species,family="quasipoisson", data=All_sp)
mod.pois1
summary(mod.pois1)

anova_pois<- anova(mod.pois1,test = "F")
anova_pois
summary(anova_pois)
