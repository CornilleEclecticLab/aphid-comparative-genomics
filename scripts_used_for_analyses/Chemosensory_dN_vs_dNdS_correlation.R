---
#title: "Chemosensory dN vs dNdS correlation test"
#original script modify by: "Sergio Olvera"
---
  
library(ggplot2)
library(viridis)
library(dplyr)
library(reshape2)

setwd("")
dN_vs_Omega <- read.table("file.txt", dec=",",sep = "\t",header = T)

p2 <- ggplot(dN_vs_Omega, aes(x=dS, y=dNdS, color=Species, fill=Species)) + 
  geom_jitter(width = 0.2) +  labs(subtitle="dNdS vs dS distribution") + scale_color_manual(values=c("Species"="#DF0101")) + theme_classic()

dN_vs_Omega$Species <- factor(dN_vs_Omega$Species, levels=c('Species'))

# Scatter plot with correlation coefficient
#:::::::::::::::::::::::::::::::::::::::::::::::::
p3 <- ggscatter(dN_vs_Omega, x = "dS", y = "dNdS", color = "Type",
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
)


# Add correlation coefficient
p3 + stat_cor(method = "pearson", label.x = 3, label.y = 30)

# Specify the number of decimal places of precision for p and r
# Using 3 decimal places for the p-value and
# 2 decimal places for the correlation coefficient (r)
p3 + stat_cor(p.accuracy = 0.001, r.accuracy = 0.01)


# Use R2 instead of R
p3 <- ggscatter(dN_vs_Omega, x = "dS", y = "dNdS", color = "Type", add = "reg.line") +
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 3
  )

