#Plots for WGS gene presence analysis
#Last updated: 12/1/2022


#Set working directory
setwd("C:/Users/lau_r/Google Drive/Penn State/Research/File for R/Thesis/WGS")

#Attach libraries
library(ggplot2)
library(readxl)
library(tidyverse)

#Import data
wgs<-read_excel("Lit review_WGS genes.xlsx", sheet=5, col_names = TRUE)

#Set Strain as factor
wgs$Strain<-as.factor(wgs$Strain)

#Subset for QAC genes plot
wgs_qac<-subset(wgs, Type=="QAC")

qac_genes<-ggplot(wgs_qac, aes(x=reorder(Gene, Function), y=Strain, fill=Presence))+
  geom_tile(color="black")+facet_grid(Family~Function, scales="free", space="free", labeller = label_wrap_gen(2))+
  theme(axis.text.x = element_text(angle=90, color="black", size=7, face='italic'), axis.text.y = element_text(color="black", size=7))+
  ylab("")+xlab("")+scale_fill_manual(values = c( "#97BC62FF", "#2C5F2D"))+
  theme(panel.background = element_rect(fill="white", color =NA), plot.background = element_rect(fill="white", color =NA))+
  theme(strip.background= element_rect(color='black', fill = "white"), strip.text.y = element_text(size=7, face='italic', angle=0),strip.text.x = element_text(size=7),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_discrete(limits=rev)
  
ggsave("QAC genes.png", plot=qac_genes, device="png", width=8, height=5, units="in", dpi=600)
ggsave("QAC genes.svg", plot=qac_genes, device="svg", width=8, height=5, units="in", dpi=600)

#Subset for BIOFILM genes plot
wgs_biofilm<-subset(wgs, Type=="Biofilm")

biofilm_genes<-ggplot(wgs_biofilm, aes(x=reorder(Gene, Function), y=Strain, fill=Presence))+
  geom_tile(color="black")+facet_grid(Family~Function, scales="free", space="free", labeller = label_wrap_gen(2))+
  theme(axis.text.x = element_text(angle=90, color="black", size=7, face='italic'), axis.text.y = element_text(color="black", size=7))+
  ylab("")+xlab("")+scale_fill_manual(values = c( "#D7C49EFF", "#755139FF"))+
  theme(panel.background = element_rect(fill="white", color =NA), plot.background = element_rect(fill="white", color =NA))+
  theme(strip.background= element_rect(color='black', fill = "white"), strip.text.y = element_text(size=3, face='italic', angle=0),strip.text.x = element_text(size=7),
        panel.border = element_rect(color="black", fill=NA))+
  scale_y_discrete(limits=rev)

ggsave("Biofilm genes.png", plot=biofilm_genes, device="png", width=8, height=5, units="in", dpi=600)
ggsave("Biofilm genes.svg", plot=biofilm_genes, device="svg", width=8, height=5, units="in", dpi=600)
