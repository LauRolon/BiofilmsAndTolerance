#Display MIC results
#USDA project

#Last updated: 03/30/2024 MLR

setwd("G:/My Drive/Penn State/Research/File for R/Thesis/03_MIC")

library(readxl)
library(ggplot2)
library(svglite)

#Import data
mic<-read_excel("MIC results.xlsx", sheet=4, col_names = TRUE)


#For PAA paper
#Subset PAA and Sterilex results
MIC_sub<-subset(mic, ATM!="BAC")

#Make barplot
MIC_plot<-ggplot(MIC_sub, aes(x=Isolate, y=MIC, fill=Organism))+
  geom_bar(stat='identity', color='black', position='dodge')+
  facet_grid(ATM~Organism, scales = "free")+
  theme(legend.text=element_text(size=8, face='italic'), legend.title= element_text(size=8), legend.position = 'bottom') +
  theme(axis.title=element_text(size=8),axis.text.x = element_text(color='black', size=7, angle=90), axis.text.y = element_text(color='black', size=8)) +
  theme(panel.grid.major.y=element_line(color='grey50', size=1), panel.grid.minor.y=element_line(color='grey50', size=0.5))+
  theme(panel.grid.major.x=element_line(color='transparent'), panel.grid.minor.x=element_line(color='transparent'))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_rect(color='black', fill="white"), strip.text.x = element_text(size=7, face='italic'),strip.text.y = element_text(size=8),
        panel.border = element_rect(color="black", fill=NA))+
  ylab("MIC (ppm)")
MIC_plot

ggsave("MIC plot PAA.png", plot=MIC_plot, device="png", width=8, height=5, units="in", dpi=600)
ggsave("MIC plot PAA.svg", plot=MIC_plot, device="svg", width=8, height=5, units="in", dpi=600)

