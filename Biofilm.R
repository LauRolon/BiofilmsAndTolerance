#Data analysis for biofilm formation assays
#MLR 
#Last updated: 01/03/2024


#Set working directory
setwd("G:/My Drive/Penn State/Research/File for R/Thesis/Biofilm") #New PC

#Attach libraries
library(ggplot2)
library(readxl)
library(dplyr)
library(psych)
library(agricolae)
library(svglite)
library(rstatix)
library(ggpubr)
library(FSA)
library(rcompanion)

#Upload data
data<-read_excel("Biofilm formation_Results.xlsx", sheet=9, col_names = TRUE)
data<-arrange(data, Trt)


#### Crystal violet ####
#ANOVA
#Check ANOVA assumptions
#Normality check using anova model residuals
cv_lm<-lm(CV ~ Trt, data = data)
ggqqplot(residuals(cv_lm)) # All values follow normality 
shapiro_test(residuals(cv_lm)) #non-significant Shapiro test for normality

#Check homogeneity of variace assumption
plot(cv_lm, 1)
data %>% levene_test(CV ~ Trt) #non-significant Levene test


#Run ANOVA model
anova_cv<-aov(CV ~ Trt, data=data)
summary(anova_cv) #p<2e-16 

#tukey test
tukey_cv<-HSD.test(anova_cv, trt="Trt") 

#Make dataframe with Tukey groups
tukey_cv_groups<-tukey_cv$groups
tukey_cv_groups$Trt<-as.character(rownames(tukey_cv_groups))
tukey_cv_groups_order<-tukey_cv_groups[with(tukey_cv_groups, order(Trt)),]

#Calculate mean and sd by treatment
cv_stat<-describeBy(data$CV, list(data$Trt), mat = TRUE) #By facility and treatment
cv_stat$Tukey<-tukey_cv_groups_order$groups
cv_stat$Label<-unique(data$Label)
cv_stat$item<-as.numeric(cv_stat$item)

#Plot
cv_plot<-ggplot(cv_stat, aes(x=reorder(Label, item), y=mean))+
  geom_bar(stat='identity', color='black', fill="#8B0AA5")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=Tukey), position=position_dodge(width=0.9), vjust=-3.5)+
  ylab("Absorvance at 570 nm") + xlab("Assemblage")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text.x = element_text(color='black', size=10, angle=90), axis.text.y = element_text(color='black', size=10),axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=10,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=10),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Crystal violet assay")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8,9))
cv_plot
ggsave("Crytal Violet Assay.png", plot=cv_plot, device="png", width=11, height=8, units="in", dpi=600)
ggsave("Crytal Violet Assay.svg", plot=cv_plot, device="svg", width=11, height=8, units="in", dpi=600)

#### Aerobic plate count ####
#ANOVA
#Check ANOVA assumptions
#Normality check using anova model residuals
apc_lm<-lm(logAPCpeg ~ Trt, data = data)
ggqqplot(residuals(apc_lm)) # Not all values follow normality
shapiro_test(residuals(apc_lm)) #significant Shapiro test for normality

#Check homogeneity of variace assumption
plot(apc_lm, 1)
data %>% levene_test(logAPCpeg ~ Trt) #significant Levene test

#ANOVA Assumptions are not met
#Run Kruskal-Wallis test

apc.kruskal <- data %>% kruskal_test(logAPCpeg ~ Trt)
apc.kruskal #p = 0.0332

#Calculate effect size
data %>% kruskal_effsize(logAPCpeg ~ Trt)

# Pairwise comparisons
DT.apc = dunnTest(logAPCpeg ~ Trt, data=data, method = "none")   

DunnLetters.apc<-cldList(P.adj ~ Comparison,
        data = DT.apc$res,
        threshold = 0.05, remove.zero = FALSE)



# #ANOVA - do not run
# anova_apc<-aov(logAPCpeg ~ Trt, data=data)
# summary(anova_apc)  #p= 1.12e-06
# 
# #tukey test
# tukey_apc<-HSD.test(anova_apc, trt="Trt") 
# 
# 
# #Make dataframe with Tukey groups
# tukey_apc_groups<-tukey_apc$groups
# tukey_apc_groups$Trt<-as.character(rownames(tukey_apc_groups))
# tukey_apc_groups_order<-tukey_apc_groups[with(tukey_apc_groups, order(Trt)),]

#Calculate mean and sd by treatment
apc_stat<-describeBy(data$logAPCpeg, list(data$Trt), mat = TRUE) #By facility and treatment
apc_stat$CLD<-DunnLetters.apc$Letter
apc_stat$Label<-unique(data$Label)
apc_stat$item<-as.numeric(apc_stat$item)

#Plot
apc_plot<-ggplot(apc_stat, aes(x=reorder(Label, item), y=mean))+
  geom_bar(stat='identity', color='black', fill="#B83289")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=CLD), position=position_dodge(width=0.9), vjust=-1)+
  ylab("log10 CFU/peg") + xlab("Assemblage")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text.x = element_text(color='black', size=10, angle=90), axis.text.y = element_text(color='black', size=10),axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=10,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,linewidth = 1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=10),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("Aerobic plate counts")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_y_continuous(limits=c(0,10),breaks = c(1,2,3,4,5,6,7,8,9,10))
apc_plot
ggsave("APC_Biofilm_Dunn.png", plot=apc_plot, device="png", width=11, height=8, units="in", dpi=600)
ggsave("APC_Biofilm_Dunn.svg", plot=apc_plot, device="svg", width=11, height=8, units="in", dpi=600)


#### L. monocytogenes MPN ####
data_lm<-subset(data, logMPNpeg!="NA")
data_lm$logMPNpeg<-as.numeric(data_lm$logMPNpeg)

#ANOVA
#Check ANOVA assumptions
#Normality check using anova model residuals
mpn_lm<-lm(logMPNpeg ~ Trt, data = data_lm)
ggqqplot(residuals(mpn_lm)) # All values follow normality
shapiro_test(residuals(mpn_lm)) #Non-significant Shapiro test for normality

#Check homogeneity of variace assumption
plot(mpn_lm, 1)
data_lm %>% levene_test(logMPNpeg ~ Trt) #non-significant Levene test



#Run ANOVA
anova_lm<-aov(logMPNswab ~ Trt, data=data_lm)
summary(anova_lm)  #p= 1.27e-10 

#tukey test
tukey_lm<-HSD.test(anova_lm, trt="Trt") 


#Make dataframe with Tukey groups
tukey_lm_groups<-tukey_lm$groups
tukey_lm_groups$Trt<-as.character(rownames(tukey_lm_groups))
tukey_lm_groups_order<-tukey_lm_groups[with(tukey_lm_groups, order(Trt)),]

#Calculate mean and sd by treatment
lm_stat<-describeBy(data_lm$logMPNswab, list(data_lm$Trt), mat = TRUE) #By facility and treatment
lm_stat$Tukey<-tukey_lm_groups_order$groups
lm_stat$Label<-unique(data_lm$Label)
lm_stat$item<-as.numeric(lm_stat$item)

#Plot
mpn_plot<-ggplot(lm_stat, aes(x=reorder(Label, item), y=mean))+
  geom_bar(stat='identity', color='black', fill="#f48849")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=Tukey), position=position_dodge(width=0.9), vjust=-1)+
  ylab("log10 MPN/peg") + xlab("Assemblage")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text.x = element_text(color='black', size=10, angle=90), axis.text.y = element_text(color='black', size=10),axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=10,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=10),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("L. monocytogenes MPN")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_y_continuous(limits = c(0,9), breaks = c(1,2,3,4,5,6,7,8,9))
mpn_plot
ggsave("MPN_Biofilm.png", plot=mpn_plot, device="png", width=11, height=8, units="in", dpi=600)
ggsave("MPN_Biofilm.svg", plot=mpn_plot, device="svg", width=11, height=8, units="in", dpi=600)


#### Increase in Lm concentration after biofilm formation ####
data_lm$logInc<-as.numeric(data_lm$logInc)

#ANOVA
#Check ANOVA assumptions
#Normality check using anova model residuals
inc_lm<-lm(logInc ~ Trt, data = data_lm)
ggqqplot(residuals(inc_lm)) # All values follow normality
shapiro_test(residuals(inc_lm)) #Non-significant Shapiro test for normality

#Check homogeneity of variace assumption
plot(inc_lm, 1)
data_lm %>% levene_test(logInc ~ Trt) #Non-significant Levene test


#Run ANOVA
anova_incr<-aov(logInc ~ Trt, data=data_lm)
summary(anova_incr)  #p= <2e-16 

#tukey test
tukey_incr<-HSD.test(anova_incr, trt="Trt") 


#Make dataframe with Tukey groups
tukey_incr_groups<-tukey_incr$groups
tukey_incr_groups$Trt<-as.character(rownames(tukey_incr_groups))
tukey_incr_groups_order<-tukey_incr_groups[with(tukey_incr_groups, order(Trt)),]

#Calculate mean and sd by treatment
incr_stat<-describeBy(data_lm$logInc, list(data_lm$Trt), mat = TRUE) #By facility and treatment
incr_stat$Tukey<-tukey_incr_groups_order$groups
incr_stat$Label<-unique(data_lm$Label)
incr_stat$item<-as.numeric(incr_stat$item)

#Plot
incr_plot<-ggplot(incr_stat, aes(x=reorder(Label, item), y=mean))+
  geom_bar(stat='identity', color='black', fill="#febd2a")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
  geom_text(aes(label=Tukey), position=position_dodge(width=0.9), vjust=-1)+
  ylab("Increase in concentration of Lm in pegs - log10 MPN/peg (Final/Initial)") + xlab("Assemblage")+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text.x = element_text(color='black', size=10, angle=90), axis.text.y = element_text(color='black', size=10),axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=10,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=10),
        panel.border = element_rect(color="black", fill=NA))+
  ggtitle("L. monocytogenes Increase")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_y_continuous(limits = c(-1,2), breaks = c(-1,0,1,22))
incr_plot
ggsave("Increase_Biofilm.png", plot=incr_plot, device="png", width=11, height=8, units="in", dpi=600)
ggsave("Increase_Biofilm.svg", plot=incr_plot, device="svg", width=11, height=8, units="in", dpi=600)


#### CLSM ####
#Upload data
clsm<-read_excel("CLSM_Results.xlsx", sheet=3, col_names = TRUE)
clsm<-arrange(clsm, Code)

#Biomass
#ANOVA
#Check ANOVA assumptions
#Normality check using anova model residuals
biomass_lm<-lm(Biomass ~ Code, data = clsm)
ggqqplot(residuals(biomass_lm)) # All values follow normality
shapiro_test(residuals(biomass_lm)) #Non-significant Shapiro test for normality

#Check homogeneity of variace assumption
plot(biomass_lm, 1)
clsm %>% levene_test(Biomass ~ Code) #significant Levene test

#Run Kruskal Wallis 
biomass.kruskal <- clsm %>% kruskal_test(Biomass ~ Code)
biomass.kruskal #p = 0.00000000575

#Calculate effect size
clsm %>% kruskal_effsize(Biomass ~ Code)

# Pairwise comparisons
DT.biomass = dunnTest(Biomass ~ Code, data=clsm, method = "none")   

DunnLetters.biomass<-cldList(P.adj ~ Comparison,
                         data = DT.biomass$res,
                         threshold = 0.05, remove.zero = FALSE)



# #ANOVA - do not use
# anova_biomass<-aov(Biomass ~ Code, data=clsm)
# summary(anova_biomass) #3.42e-14 *** 
# 
# #tukey test
# tukey_biomass<-HSD.test(anova_biomass, trt="Code") 
# 
# #Make dataframe with Tukey groups
# tukey_biomass_groups<-tukey_biomass$groups
# tukey_biomass_groups$Code<-as.character(rownames(tukey_biomass_groups))
# tukey_biomass_groups_order<-tukey_biomass_groups[with(tukey_biomass_groups, order(Code)),]
# 
# #Calculate mean and sd by treatment
# biomass_stat<-describeBy(clsm$Biomass, list(clsm$Code), mat = TRUE) #By facility and treatment
# biomass_stat$Tukey<-tukey_biomass_groups_order$groups
# biomass_stat$Label<-unique(clsm$Label)
# biomass_stat$item<-as.numeric(biomass_stat$item)
# 
# #Plot
# biomass_plot<-ggplot(biomass_stat, aes(x=reorder(Label,item), y=mean))+
#   geom_bar(stat='identity', color='black', fill="#8B0AA5")+
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
#   #geom_text(aes(label=Tukey), position=position_dodge(width=0.9), vjust=-3.5)+
#   ylab("Biomass (um3/um2)") + xlab("Assemblage")+
#   theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
#   theme(axis.text.x = element_text(color='black', size=10, angle=90), axis.text.y = element_text(color='black', size=10),axis.ticks = element_line(color='black')) +
#   theme(axis.title = element_text(size=10,color='black')) +
#   theme(panel.background = element_rect(fill='white', color = NA),
#         plot.background = element_rect(fill = 'white',color = NA), 
#         panel.border = element_rect(color='black',fill = NA,size=1))+
#   #theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
#   theme(strip.background= element_blank(), strip.text = element_text(size=10),
#         panel.border = element_rect(color="black", fill=NA))+
#   ggtitle("CLSM - Biomass")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
#   scale_y_continuous(breaks = c(1,2,3,4))
# biomass_plot
# ggsave("CLSM-Biomass.png", plot=biomass_plot, device="png", width=8, height=6, units="in", dpi=600)
# ggsave("CLSM-Biomass.svg", plot=biomass_plot, device="svg", width=8, height=6, units="in", dpi=600)
# 

#Roughness
#ANOVA
roughness_lm<-lm(Roughness ~ Code, data = clsm)
ggqqplot(residuals(roughness_lm)) # All values follow normality
shapiro_test(residuals(roughness_lm)) #Non-significant Shapiro test for normality

#Check homogeneity of variace assumption
plot(roughness_lm, 1)
clsm %>% levene_test(Roughness ~ Code) #significant Levene test

# Run Kruskal Wallis test
roughness.kruskal <- clsm %>% kruskal_test(Roughness ~ Code)
roughness.kruskal #p = 9.38*10^-10

#Calculate effect size
clsm %>% kruskal_effsize(Roughness ~ Code)

# Pairwise comparisons
DT.roughness = dunnTest(Roughness ~ Code, data=clsm, method = "none")   

DunnLetters.roughness<-cldList(P.adj ~ Comparison,
                             data = DT.roughness$res,
                             threshold = 0.05, remove.zero = FALSE)


# #ANOVA - do not use
# anova_roughness<-aov(Roughness ~ Code, data=clsm)
# summary(anova_roughness) #<2e-16 ***
# 
# #tukey test
# tukey_roughness<-HSD.test(anova_roughness, trt="Code") 
# 
# #Make dataframe with Tukey groups
# tukey_roughness_groups<-tukey_roughness$groups
# tukey_roughness_groups$Code<-as.character(rownames(tukey_roughness_groups))
# tukey_roughness_groups_order<-tukey_roughness_groups[with(tukey_roughness_groups, order(Code)),]
# 
# #Calculate mean and sd by treatment
# roughness_stat<-describeBy(clsm$Roughness, list(clsm$Code), mat = TRUE) #By facility and treatment
# roughness_stat$Tukey<-tukey_roughness_groups_order$groups
# roughness_stat$Label<-unique(clsm$Label)
# roughness_stat$item<-as.numeric(roughness_stat$item)
# 
# #Plot
# roughness_plot<-ggplot(roughness_stat, aes(x=reorder(Label,item), y=mean))+
#   geom_bar(stat='identity', color='black', fill="#8B0AA5")+
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9))+
#   geom_text(aes(label=Tukey), position=position_dodge(width=0.9), vjust=-3.5)+
#   ylab("Roughness") + xlab("Assemblage")+
#   theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
#   theme(axis.text.x = element_text(color='black', size=10, angle=90), axis.text.y = element_text(color='black', size=10),axis.ticks = element_line(color='black')) +
#   theme(axis.title = element_text(size=10,color='black')) +
#   theme(panel.background = element_rect(fill='white', color = NA),
#         plot.background = element_rect(fill = 'white',color = NA), 
#         panel.border = element_rect(color='black',fill = NA,size=1))+
#   #theme(panel.grid.major.y = element_line(color = "grey80"), panel.grid.minor.y = element_line(color = "grey80"))+
#   theme(strip.background= element_blank(), strip.text = element_text(size=10),
#         panel.border = element_rect(color="black", fill=NA))+
#   ggtitle("CLSM - Roughness")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
#   #scale_y_continuous(breaks = c(0,1,2,3,4))
# roughness_plot
# ggsave("CLSM-Roughness.png", plot=roughness_plot, device="png", width=8, height=6, units="in", dpi=600)
# ggsave("CLSM-Roughness.svg", plot=roughness_plot, device="svg", width=8, height=6, units="in", dpi=600)

##Make dataframe for CLSM results
#Calculate mean and sd by treatment
biomass_stat<-describeBy(clsm$Biomass, list(clsm$Code), mat = TRUE) #By facility and treatment
roughness_stat<-describeBy(clsm$Roughness, list(clsm$Code), mat = TRUE)

#Add Dunn's letters
biomass_stat$CLR<-DunnLetters.biomass$Letter
roughness_stat$CLR<-DunnLetters.roughness$Letter

write.csv(biomass_stat, "CLSM_biomass_Dunn.csv")
write.csv(roughness_stat, "CLSM_roughness_Dunn.csv")

#### Growth curves ####
setwd("G:/My Drive/Penn State/Research/File for R/Thesis/Tolerance") #New PC
growth_curves<-read_excel("Tolerance_BAC_Results_NEW.xlsx", sheet=5, col_names = TRUE)

growth_plot<-ggplot(growth_curves, aes(x=Time, y=Average, group=Strain, color=Group))+
  geom_point()+geom_line()+facet_wrap(vars(Group))+
  geom_text(data = subset(growth_curves, Time == 69.1), aes(label = Strain, colour = Group, x = 80, y = Average), hjust = -.1, size=3)+
  scale_x_continuous(limits = c(0,100), breaks = c(0,20,40,60,80,100))+
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1))+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=13), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  #theme(panel.grid.major = element_line(color = "grey50"), panel.grid.minor.y = element_line(color = "grey50"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position = "bottom")+
  ylab("OD600")+xlab("Time (hours)")+
  ggtitle("Growth curves - Isolares")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
growth_plot
ggsave("Preliminary growth curve.png", plot=growth_plot, device="png", width=10, height=8, units="in", dpi=600)
ggsave("Preliminary growth curve.svg", plot=growth_plot, device="svg", width=10, height=8, units="in", dpi=600)
