#USDA project
#Tolerance data analyses

#Last updated: 07/24/2023 MLR

#Set working directory
setwd("G:/My Drive/Penn State/Research/File for R/Thesis/Tolerance") 

#Attach libraries
library(readxl)
library(ggplot2)
library(svglite)
library(psych)
library(agricolae)

#### Biofilm only ####
#Input data
bac<-read_excel("Tolerance_BAC_Results_NEW.xlsx", sheet=1, col_names = TRUE)


#APC
bac2<-subset(bac, logAPC!="NA")
bac2$logAPC<-as.numeric(bac2$logAPC)


#Calculate summary statistics
stat_bac2<-describeBy(bac2$logAPC, list(bac2$Trt, bac2$Timepoint), mat = TRUE) #Averages all reps by timepoint and treatment
stat_bac2$group2<-as.numeric(stat_bac2$group2)
stat_bac2 <- stat_bac2[order(stat_bac2$group1),]
stat_bac2$Number<-c(rep(1,5),rep(2,10),rep(3,5), rep(2,5), rep(3,10), rep(4,5), rep(3,5),rep(4,10), rep(5,5),rep(2,5), rep(3,10), rep(4,5))
stat_bac2<-subset(stat_bac2, vars != "NA")

#Plot
apc<-ggplot(stat_bac2, aes(x=group2, y=mean, group=group1, color=group1))+
  geom_point()+geom_line()+facet_grid(.~Number)+
  geom_text(data = subset(stat_bac2, group2 == 4), aes(label = group1, colour = group1, x = 4.05, y = mean), hjust = -.1, size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2)+
  geom_hline(yintercept = 0.6, color="grey80", linetype=2)+
  scale_x_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5))+
  scale_y_continuous(limits = c(0,8), breaks = c(0,1,2,3,4,5,6,7,8))+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=13), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=13,color='black'), axis.text.x = element_text(angle=90)) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA))+
  #theme(panel.grid.major = element_line(color = "grey50"), panel.grid.minor.y = element_line(color = "grey50"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position = "bottom")+
  ylab("log10 APC/peg")+xlab("Time (hours)")+
  ggtitle("Aerobic plate count - Tolerance")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
apc
ggsave("Tolerance_Biofilm_APC_average_2h.png", plot=apc, device="png", width=10, height=8, units="in", dpi=600)
ggsave("Tolerance_Biofilm_APC_average_2h.svg", plot=apc, device="svg", width=10, height=8, units="in", dpi=600)



#MPN
#Calculate summary statistics
stat_mpn<-describeBy(bac$logMPN, list(bac$Trt, bac$Timepoint), mat = TRUE) #By facility and treatment
stat_mpn$group2<-as.numeric(stat_mpn$group2)
stat_mpn <- stat_mpn[order(stat_mpn$group1),]
stat_mpn$Number<-c(rep(1,5),rep(2,10),rep(3,5), rep(2,5), rep(3,10), rep(4,5), rep(3,5),rep(4,10), rep(5,5),rep(2,5), rep(3,10), rep(4,5))

mpn<-ggplot(stat_mpn, aes(x=group2, y=mean, group=group1, color=group1))+
  geom_point()+geom_line()+facet_grid(.~Number)+
  geom_text(data = subset(stat_mpn, group2 == 4), aes(label = group1, colour = group1, x = 4.05, y = mean), hjust = -.1, size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2)+
  geom_hline(yintercept = 0.9, color="grey80", linetype=2)+
  scale_x_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5))+
  scale_y_continuous(limits = c(0,8), breaks = c(0,1,2,3,4,5,6,7,8))+
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
  ylab("log10 MPN/peg")+xlab("Time (hours)")+
  ggtitle("L. monocytogenes MPN - Tolerance")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
mpn
ggsave("Tolerance_Biofilm_MPN_average_2h.png", plot=mpn, device="png", width=10, height=8, units="in", dpi=600)
ggsave("Tolerance_Biofilm_MPN_average_2h.svg", plot=mpn, device="svg", width=10, height=8, units="in", dpi=600)



#### Planktonic only ####
#Input data
bac_plank<-read_excel("Tolerance_BAC_Results_NEW.xlsx", sheet=2, col_names = TRUE)


#APC

#Calculate summary statistics
stat_bac_plank<-describeBy(bac_plank$logAPC, list(bac_plank$Trt, bac_plank$Timepoint), mat = TRUE) #By facility and treatment
stat_bac_plank$group2<-as.numeric(stat_bac_plank$group2)
stat_bac_plank <- stat_bac_plank[order(stat_bac_plank$group1),]
stat_bac_plank$Number<-c(rep(1,5),rep(2,10),rep(3,5), rep(2,5), rep(3,10), rep(4,5), rep(3,5),rep(4,10), rep(5,5),rep(2,5), rep(3,10), rep(4,5))


#Plot
apc_bac_plank<-ggplot(stat_bac_plank, aes(x=group2, y=mean, group=group1, color=group1))+
  geom_point()+geom_line()+facet_grid(.~Number)+
  geom_text(data = subset(stat_bac_plank, group2 == 2), aes(label = group1, colour = group1, x = 2.05, y = mean), hjust = -.1, size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2)+
  geom_hline(yintercept = 0.6, color="grey80", linetype=2)+
  scale_x_continuous(limits = c(0,3), breaks = c(0,1,2,3))+
  scale_y_continuous(limits = c(0,8), breaks = c(0,1,2,3,4,5,6,7,8))+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=13), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=13,color='black'), axis.text.x = element_text(angle=90)) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  #theme(panel.grid.major = element_line(color = "grey50"), panel.grid.minor.y = element_line(color = "grey50"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position = "bottom")+
  ylab("log10 APC/ml")+xlab("Time (hours)")+
  ggtitle("Aerobic plate count - Tolerance Planktonic")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
apc_bac_plank
ggsave("ToleranceAPC_Planktonic.png", plot=apc_bac_plank, device="png", width=10, height=8, units="in", dpi=600)
ggsave("ToleranceAPC_Planktonic.svg", plot=apc_bac_plank, device="svg", width=10, height=8, units="in", dpi=600)



#MPN
#Calculate summary statistics
stat_mpn_bac_plank<-describeBy(bac_plank$logMPN, list(bac_plank$Trt, bac_plank$Timepoint), mat = TRUE) #By facility and treatment
stat_mpn_bac_plank$group2<-as.numeric(stat_mpn_bac_plank$group2)
stat_mpn_bac_plank <- stat_mpn_bac_plank[order(stat_mpn_bac_plank$group1),]
stat_mpn_bac_plank$Number<-c(rep(1,5),rep(2,10),rep(3,5), rep(2,5), rep(3,10), rep(4,5), rep(3,5),rep(4,10), rep(5,5),rep(2,5), rep(3,10), rep(4,5))


mpn_bac_plank<-ggplot(stat_mpn_bac_plank, aes(x=group2, y=mean, group=group1, color=group1))+
  geom_point()+geom_line()+facet_grid(.~Number)+
  geom_text(data = subset(stat_mpn_bac_plank, group2 == 2), aes(label = group1, colour = group1, x = 2.05, y = mean), hjust = -.1, size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2)+
  geom_hline(yintercept = 0.9, color="grey80", linetype=2)+
  scale_x_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5))+
  scale_y_continuous(limits = c(0,8), breaks = c(0,1,2,3,4,5,6,7,8))+
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
  ylab("log10 MPN/ml")+xlab("Time (hours)")+
  ggtitle("L. monocytogenes MPN - Tolerance Planktonic")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
mpn_bac_plank
ggsave("ToleranceMPN__Planktonic.png", plot=mpn_bac_plank, device="png", width=10, height=8, units="in", dpi=600)
ggsave("ToleranceMPN__Planktonic.svg", plot=mpn_bac_plank, device="svg", width=10, height=8, units="in", dpi=600)



#### Biofilm and Planktonic together ####
#Input data
bac_both<-read_excel("Tolerance_BAC_Results_NEW.xlsx", sheet=3, col_names = TRUE)

#Add column specifying combinations of temperature and concentration
bac_both$factorABC <- with(bac_both, interaction(Trt, Timepoint, Form))


#APC
bac_both_2<-subset(bac_both, logAPC!="NA")
bac_both_2$logAPC<-as.numeric(bac_both_2$logAPC)

#ANOVA for the interaction effect (since I care about each independent variable)
anova_apc<-aov(logAPC ~ factorABC, data=bac_both_2)
summary(anova_apc)


#tukey test
tukey_apc<-HSD.test(anova_apc, trt="factorABC") 
tukey_apc


#Calculate summary statistics
stat_bac_both_2<-describeBy(bac_both_2$logAPC, list(bac_both_2$Trt, bac_both_2$Timepoint, bac_both_2$Form), mat = TRUE) #By facility and treatment
stat_bac_both_2$group2<-as.numeric(stat_bac_both_2$group2)
stat_bac_both_2 <- stat_bac_both_2[order(stat_bac_both_2$group1),]
stat_bac_both_2$Number<-c(rep(1,10),rep(2,20),rep(3,10), rep(2,10), rep(3,20), rep(4,10), rep(3,10),rep(4,20), rep(5,10),rep(2,10), rep(3,20), rep(4,10))
stat_bac_both_2<-subset(stat_bac_both_2, vars !="NA")

#Plot
apc_bac_both<-ggplot(stat_bac_both_2, aes(x=group2, y=mean, group=interaction(group1,group3), color=group1, linetype=group3, shape=group3))+
  geom_point()+geom_line()+facet_wrap(vars(group1))+
  geom_text(data = subset(stat_bac_both_2, group2 == 2), aes(label = group1, colour = group1, x = 2.05, y = mean), hjust = -.1, size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, linetype=1)+
  geom_hline(yintercept = 0.6, color="grey80", linetype=2)+
  scale_x_continuous(limits = c(0,3), breaks = c(0,1,2,3))+
  scale_y_continuous(limits = c(0,8), breaks = c(0,1,2,3,4,5,6,7,8))+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_text(color='black', size=13), axis.ticks = element_line(color='black')) +
  theme(axis.title = element_text(size=13,color='black'), axis.text.x = element_text(angle=90)) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  #theme(panel.grid.major = element_line(color = "grey50"), panel.grid.minor.y = element_line(color = "grey50"))+
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.position = "bottom")+
  ylab("log10 APC/peg biofilms or log 10 APC/ml planktonic")+xlab("Time (hours)")+
  ggtitle("Aerobic plate count - Tolerance Both")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
apc_bac_both
ggsave("Tolerance_APC_Both.png", plot=apc_bac_both, device="png", width=10, height=8, units="in", dpi=600)
ggsave("Tolerance_APC_Both.svg", plot=apc_bac_both, device="svg", width=10, height=8, units="in", dpi=600)

#Anova by treatment
#Subset by treatment
T1<-subset(bac_both_2, Code=="T1" )
T6<-subset(bac_both_2, Code=="T6" )
T7<-subset(bac_both_2, Code=="T7" )
T8<-subset(bac_both_2, Code=="T8" )
T9<-subset(bac_both_2, Code=="T9" )
T10<-subset(bac_both_2, Code=="T10" )
T11<-subset(bac_both_2, Code=="T11" )
T12<-subset(bac_both_2, Code=="T12" )
T13<-subset(bac_both_2, Code=="T13" )
T14<-subset(bac_both_2, Code=="T14" )
T15<-subset(bac_both_2, Code=="T15" )
T16<-subset(bac_both_2, Code=="T16" )
T17<-subset(bac_both_2, Code=="T17" )
T18<-subset(bac_both_2, Code=="T18" )
T19<-subset(bac_both_2, Code=="T19" )
T20<-subset(bac_both_2, Code=="T20" )


#ANOVA for the interaction effect (since I care about each independent variable)
anova_apc_T1<-aov(logAPC ~ factorABC, data=T1)
summary(anova_apc_T1)

anova_apc_T6<-aov(logAPC ~ factorABC, data=T6)
summary(anova_apc_T6)

anova_apc_T7<-aov(logAPC ~ factorABC, data=T7)
summary(anova_apc_T7)

anova_apc_T8<-aov(logAPC ~ factorABC, data=T8)
summary(anova_apc_T8)

anova_apc_T9<-aov(logAPC ~ factorABC, data=T9)
summary(anova_apc_T9)

anova_apc_T10<-aov(logAPC ~ factorABC, data=T10)
summary(anova_apc_T10)

anova_apc_T11<-aov(logAPC ~ factorABC, data=T11)
summary(anova_apc_T11)

anova_apc_T12<-aov(logAPC ~ factorABC, data=T12)
summary(anova_apc_T12)

anova_apc_T13<-aov(logAPC ~ factorABC, data=T13)
summary(anova_apc_T13)

anova_apc_T14<-aov(logAPC ~ factorABC, data=T14)
summary(anova_apc_T14)

anova_apc_T15<-aov(logAPC ~ factorABC, data=T15)
summary(anova_apc_T15)

anova_apc_T16<-aov(logAPC ~ factorABC, data=T16)
summary(anova_apc_T16)

anova_apc_T17<-aov(logAPC ~ factorABC, data=T17)
summary(anova_apc_T17)

anova_apc_T18<-aov(logAPC ~ factorABC, data=T18)
summary(anova_apc_T18)

anova_apc_T19<-aov(logAPC ~ factorABC, data=T19)
summary(anova_apc_T19)

anova_apc_T20<-aov(logAPC ~ factorABC, data=T20)
summary(anova_apc_T20)

#tukey test
tukey_apc_T1<-HSD.test(anova_apc_T1, trt="factorABC") 
tukey_apc_T1

tukey_apc_T6<-HSD.test(anova_apc_T6, trt="factorABC") 
tukey_apc_T6

tukey_apc_T7<-HSD.test(anova_apc_T7, trt="factorABC") 
tukey_apc_T7

tukey_apc_T8<-HSD.test(anova_apc_T8, trt="factorABC") 
tukey_apc_T8

tukey_apc_T9<-HSD.test(anova_apc_T9, trt="factorABC") 
tukey_apc_T9

tukey_apc_T1<-HSD.test(anova_apc_T1, trt="factorABC") 
tukey_apc_T1

tukey_apc_T10<-HSD.test(anova_apc_T10, trt="factorABC") 
tukey_apc_T10

tukey_apc_T11<-HSD.test(anova_apc_T11, trt="factorABC") 
tukey_apc_T11

tukey_apc_T12<-HSD.test(anova_apc_T12, trt="factorABC") 
tukey_apc_T12

tukey_apc_T13<-HSD.test(anova_apc_T13, trt="factorABC") 
tukey_apc_T13

tukey_apc_T14<-HSD.test(anova_apc_T14, trt="factorABC") 
tukey_apc_T14

tukey_apc_T15<-HSD.test(anova_apc_T15, trt="factorABC") 
tukey_apc_T15

tukey_apc_T16<-HSD.test(anova_apc_T16, trt="factorABC") 
tukey_apc_T16

tukey_apc_T17<-HSD.test(anova_apc_T17, trt="factorABC") 
tukey_apc_T17

tukey_apc_T18<-HSD.test(anova_apc_T18, trt="factorABC") 
tukey_apc_T18

tukey_apc_T19<-HSD.test(anova_apc_T19, trt="factorABC") 
tukey_apc_T19

tukey_apc_T20<-HSD.test(anova_apc_T20, trt="factorABC") 
tukey_apc_T20

tukey_groups_apc<-bind_rows(tukey_apc_T1$groups, tukey_apc_T6$groups,tukey_apc_T7$groups,tukey_apc_T8$groups,
                            tukey_apc_T9$groups,tukey_apc_T10$groups,tukey_apc_T11$groups,tukey_apc_T12$groups,
                            tukey_apc_T13$groups,tukey_apc_T14$groups,tukey_apc_T15$groups,tukey_apc_T16$groups,
                            tukey_apc_T17$groups,tukey_apc_T18$groups,tukey_apc_T19$groups,tukey_apc_T20$groups)

tukey_groups_apc<-tukey_groups_apc[ order(row.names(tukey_groups_apc)), ]
tukey_groups_apc$Code<-c(rep("L+F",9),rep("L+M+F",9),rep("L+M",9),rep("L+P",9),rep("L+P+F",9),rep("L+P+M",9),
                         rep("L+P+M+F",9),rep("L+P+X",9),rep("L+P+X+F",9),rep("L+P+X+M+F",9),
                         rep("L+P+X+M",9),rep("L+X+F",9),rep("L+X+M",9),rep("L+X+M+F",9),rep("L+X",9),rep("L",9))
tukey_groups_apc$Timepoint<-rep(c(0.25,0.5,0.5,0,0,1,1,2,2),16)
tukey_groups_apc$Form<-rep(c(rep(c("Planktonic","Biofilm"),4),"Planktonic"),16)


tukey_means_apc<-bind_rows(tukey_apc_T9$means,tukey_apc_T15$means,tukey_apc_T8$means,tukey_apc_T6$means,
                           tukey_apc_T12$means,tukey_apc_T11$means,tukey_apc_T18$means,tukey_apc_T10$means,
                           tukey_apc_T17$means,tukey_apc_T20$means,tukey_apc_T16$means,tukey_apc_T14$means,
                           tukey_apc_T13$means,tukey_apc_T19$means,tukey_apc_T7$means,tukey_apc_T1$means)

tukey_apc_sum<-bind_cols(tukey_groups_apc,tukey_means_apc)


#plot apc results with tukey letters
apc_bac_tukey<-ggplot(tukey_apc_sum, aes(x=Timepoint, y=`logAPC...1`, group=interaction(Form,Code), color=Code, linetype=Form, shape=Form))+
  geom_point()+geom_line()+facet_wrap(vars(Code))+
  geom_text(aes(label=groups), position=position_dodge(width=0), vjust=-1)+
  geom_errorbar(aes(ymin=`logAPC...1`-se, ymax=`logAPC...1`+se), width=.07, linetype=1)+
  geom_hline(yintercept = 0.9, color="grey80", linetype=2)+
  scale_x_continuous(limits = c(-0.15,2.5), breaks = c(0,0.5,1,1.5,2,2.5))+
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10), minor_breaks = c(1,3,5,7,9))+
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
  ylab("log10 CFU/peg biofilms or log 10 CFU/ml planktonic")+xlab("Time (hours)")+
  ggtitle("Total count - Tolerance Both")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
apc_bac_tukey
ggsave("Tolerance_apc_Both_Tukey.png", plot=apc_bac_tukey, device="png", width=10, height=8, units="in", dpi=600)
ggsave("Tolerance_apc_Both_Tukey.svg", plot=apc_bac_tukey, device="svg", width=10, height=8, units="in", dpi=600)





#MPN
#ANOVA for the interaction effect (since I care about each independent variable)
anova_mpn<-aov(logMPN ~ factorABC, data=bac_both)
summary(anova_mpn)

#tukey test
tukey_mpn<-HSD.test(anova_mpn, trt="factorABC") 
tukey_mpn

#Calculate summary statistics
stat_mpn_bac_both<-describeBy(bac_both$logMPN, list(bac_both$Trt, bac_both$Timepoint, bac_both$Form), mat = TRUE) #By facility and treatment
stat_mpn_bac_both$group2<-as.numeric(stat_mpn_bac_both$group2)
stat_mpn_bac_both <- stat_mpn_bac_both[order(stat_mpn_bac_both$group1),]
stat_mpn_bac_both$Number<-c(rep(1,10),rep(2,20),rep(3,10), rep(2,10), rep(3,20), rep(4,10), rep(3,10),rep(4,20), rep(5,10),rep(2,10), rep(3,20), rep(4,10))
stat_mpn_bac_both<-subset(stat_mpn_bac_both, vars !="NA")


mpn_bac_both<-ggplot(stat_mpn_bac_both, aes(x=group2, y=mean, group=interaction(group1,group3), color=group1, linetype=group3, shape=group3))+
  geom_point()+geom_line()+facet_wrap(vars(group1))+
  geom_text(data = subset(stat_mpn_bac_both, group2 == 2), aes(label = group1, colour = group1, x = 2.05, y = mean), hjust = -.1, size=3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2)+
  geom_hline(yintercept = 0.9, color="grey80", linetype=2)+
  scale_x_continuous(limits = c(0,3), breaks = c(0,1,2,3))+
  scale_y_continuous(limits = c(0,8), breaks = c(0,1,2,3,4,5,6,7,8))+
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
  ylab("log10 MPN/peg biofilms or log 10 MPN/ml planktonic")+xlab("Time (hours)")+
  ggtitle("L. monocytogenes MPN - Tolerance Both")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
mpn_bac_both
ggsave("Tolerance_MPN_Both.png", plot=mpn_bac_both, device="png", width=10, height=8, units="in", dpi=600)
ggsave("Tolerance_MPN_Both.svg", plot=mpn_bac_both, device="svg", width=10, height=8, units="in", dpi=600)

#Anova by treatment
#Subset by treatment
T1<-subset(bac_both, Code=="T1" )
T6<-subset(bac_both, Code=="T6" )
T7<-subset(bac_both, Code=="T7" )
T8<-subset(bac_both, Code=="T8" )
T9<-subset(bac_both, Code=="T9" )
T10<-subset(bac_both, Code=="T10" )
T11<-subset(bac_both, Code=="T11" )
T12<-subset(bac_both, Code=="T12" )
T13<-subset(bac_both, Code=="T13" )
T14<-subset(bac_both, Code=="T14" )
T15<-subset(bac_both, Code=="T15" )
T16<-subset(bac_both, Code=="T16" )
T17<-subset(bac_both, Code=="T17" )
T18<-subset(bac_both, Code=="T18" )
T19<-subset(bac_both, Code=="T19" )
T20<-subset(bac_both, Code=="T20" )


#ANOVA for the interaction effect (since I care about each independent variable)
anova_mpn_T1<-aov(logMPN ~ factorABC, data=T1)
summary(anova_mpn_T1)

anova_mpn_T6<-aov(logMPN ~ factorABC, data=T6)
summary(anova_mpn_T6)

anova_mpn_T7<-aov(logMPN ~ factorABC, data=T7)
summary(anova_mpn_T7)

anova_mpn_T8<-aov(logMPN ~ factorABC, data=T8)
summary(anova_mpn_T8)

anova_mpn_T9<-aov(logMPN ~ factorABC, data=T9)
summary(anova_mpn_T9)

anova_mpn_T10<-aov(logMPN ~ factorABC, data=T10)
summary(anova_mpn_T10)

anova_mpn_T11<-aov(logMPN ~ factorABC, data=T11)
summary(anova_mpn_T11)

anova_mpn_T12<-aov(logMPN ~ factorABC, data=T12)
summary(anova_mpn_T12)

anova_mpn_T13<-aov(logMPN ~ factorABC, data=T13)
summary(anova_mpn_T13)

anova_mpn_T14<-aov(logMPN ~ factorABC, data=T14)
summary(anova_mpn_T14)

anova_mpn_T15<-aov(logMPN ~ factorABC, data=T15)
summary(anova_mpn_T15)

anova_mpn_T16<-aov(logMPN ~ factorABC, data=T16)
summary(anova_mpn_T16)

anova_mpn_T17<-aov(logMPN ~ factorABC, data=T17)
summary(anova_mpn_T17)

anova_mpn_T18<-aov(logMPN ~ factorABC, data=T18)
summary(anova_mpn_T18)

anova_mpn_T19<-aov(logMPN ~ factorABC, data=T19)
summary(anova_mpn_T19)

anova_mpn_T20<-aov(logMPN ~ factorABC, data=T20)
summary(anova_mpn_T20)

#tukey test
tukey_mpn_T1<-HSD.test(anova_mpn_T1, trt="factorABC") 
tukey_mpn_T1

tukey_mpn_T6<-HSD.test(anova_mpn_T6, trt="factorABC") 
tukey_mpn_T6

tukey_mpn_T7<-HSD.test(anova_mpn_T7, trt="factorABC") 
tukey_mpn_T7

tukey_mpn_T8<-HSD.test(anova_mpn_T8, trt="factorABC") 
tukey_mpn_T8

tukey_mpn_T9<-HSD.test(anova_mpn_T9, trt="factorABC") 
tukey_mpn_T9

tukey_mpn_T1<-HSD.test(anova_mpn_T1, trt="factorABC") 
tukey_mpn_T1

tukey_mpn_T10<-HSD.test(anova_mpn_T10, trt="factorABC") 
tukey_mpn_T10

tukey_mpn_T11<-HSD.test(anova_mpn_T11, trt="factorABC") 
tukey_mpn_T11

tukey_mpn_T12<-HSD.test(anova_mpn_T12, trt="factorABC") 
tukey_mpn_T12

tukey_mpn_T13<-HSD.test(anova_mpn_T13, trt="factorABC") 
tukey_mpn_T13

tukey_mpn_T14<-HSD.test(anova_mpn_T14, trt="factorABC") 
tukey_mpn_T14

tukey_mpn_T15<-HSD.test(anova_mpn_T15, trt="factorABC") 
tukey_mpn_T15

tukey_mpn_T16<-HSD.test(anova_mpn_T16, trt="factorABC") 
tukey_mpn_T16

tukey_mpn_T17<-HSD.test(anova_mpn_T17, trt="factorABC") 
tukey_mpn_T17

tukey_mpn_T18<-HSD.test(anova_mpn_T18, trt="factorABC") 
tukey_mpn_T18

tukey_mpn_T19<-HSD.test(anova_mpn_T19, trt="factorABC") 
tukey_mpn_T19

tukey_mpn_T20<-HSD.test(anova_mpn_T20, trt="factorABC") 
tukey_mpn_T20

tukey_groups_mpn<-bind_rows(tukey_mpn_T1$groups, tukey_mpn_T6$groups,tukey_mpn_T7$groups,tukey_mpn_T8$groups,
                            tukey_mpn_T9$groups,tukey_mpn_T10$groups,tukey_mpn_T11$groups,tukey_mpn_T12$groups,
                            tukey_mpn_T13$groups,tukey_mpn_T14$groups,tukey_mpn_T15$groups,tukey_mpn_T16$groups,
                            tukey_mpn_T17$groups,tukey_mpn_T18$groups,tukey_mpn_T19$groups,tukey_mpn_T20$groups)

tukey_groups_mpn<-tukey_groups_mpn[ order(row.names(tukey_groups_mpn)), ]
tukey_groups_mpn$Code<-c(rep("L+F",9),rep("L+M+F",9),rep("L+M",9),rep("L+P",9),rep("L+P+F",9),rep("L+P+M",9),
                         rep("L+P+M+F",9),rep("L+P+X",9),rep("L+P+X+F",9),rep("L+P+X+M+F",9),
                         rep("L+P+X+M",9),rep("L+X+F",9),rep("L+X+M",9),rep("L+X+M+F",9),rep("L+X",9),rep("L",9))
tukey_groups_mpn$Timepoint<-rep(c(0.25,0.5,0.5,0,0,1,1,2,2),16)
tukey_groups_mpn$Form<-rep(c(rep(c("Planktonic","Biofilm"),4),"Planktonic"),16)


tukey_means_mpn<-bind_rows(tukey_mpn_T9$means,tukey_mpn_T15$means,tukey_mpn_T8$means,tukey_mpn_T6$means,
                                tukey_mpn_T12$means,tukey_mpn_T11$means,tukey_mpn_T18$means,tukey_mpn_T10$means,
                                tukey_mpn_T17$means,tukey_mpn_T20$means,tukey_mpn_T16$means,tukey_mpn_T14$means,
                                tukey_mpn_T13$means,tukey_mpn_T19$means,tukey_mpn_T7$means,tukey_mpn_T1$means)

tukey_mpn_sum<-bind_cols(tukey_groups_mpn,tukey_means_mpn)


#plot mpn results with tukey letters
mpn_bac_tukey<-ggplot(tukey_mpn_sum, aes(x=Timepoint, y=`logMPN...1`, group=interaction(Form,Code), color=Code, linetype=Form, shape=Form))+
  geom_point()+geom_line()+facet_wrap(vars(Code))+
  geom_text(aes(label=groups), position=position_dodge(width=0), vjust=-1)+
  geom_errorbar(aes(ymin=`logMPN...1`-se, ymax=`logMPN...1`+se), width=.07, linetype=1)+
  geom_hline(yintercept = 0.9, color="grey80", linetype=2)+
  scale_x_continuous(limits = c(-0.15,2.5), breaks = c(0,0.5,1,1.5,2,2.5))+
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10), minor_breaks = c(1,3,5,7,9))+
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
  ylab("log10 MPN/peg biofilms or log 10 MPN/ml planktonic")+xlab("Time (hours)")+
  ggtitle("L. monocytogenes MPN - Tolerance Both")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
mpn_bac_tukey
ggsave("Tolerance_MPN_Both_Tukey.png", plot=mpn_bac_tukey, device="png", width=10, height=8, units="in", dpi=600)
ggsave("Tolerance_MPN_Both_Tukey.svg", plot=mpn_bac_tukey, device="svg", width=10, height=8, units="in", dpi=600)


#### Model death curves ####

#Import data
bac_both_diff<-read_excel("Tolerance_BAC_Results_NEW.xlsx", sheet=4, col_names = TRUE)
bac_both_diff$logdiffAPC<-as.numeric(bac_both_diff$logdiffAPC)

#Subset data
#Planktonic
T1_plank<-subset(bac_both_diff, Code=="T1" & Form=="Planktonic")
T6_plank<-subset(bac_both_diff, Code=="T6" & Form=="Planktonic")
T7_plank<-subset(bac_both_diff, Code=="T7" & Form=="Planktonic")
T8_plank<-subset(bac_both_diff, Code=="T8" & Form=="Planktonic")
T9_plank<-subset(bac_both_diff, Code=="T9" & Form=="Planktonic")
T10_plank<-subset(bac_both_diff, Code=="T10" & Form=="Planktonic")
T11_plank<-subset(bac_both_diff, Code=="T11" & Form=="Planktonic")
T12_plank<-subset(bac_both_diff, Code=="T12" & Form=="Planktonic")
T13_plank<-subset(bac_both_diff, Code=="T13" & Form=="Planktonic")
T14_plank<-subset(bac_both_diff, Code=="T14" & Form=="Planktonic")
T15_plank<-subset(bac_both_diff, Code=="T15" & Form=="Planktonic")
T16_plank<-subset(bac_both_diff, Code=="T16" & Form=="Planktonic")
T17_plank<-subset(bac_both_diff, Code=="T17" & Form=="Planktonic")
T18_plank<-subset(bac_both_diff, Code=="T18" & Form=="Planktonic")
T19_plank<-subset(bac_both_diff, Code=="T19" & Form=="Planktonic")
T20_plank<-subset(bac_both_diff, Code=="T20" & Form=="Planktonic")

#Biofilm
T1_biof<-subset(bac_both_diff, Code=="T1" & Form=="Biofilm")
T6_biof<-subset(bac_both_diff, Code=="T6" & Form=="Biofilm")
T7_biof<-subset(bac_both_diff, Code=="T7" & Form=="Biofilm")
T8_biof<-subset(bac_both_diff, Code=="T8" & Form=="Biofilm")
T9_biof<-subset(bac_both_diff, Code=="T9" & Form=="Biofilm")
T10_biof<-subset(bac_both_diff, Code=="T10" & Form=="Biofilm")
T11_biof<-subset(bac_both_diff, Code=="T11" & Form=="Biofilm")
T12_biof<-subset(bac_both_diff, Code=="T12" & Form=="Biofilm")
T13_biof<-subset(bac_both_diff, Code=="T13" & Form=="Biofilm")
T14_biof<-subset(bac_both_diff, Code=="T14" & Form=="Biofilm")
T15_biof<-subset(bac_both_diff, Code=="T15" & Form=="Biofilm")
T16_biof<-subset(bac_both_diff, Code=="T16" & Form=="Biofilm")
T17_biof<-subset(bac_both_diff, Code=="T17" & Form=="Biofilm")
T18_biof<-subset(bac_both_diff, Code=="T18" & Form=="Biofilm")
T19_biof<-subset(bac_both_diff, Code=="T19" & Form=="Biofilm")
T20_biof<-subset(bac_both_diff, Code=="T20" & Form=="Biofilm")


#Calculate the d values at different temperatures

#Fit log-linear
#Planktonic
T1_plank_L <- lm(logDiffMPN ~Timepoint, T1_plank)
summary(T1_plank_L)
rsquare.T1_plank_L <- cor(T1_plank$logDiffMPN, predict(T1_plank_L ))^2

T6_plank_L <- lm(logDiffMPN ~Timepoint, T6_plank)
summary(T6_plank_L)
rsquare.T6_plank_L <- cor(T6_plank$logDiffMPN, predict(T6_plank_L ))^2

T7_plank_L <- lm(logDiffMPN ~Timepoint, T7_plank)
summary(T7_plank_L)
rsquare.T7_plank_L <- cor(T7_plank$logDiffMPN, predict(T7_plank_L ))^2

T8_plank_L <- lm(logDiffMPN ~Timepoint, T8_plank)
summary(T8_plank_L)
rsquare.T8_plank_L <- cor(T8_plank$logDiffMPN, predict(T8_plank_L ))^2

T9_plank_L <- lm(logDiffMPN ~Timepoint, T9_plank)
summary(T9_plank_L)
rsquare.T9_plank_L <- cor(T9_plank$logDiffMPN, predict(T9_plank_L ))^2

T10_plank_L <- lm(logDiffMPN ~Timepoint, T10_plank)
summary(T10_plank_L)
rsquare.T10_plank_L <- cor(T10_plank$logDiffMPN, predict(T10_plank_L ))^2

T11_plank_L <- lm(logDiffMPN ~Timepoint, T11_plank)
summary(T11_plank_L)
rsquare.T11_plank_L <- cor(T11_plank$logDiffMPN, predict(T11_plank_L ))^2

T12_plank_L <- lm(logDiffMPN ~Timepoint, T12_plank)
summary(T12_plank_L)
rsquare.T12_plank_L <- cor(T12_plank$logDiffMPN, predict(T12_plank_L ))^2

T13_plank_L <- lm(logDiffMPN ~Timepoint, T13_plank)
summary(T13_plank_L)
rsquare.T13_plank_L <- cor(T13_plank$logDiffMPN, predict(T13_plank_L ))^2

T14_plank_L <- lm(logDiffMPN ~Timepoint, T14_plank)
summary(T14_plank_L)
rsquare.T14_plank_L <- cor(T14_plank$logDiffMPN, predict(T14_plank_L ))^2

T15_plank_L <- lm(logDiffMPN ~Timepoint, T15_plank)
summary(T15_plank_L)
rsquare.T15_plank_L <- cor(T15_plank$logDiffMPN, predict(T15_plank_L ))^2

T16_plank_L <- lm(logDiffMPN ~Timepoint, T16_plank)
summary(T16_plank_L)
rsquare.T16_plank_L <- cor(T16_plank$logDiffMPN, predict(T16_plank_L ))^2

T17_plank_L <- lm(logDiffMPN ~Timepoint, T17_plank)
summary(T17_plank_L)
rsquare.T17_plank_L <- cor(T17_plank$logDiffMPN, predict(T17_plank_L ))^2

T18_plank_L <- lm(logDiffMPN ~Timepoint, T18_plank)
summary(T18_plank_L)
rsquare.T18_plank_L <- cor(T18_plank$logDiffMPN, predict(T18_plank_L ))^2

T19_plank_L <- lm(logDiffMPN ~Timepoint, T19_plank)
summary(T19_plank_L)
rsquare.T19_plank_L <- cor(T19_plank$logDiffMPN, predict(T19_plank_L ))^2

T20_plank_L <- lm(logDiffMPN ~Timepoint, T20_plank)
summary(T20_plank_L)
rsquare.T20_plank_L <- cor(T20_plank$logDiffMPN, predict(T20_plank_L ))^2


#Biofilm
T1_biof_L <- lm(logDiffMPN ~Timepoint, T1_biof)
summary(T1_biof_L)
rsquare.T1_biof_L <- cor(T1_biof$logDiffMPN, predict(T1_biof_L ))^2

T6_biof_L <- lm(logDiffMPN ~Timepoint, T6_biof)
summary(T6_biof_L)
rsquare.T6_biof_L <- cor(T6_biof$logDiffMPN, predict(T6_biof_L ))^2

T7_biof_L <- lm(logDiffMPN ~Timepoint, T7_biof)
summary(T7_biof_L)
rsquare.T7_biof_L <- cor(T7_biof$logDiffMPN, predict(T7_biof_L ))^2

T8_biof_L <- lm(logDiffMPN ~Timepoint, T8_biof)
summary(T8_biof_L)
rsquare.T8_biof_L <- cor(T8_biof$logDiffMPN, predict(T8_biof_L ))^2

T9_biof_L <- lm(logDiffMPN ~Timepoint, T9_biof)
summary(T9_biof_L)
rsquare.T9_biof_L <- cor(T9_biof$logDiffMPN, predict(T9_biof_L ))^2

T10_biof_L <- lm(logDiffMPN ~Timepoint, T10_biof)
summary(T10_biof_L)
rsquare.T10_biof_L <- cor(T10_biof$logDiffMPN, predict(T10_biof_L ))^2

T11_biof_L <- lm(logDiffMPN ~Timepoint, T11_biof)
summary(T11_biof_L)
rsquare.T11_biof_L <- cor(T11_biof$logDiffMPN, predict(T11_biof_L ))^2

T12_biof_L <- lm(logDiffMPN ~Timepoint, T12_biof)
summary(T12_biof_L)
rsquare.T12_biof_L <- cor(T12_biof$logDiffMPN, predict(T12_biof_L ))^2

T13_biof_L <- lm(logDiffMPN ~Timepoint, T13_biof)
summary(T13_biof_L)
rsquare.T13_biof_L <- cor(T13_biof$logDiffMPN, predict(T13_biof_L ))^2

T14_biof_L <- lm(logDiffMPN ~Timepoint, T14_biof)
summary(T14_biof_L)
rsquare.T14_biof_L <- cor(T14_biof$logDiffMPN, predict(T14_biof_L ))^2

T15_biof_L <- lm(logDiffMPN ~Timepoint, T15_biof)
summary(T15_biof_L)
rsquare.T15_biof_L <- cor(T15_biof$logDiffMPN, predict(T15_biof_L ))^2

T16_biof_L <- lm(logDiffMPN ~Timepoint, T16_biof)
summary(T16_biof_L)
rsquare.T16_biof_L <- cor(T16_biof$logDiffMPN, predict(T16_biof_L ))^2

T17_biof_L <- lm(logDiffMPN ~Timepoint, T17_biof)
summary(T17_biof_L)
rsquare.T17_biof_L <- cor(T17_biof$logDiffMPN, predict(T17_biof_L ))^2

T18_biof_L <- lm(logDiffMPN ~Timepoint, T18_biof)
summary(T18_biof_L)
rsquare.T18_biof_L <- cor(T18_biof$logDiffMPN, predict(T18_biof_L ))^2

T19_biof_L <- lm(logDiffMPN ~Timepoint, T19_biof)
summary(T19_biof_L)
rsquare.T19_biof_L <- cor(T19_biof$logDiffMPN, predict(T19_biof_L ))^2

T20_biof_L <- lm(logDiffMPN ~Timepoint, T20_biof)
summary(T20_biof_L)
rsquare.T20_biof_L <- cor(T20_biof$logDiffMPN, predict(T20_biof_L ))^2



library(dplyr)
library(broom)

loglin_diff<-bind_rows(T1_plank_L$coefficients,T6_plank_L$coefficients,T7_plank_L$coefficients,T8_plank_L$coefficients,
                  T9_plank_L$coefficients,T10_plank_L$coefficients,T11_plank_L$coefficients,T12_plank_L$coefficients,
                  T13_plank_L$coefficients,T14_plank_L$coefficients,T15_plank_L$coefficients,T16_plank_L$coefficients,
                  T17_plank_L$coefficients,T18_plank_L$coefficients,T19_plank_L$coefficients,T20_plank_L$coefficients,
                  T1_biof_L$coefficients,T6_biof_L$coefficients,T7_biof_L$coefficients,T8_biof_L$coefficients,
                  T9_biof_L$coefficients,T10_biof_L$coefficients,T11_biof_L$coefficients,T12_biof_L$coefficients,
                  T13_biof_L$coefficients,T14_biof_L$coefficients,T15_biof_L$coefficients,T16_biof_L$coefficients,
                  T17_biof_L$coefficients,T18_biof_L$coefficients,T19_biof_L$coefficients,T20_biof_L$coefficients)



loglin_stats<-bind_rows(glance(T1_plank_L),glance(T6_plank_L),glance(T7_plank_L),glance(T8_plank_L),
                  glance(T9_plank_L),glance(T10_plank_L),glance(T11_plank_L),glance(T12_plank_L),
                  glance(T13_plank_L),glance(T14_plank_L),glance(T15_plank_L),glance(T16_plank_L),
                  glance(T17_plank_L),glance(T18_plank_L),glance(T19_plank_L),glance(T20_plank_L),
                  glance(T1_biof_L),glance(T6_biof_L),glance(T7_biof_L),glance(T8_biof_L),
                  glance(T9_biof_L),glance(T10_biof_L),glance(T11_biof_L),glance(T12_biof_L),
                  glance(T13_biof_L),glance(T14_biof_L),glance(T15_biof_L),glance(T16_biof_L),
                  glance(T17_biof_L),glance(T18_biof_L),glance(T19_biof_L),glance(T20_biof_L))


model_loglin<-bind_cols(loglin_diff,loglin_stats)
model_loglin$Code<-rep(c("T1","T6","T7","T8","T9","T10","T11","T12","T13","T14","T15","T16","T17","T18","T19","T20"),2)

model_loglin$Form<-c(rep("Planktonic",16),rep("Biofilm",16))

write.csv(model_loglin, file="loglineal models tolerance 20230731.csv")

