#USDA biofilm project microbiome data analysis 
#Got ASVs from DADA2 pipeline using default parameters

#Last updated: MLR 8/8/2023

#Set working directory
setwd("G:/My Drive/Penn State/Research/File for R/Thesis/Amplicon")

#Attach libraries
library(ggplot2)
library(dplyr)
library(phyloseq)
library(zCompositions)
library(compositions)
library(viridis)
library(svglite)
library(pairwiseAdonis)
library(decontam)

#### Import data ####
asvs_all<-read.csv('ASV_16s.csv', header = TRUE, row.names = 1)
taxon_all<-as.data.frame(read.csv('Taxon_16s.csv', header = TRUE, row.names = 1))
metadata_all<-read.csv('Metadata_biofilm.csv', header = TRUE, row.names = 1)

#Add '_unclassified' marker to NAs in taxon table
taxon_all$Phylum<-ifelse(is.na(taxon_all$Phylum), paste(taxon_all$Kingdom, "unclassified", sep = '_'), taxon_all$Phylum)
taxon_all$Class<-ifelse(is.na(taxon_all$Class), paste(taxon_all$Phylum, "unclassified", sep = '_'), taxon_all$Class)
taxon_all$Order<-ifelse(is.na(taxon_all$Order), paste(taxon_all$Class, "unclassified", sep = '_'), taxon_all$Order)
taxon_all$Family<-ifelse(is.na(taxon_all$Family), paste(taxon_all$Order, "unclassified", sep = '_'), taxon_all$Family)
taxon_all$Genus<-ifelse(is.na(taxon_all$Genus), paste(taxon_all$Family, "unclassified", sep = '_'), taxon_all$Genus)
taxon_all$Species<-ifelse(is.na(taxon_all$Species), paste(taxon_all$Genus, "unclassified", sep = '_'), taxon_all$Species)

#Remove extra _unclassified
taxon_all$Class<-gsub("_unclassified_unclassified", '_unclassified', taxon_all$Class)
taxon_all$Order<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Order)
taxon_all$Order<-gsub("_unclassified_unclassified", '_unclassified', taxon_all$Order)
taxon_all$Family<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Family)
taxon_all$Family<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Family)
taxon_all$Family<-gsub("_unclassified_unclassified", '_unclassified', taxon_all$Family)
taxon_all$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Genus)
taxon_all$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Genus)
taxon_all$Genus<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Genus)
taxon_all$Genus<-gsub("_unclassified_unclassified", '_unclassified', taxon_all$Genus)
taxon_all$Species<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Species)
taxon_all$Species<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Species)
taxon_all$Species<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Species)
taxon_all$Species<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_all$Species)
taxon_all$Species<-gsub("_unclassified_unclassified", '_unclassified', taxon_all$Species)

#Convert asv and taxon tables to matrix
asvs_all<-as.matrix(asvs_all)
taxon_all<-as.matrix(taxon_all)

#Make phyloseq object
phyloseq_16s<-phyloseq(otu_table(asvs_all, taxa_are_rows = FALSE), tax_table(taxon_all), sample_data(metadata_all))

#Remove Chloroplast and Mitochondria reads from ASV table
physeq_16s <- subset_taxa(phyloseq_16s, Order!="Chloroplast" )
physeq_16s <- physeq_16s %>% subset_taxa( Family!="Mitochondria" )


#Get ASV table from phyloseq object
asv.16s<-as.data.frame(t(otu_table(physeq_16s)))
tail(rowSums(asv.16s))

#Remove ASVs with zero counts in all samples
asv.16s<-asv.16s[ which(rowSums(asv.16s)>0),]
asv.16s<-t(asv.16s)

taxon_16s_clean<-as.matrix(tax_table(physeq_16s))

####  Identify ASVs that match positive control strains ####
#Need to do this step before decontaminating data because most of the ASVs in the positive controls are not found in dataset and decontam will consider them as contaminants

#Obtain positive control declared composition and sequences from Zymo website: https://s3.amazonaws.com/zymofiles/BioPool/ZymoBIOMICS.STD.refseq.v2.zip.  
#Match all ASVs of Pseudomonas, E. coli, Salmonella, Lactobacillus, Enterococcus, Staphylococcus, Listeria and Bacillus to Zymo reference genomes using Mega10 to identify ASVs that are organisms in positive control
#Select ASVs that match with 0 SNPs to reference seqs from Zymo
#ASV25 - Bacillus
#ASV48 - Enterococcus
#ASV31 - Escherichia
#ASV1 - Listeria
#ASV9 - Pseudomonas
#ASV29 Salmonella
#ASV33 - Staphylococcus
#ASV19 - Limosilactobacillus / formerly Lactobacillus

#Make Phyloseq with cleaned ASV
phyloseq_16s<-phyloseq(otu_table(asv.16s, taxa_are_rows = FALSE), tax_table(taxon_16s_clean), sample_data(metadata_all))

phyloseq_PC<-subset_samples(phyloseq_16s, Trt =="PosControl_Mock")
phyloseq_PCDNA<-subset_samples(phyloseq_16s, Trt =="PosControl_DNA")

phyloseq_PC_clean<-subset_taxa(phyloseq_PC, rownames(tax_table(phyloseq_PC))=="ASV25"|rownames(tax_table(phyloseq_PC))=="ASV48"|rownames(tax_table(phyloseq_PC))=="ASV31"|rownames(tax_table(phyloseq_PC))=="ASV1"|rownames(tax_table(phyloseq_PC))=="ASV9"|rownames(tax_table(phyloseq_PC))=="ASV29"|rownames(tax_table(phyloseq_PC))=="ASV33"|rownames(tax_table(phyloseq_PC))=="ASV19")
phyloseq_PCDNA_clean<-subset_taxa(phyloseq_PCDNA, rownames(tax_table(phyloseq_PCDNA))=="ASV25"|rownames(tax_table(phyloseq_PC))=="ASV48"|rownames(tax_table(phyloseq_PC))=="ASV31"|rownames(tax_table(phyloseq_PC))=="ASV1"|rownames(tax_table(phyloseq_PC))=="ASV9"|rownames(tax_table(phyloseq_PC))=="ASV29"|rownames(tax_table(phyloseq_PC))=="ASV33"|rownames(tax_table(phyloseq_PC))=="ASV19")

#See positive extraction control by rep
PC_RA_rep<-phyloseq_PC_clean%>%
  transform_sample_counts(function(x) x/sum(x)*100)%>%
  psmelt()

barplot_PC_reps<-ggplot(PC_RA_rep, aes(x=SampleID,y=Abundance, fill=Genus))+
  geom_bar(stat='identity', color='black')+
  geom_text(aes(label=Abundance), position = position_stack(vjust=0.5), size=4)+  ylab("Relative abundance (%)")+
  ggtitle("PC barplot by rep")+xlab("")+
  theme(axis.text.y=element_text(color='black', size=13), axis.ticks=element_line(color='black'),
        axis.text.x = element_text(angle=90, size=13, color='black'), axis.title.y = element_text(color='black', size=13))+
  theme(panel.background = element_rect(fill="grey99", color =NA),plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=10, angle=0, face = 'italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))

PCDNA_RA<-phyloseq_PCDNA_clean%>%
  transform_sample_counts(function(x) x/sum(x)*100)%>%
  psmelt()

PC_Zymo<-data.frame(OTU=c("ASV25","ASV48","ASV31","ASV1","ASV9","ASV29","ASV33","ASV19"),
                    Sample=rep("Zymo",8),
                    Genus=c("Bacillus","Enterococcus","Escherichia-Shigella","Listeria","Pseudomonas","Salmonella","Staphylococcus","Limosilactobacillus"),
                    Abundance=c(17.4,9.9,10.1,14.1,4.2,10.4,15.5,18.4))

PC_data_all<-bind_rows(PC_RA_rep,PCDNA_RA,PC_Zymo)
PC_data_all$Abundance<-round(PC_data_all$Abundance, digits = 1)
PC_data_all<-PC_data_all[order(PC_data_all$Sample),]
PC_data_all$SampleOrder<-c(rep(2,8),rep(3,8),rep(4,8),rep(5,8),rep(1,8))

barplot_PC_reps<-ggplot(PC_data_all, aes(x=reorder(Sample,SampleOrder),y=Abundance, fill=Genus))+
  geom_bar(stat='identity', color='black')+
  geom_text(aes(label=Abundance), position = position_stack(vjust=0.5), size=4)+  ylab("Relative abundance (%)")+
  ggtitle("PC barplot")+xlab("")+
  theme(axis.text.y=element_text(color='black', size=13), axis.ticks=element_line(color='black'),
        axis.text.x = element_text(angle=90, size=13, color='black'), axis.title.y = element_text(color='black', size=13))+
  theme(panel.background = element_rect(fill="grey99", color =NA),plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=10, angle=0, face = 'italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Positive controls_reps.png", plot =barplot_PC_reps, device="png", width=20, height=10, units="in",dpi=600)
ggsave("Positive controls_reps.svg", plot =barplot_PC_reps, device="svg", width=20, height=10, units="in",dpi=600)

#Calculate relative abundance and average across replicates
PC_RA<-phyloseq_PC_clean%>%
  transform_sample_counts(function(x) x/sum(x)*100)%>%
  psmelt()%>%
  group_by(OTU, Genus, Trt)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

PCDNA_RA<-phyloseq_PCDNA_clean%>%
  transform_sample_counts(function(x) x/sum(x)*100)%>%
  psmelt()%>%
  group_by(OTU, Genus, Trt)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

PC_Zymo<-data.frame(OTU=c("ASV25","ASV48","ASV31","ASV1","ASV9","ASV29","ASV33","ASV19"),
                    Genus=c("Bacillus","Enterococcus","Escherichia-Shigella","Listeria","Pseudomonas","Salmonella","Staphylococcus","Limosilactobacillus"),
                    Trt=rep("Zymo",8),
                    Mean=c(17.4,9.9,10.1,14.1,4.2,10.4,15.5,18.4))
                    
PC_data<-bind_rows(PC_RA,PCDNA_RA,PC_Zymo)
PC_data$Order<-c(rep(2,8),rep(3,8),rep(1,8))

#Plot with average 
barplot_PC_all<-ggplot(PC_data, aes(x=reorder(Trt,Order),y=Mean, fill=Genus))+
  geom_bar(stat='identity', color='black')+
  geom_text(aes(label=Mean), position = position_stack(vjust=0.5), size=4)+  ylab("Relative abundance (%)")+
  ggtitle("PC barplot")+xlab("")+
  theme(axis.text.y=element_text(color='black', size=13), axis.ticks=element_line(color='black'),
        axis.text.x = element_text(angle=90, size=13, color='black'), axis.title.y = element_text(color='black', size=13))+
  theme(panel.background = element_rect(fill="grey99", color =NA),plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=10, angle=0, face = 'italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Positive controls.png", plot =barplot_PC_all, device="png", width=20, height=10, units="in",dpi=600)
ggsave("Positive controls.svg", plot =barplot_PC_all, device="svg", width=20, height=10, units="in",dpi=600)




####Check negative control and decontaminate ASVs ####

phyloseq_16s<-phyloseq(otu_table(asv.16s, taxa_are_rows = FALSE), tax_table(taxon_16s_clean), sample_data(metadata_all))

#Subset negative controls
phyloseq_NC<-subset_samples(phyloseq_16s, Code =="NegControl")

#Make long table
phyloseq.NC<-psmelt(phyloseq_NC)

#Remove rows with less than 100 reads
phyloseq.NC<-subset(phyloseq.NC, Abundance>100)

#Plot reads of control by  NC
barplot_NC<-ggplot(phyloseq.NC, aes(x=reorder(OTU, desc(Abundance)),y=Abundance, fill=SampleID))+
  geom_bar(stat='identity', color='black')+ 
  geom_text(aes(label=Genus, angle=90, hjust=0, size=13))+  ylab("Reads (x10,000)")+
  scale_y_continuous(breaks= seq(0,75000, 5000), labels = function(x){x/10000}, limits=c(0,75000)) + 
  ggtitle("Negative control read counts")+xlab("")+
  theme(axis.text.y=element_text(color='black', size=13), axis.ticks=element_line(color='black'),
        axis.text.x = element_text(angle=90, size=13, color='black'), axis.title.y = element_text(color='black', size=13))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 1)) +
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill="grey99", color =NA),plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Controls_NC.png", plot =barplot_NC, device="png", width=20, height=10, units="in",dpi=600)


#Decontaminate sequencing data

#Detect contaminants with prevalence method (threshold=0.1)
sample_data(phyloseq_16s)$is.neg <- sample_data(phyloseq_16s)$Code == "NegControl"
contamdf.prev <- isContaminant(phyloseq_16s, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev$contaminant)
which(contamdf.prev$contaminant) 
#Note decontam with threshold 0.5 is removing a most of the ASVs that match the strains in the experiment

#Remove identified contaminants (prevalence method) from ASV table
phyloseq.noncontam <- prune_taxa(!contamdf.prev$contaminant, phyloseq_16s)
phyloseq.noncontam

#Remove negative and positive controls from ASV table
phyloseq_clean<-subset_samples(phyloseq.noncontam, Type!="Control" )

#Get metadata 
metadata_clean<-subset(metadata_all, Type != "Control")


#Get ASV table from phyloseq object
asv_clean<-t(as.data.frame(otu_table(phyloseq_clean)))
tail(rowSums(asv_clean))

#Remove ASVs with zero counts in all samples
asv_clean<-asv_clean[ which(rowSums(asv_clean)>0),]




#Compositional analysis of microbiome data -based on Microbiome Analysis in R. Chap 10. #
#Step 1: Convert ASV table to appropriate format. Following step requires samples on rows and ASVs in columns
head(t(asv_clean))

#Step 2: Replace zero values before clr transformation. Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
asv.n0<-t(cmultRepl(t(asv_clean), label=0, method="CZM", output="p-counts")) #

head(asv.n0) #output table needs to have samples in columns and ASVs in rows

#Note: Check the output to make sure there are no negative numbers. If samples or ASV are sparse, the CZM method will 
#add a negative number that interferes with the log normalization. If necessary use function below to convert negative values
#into positives
asv_n0<-ifelse(asv.n0 < 0, asv.n0*(-1), asv.n0)


#Step 3: Convert data to proportions
asv.n0_prop<-apply(asv_n0, 2, function(x) {x/sum(x)})

#Step 4: Perform abundance and sample filtering and deal sparsity
#Filter the data to remove all taxa with less than 0.00001% abundance in any sample
asv.n0_prop_f<-asv.n0[apply(asv.n0_prop, 1, min) > 0.00000000000000000001, ]
head(asv.n0_prop_f) #Check that samples are on columns and asvs in rows


#Step 5: perform CLR transformation
asv.n0.clr<-t(apply(asv.n0_prop, 2, function(x){log(x)-mean(log(x))}))
head(asv.n0.clr) #Check output table. Samples should be in rows and asvs in columns

#Step 6: Perform Singular Value Decomposition (PCA)
pc.clr<-prcomp(asv.n0.clr)

png("Screeplot - PCA.png")
par(mar=c(2,2,2,2))
screeplot(pc.clr, type='lines', main="PCA")
dev.off()

#Calculate total variance of the data
mvar.clr<-mvar(asv.n0.clr)

#Display results
row<-rownames(asv.n0.clr) #Make vector with sample names
pc_out<-as.data.frame(pc.clr$x[,1:2]) #Get PC1 and PC2
pc_out_meta<-as.data.frame(bind_cols(pc_out,metadata_clean)) #Add metadata information
row.names(pc_out_meta)<-row #Add rownames to dataframe
pc_out_meta$Rep<-as.factor(pc_out_meta$Rep)

# Make PCA plot
PCA <- ggplot(pc_out_meta, aes(x=PC1,y=PC2, color=Trt, shape=Rep))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13)) +
  theme(legend.position = 'bottom')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  theme(panel.background = element_rect(fill='transparent', color = NA),
        plot.background = element_rect(fill = 'transparent',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr$sdev[1]^2/mvar.clr*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr$sdev[2]^2/mvar.clr*100, digits=1), "%", sep="")) +
  ggtitle("PCA - Obj 2 by facility")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA
ggsave("PCA_Obj2.png", plot=PCA_Obj2, device="png", width=8, height=8, units="in", dpi=600)
ggsave("PCA_Obj2.svg", plot=PCA_Obj2, device="svg", width=8, height=8, units="in", dpi=600)

# # PERMANOVA
# #Calculate Aitchinson distance
# dist_obj2<-dist(asv.n0.clr_obj2, method='euclidean')
# 
# #By Facility
# permanova_obj2_fac<-pairwise.adonis2(dist_obj2~Facility, data=metadata_obj2_clean, perm = 999, p.adjust.m = 'bonferroni')
# permanova_obj2_fac
# 
# 
# #By Treatment
# permanova_obj2_trt<-pairwise.adonis2(dist_obj2~Treatment, data=metadata_obj2_clean, perm = 999, p.adjust.m = 'bonferroni')
# permanova_obj2_trt


# STACKED BARPLOT 
#Note: used compositional approach to transform the sample counts to compositions. 


#Transform sample counts into compositions
asv.n0.acomp<-as.data.frame(acomp(t(asv_n0)), total=1)
rowSums(asv.n0.acomp) #Verify rows sum to 1

#Make Phyloseq object
phyloseq_RA <- phyloseq(otu_table(asv.n0.acomp, taxa_are_rows = FALSE), tax_table(taxon_all), sample_data(metadata_clean))

#Make long format table from Phyloseq object
asv_long <- phyloseq_RA %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance as percentage
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

# Calculate mean relative abundance by Treatment for each ASV
asv_mean<-asv_long %>%
  group_by(OTU, Trt, PMA, Type, Family, Genus, SampleOrder)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

#Filter table to obtain only ASVs with over 2% in at least one sample
asv_over2abund <- filter(asv_mean, Mean>8)

#Stacked barplot by treatment at the ASV level
barplot_biofilms<-ggplot(asv_over2abund, aes(x=reorder(Trt,SampleOrder), y=Mean, fill=Family))+
  geom_bar(stat='identity', color='black')+facet_grid(Type~PMA, scales = "free", space = 'free')+
  geom_text(aes(label=OTU), position = position_stack(vjust=0.5), size=3)+
  ylim(0,100)+
  theme(axis.title = element_text(color='black'), axis.text.x=element_text(color='black', size=10), 
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 4)) +
  theme(legend.position = "bottom")+ylab("Relative Abundance (%)")+xlab("Treatment")+
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Microbiota of biofilms - Obj 2")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Barplots_biofilm ai.png", plot=barplot_biofilms, device="png", width=13, height=10, units="in", dpi=600)
ggsave("Barplots_biofilm ai.svg", plot=barplot_biofilms, device="svg", width=13, height=10, units="in", dpi=600)


#### Alternative pipeline ####

#Get metadata 
metadata_clean<-subset(metadata_all, Type != "Control")

#Subset only ASVs that match the strains used in the study
phyloseq_trt<-subset_samples(phyloseq_16s, Type !="Control")

phyloseq_trt_clean<-subset_taxa(phyloseq_trt, rownames(tax_table(phyloseq_trt))=="ASV5"|rownames(tax_table(phyloseq_trt))=="ASV20"|rownames(tax_table(phyloseq_trt))=="ASV1"|rownames(tax_table(phyloseq_trt))=="ASV15"|
                                  rownames(tax_table(phyloseq_trt))=="ASV58"|rownames(tax_table(phyloseq_trt))=="ASV32"|rownames(tax_table(phyloseq_trt))=="ASV17"|rownames(tax_table(phyloseq_trt))=="ASV42"
                                |rownames(tax_table(phyloseq_trt))=="ASV36"|rownames(tax_table(phyloseq_trt))=="ASV6"|rownames(tax_table(phyloseq_trt))=="ASV3"|rownames(tax_table(phyloseq_trt))=="ASV35"
                                |rownames(tax_table(phyloseq_trt))=="ASV4"|rownames(tax_table(phyloseq_trt))=="ASV308"|rownames(tax_table(phyloseq_trt))=="ASV45"|rownames(tax_table(phyloseq_trt))=="ASV30"|rownames(tax_table(phyloseq_trt))=="ASV18"|rownames(tax_table(phyloseq_trt))=="ASV2"
                                |rownames(tax_table(phyloseq_trt))=="ASV26"|rownames(tax_table(phyloseq_trt))=="ASV7"|rownames(tax_table(phyloseq_trt))=="ASV21"|rownames(tax_table(phyloseq_trt))=="ASV44"|rownames(tax_table(phyloseq_trt))=="ASV13"|rownames(tax_table(phyloseq_trt))=="ASV11")

#Subset samples by treatment and keep taxa that were added to each treatment
phyloseq_T1<-phyloseq_trt_clean %>%
  subset_samples(Code =="T1")%>%
  subset_taxa(Family =="Listeriaceae")

phyloseq_T2<-phyloseq_trt_clean %>%
  subset_samples(Code =="T2")%>%
  subset_taxa(Family =="Pseudomonadaceae")

phyloseq_T3<-phyloseq_trt_clean %>%
  subset_samples(Code =="T3")%>%
  subset_taxa(Family =="Xanthomonadaceae")

phyloseq_T4<-phyloseq_trt_clean %>%
  subset_samples(Code =="T4")%>%
  subset_taxa(Family =="Microbacteriaceae")

phyloseq_T5<-phyloseq_trt_clean %>%
  subset_samples(Code =="T5")%>%
  subset_taxa(Family =="Flavobacteriaceae" | Family=="Sphingobacteriaceae")

phyloseq_T6<-phyloseq_trt_clean %>%
  subset_samples(Code =="T6")%>%
  subset_taxa(Family =="Listeriaceae"| Family =="Pseudomonadaceae")

phyloseq_T7<-phyloseq_trt_clean %>%
  subset_samples(Code =="T7")%>%
  subset_taxa(Family =="Listeriaceae"| Family =="Xanthomonadaceae")

phyloseq_T8<-phyloseq_trt_clean %>%
  subset_samples(Code =="T8")%>%
  subset_taxa(Family =="Listeriaceae"| Family =="Microbacteriaceae")

phyloseq_T9<-phyloseq_trt_clean %>%
  subset_samples(Code =="T9")%>%
  subset_taxa(Family =="Listeriaceae"| Family =="Flavobacteriaceae" | Family=="Sphingobacteriaceae")

phyloseq_T10<-phyloseq_trt_clean %>%
  subset_samples(Code =="T10")%>%
  subset_taxa(Family =="Listeriaceae"| Family =="Pseudomonadaceae"| Family =="Xanthomonadaceae")

phyloseq_T11<-phyloseq_trt_clean %>%
  subset_samples(Code =="T11")%>%
  subset_taxa(Family =="Listeriaceae"| Family =="Pseudomonadaceae"| Family =="Microbacteriaceae")

phyloseq_T12<-phyloseq_trt_clean %>%
  subset_samples(Code =="T12")%>%
  subset_taxa(Family =="Listeriaceae"| Family =="Pseudomonadaceae"| Family =="Flavobacteriaceae"| Family=="Sphingobacteriaceae")

phyloseq_T13<-phyloseq_trt_clean %>%
  subset_samples(Code =="T13")%>%
  subset_taxa(Family =="Listeriaceae"| Family =="Xanthomonadaceae"| Family =="Microbacteriaceae")

phyloseq_T14<-phyloseq_trt_clean %>%
  subset_samples(Code =="T14")%>%
  subset_taxa(Family =="Listeriaceae"| Family =="Xanthomonadaceae"| Family =="Flavobacteriaceae" | Family=="Sphingobacteriaceae")

phyloseq_T15<-phyloseq_trt_clean %>%
  subset_samples(Code =="T15")%>%
  subset_taxa(Family =="Listeriaceae"| Family =="Microbacteriaceae"| Family =="Flavobacteriaceae" | Family=="Sphingobacteriaceae")

phyloseq_T16<-phyloseq_trt_clean %>%
  subset_samples(Code =="T16")%>%
  subset_taxa(Family =="Listeriaceae"| Family =="Pseudomonadaceae"| Family =="Xanthomonadaceae"| Family =="Microbacteriaceae")

phyloseq_T17<-phyloseq_trt_clean %>%
  subset_samples(Code =="T17")%>%
  subset_taxa(Family =="Listeriaceae"| Family =="Pseudomonadaceae"| Family =="Xanthomonadaceae"| Family =="Flavobacteriaceae" | Family=="Sphingobacteriaceae")

phyloseq_T18<-phyloseq_trt_clean %>%
  subset_samples(Code =="T18")%>%
  subset_taxa(Family =="Listeriaceae"| Family =="Pseudomonadaceae"| Family =="Microbacteriaceae"|Family =="Flavobacteriaceae" | Family=="Sphingobacteriaceae")

phyloseq_T19<-phyloseq_trt_clean %>%
  subset_samples(Code =="T19")%>%
  subset_taxa(Family =="Listeriaceae"| Family =="Xanthomonadaceae"| Family =="Microbacteriaceae"|Family =="Flavobacteriaceae" | Family=="Sphingobacteriaceae")

phyloseq_T20<-phyloseq_trt_clean %>%
  subset_samples(Code =="T20")%>%
  subset_taxa(Family =="Listeriaceae"| Family =="Pseudomonadaceae"| Family =="Xanthomonadaceae"| Family =="Microbacteriaceae"|Family =="Flavobacteriaceae" | Family=="Sphingobacteriaceae")

#Merge Phyloseqs together
physeq_merge<-merge_phyloseq(phyloseq_T1, phyloseq_T2,phyloseq_T3,phyloseq_T4,phyloseq_T5,
                             phyloseq_T6,phyloseq_T7,phyloseq_T8,phyloseq_T9,phyloseq_T10,
                             phyloseq_T11,phyloseq_T12,phyloseq_T13,phyloseq_T14,phyloseq_T15,
                             phyloseq_T16,phyloseq_T17,phyloseq_T18,phyloseq_T19,phyloseq_T20)



#Get ASV table from phyloseq object
asv_trt_clean<-t(as.data.frame(otu_table(physeq_merge)))
tail(rowSums(asv_trt_clean))

#Compositional analysis of microbiome data -based on Microbiome Analysis in R. Chap 10. #
#Step 1: Convert ASV table to appropriate format. Following step requires samples on rows and ASVs in columns
head(t(asv_trt_clean))

#Step 2: Replace zero values before clr transformation. Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
asv.n0_trt<-t(cmultRepl(t(asv_trt_clean), label=0, method="CZM", output="p-counts")) #

head(asv.n0_trt) #output table needs to have samples in columns and ASVs in rows

#Note: Check the output to make sure there are no negative numbers. If samples or ASV are sparse, the CZM method will 
#add a negative number that interferes with the log normalization. If necessary use function below to convert negative values
#into positives
asv_n0_trt<-ifelse(asv.n0_trt < 0, asv.n0_trt*(-1), asv.n0_trt)

#Step 3: Convert data to proportions
asv.n0_trt_prop<-apply(asv_n0_trt, 2, function(x) {x/sum(x)})

#Step 4: Perform abundance and sample filtering and deal sparsity
#Filter the data to remove all taxa with less than 0.00001% abundance in any sample
asv.n0_trt_prop_f<-asv.n0_trt[apply(asv.n0_trt_prop, 1, min) > 0.00000000000000000001, ]
head(asv.n0_trt_prop_f) #Check that samples are on columns and asvs in rows


#Step 5: perform CLR transformation
asv.n0.clr_trt<-t(apply(asv.n0_trt_prop, 2, function(x){log(x)-mean(log(x))}))
head(asv.n0.clr_trt) #Check output table. Samples should be in rows and asvs in columns

#Step 6: Perform Singular Value Decomposition (PCA)
pc.clr_trt<-prcomp(asv.n0.clr_trt)

png("Screeplot - PCA.png")
par(mar=c(2,2,2,2))
screeplot(pc.clr_trt, type='lines', main="PCA")
dev.off()

#Calculate total variance of the data
mvar.clr_trt<-mvar(asv.n0.clr_trt)

#Display results
row_trt<-rownames(asv.n0.clr_trt) #Make vector with sample names
pc_out_trt<-as.data.frame(pc.clr_trt$x[,1:2]) #Get PC1 and PC2
pc_out_meta_trt<-as.data.frame(bind_cols(pc_out_trt,metadata_clean)) #Add metadata information
row.names(pc_out_meta_trt)<-row_trt #Add rownames to dataframe
pc_out_meta_trt$Rep<-as.factor(pc_out_meta_trt$Rep)

# Make PCA plot
PCA <- ggplot(pc_out_meta_trt, aes(x=PC1,y=PC2, color=Trt, shape=Type))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=13,color='black'), legend.title= element_text(size=13)) +
  theme(legend.position = 'bottom')+
  theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=13,color='black')) +
  # theme(panel.background = element_rect(fill='transparent', color = NA),
  #       plot.background = element_rect(fill = 'transparent',color = NA),
  #       panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_trt$sdev[1]^2/mvar.clr_trt*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_trt$sdev[2]^2/mvar.clr_trt*100, digits=1), "%", sep="")) +
  ggtitle("PCA - Biofilms and Cocktails")+theme(plot.title = element_text(hjust = 0.5, face = 'bold'))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA
ggsave("PCA_Biofilms.png", plot=PCA, device="png", width=8, height=8, units="in", dpi=600)
ggsave("PCA_Biofilms.svg", plot=PCA, device="svg", width=8, height=8, units="in", dpi=600)


#Transform sample counts into compositions
asv.n0.acomp_trt<-as.data.frame(acomp(t(asv_trt_clean)), total=1)
rowSums(asv.n0.acomp_trt) #Verify rows sum to 1

#Make Phyloseq object
phyloseq_RA_trt <- phyloseq(otu_table(asv.n0.acomp_trt, taxa_are_rows = FALSE), tax_table(taxon_all), sample_data(metadata_clean))

#Make long format table from Phyloseq object
asv_long_trt <- phyloseq_RA_trt %>%  
  transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance as percentage
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

# Calculate mean relative abundance by Treatment for each ASV
asv_mean_trt<-asv_long_trt %>%
  group_by(OTU, Trt, PMA, Type, Family, SampleOrder)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))


#Stacked barplot by treatment at the ASV level
barplot_biofilms_trt<-ggplot(asv_mean_trt, aes(x=reorder(Trt,SampleOrder), y=Mean, fill=Family))+
  geom_bar(stat='identity', color='black')+facet_grid(Type~PMA, scales = "free", space = 'free')+
  geom_text(aes(label=OTU), position = position_stack(vjust=0.5), size=2)+
  ylim(0,100)+
  theme(axis.title = element_text(color='black'), axis.text.x=element_text(color='black', size=10, angle=90), 
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) + 
  ylab("Relative Abundance (%)")+xlab("Treatment")+
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black',face="italic"), legend.title= element_text(size=10)) +
  ggtitle("Biofilms composition with strain ASV")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Barplots_biofilm by strain.png", plot=barplot_biofilms_trt, device="png", width=13, height=10, units="in", dpi=600)
ggsave("Barplots_biofilm by strain.svg", plot=barplot_biofilms_trt, device="svg", width=13, height=10, units="in", dpi=600)

#Subset data by type and print each barplot separately
asv_mean_trt_PMA<-subset(asv_mean_trt, PMA=="Y" & Mean >0)
asv_mean_trt_noPMA<-subset(asv_mean_trt, PMA=="N" & Type=="Biofilm" & Mean >0)
asv_mean_trt_cocktail<-subset(asv_mean_trt,  Type=="Cocktail"& Mean >0)

#Make barplots
barplot_biofilms_pma<-ggplot(asv_mean_trt_PMA, aes(x=reorder(Trt,SampleOrder), y=Mean, fill=Family))+
  geom_bar(stat='identity', color='black')+
  geom_text(aes(label=OTU), position = position_stack(vjust=0.5), size=2)+
  #ylim(0,100)+
  theme(axis.title = element_text(color='black'), axis.text.x=element_text(color='black', size=10, angle=90), 
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) + 
  ylab("Relative Abundance (%)")+xlab("Treatment")+
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black',face="italic"), legend.title= element_text(size=10)) +
  ggtitle("Biofilms composition with strain ASV - PMA")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Barplots_biofilm by strain PMA.png", plot=barplot_biofilms_pma, device="png", width=13, height=10, units="in", dpi=600)
ggsave("Barplots_biofilm by strain PMA.svg", plot=barplot_biofilms_pma, device="svg", width=13, height=10, units="in", dpi=600)

barplot_biofilms_nopma<-ggplot(asv_mean_trt_noPMA, aes(x=reorder(Trt,SampleOrder), y=Mean, fill=Family))+
  geom_bar(stat='identity', color='black')+
  geom_text(aes(label=OTU), position = position_stack(vjust=0.5), size=2)+
  ylim(0,100)+
  theme(axis.title = element_text(color='black'), axis.text.x=element_text(color='black', size=10, angle=90), 
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) + 
  ylab("Relative Abundance (%)")+xlab("Treatment")+
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black',face="italic"), legend.title= element_text(size=10)) +
  ggtitle("Biofilms composition with strain ASV - no PMA")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Barplots_biofilm by strain noPMA.png", plot=barplot_biofilms_nopma, device="png", width=13, height=10, units="in", dpi=600)
ggsave("Barplots_biofilm by strain noPMA.svg", plot=barplot_biofilms_nopma, device="svg", width=13, height=10, units="in", dpi=600)


barplot_biofilms_cocktail<-ggplot(asv_mean_trt_cocktail, aes(x=reorder(Trt,SampleOrder), y=Mean, fill=Family))+
  geom_bar(stat='identity', color='black')+
  geom_text(aes(label=OTU), position = position_stack(vjust=0.5), size=2)+
  ylim(0,100)+
  theme(axis.title = element_text(color='black'), axis.text.x=element_text(color='black', size=10, angle=90), 
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) + 
  ylab("Relative Abundance (%)")+xlab("Treatment")+
  theme(panel.background = element_rect(fill="white", color =NA),
        plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black',face="italic"), legend.title= element_text(size=10)) +
  ggtitle("Biofilms composition with strain ASV - Cocktail")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Barplots_ by strain Cocktail.png", plot=barplot_biofilms_cocktail, device="png", width=13, height=10, units="in", dpi=600)
ggsave("Barplots_ by strain Cocktail.svg", plot=barplot_biofilms_cocktail, device="svg", width=13, height=10, units="in", dpi=600)
