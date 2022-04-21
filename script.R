##1 Sequences processing #####
library(dada2) #packageVersion("dada2") ‘1.16.0’
library(phyloseq) #packageVersion("phyloseq") ‘1.32.0’

rm(list=ls())
samples <- scan("samples", what="character")

forward_reads <- paste0(samples, "_S1_L001_R1_001.fastq")
reverse_reads <- paste0(samples, "_S1_L001_R2_001.fastq")

filtered_forward_reads <- paste0(samples, "_R1_filtered.fastq.gz")
filtered_reverse_reads <- paste0(samples, "_R2_filtered.fastq.gz")

#_1.1 Quality filtering #####
plotQualityProfile(forward_reads[1:5])
plotQualityProfile(reverse_reads[1:5])

filtered_out<-filterAndTrim(forward_reads, filtered_forward_reads,
                            reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
                            rm.phix=TRUE, minLen=c(250,210), truncLen=c(250,210))

plotQualityProfile(filtered_forward_reads[1:5])
plotQualityProfile(filtered_reverse_reads[1:5])

err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE)

plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)

#_1.2 Dereplication #####
derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples

#_1.3 Infering ASVs #####
dada_forward <- dada(derep_forward, err=err_forward_reads, pool=TRUE, multithread=TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool=TRUE, multithread=TRUE)

#_1.4 Merging paired-end reads #####
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse, returnRejects=F,
                               propagateCol="sequence", derep_reverse, maxMismatch=0, minOverlap=12)

#_1.5 Distribution seqs length #####
seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab); dim(seqtab)

table(nchar(getSequences(seqtab)))

seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 398:448]
table(nchar(getSequences(seqtab)))

#_1.6 Chimera identification  #####
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T, multithread=T, minFoldParentOverAbundance=1)
sum(seqtab)
sum(seqtab.nochim)/sum(seqtab)

#_1.7 Cleaning table #####
getN <- function(x) sum(getUniques(x))
cleaning_table=data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          seqtab=rowSums(seqtab), nonchim=rowSums(seqtab.nochim),
                          perc_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))
cleaning_table
write.table(cleaning_table, "Filtering_procedure_table.csv", sep="\t", quote=F, col.names=NA)

length(as.data.frame(seqtab.nochim)) #number of asv identified
sum(colSums(as.data.frame(seqtab.nochim))) #number of total seqs
#save.image(file="output.poisson_pahs.dada2.R")

#_1.8 Assigning taxonomy #####
taxo <- assignTaxonomy(seqtab.nochim, "D:/silva_nr_v138_train_set.fa.gz", multithread=TRUE, tryRC=F)

#_1.9 Phyloseq object  #####
metadata <- read.table("02.samplesID_metabarcoding.csv", header=T, row.names=1, check.names=F, sep=";")

physeq=phyloseq(otu_table(seqtab.nochim, taxa_are_rows=F),
                tax_table(taxo),
                sample_data(metadata))

#_1.10 Getting ASV sequences #####
dna <- Biostrings::DNAStringSet(taxa_names(physeq))
names(dna) <- taxa_names(physeq)
physeq <- merge_phyloseq(physeq, dna)
taxa_names(physeq) <- paste0("ASV", seq(ntaxa(physeq)))
refseq(physeq)

sample_sums(physeq)

#_1.11 Calculating relative abundance #####
physeq.rel <- transform_sample_counts(physeq, function(OTU) (OTU/sum(OTU)*100))
sample_sums(physeq.rel)

#_1.12 Saving phyloseq object  #####
save(physeq, physeq.rel, file = "chub_pahs_parasite_phyloseq.RData")

##2 alpha diversity #####
library(ggplot2); packageVersion("ggplot2") #‘3.3.5’
library(phyloseq); packageVersion("phyloseq") #‘1.32.0’
library(ggpubr); packageVersion("ggpubr") #‘0.4.0’

rm(list=ls())
load("chub_pahs_parasite_phyloseq.RData")
physeq@sam_data$PAH <- gsub("0.1X", "CTRL", physeq@sam_data$PAH, fixed=TRUE)
physeq@sam_data$PAH <- gsub("10X", "PAHs", physeq@sam_data$PAH, fixed=TRUE)
physeq@sam_data$Parasite_PAH <- paste0(physeq@sam_data$PAH,"_", physeq@sam_data$Parasite)

#_2.1 rarefaction curves #####
vegan::rarecurve((otu_table(physeq)), step=500, cex=0.5,
                 xlab="number of reads", ylab="number of ASVs")

#_2.2 alpha diversity indexes #####
physeq = rarefy_even_depth(physeq, sample.size = min(sample_sums(physeq)), replace = F)
richness<-estimate_richness(physeq, split = TRUE, measures = c("Observed", "Shannon", "Simpson"))

richness=richness[match(rownames(sample_data(physeq)), rownames(richness)),]
sample_data(physeq)=cbind(sample_data(physeq), richness)
sample_data(physeq)

richness$PAH <- factor(richness$PAH, levels = c("CTRL", "PAHs"))
richness$Parasite <- factor(richness$Parasite, levels = c("Uninfected", "Infected"))
richness$Parasite_PAH <- factor(richness$Parasite_PAH,
                                levels = c("CTRL_Uninfected", "CTRL_Infected", "PAHs_Uninfected", "PAHs_Infected"))

library(phyloseq)
library(dplyr)
library(stringr)
library(tidyverse)

sample_data(physeq) %>% 
  as("matrix") %>%
  as_tibble(rownames = "OTU") %>%
  mutate(Observed = as.numeric(Observed),
         Shannon = as.numeric(Shannon),
         Simpson = as.numeric(Simpson))%>%
  select(-Parasite, -PAH) %>%
  group_by(Parasite_PAH)%>%
  summarize(mean_obs=mean(Observed), sd_obs=sd(Observed),
            mean_H=mean(Shannon), sd_H=sd(Shannon),
            mean_S=mean(Simpson),sd_S=sd(Simpson))
  
sample_data(physeq) %>% 
  as("matrix") %>%
  as_tibble(rownames = "OTU") %>%
  mutate(Observed = as.numeric(Observed),
         Shannon = as.numeric(Shannon),
         Simpson = as.numeric(Simpson))%>%
  select(-Parasite_PAH, -PAH) %>%
  group_by(Parasite)%>%
  summarize(mean_obs=mean(Observed), sd_obs=sd(Observed),
            mean_H=mean(Shannon), sd_H=sd(Shannon),
            mean_S=mean(Simpson),sd_S=sd(Simpson))

sample_data(physeq) %>% 
  as("matrix") %>%
  as_tibble(rownames = "OTU") %>%
  mutate(Observed = as.numeric(Observed),
         Shannon = as.numeric(Shannon),
         Simpson = as.numeric(Simpson))%>%
  select(-Parasite_PAH, -Parasite) %>%
  group_by(PAH)%>%
  summarize(mean_obs=mean(Observed), sd_obs=sd(Observed),
            mean_H=mean(Shannon), sd_H=sd(Shannon),
            mean_S=mean(Simpson),sd_S=sd(Simpson))

#_2.3 plotting indexes (FIGURE 1) #####
my_comparisons <- list( c("CTRL_Uninfected", "PAHs_Uninfected"),
                        c("CTRL_Infected", "PAHs_Infected"),
                        c("CTRL_Uninfected", "CTRL_Infected"),
                        c("PAHs_Uninfected", "PAHs_Infected"),
                        c("PAHs", "CTRL"),
                        c("Uninfected", "Infected"))
richness=sample_data(physeq)
p0<-ggplot(richness, aes (x = PAH, y = Observed),
           color="black", fill= "grey90", alpha=0.9)+
  geom_violin(trim=T, fill="gray90")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  stat_summary(fun=mean, geom="point", shape=18, size=4,
               color="cadetblue", alpha=0.8)+
  xlab("")+ ylab("Richness (Sobs)")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
  stat_compare_means(label.y = 1.7, size=3)+
  labs(title = "")+
  theme_bw()+
  theme(legend.position = "none", axis.title.y = element_text(size=9),
        plot.title = element_text(hjust = 0.5, size=10))

p1<-ggplot(richness, aes(x = Parasite, y = Observed),
           color="black", fill= "grey90", alpha=0.9)+
  geom_violin(trim=T, fill="gray90")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  stat_summary(fun=mean, geom="point", shape=18, size=4,
               color="cadetblue", alpha=0.8)+
  xlab("")+ ylab("Richness (Sobs)")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
  stat_compare_means(label.y = 1.7, size=3)+
  labs(title = "")+
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_text(size=9),
        plot.title = element_text(hjust = 0.5, size=10))

p2<-ggplot(richness, aes(x =PAH, y = Shannon),
           color="black", fill= "grey90", alpha=0.9)+
  geom_violin(trim=T, fill="gray90")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  stat_summary(fun=mean, geom="point", shape=18, size=4,
               color="cadetblue", alpha=0.8)+
  xlab("")+ ylab("Shannon (H')")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
  stat_compare_means(label.y = 1, size=3)+
  labs(title = "")+
  theme_bw()+
  theme(legend.position = "none", axis.title.y = element_text(size=9),
        plot.title = element_text(hjust = 0.5, size=10))


p3<-ggplot(richness, aes(x = Parasite, y = Shannon),
           color="black", fill= "grey90", alpha=0.9)+
  geom_violin(trim=T, fill="gray90")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  stat_summary(fun=mean, geom="point", shape=18, size=4,
               color="cadetblue", alpha=0.8)+
  xlab("")+ ylab("Shannon (H')")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
  stat_compare_means(label.y = 1, size=3)+
  labs(title = "")+
  theme_bw() +
  theme(legend.position = "none", axis.title.y = element_text(size=9),
        plot.title = element_text(hjust = 0.5, size=10))

p4<-ggplot(richness, aes(x=PAH, y=Observed),
           color="black", fill= "grey90", alpha=0.9)+
  geom_violin(trim=T, fill="gray90")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  stat_summary(fun=mean, geom="point", shape=18, size=4,
               color="cadetblue", alpha=0.8)+
  xlab("")+ ylab("Richness (Sobs)")+
  facet_wrap(facets="Parasite")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
  stat_compare_means(label.y = 1, size=3)+
  theme_bw()+
  theme(legend.position = "none", strip.background = element_blank(),
        axis.title.y = element_text(size=9))

p5<-ggplot(richness, aes(x = PAH, y = Shannon),
           color="black", fill= "grey90", alpha=0.9)+
  geom_violin(trim=T, fill="gray90")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1)+
  stat_summary(fun=mean, geom="point", shape=18, size=4,
               color="cadetblue", alpha=0.8)+
  xlab("")+ylab("Shannon (H')")+
  facet_wrap(facets="Parasite")+
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")+
  stat_compare_means(label.y = 1, size=3)+
  theme_bw() +
  theme(legend.position="right", legend.title=element_blank(),
        axis.title.y = element_text(size=9),
        strip.background = element_blank())

ggpubr::ggarrange(p1, p0, p4, p3, p2, p5, nrow=2, ncol=3, widths = c(1,1,1.7),
                  labels=c("A", "B", "C", "D", "E", "F"))
sample_sums(physeq)

##3 beta diversity #####
library(ggplot2); packageVersion("ggplot2") #‘3.3.5’
library(vegan); packageVersion("vegan") #‘2.5.7’

rm(list=ls())
load("chub_pahs_parasite_phyloseq.RData")
physeq@sam_data$PAH <- gsub("0.1X", "CTRL", physeq@sam_data$PAH, fixed=TRUE)
physeq@sam_data$PAH <- gsub("10X", "PAHs", physeq@sam_data$PAH, fixed=TRUE)
physeq@sam_data

#_3.1 PCoA plot (FIGURE 2) #####
asv<-as.data.frame(otu_table(physeq))
dist<-dist(decostand(asv, method="hellinger"))

#PCoA analysis
pcoa <- cmdscale (dist, eig = TRUE)
pcoa$eig[1]/sum(pcoa$eig)*100
pcoa$eig[2]/sum(pcoa$eig)*100

df<-data.frame(pcoa$points[,1:2])
metadata <- as(sample_data(physeq), "data.frame")
df<-cbind(df, metadata)

ggplot(df, aes(x=X1 , y=X2))+
  geom_point(aes(shape=PAH, color=Parasite), size=5, alpha=0.8, show.legend=T)+
  scale_color_manual(values=c("#2166AC", "#B2182B"))+
  labs(x = "PCoA1 (30.7%)", y = "PCoA2 (12.3%)") +
  theme_bw()

#_3.2 PERMANOVA tests #####
ps.adonis.parasite<-adonis(dist ~ Parasite, data=metadata, perm=9999)
ps.adonis.parasite

ps.adonis.pah<-adonis(dist ~ PAH, data=metadata, perm=9999)
ps.adonis.pah

#_3.3 Variance partitioning #####
var.parasite<-subset(metadata, select=c(Parasite))
var.PAH<-subset(metadata, select=c(PAH))

asv.part<-varpart(decostand(asv, method="hellinger"),
                  var.parasite, var.PAH)

## 4 Phylum_composition barplot (FIGURE S2 & TABLE S4) ####
library(phyloseq); packageVersion("phyloseq") #‘1.32.0’
library(tidyverse); packageVersion("tidyverse") #‘1.3.1’
library(tidyr); packageVersion("tidyr") #‘1.1.3’
library(dplyr); packageVersion("dplyr") #‘1.0.7’
library(stringr); packageVersion("stringr") #‘1.4.0’
library(broom); packageVersion("broom") #‘0.7.6’
library(ggtext); packageVersion("ggtext") #‘0.1.1’

rm(list=ls())
load("chub_pahs_parasite_phyloseq.RData")

## asv table
asv_table<-otu_table(physeq) %>%  # Note, need taxa_as_rows = TRUE
  as("matrix") %>%
  as_tibble(rownames="OTU") %>%
  rename(samples = OTU) %>%
  pivot_longer(-samples, names_to="asv_id", values_to="count")

## tax table
taxo_table<-tax_table(physeq) %>%
  as("matrix") %>%
  as_tibble(rownames = "OTU") %>%
  rename_all(tolower) %>%
  rename(asv_id= otu) %>%
  mutate(phylum= ifelse(phylum == "Proteobacteria", class, phylum)) %>%
  mutate(phylum=str_replace(phylum, "(.*)", "*\\1*")) %>%
  select(asv_id, phylum)

## metadata table
metadata<-sample_data(physeq) %>%
  as("matrix") %>%
  as_tibble(rownames = "samples") %>%
  mutate(PAH=str_replace(PAH, "0.1X", "CTRL"),
         PAH=str_replace(PAH, "10X", "PAHs")) %>%
  mutate(Parasite_PAH=paste0(PAH,"_",Parasite)) %>%
  select(samples, PAH, Parasite, Parasite_PAH)

## composite
composite<-inner_join(asv_table, taxo_table, by="asv_id") %>%
  group_by(samples, phylum) %>%
  summarize(count = sum(count), .groups="drop") %>%
  group_by(samples) %>%
  mutate(rel_abund=count/sum(count)*100) %>%
  #summarize(total=sum(rel_abund))
  ungroup() %>%
  group_by(phylum) %>%
  filter(mean(rel_abund) > 0.5) %>%
  ungroup() %>%
  select(-count) %>%
  inner_join(., metadata, by="samples") %>%
  group_by(Parasite, PAH, phylum) %>%
  summarize(mean_rel_abund=mean(rel_abund), .groups="drop") %>%
  spread(phylum, mean_rel_abund) %>%
  mutate(Others = 100-rowSums(.[-c(1:2)])) %>%
  pivot_longer(!c(Parasite, PAH), names_to = "phylum", values_to = "mean_rel_abund")

##_4.1 barplot ####
library(RColorBrewer)
colourCount<-length(unique(composite$phylum))
getPalette<-colorRampPalette(brewer.pal(12, "Paired"))

composite$Parasite<-factor(composite$Parasite,
                           levels = c("Uninfected", "Infected"))

ggplot(composite, aes (x=PAH, y=mean_rel_abund, fill = phylum)) +
  geom_bar(position = "fill",stat = "identity")+ theme_classic()+
  theme(axis.text.x = element_text(angle=0), legend.position = "right")+
  scale_fill_manual(values = getPalette(colourCount))+
  xlab("")+ylab("Relative Abundance")+
  facet_wrap(facets="Parasite")+
  theme(legend.text= element_markdown(), legend.title=element_blank())

##_4.2 Table abundance main phyla####
table<-inner_join(asv_table, taxo_table, by="asv_id") %>%
  group_by(samples, phylum) %>%
  summarize(count = sum(count), .groups="drop") %>%
  group_by(samples) %>%
  mutate(rel_abund=count/sum(count)*100) %>%
  ungroup() %>%
  group_by(phylum) %>%
  filter(mean(rel_abund) > 0.5) %>%
  ungroup() %>%
  select(-count) %>%
  inner_join(., metadata, by="samples") %>%
  group_by(Parasite_PAH, phylum) %>%
  summarize(mean_rel_abund=round(mean(rel_abund),2),
            sd_rel_abund=round(sd(rel_abund),2), .groups="drop") %>%
  unite(mean_sd, mean_rel_abund, sd_rel_abund, sep=" ± ") %>%
  spread(phylum, mean_sd)%>%
  t

#5 plot phyla abundances & wilcox.test (FIGURE 3 and TABLE S5) #####
library(phyloseq); packageVersion("phyloseq") #‘1.32.0’
library(tidyverse); packageVersion("tidyverse") #‘1.3.1’
library(tidyr); packageVersion("tidyr") #‘1.1.3’
library(dplyr); packageVersion("dplyr") #‘1.0.7’
library(broom); packageVersion("broom") #‘0.7.6’
library(ggtext); packageVersion("ggtext") #‘0.1.1’

rm(list=ls())
load("chub_pahs_parasite_phyloseq.RData")

## asv table
asv_table<-otu_table(physeq) %>%  # Note, need taxa_as_rows = TRUE
  as("matrix") %>%
  as_tibble(rownames="OTU") %>%
  rename(samples = OTU) %>%
  pivot_longer(-samples, names_to="asv_id", values_to="count")

## tax table
taxo_table<-tax_table(physeq) %>%
  as("matrix") %>%
  as_tibble(rownames = "OTU") %>%
  rename_all(tolower) %>%
  rename(asv_id= otu) %>%
  mutate(phylum= ifelse(phylum == "Proteobacteria", class, phylum)) %>%        
  select(asv_id, phylum)

## metadata table
metadata<-sample_data(physeq) %>%
  as("matrix") %>%
  as_tibble(rownames = "samples") %>%
  mutate(PAH=str_replace(PAH, "0.1X", "CTRL"),
         PAH=str_replace(PAH, "10X", "PAHs")) %>%
  mutate(Parasite_PAH=paste0(PAH,"_",Parasite)) %>%
  select(samples, PAH, Parasite, Parasite_PAH)

## composite
composite<-inner_join(asv_table, taxo_table, by="asv_id") %>%
  group_by(samples, phylum) %>%
  summarize(count = sum(count), .groups="drop") %>%
  group_by(samples) %>%
  mutate(rel_abund=count/sum(count)*100) %>%
  #summarize(total=sum(rel_abund))
  ungroup() %>%
  group_by(phylum) %>%
  filter(mean(rel_abund) > 0.5) %>%
  ungroup() %>%
  select(-count) %>%
  inner_join(., metadata, by="samples")

vect_main_phylum<-unique(composite$phylum)
#save(vect_main_phylum, file="vect_main_phylum")

#_5.1 wilcox.test infected (n=10) vs uninfected (n=10) #####
wilcox_test_parasite_phylum<-composite %>%
  nest(data=-phylum) %>%
  mutate(test.res=map(.x=data, ~wilcox.test(rel_abund ~ Parasite,
                                            data=.x) %>% tidy)) %>%
  unnest(test.res) %>%
  mutate(p.adj=p.adjust(p.value, method="BH")) %>%
  select(phylum, p.value, p.adj)

#_5.2 wilcox.test ctrl (n=10) vs PAHs (n=10) #####
wilcox_test_PAH_phylum=composite %>%
  nest(data=-phylum) %>%
  mutate(test.res=map(.x=data, ~wilcox.test(rel_abund ~ PAH,
                                            data=.x) %>% tidy)) %>%
  unnest(test.res) %>%
  mutate(p.adj=p.adjust(p.value, method="BH")) %>%
  select(phylum, p.value, p.adj)

#_5.3 plots #####
p1=composite %>%
  mutate(rel_abund=rel_abund+1/20000,
         phylum=str_replace(phylum, "(.*)", "*\\1*")) %>%
  ggplot(aes(rel_abund, y=phylum, fill=Parasite))+
  geom_boxplot(outlier.size = 0.5, alpha=0.4, lwd=0.2)+
  labs(x="Log10 (Relative Abundances)", y=NULL)+
  scale_fill_manual(values=c("#2166AC", "black"))+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text.y = element_markdown(),
        axis.title.x = element_blank(),
        legend.title=element_blank(),
        legend.position="top", legend.margin=margin(10,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

p2<-composite %>%
  mutate(rel_abund=rel_abund+1/20000,
         phylum=str_replace(phylum, "(.*)", "*\\1*")) %>%
  ggplot(aes(rel_abund, y=phylum, fill=PAH))+
  geom_boxplot(outlier.size = 0.5, alpha=0.4, lwd=0.2)+
  labs(x="Log10 (Relative Abundance)", y=NULL)+
  scale_fill_manual(values=c("#2166AC", "black"))+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text.y = element_markdown(),
        axis.title.x = element_text(size=9),
        legend.title=element_blank(),
        legend.position="top", legend.margin=margin(10,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

ggpubr::ggarrange(p1, p2, nrow=2, heights=c(1,1),
                  labels=c("A", "B"))

#6 plot genus abundances & wilcox.test (FIGURE 4 and TABLE S6) #####
library(phyloseq); packageVersion("phyloseq") #‘1.32.0’
library(tidyverse); packageVersion("tidyverse") #‘1.3.1’
library(tidyr); packageVersion("tidyr") #‘1.1.3’
library(dplyr); packageVersion("dplyr") #‘1.0.7’
library(broom); packageVersion("broom") #‘0.7.6’
library(ggtext); packageVersion("ggtext") #‘0.1.1’

rm(list=ls())
load("chub_pahs_parasite_phyloseq.RData")

## asv table
asv_table<-otu_table(physeq) %>%  # Note, need taxa_as_rows = TRUE
  as("matrix") %>%
  as_tibble(rownames="OTU") %>%
  rename(samples = OTU) %>%
  pivot_longer(-samples, names_to="asv_id", values_to="count")

## tax table
taxo_table<-tax_table(physeq) %>%
  as("matrix") %>%
  as_tibble(rownames = "OTU") %>%
  rename_all(tolower) %>%
  rename(asv_id= otu) %>%
  mutate(phylum= ifelse(phylum == "Proteobacteria", class, phylum)) %>%        
  select(asv_id, phylum, genus) %>%
  mutate(genus=replace_na(genus, "unclass.")) %>%
  mutate(genus=paste0("(",phylum,") ", genus)) %>%
  filter(across(genus, ~ !grepl('unclass.', .))) %>%
  select(asv_id, genus)

## metadata table
metadata<-sample_data(physeq) %>%
  as("matrix") %>%
  as_tibble(rownames = "samples") %>%
  mutate(PAH=str_replace(PAH, "0.1X", "CTRL"),
         PAH=str_replace(PAH, "10X", "PAH")) %>%
  mutate(Parasite_PAH=paste0(PAH,"_",Parasite)) %>%
  select(samples, PAH, Parasite, Parasite_PAH)

## composite
composite<-inner_join(asv_table, taxo_table, by="asv_id") %>%
  group_by(samples, genus) %>%
  summarize(count = sum(count), .groups="drop") %>%
  group_by(samples) %>%
  mutate(rel_abund=count/sum(count)*100) %>%
  #summarize(total=sum(rel_abund))
  ungroup() %>%
  group_by(genus) %>%
  filter(mean(rel_abund) > 1) %>%
  ungroup() %>%
  select(-count) %>%
  inner_join(., metadata, by="samples")

vect_main_genus<-unique(composite$genus)
#save(vect_main_genus, file="vect_main_genus")

## wilcox.test infected (n=10) vs uninfected (n=10)
wilcox_test_parasite_genus<-composite %>%
  nest(data=-genus) %>%
  mutate(test.res=map(.x=data, ~wilcox.test(rel_abund ~ Parasite,
                                            data=.x) %>% tidy)) %>%
  unnest(test.res) %>%
  mutate(p.adj=p.adjust(p.value, method="BH"))  %>%
  #filter(p.adj<0.05) %>%
  select(genus, p.value, p.adj)

## wilcox.test ctrl (n=10) vs pah (n=10)
wilcox_test_pah_genus<-composite %>%
  nest(data=-genus) %>%
  mutate(test.res=map(.x=data, ~wilcox.test(rel_abund ~ PAH,
                                            data=.x) %>% tidy)) %>%
  unnest(test.res) %>%
  mutate(p.adj=p.adjust(p.value, method="BH")) %>%
  #filter(p.adj<0.05) %>%
  select(genus, p.value, p.adj)

## plot genus
p3<-composite %>%
  mutate(rel_abund=rel_abund+1/20000,
         genus=str_replace(genus, "(.*)", "*\\1*")) %>%
  ggplot(aes(rel_abund, y=genus, fill=Parasite))+
  geom_boxplot(outlier.size = 0.5, alpha=0.4, lwd=0.2)+
  labs(x="Log10 (Relative Abundances)", y=NULL)+
  scale_fill_manual(values=c("#2166AC", "black"))+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text.y = element_markdown(),
        axis.title.x = element_blank(),
        legend.title=element_blank(),
        legend.position="top", legend.margin=margin(10,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))+
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

p4<-composite %>%
  mutate(rel_abund=rel_abund+1/20000,
         genus=str_replace(genus, "(.*)", "*\\1*")) %>%
  ggplot(aes(rel_abund, y=genus, fill=PAH))+
  geom_boxplot(outlier.size = 0.5, alpha=0.4, lwd=0.2)+
  labs(x="Log10 (Relative Abundance)", y=NULL)+
  scale_fill_manual(values=c("#2166AC", "black"))+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text.y = element_markdown(),
        axis.title.x = element_text(size=9),
        legend.title=element_blank(),
        legend.position="top", legend.margin=margin(10,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))+
  guides(color = guide_legend(nrow = 1, byrow = TRUE))

ggpubr::ggarrange(p3, p4, nrow=2, heights=c(1,1),
                  labels=c("A", "B"))

#7 plot phyla abundances and associated wilcox.test independently for infected and uninfected chub (FIGURE S3 and TABLE S7) #####
library(phyloseq); packageVersion("phyloseq") #‘1.32.0’
library(tidyverse); packageVersion("tidyverse") #‘1.3.1’
library(tidyr); packageVersion("tidyr") #‘1.1.3’
library(dplyr); packageVersion("dplyr") #‘1.0.7’
library(broom); packageVersion("broom") #‘0.7.6’
library(ggtext); packageVersion("ggtext") #‘0.1.1’

rm(list=ls())
load("chub_pahs_parasite_phyloseq.RData")
load("vect_main_phylum")

#_7.1 uninfected chub  (n=10 samples) ####
physeq.uninfected<-subset_samples(physeq, Parasite=="Uninfected")
physeq.uninfected<-prune_taxa(taxa_sums(physeq.uninfected)>0, physeq.uninfected)

## asv table
asv_table<-otu_table(physeq.uninfected) %>%
  as("matrix") %>%
  as_tibble(rownames="OTU") %>%
  rename(samples = OTU) %>%
  pivot_longer(-samples, names_to="asv_id", values_to="count")

## tax table
taxo_table<-tax_table(physeq.uninfected) %>%
  as("matrix") %>%
  as_tibble(rownames = "OTU") %>%
  rename_all(tolower) %>%
  rename(asv_id= otu) %>%
  mutate(phylum= ifelse(phylum == "Proteobacteria", class, phylum)) %>%        
  select(asv_id, phylum)

## metadata table
metadata<-sample_data(physeq.uninfected) %>%
  as("matrix") %>%
  as_tibble(rownames = "samples") %>%
  mutate(PAH=str_replace(PAH, "0.1X", "CTRL"),
         PAH=str_replace(PAH, "10X", "PAHs")) %>%
  mutate(Parasite_PAH=paste0(PAH,"_",Parasite)) %>%
  select(samples, PAH, Parasite, Parasite_PAH)

## composite
composite<-inner_join(asv_table, taxo_table, by="asv_id") %>%
  group_by(samples, phylum) %>%
  summarize(count = sum(count), .groups="drop") %>%
  group_by(samples) %>%
  mutate(rel_abund=count/sum(count)*100) %>%
  ungroup() %>%
  group_by(phylum) %>%
  filter(phylum %in% vect_main_phylum) %>%
  ungroup() %>%
  select(-count) %>%
  inner_join(., metadata, by="samples")

## wilcox.test ctrl (n=5) vs pah (n=5)
wilcox_test_pah_phylum<-composite %>%
  nest(data=-phylum) %>%
  mutate(test.res=map(.x=data, ~wilcox.test(rel_abund ~ PAH,
                                            data=.x) %>% tidy)) %>%
  unnest(test.res) %>%
  mutate(p.adj=p.adjust(p.value, method="BH")) %>%
  select(phylum, p.value, p.adj)

## plot phylum
p1<-composite %>%
  mutate(rel_abund=rel_abund+1/20000,
         phylum=str_replace(phylum, "(.*)", "*\\1*")) %>%
  ggplot(aes(rel_abund, y=phylum, fill=PAH))+
  geom_boxplot(outlier.size = 0.5, alpha=0.4, lwd=0.2)+
  labs(x="Log10 (Relative Abundances)", y=NULL)+
  scale_fill_manual(values=c("#2166AC", "black"))+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text.y = element_markdown(),
        axis.title.x = element_blank(),
        legend.title=element_blank(),
        legend.position="top", legend.margin=margin(10,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

#_7.2 infected chub (n=10 samples) ####
physeq.infected<-subset_samples(physeq, Parasite=="Infected")
physeq.infected<-prune_taxa(taxa_sums(physeq.infected)>0, physeq.infected)

## asv table
asv_table<-otu_table(physeq.infected) %>%
  as("matrix") %>%
  as_tibble(rownames="OTU") %>%
  rename(samples = OTU) %>%
  pivot_longer(-samples, names_to="asv_id", values_to="count")

## tax table
taxo_table<-tax_table(physeq.infected) %>%
  as("matrix") %>%
  as_tibble(rownames = "OTU") %>%
  rename_all(tolower) %>%
  rename(asv_id= otu) %>%
  mutate(phylum= ifelse(phylum == "Proteobacteria", class, phylum)) %>%        
  select(asv_id, phylum)

## metadata table
metadata<-sample_data(physeq.infected) %>%
  as("matrix") %>%
  as_tibble(rownames = "samples") %>%
  mutate(PAH=str_replace(PAH, "0.1X", "CTRL"),
         PAH=str_replace(PAH, "10X", "PAHs")) %>%
  mutate(Parasite_PAH=paste0(PAH,"_",Parasite)) %>%
  select(samples, PAH, Parasite, Parasite_PAH)

## composite
composite<-inner_join(asv_table, taxo_table, by="asv_id") %>%
  group_by(samples, phylum) %>%
  summarize(count = sum(count), .groups="drop") %>%
  group_by(samples) %>%
  mutate(rel_abund=count/sum(count)*100) %>%
  ungroup() %>%
  group_by(phylum) %>%
  filter(phylum %in% vect_main_phylum) %>%
  ungroup() %>%
  select(-count) %>%
  inner_join(., metadata, by="samples")

## wilcox.test ctrl (n=5) vs pah (n=5)
wilcox_test_pah_phylum<-composite %>%
  nest(data=-phylum) %>%
  mutate(test.res=map(.x=data, ~wilcox.test(rel_abund ~ PAH,
                                            data=.x) %>% tidy)) %>%
  unnest(test.res) %>%
  mutate(p.adj=p.adjust(p.value, method="BH")) %>%
  select(phylum, p.value, p.adj)

## plot phylum
p2<-composite %>%
  mutate(rel_abund=rel_abund+1/20000,
         phylum=str_replace(phylum, "(.*)", "*\\1*")) %>%
  ggplot(aes(rel_abund, y=phylum, fill=PAH))+
  geom_boxplot(outlier.size = 0.5, alpha=0.4, lwd=0.2)+
  labs(x="Log10 (Relative Abundances)", y=NULL)+
  scale_fill_manual(values=c("#2166AC", "black"))+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text.y = element_markdown(),
        axis.title.x = element_text(size=9),
        legend.title=element_blank(),
        legend.position="top", legend.margin=margin(10,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

ggpubr::ggarrange(p1, p2, nrow=2, heights=c(1,1),
                  labels=c("A - Uninfected", "B - Infected"))

#8 plot genera abundances and associated wilcox.test independently for infected and uninfected chub (FIGURE S4 and TABLE S8) #####
library(phyloseq); packageVersion("phyloseq") #‘1.32.0’
library(tidyverse); packageVersion("tidyverse") #‘1.3.1’
library(tidyr); packageVersion("tidyr") #‘1.1.3’
library(dplyr); packageVersion("dplyr") #‘1.0.7’
library(broom); packageVersion("broom") #‘0.7.6’
library(ggtext); packageVersion("ggtext") #‘0.1.1’

rm(list=ls())
load("chub_pahs_parasite_phyloseq.RData")
load("vect_main_genus")

#_8.1 uninfected chub  (n=10 samples) ####
physeq.uninfected<-subset_samples(physeq, Parasite=="Uninfected")
physeq.uninfected<-prune_taxa(taxa_sums(physeq.uninfected)>0, physeq.uninfected)

## asv table
asv_table<-otu_table(physeq.uninfected) %>%
  as("matrix") %>%
  as_tibble(rownames="OTU") %>%
  rename(samples = OTU) %>%
  pivot_longer(-samples, names_to="asv_id", values_to="count")

## tax table
taxo_table<-tax_table(physeq.uninfected) %>%
  as("matrix") %>%
  as_tibble(rownames = "OTU") %>%
  rename_all(tolower) %>%
  rename(asv_id= otu) %>%
  mutate(phylum= ifelse(phylum == "Proteobacteria", class, phylum)) %>%        
  select(asv_id, phylum, genus) %>%
  mutate(genus=replace_na(genus, "unclass.")) %>%
  mutate(genus=paste0("(",phylum,") ", genus)) %>%
  filter(across(genus, ~ !grepl('unclass.', .))) %>%
  select(asv_id, genus)

## metadata table
metadata<-sample_data(physeq.uninfected) %>%
  as("matrix") %>%
  as_tibble(rownames = "samples") %>%
  mutate(PAH=str_replace(PAH, "0.1X", "CTRL"),
         PAH=str_replace(PAH, "10X", "PAHs")) %>%
  mutate(Parasite_PAH=paste0(PAH,"_",Parasite)) %>%
  select(samples, PAH, Parasite, Parasite_PAH)

## composite
composite<-inner_join(asv_table, taxo_table, by="asv_id") %>%
  group_by(samples, genus) %>%
  summarize(count = sum(count), .groups="drop") %>%
  group_by(samples) %>%
  mutate(rel_abund=count/sum(count)*100) %>%
  ungroup() %>%
  group_by(genus) %>%
  filter(genus %in% vect_main_genus) %>%
  ungroup() %>%
  select(-count) %>%
  inner_join(., metadata, by="samples")

## wilcox.test ctrl (n=5) vs pah (n=5)
wilcox_test_pah_genus<-composite %>%
  nest(data=-genus) %>%
  mutate(test.res=map(.x=data, ~wilcox.test(rel_abund ~ PAH,
                                            data=.x) %>% tidy)) %>%
  unnest(test.res) %>%
  mutate(p.adj=p.adjust(p.value, method="BH")) %>%
  select(genus, p.value, p.adj)

p1<-composite %>%
  mutate(rel_abund=rel_abund+1/20000,
         genus=str_replace(genus, "(.*)", "*\\1*")) %>%
  ggplot(aes(rel_abund, y=genus, fill=PAH))+
  geom_boxplot(outlier.size = 0.5, alpha=0.4, lwd=0.2)+
  labs(x="Log10 (Relative Abundances)", y=NULL)+
  scale_fill_manual(values=c("#2166AC", "black"))+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text.y = element_markdown(),
        axis.title.x = element_blank(),
        legend.title=element_blank(),
        legend.position="top", legend.margin=margin(10,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

#_8.2 infected chub  (n=10 samples) ####
physeq.infected<-subset_samples(physeq, Parasite=="Infected")
physeq.infected<-prune_taxa(taxa_sums(physeq.infected)>0, physeq.infected)

## asv table
asv_table<-otu_table(physeq.infected) %>%
  as("matrix") %>%
  as_tibble(rownames="OTU") %>%
  rename(samples = OTU) %>%
  pivot_longer(-samples, names_to="asv_id", values_to="count")

## tax table
taxo_table<-tax_table(physeq.infected) %>%
  as("matrix") %>%
  as_tibble(rownames = "OTU") %>%
  rename_all(tolower) %>%
  rename(asv_id= otu) %>%
  mutate(phylum= ifelse(phylum == "Proteobacteria", class, phylum)) %>%        
  select(asv_id, phylum, genus) %>%
  mutate(genus=replace_na(genus, "unclass.")) %>%
  mutate(genus=paste0("(",phylum,") ", genus)) %>%
  filter(across(genus, ~ !grepl('unclass.', .))) %>%
  select(asv_id, genus)

## metadata table
metadata<-sample_data(physeq.infected) %>%
  as("matrix") %>%
  as_tibble(rownames = "samples") %>%
  mutate(PAH=str_replace(PAH, "0.1X", "CTRL"),
         PAH=str_replace(PAH, "10X", "PAHs")) %>%
  mutate(Parasite_PAH=paste0(PAH,"_",Parasite)) %>%
  select(samples, PAH, Parasite, Parasite_PAH)

## composite
composite<-inner_join(asv_table, taxo_table, by="asv_id") %>%
  group_by(samples, genus) %>%
  summarize(count = sum(count), .groups="drop") %>%
  group_by(samples) %>%
  mutate(rel_abund=count/sum(count)*100) %>%
  ungroup() %>%
  group_by(genus) %>%
  filter(genus %in% vect_main_genus) %>%
  ungroup() %>%
  select(-count) %>%
  inner_join(., metadata, by="samples")

## wilcox.test ctrl (n=5) vs pah (n=5)
wilcox_test_pah_genus<-composite %>%
  nest(data=-genus) %>%
  mutate(test.res=map(.x=data, ~wilcox.test(rel_abund ~ PAH,
                                            data=.x) %>% tidy)) %>%
  unnest(test.res) %>%
  mutate(p.adj=p.adjust(p.value, method="BH")) %>%
  select(genus, p.value, p.adj)

p2<-composite %>%
  mutate(rel_abund=rel_abund+1/20000,
         genus=str_replace(genus, "(.*)", "*\\1*")) %>%
  ggplot(aes(rel_abund, y=genus, fill=PAH))+
  geom_boxplot(outlier.size = 0.5, alpha=0.4, lwd=0.2)+
  labs(x="Log10 (Relative Abundances)", y=NULL)+
  scale_fill_manual(values=c("#2166AC", "black"))+
  scale_x_log10()+
  theme_bw()+
  theme(axis.text.y = element_markdown(),
        axis.title.x = element_text(size=9),
        legend.title=element_blank(),
        legend.position="top", legend.margin=margin(10,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

ggpubr::ggarrange(p1, p2, nrow=2, heights=c(1,1),
                  labels=c("A - Uninfected", "B - Infected"))
