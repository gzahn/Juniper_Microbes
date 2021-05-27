# Community analyses

# packages ####
library(phyloseq)
library(tidyverse)
library(vegan)
library(plotly)
library(ggbiplot)
library(patchwork)
library(microbiome)
library(ade4)
library(janitor)
library(ecodist)

# functions and palettes ####
source("./R/palettes.R")
source("./R/plot_bar2.R")

# load data ####
fung <- readRDS("./output/juniper_phyloseq_object_clean_ITS2.RDS")
bact <- readRDS("./output/juniper_phyloseq_object_clean_16S.RDS") # may need to go back and build a phylogeny?
meta <- read_csv("./metadata/cleaned_metadata.csv")
bact
fung
# chemical compound names
compound_names <- c("alpha.pinene","para.cymene","alpha.terpineol","cedr.9.ene","alpha.cedrene","beta.cedrene",
                    "cis.thujopsene","alpha.himachalene","beta.chamigrene","cuparene","compound.1","alpha.chamigrene",
                    "widdrol","cedrol","beta.acorenol","alpha.acorenol","gamma.eudesmol","beta.eudesmol",
                    "alpha.eudesmol","cedr.8.en.13.ol","cedr.8.en.15.ol","compound.2","thujopsenal")

# Clean bacterial taxonomy
bact <- subset_taxa(bact, Kingdom == "Bacteria")
bact <- subset_taxa(bact,Class != "Chloroplast")

# remove empty samples/taxa
bact <- subset_taxa(bact, taxa_sums(bact) > 0)
bact <- subset_samples(bact, sample_sums(bact) > 0)
fung <- subset_taxa(fung, taxa_sums(fung) > 0)
fung <- subset_samples(fung, sample_sums(fung) > 0)

# merge ps objects
full <- merge_phyloseq(fung,bact)

# convert to relabund ####
fung_ra <- transform_sample_counts(fung,function(x){x/sum(x)})
bact_ra <- transform_sample_counts(bact,function(x){x/sum(x)})
full_ra <- transform_sample_counts(full,function(x){x/sum(x)})


# quick barplots ####
p1 <- plot_bar2(full_ra,fill="Kingdom") + scale_fill_manual(values = pal.discrete)
p1
ggsave("./output/figs/Fig_Kingdom_by_sample.png",dpi=300)


p2 <- plot_bar2(fung_ra,fill="Phylum") + scale_fill_manual(values = pal.discrete) + labs(y="Relative abundance",title = "Fungi") + theme(panel.background = element_blank())
p3 <- plot_bar2(bact_ra,fill="Phylum") + scale_fill_manual(values = pal.discrete) + labs(y="Relative abundance",title = "Bacteria")+ theme(panel.background = element_blank())

p2 / p3
ggsave("./output/figs/phylum_barplot.png",height = 8,width = 6)


sampleid <- sample_data(full)$Sample.Name. %>% str_split("-") %>% map_chr(3)

fullmerged <- merge_samples(full,sampleid)

fullmerged %>% sample_data()
fullmerged %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>% 
  plot_bar2(fill="Phylum") + 
  scale_fill_manual(values = pal.discrete[2:3],labels=c("Bacteria","Fungi")) + 
  labs(y="Relative abundance",caption="Fig A") + 
  theme(panel.background = element_blank()) 
ggsave("./output/figs/bacteria_vs_fungi_barplot.png")
  
  
plot_bar2(full_ra,fill="Phylum") + scale_fill_manual(values = pal.discrete) + labs(y="Relative abundance") + theme(panel.background = element_blank())

# estimate alpha diversity ####
fung_ra@sam_data$richness <- estimate_richness(fung)$Observed
fung_ra@sam_data$shannon <- estimate_richness(fung)$Shannon
bact_ra@sam_data$richness <- estimate_richness(bact)$Observed
bact_ra@sam_data$shannon <- estimate_richness(bact)$Shannon


sample_data(fung_ra) %>% names


# PCoA ordination
fung_PCoA <- ordinate(fung_ra,method="PCoA",na.rm = TRUE)
bact_PCoA <- ordinate(bact_ra,method="PCoA",na.rm = TRUE)

fung_ra %>% sample_data()
# plot ordinations
p2 <- plot_ordination(fung_ra, fung_PCoA,color="YearsSinceBurn",title = "Fungi") +
  theme_classic() + theme(legend.position = "none")
p3 <- plot_ordination(bact_ra, bact_PCoA,color="YearsSinceBurn",title="Bacteria") +
  theme_classic() + labs(color="Years since\nburn")
p2 + p3



# NMDS Ordinations ####

# Get data separated
fung_asv <- otu_table(fung_ra) %>% as("matrix") %>% as.data.frame
bact_asv <- otu_table(bact_ra) %>% as("matrix") %>% as.data.frame

fung_compounds <- sample_data(fung_ra) %>% 
  meta() %>% 
  select(compound_names,YearsSinceBurn)

bact_compounds <- sample_data(bact_ra) %>% 
  meta() %>% 
  select(compound_names,YearsSinceBurn)

# NMDS
fung_NMDS <- metaMDS(fung_asv, k = 2, trymax = 100, trace = F, autotransform = TRUE, distance="bray")
fung_ef <- envfit(fung_NMDS,fung_compounds,permutations = 999)

bact_NMDS <- metaMDS(bact_asv, k = 2, trymax = 100, trace = F, autotransform = TRUE, distance="bray")
bact_ef <- envfit(bact_NMDS,bact_compounds,permutations = 999)

# plot NMDS w/ compounds
png("./output/figs/fungal_biplot_NMDS.png")
plot(fung_NMDS, type = "p", display = "sites",main="Fungal NMDS")
plot(fung_ef, p.max = 0.05)
dev.off()

png("./output/figs/bacterial_biplot_NMDS.png")
plot(bact_NMDS, type = "p", display = "sites",main="Bacterial NMDS")
plot(bact_ef, p.max = 0.05)
dev.off()

sample_data(fung_ra) %>% 
  meta()


# Mantel tests ####
fung_spatial.dist = dist(cbind(fung_ra@sam_data$Latitude, fung_ra@sam_data$Longitude))
fung_comm.dist =   vegdist(fung_asv)
bact_spatial.dist = dist(cbind(bact_ra@sam_data$Latitude, bact_ra@sam_data$Longitude))
bact_comm.dist =   vegdist(bact_asv)

fung_mantel.test = mantel.rtest(fung_spatial.dist, fung_comm.dist, nrepet = 9999)
bact_mantel.test = mantel.rtest(bact_spatial.dist, bact_comm.dist, nrepet = 9999)

sink("./output/mantel_test_results.txt")
fung_mantel.test
bact_mantel.test
sink(NULL)





# yield percent effect on alpha diversity

sample_data(fung_ra) %>% 
  meta() %>% 
  glm(formula = shannon ~ Yield_percent,
      data = .) %>% 
  summary()

sample_data(bact_ra) %>% 
  meta() %>% 
  glm(formula = shannon ~ Yield_percent,
      data = .) %>% 
  summary()


  

# sig compound model
fung_meta <- sample_data(fung_ra) %>% 
  meta() 
names(fung_meta) <- fung_meta %>% names() %>% make_clean_names()

bact_meta <- sample_data(bact_ra) %>% 
  meta() 
names(bact_meta) <- bact_meta %>% names() %>% make_clean_names()

sink("./output/beta_alpha_acorenal_effect_on_fungal_diversity.txt")
fung_meta %>% 
  glm(formula = shannon ~ beta_acorenol * alpha_acorenol + years_since_burn,
      data = .) %>% 
  summary()
sink(NULL)

fung_meta %>% 
  glm(formula = shannon ~ beta_acorenol * alpha_acorenol + years_since_burn,
      data = .) %>% summary()
  report::report(.)

bact_meta %>% 
  glm(formula = shannon ~ beta_acorenol * alpha_acorenol + years_since_burn,
      data = .) %>% summary()
  report::report(.)


fung_meta %>% 
  ggplot(aes(x=alpha_acorenol,y=beta_acorenol,color=shannon)) + 
  geom_jitter(size=6) +
  theme_classic() +
  scale_color_viridis_c() +
  labs(color="Shannon\ndiversity") +
  theme(axis.title = element_text(face="bold",size=12))
ggsave("./output/figs/fungal_diversity_alpha-beta-acorenal.png",dpi=300)



fung_meta %>% 
  select(shannon,alpha_acorenol,beta_acorenol) %>% 
  mutate(total_acorenal=alpha_acorenol+beta_acorenol) %>% 
  ggplot(aes(x=total_acorenal,y=shannon)) + 
  geom_point() + geom_smooth(method='lm')
  theme_classic() +
  scale_color_viridis_c() +
  labs(color="Shannon\ndiversity") +
  theme(axis.title = element_text(face="bold",size=12))


# what compounds important


# permanova ####
sink("./output/permanova_table_fungi-alpha-beta-acorenal_yearssinceburn.txt")
adonis(fung_asv ~ fung_meta$years_since_burn + fung_meta$alpha_acorenol * fung_meta$beta_acorenol)
sink(NULL)

sink("./output/permanova_table_bacteria-alpha-beta-acorenal_yearssinceburn.txt")
adonis(bact_asv ~ bact_meta$years_since_burn + bact_meta$alpha_acorenol * bact_meta$beta_acorenol)
sink(NULL)


# multiple regression on matrices ####

# new vector of compounds with clean names
compound_names <- c("alpha_pinene","para_cymene","alpha_terpineol","cedr_9_ene","alpha_cedrene","beta_cedrene",
                    "cis_thujopsene","alpha_himachalene","beta_chamigrene","cuparene","compound_1","alpha_chamigrene",
                    "widdrol","cedrol","beta_acorenol","alpha_acorenol","gamma_eudesmol","beta_eudesmol",
                    "alpha_eudesmol","cedr_8_en_13_ol","cedr_8_en_15_ol","compound_2","thujopsenal")


fung_asv_pa <- fung_asv
bact_asv_pa <- bact_asv
fung_asv_pa[fung_asv_pa > 0] <- 1
bact_asv_pa[bact_asv_pa > 0] <- 1


fung_compound_dist <- fung_meta %>% 
  select(compound_names) %>% 
  vegdist()

bact_compound_dist <- bact_meta %>% 
  select(compound_names) %>% 
  vegdist()


# what compounds important?
fung_meta %>% 
  pivot_longer(compound_names,names_to="compound") %>% 
  glm(data = .,
      formula = shannon ~ compound * value) %>% 
  summary()



# using relabund
fung_dist_MRM <- MRM(fung_comm.dist ~ fung_compound_dist, nperm = 9999)
bact_dist_MRM <- MRM(bact_comm.dist ~ bact_compound_dist, nperm = 9999)


# using presence/absence
fung_dist_MRM_logistic <- MRM(vegdist(fung_asv_pa) ~ fung_compound_dist, nperm = 999,method = "logistic")
bact_dist_MRM_logistic <- MRM(vegdist(bact_asv_pa) ~ bact_compound_dist, nperm = 999,method = "logistic")



# stacked barcharts ####
p5 <- fung_ra %>% 
  plot_bar2(fill = "Class") +
  scale_fill_manual(values = pal.discrete) +
  theme(panel.background = element_blank(),
        axis.title = element_text(size=12,face="bold"),
        legend.title = element_text(size=12,face="bold")) +
  labs(y="Relative abundance",fill="Fungal class")


p6 <- bact_ra %>% 
  plot_bar2(fill = "Class") +
  scale_fill_manual(values = pal.discrete) +
  theme(panel.background = element_blank(),
        axis.title = element_text(size=12,face="bold"),
        legend.title = element_text(size=12,face="bold")) +
  labs(y="",fill="Bacterial class")

p5 + p6
ggsave("./output/figs/fungal-and-bacterial_class_barchart.png",dpi=300,width = 10,height = 6)



# plots against chemical components ####

fung_meta %>% 
  pivot_longer(compound_names,names_to= "compound",values_to="amount") %>% 
ggplot(aes(x=amount,y=shannon,color=compound)) + geom_point() + theme_minimal() +
  geom_smooth(method = "lm") + labs(y="Fungal α-Diversity (Shannon)",x="Compound amount (unit?)",title="Fungal α-Diversity") + 
  facet_wrap(~compound,scales="free_x") +
  scale_color_manual(values=pal.discrete) + theme(strip.text = element_text(face="bold",size=16),
                                                  axis.title = element_text(face="bold",size=18),
                                                  legend.position = "none")
ggsave("./output/figs/Fungal_α-Diversity_over_chem-compound.png",dpi=300,width = 12,height = 10)

#bacteria
bact_meta %>% 
  pivot_longer(compound_names,names_to= "compound",values_to="amount") %>% 
  ggplot(aes(x=amount,y=shannon,color=compound)) + geom_point() + theme_minimal() +
  geom_smooth(method = "lm") + labs(y="Bacterial α-Diversity (Shannon)",x="Compound amount (unit?)",title="Bacterial α-Diversity") + 
  facet_wrap(~compound,scales="free_x") +
  scale_color_manual(values=pal.discrete) + theme(strip.text = element_text(face="bold",size=16),
                                                  axis.title = element_text(face="bold",size=18),
                                                  legend.position = "none")
ggsave("./output/figs/Bacterial_α-Diversity_over_chem-compound.png",dpi=300,width = 12,height = 10)

