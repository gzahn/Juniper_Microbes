# Community analyses

# packages ####
library(phyloseq)
library(tidyverse)
library(vegan)
library(plotly)


# functions and palettes ####
source("./R/palettes.R")
source("./R/plot_bar2.R")

# load data ####
fung <- readRDS("./output/juniper_phyloseq_object_clean_ITS2.RDS")
bact <- readRDS("./output/juniper_phyloseq_object_clean_16S.RDS") # may need to go back and build a phylogeny?
meta <- read_csv("./metadata/cleaned_metadata.csv")

sample_names(fung) <- meta$Sample_ID
sample_names(bact) <- meta$Sample_ID

# convert to relabund ####
fung_ra <- transform_sample_counts(fung, function(x) x/sum(x))
bact_ra <- transform_sample_counts(bact, function(x) x/sum(x))


# quick barplots ####
p1 <- plot_bar2(fung_ra,fill="Phylum") + scale_fill_manual(values = pal.discrete)
p1
ggsave("./output/figs/Fig9_fungal_phyla_by_sample.png",dpi=300)

p2 <- plot_bar2(fung_ra,fill="Class") + scale_fill_manual(values = pal.discrete)
p2
ggsave("./output/figs/Fig10_fungal_class_by_sample.png",dpi=300)

p3 <- plot_bar2(fung_ra,fill="Order") + scale_fill_manual(values = pal.continuous[round(seq(from=1,to=180,length.out = 41))])


ggplotly(p3)


plot_bar2(bact_ra,fill="Phylum")



# estimate alpha diversity ####   ?????
                # fung.alphadiv <- estimate_richness(fung_ra)
                # bact.alphadiv <- estimate_richness(bact_ra)
                # 
                # alpha.bact <- merge(meta,fung.alphadiv)
                # alpha.fung <- merge(meta,fung.alphadiv)

meta$fung_alphadiv <- diversity(otu_table(fung_ra))
meta$bact_alphadiv <- diversity(otu_table(bact_ra))


# simple plots against burn year
ggplot(meta, aes(x=YearsSinceBurn,y=fung_alphadiv)) + geom_point() + theme_minimal() +
  geom_smooth(method = "lm") + labs(y="Fungal α-Diversity (Shannon)",x="Years Since Burn")
ggsave("./output/figs/Fig7_Fungal_alphdiv_over_burnyear.png",dpi=300)



ggplot(meta, aes(x=YearsSinceBurn,y=bact_alphadiv)) + geom_point() + theme_minimal() +
  geom_smooth(method = "lm") + labs(y="Bacterial α-Diversity (Shannon)",x="Years Since Burn")
ggsave("./output/figs/Fig8_Bacterial_alphdiv_over_burnyear.png",dpi=300)


# make long-form of chemical componenets
chemcols <- names(meta)[9:31]
meta_long <- gather(meta, "Compound", "Amount",chemcols)


# plots against chemical components

#fungi
ggplot(meta_long, aes(x=Amount,y=fung_alphadiv,color=Compound)) + geom_point() + theme_minimal() +
  geom_smooth(method = "lm") + labs(y="Fungal α-Diversity (Shannon)",x="Compound Amount (unit?)",title="Fungal α-Diversity") + facet_wrap(~Compound,scales="free_x") +
  scale_color_manual(values=pal.discrete) + theme(strip.text = element_text(face="bold",size=16),
                                                  axis.title = element_text(face="bold",size=18),
                                                  legend.position = "none")
ggsave("./output/figs/Fig5_Fungal_α-Diversity_over_chem-compound.png",dpi=300,width = 12,height = 10)

#bacteria
ggplot(meta_long, aes(x=Amount,y=bact_alphadiv,color=Compound)) + geom_point() + theme_minimal() +
  geom_smooth(method = "lm") + labs(y="Bacterial α-Diversity (Shannon)",x="Compound Amount (unit?)",title="Bacterial α-Diversity") + facet_wrap(~Compound,scales="free_x") +
  scale_color_manual(values=pal.discrete) + theme(strip.text = element_text(face="bold",size=16),
                                                  axis.title = element_text(face="bold",size=18),
                                                  legend.position = "none")
ggsave("./output/figs/Fig6_Bacterial_α-Diversity_over_chem-compound.png",dpi=300,width = 12,height = 10)





# to-do list:
# find taxa that are abundant
  # compare each fungal class over years since burn and chemical abundance
  # which fungal groups are indicitave of years since burn or chemical analytics (corncob!?)





# model class relabundance against years since burn ####
mod1 <- 


