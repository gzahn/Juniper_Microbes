#####################################################
#          Explore and visualize metadata           #
# This includes chemical and entomological data,    #
# as well as location and age                       #
# Author: Geoff Zahn                                #
#####################################################


# Packages
library(tidyverse)
library(readxl)
library(skimr)
library(patchwork)

# Custom color palettes
source("./R/palettes.R")

# Load metadata
meta <- read_xlsx("./metadata/metadata.xlsx")

# First look
glimpse(meta)
skim(meta)

# Extra tidying and prep ####

# Fix column classes
meta$Living_Larvae <- as.numeric(meta$Living_Larvae)
meta$Raw_Exit_Holes <- as.numeric(meta$Raw_Exit_Holes)

# Find analytical chemistry columns / rename compounds
chem_cols <- grep(pattern = "CHEM_",names(meta))
names(meta)[chem_cols] <- str_remove(names(meta)[chem_cols],"CHEM_")

chem_cols <- names(meta)[chem_cols]

# Add total chemical (relative area %) and mean
meta$ChemTotal <- rowSums(meta[,chem_cols])
meta$ChemMean <- meta$ChemTotal / length(chem_cols)

# Find entomology columns
names(meta)
ento_cols <- c("Bolt_Surface_Area_cm2","Raw_Exit_Holes_per_cm2","Raw_Exit_Holes","Living_Larvae")


# Rearrange rows by burn date
meta <- arrange(meta,desc(BurnYear))

# convert Sample_ID to factor
meta$Sample_ID <- factor(meta$Sample_ID,levels = meta$Sample_ID)

# Years since burn
meta$YearsSinceBurn <- 2020 - meta$BurnYear

# Export cleaned and prepped metadata ####
write_csv(meta,"./metadata/cleaned_metadata.csv")

# gather long version of meta with "Chem" column
longmeta <- gather(meta,key = Compound,value = Area,chem_cols)

# Plots of chemistry ####
ggplot(meta,aes(x=Sample_ID,y=ChemTotal)) + geom_col()

# chem over time (by compound)
chem_time_plot <- ggplot(longmeta,aes(x=YearsSinceBurn,y=Area,colour=Compound)) + 
  geom_point(size=4,alpha=.5) + geom_smooth(se=FALSE,method="lm",color="Black",linetype=2) +
  facet_wrap(~Compound) + 
  scale_color_viridis_d() + 
  labs(x="Years since tree death",y="Relative Area %") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle=60,hjust=1,face="bold"),
        strip.text = element_text(size=12,face="bold"),
        axis.title = element_text(size=14,face="bold")) 

chem_time_plot  
ggsave("./output/figs/Fig1_Chem_Area_over_Time.png",dpi=300,width = 12,height = 8)

# chem over time (mean)
chem_total_plot <- ggplot(meta,aes(x=YearsSinceBurn,y=ChemMean)) + 
  geom_point() + 
  geom_smooth(se=TRUE,method="lm",color="Black",linetype=2) + 
  theme_bw() +
  labs(x="Years since tree death",y="Total measured chem relative area %")
chem_total_plot

chem_mean_plot <- ggplot(longmeta,aes(x=YearsSinceBurn,y=Area,color=Raw_Exit_Holes)) + 
  geom_jitter(size=3,alpha=.5,height = 0,width = .2) +
  geom_smooth(se=TRUE,method="lm",color="Black",linetype=2) + 
  labs(x="Years since tree death",y="Relative area %",color="Number of insect\nexit holes",
       caption = "Relative area % shown for all chemical compounds") +
  theme_bw() +
  theme(axis.title = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"))
chem_mean_plot 
ggsave("./output/figs/Fig2_Chem_Totals_over_Time.png",dpi=300,width = 4,height = 4)

# Yield percent over time
yield_time_plot <- ggplot(longmeta,aes(x=YearsSinceBurn,y=Yield_percent,color=Raw_Exit_Holes)) + 
  geom_point(size=4) + geom_smooth(se=TRUE,method="lm",color="Black",linetype=2) +
  theme_bw() +
  labs(x="Years since tree death",y="Yield (% w/w)",color= "Number of insect\nexit holes") +
  theme(axis.title = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"))
yield_time_plot
ggsave("./output/figs/Fig3_Yield_over_Time.png",dpi=300,width = 6,height = 6)

names(longmeta)


#insect holes over time
insect1 <- ggplot(meta,aes(x=YearsSinceBurn,y=Raw_Exit_Holes)) + geom_point(size=2) + theme_bw() + 
  geom_smooth(method = "lm",color="Black",linetype=2) +
  labs(x="",y="Insect exit hole count") + 
  theme(axis.title = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"))

insect2 <- ggplot(meta,aes(x=YearsSinceBurn,y=Raw_Exit_Holes_per_cm2)) + geom_point(size=2) + theme_bw() +
  geom_smooth(method = "lm",color="Black",linetype=2) +
  labs(x="Years since tree death",y="Insect exit holes per cm2")  + 
  theme(axis.title = element_text(size=14,face="bold"),
        axis.text = element_text(size=12,face="bold"))

insect1 / insect2
ggsave("./output/figs/Fig4_Insect_bore-holes_over_Time.png",dpi=300,height = 4,width = 6)


# Models ####
chem_mod <- glm(data=longmeta,Area ~ YearsSinceBurn * Compound)

sink("./output/stats/glm_compounds_over_time.txt")
summary(chem_mod)
sink(NULL)

