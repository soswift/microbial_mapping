# This script is for generating stacked barplots showing spread across trophic levels and EMPO for top Waimea ASVs
# Sean Swift 10/20

library(ggplot2)
library(phyloseq)
library(data.table)

# Set up environment
# get standardized colors for plotting ecological variables
source("../src/get_colors.R")
waimea_colors <- get_waimea_colors(return_all = T)

# set ggplot theme
theme_set(theme_minimal())

# Read in phyloseq object
MM_16S_phy <- readRDS("../data/processed/cleaned/MM_16S_phy.rds")

# Pull out top 200 ASVs
top_asv_names <- names(sort(taxa_sums(MM_16S_phy), decreasing = TRUE)[1:200])
top_phy       <- prune_taxa(top_asv_names, MM_16S_phy)

# pull out abundance counts and sample data as data.tables
top_asv_table  <- setnames( data.table(as(otu_table(top_phy),"matrix"),
                                       keep.rownames = T)
                            ,"rn", "OTU_ID")
top_asv_sample <- setnames(data.table(as(sample_data(top_phy), "matrix"),
                                      keep.rownames = T),
                           "rn", "sequencing_id")

# reshape into long format for ggplot
long_asvs <- melt.data.table(top_asv_table,
                id.vars = "OTU_ID",
                variable.name = "sequencing_id",
                value.name = "rel_abund")

# add in relevant sample data for filling stacked barplots
eco_vars <-c("sequencing_id", "habitat", "trophic", "empo_1", "empo_2", "empo_3")


long_asvs <- merge(long_asvs,
                   top_asv_sample[ , ..eco_vars],
                   by = "sequencing_id")


# plot out stacked barcharts differentially filled with ecologically informative variables (habitat type, trophic leve, and EMPO)


eco_plots <- lapply(
                  eco_vars[-1],
                  FUN = function(an_eco_var) {
                    p <- ggplot(long_asvs, aes_(quote(OTU_ID),quote(rel_abund), fill = as.name(an_eco_var))) +
                         geom_bar(position = "fill", stat = "identity") +
                         labs(x = "Top ASVs", y = "Percent Of Total Relative Abundance")+
                         # Standardized color scheme
                         scale_fill_manual(values = waimea_colors, drop = T) +
                         theme(axis.text = element_blank() ,
                               axis.ticks = element_blank())
                         
                    
                    ggsave(
                      filename = paste(
                        "outputs/stacked_barplots/", an_eco_var, "_stacked_barplot.png"),
                      dpi = 300,
                      plot = p )
                    
                    return(p)
                  })
