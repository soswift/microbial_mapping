# This script is for generating stacked barplots showing spread across trophic levels and EMPO for top Waimea ASVs
# Sean Swift 10/20

library(ggplot2)
library(phyloseq)
library(data.table)
library(ggpubr)

# Set up environment
# get standardized colors for plotting ecological variables
source("../src/get_colors.R")
waimea_colors <- get_waimea_colors()

# set ggplot theme
theme_set(theme_minimal())

# Read in phyloseq object
MM_16S_phy <- readRDS("../data/processed/cleaned/MM_16S_phy.rds")


# simplify trophic levels to Environmental, PrimaryProducer, and Consumer
sample_data(MM_16S_phy)$trophic[!(sample_data(MM_16S_phy)$trophic %in%
                                    c("Environmental", "PrimaryProducer"))] <- "Consumer"

# Pull out the most abundant 200 ASVs from the whole study
top_asv_names <- names(sort(taxa_sums(MM_16S_phy), decreasing = TRUE))#[1:5000])
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
eco_vars <-c("sequencing_id", "habitat", "trophic", "empo_1", "empo_2", "empo_3", "site_name")


long_asvs <- merge(long_asvs,
                   top_asv_sample[ , ..eco_vars],
                   by = "sequencing_id")


# plot out stacked barcharts filled with ecologically informative variables (habitat type, trophic leve, and EMPO)


eco_plots <- lapply(
                  eco_vars[-1],
                  FUN = function(an_eco_var) {
                    
                    dat <- copy(long_asvs[ , .(rel_abund = sum(rel_abund)), by = c("OTU_ID",an_eco_var)])
                    
                    # Order the bars by distribution across the ecological variable (distance matrix + hclust)
                    
                    cast_f <- as.formula(paste("OTU_ID ~", an_eco_var))
                    wide <- dcast(dat, cast_f, value.var = "rel_abund")
                    wide_mat <- as.matrix(wide[ , -1])
                    row.names(wide_mat) <- wide[ , OTU_ID]

                    bar_order <- hclust(d = dist(wide_mat), method = "ward.D")
                    dat[ , OTU_ID := factor(OTU_ID, levels = bar_order$labels[bar_order$order])]
                    
                    # Plot the stacked bars
                    p <- ggplot(dat, aes_(quote(OTU_ID),quote(rel_abund), fill = as.name(an_eco_var)), wide_mat) +
                         geom_bar(position = "fill", stat = "identity") +
                         labs(x = "Top ASVs", y = "Percent Of Total Relative Abundance")+
                         # Standardized color scheme
                         scale_fill_manual(values = waimea_colors, drop = T) +
                         theme_pubclean()+
                         theme(axis.text.x = element_blank() ,
                               axis.ticks.x = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.background = element_blank(),
                               legend.position = "right")
                    
                    ggsave(
                      filename = paste(
                        "outputs/stacked_barplots/", an_eco_var, "_stacked_barplot.png"),
                      dpi = 300,
                      plot = p )
                    
                    return(p)
                  })

# convert plot to grob
stack_grobs <-lapply(eco_plots, ggplotGrob)

# set a standard width for all grobs
std_width <- stack_grobs[[1]]$widths
stack_std <- lapply(stack_grobs, function(x) {
  x$widths <- std_width
  return(x)})
# arrange ordinations
g <- ggarrange(plotlist=stack_std,
               nrow = 3, ncol = 2)


ggsave(paste0("outputs/stacked_barplots/facet_stacked_barplot.pdf"),
       plot = g, 
       width = 10,
       height = 15)
