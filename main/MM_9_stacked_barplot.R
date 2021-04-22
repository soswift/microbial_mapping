# This script is for generating stacked barplots showing spread across trophic levels and EMPO for Waimea ASVs
# Sean Swift 4/21

library(ggplot2)
library(data.table)
library(ggpubr)
library(parallel)

# Set up environment -----------------
# get standardized colors for plotting ecological variables
source("../src/get_colors.R")
waimea_colors <- get_waimea_colors()

# params
set.seed(2020)
options(mc.cores = parallel::detectCores() -1)

# set ggplot theme
theme_set(theme_minimal())

# Define functions ---------------------
# read_hpc_abund() is a function for reading and cleaning C-MAIKI pipeline OTU table
read_hpc_abund <- function(abundance_file){

  # read in
  abundance_table <- fread(abundance_file)

  # drop extra columns
  abundance_table[, c('label','numOtus') :=NULL]

  # transpose, keep OTU names
  abundance_table <- transpose(abundance_table,
                               keep.names = "OTU_ID",
                               make.names = "Group")

  setkey(abundance_table,"OTU_ID")

  return(abundance_table)
}


# Identify data files --------------------
#HPC file locations
 otu_file <-   "/home/seanos/cmaiki-lts/seanos/waimea/pipeline/waimea_16S_v3/abundance_table_100.shared"
 meta_file  <- "/home/seanos/cmaiki-lts/seanos/waimea/pipeline/waimea_16S_v3/combined_prep_sample_metadata_16s.csv"

asv_table <- read_hpc_abund(otu_file)

# read in metadata
sam_dat <- fread(meta_file, colClasses = "character")

# re-classify trophic level
sam_dat[!(trophic %in% c("Environmental", "PrimaryProducer")), trophic := "Consumer"]

# make sure metadata plays nice with asv table
sam_dat <- sam_dat[sample_type != "" & sequencing_id %in% names(asv_table)]

# subset ASVs
keep_cols <- c("OTU_ID", sam_dat$sequencing_id)
asv_table <- asv_table[ , ..keep_cols]
asv_table <- asv_table[rowSums(asv_table[ , -1]) > 0 ]

some_asvs <- sample(1:nrow(asv_table), 5000)
asv_table <- asv_table[ some_asvs ,]



# reshape data into long format for ggplot
long_asvs <- melt.data.table(asv_table,
                             id.vars = "OTU_ID",
                             variable.name = "sequencing_id",
                             value.name = "rel_abund")



# identify relevant sample data for filling stacked barplots
eco_vars <-c("sequencing_id",
             "habitat",
             "trophic",
             "empo_1",
             "empo_2",
             "empo_3",
             "site_name")

asvs_w_meta <- merge(long_asvs, sam_dat[ , ..eco_vars])

fwrite(asvs_w_meta, "random_5k_asvs.csv")

# Plot out stacked barcharts ------------
# filled barcharts with ecologically informative variables (habitat type, trophic level, and EMPO)


# TODO: Change to pa instead of abund  ************************

# Restart by reading in the intermediate file if you need to alter the plots
asvs_w_meta <- fread("../data/processed/stacked_bars/random_5k_asvs.csv",
                     header = T)
# add a column for sample occupancy rather than read abundance
asvs_w_meta[ , pa := rel_abund][ pa > 0 , pa := 1]

eco_plots <-  lapply(
                  X = eco_vars[-1],
                  FUN = function(an_eco_var) {
                      
                    print(an_eco_var)
                    dat <- copy(asvs_w_meta)
                    
                    # sum by abundance by the ecological variable
                    dat <- dat[ ,.(pa = sum(pa)), by = c(an_eco_var, "OTU_ID")]
                    
                    
                    # order the bars using clustering (distance matrix + hclust)
                    cast_f <- as.formula(paste("OTU_ID ~", an_eco_var))
                    wide <- dcast(dat,
                                  cast_f,
                                  value.var = "pa")
                    
                    wide_mat <- as.matrix(wide[ , -1])
                    row.names(wide_mat) <- wide[ , OTU_ID]
                    print("clustering")
                    
                    bar_order <- hclust(d = dist(wide_mat),
                                        method = "ward.D2")
                    
                    dat[ , OTU_ID := factor(OTU_ID, levels = bar_order$labels[bar_order$order])]
                    print("plotting")
                    
                    # Plot the stacked bars
                    p <- ggplot(dat, aes_(x = quote(OTU_ID), 
                                          y = quote(pa),
                                          fill = as.name(an_eco_var))) +
                         geom_bar(position = "fill",
                                  stat = "identity") +
                         labs(x = "ASVs", y = "Proportion of total occurances")+
                         # Standardized color scheme
                         scale_fill_manual(values = waimea_colors, drop = T) +
                         theme_pubclean()+
                         theme(axis.text.x = element_blank() ,
                               axis.ticks.x = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               panel.background = element_blank(),
                               legend.position = "right")
                    print("saving")
                    ggsave(
                      filename = paste0("outputs/stacked_barplots/full_data_subset/", an_eco_var, "_stacked_barplot.pdf"),
                      width = 10,
                      heigh = 7,
                      plot = p )
                    
                    return(p)
                  }
                  )
saveRDS(eco_plots, "../data/processed/stacked_bars/plots.rds")


# convert plot to grob
stack_grobs <-mclapply(eco_plots, ggplotGrob)

# set a standard width for all grobs
std_width <- stack_grobs[[1]]$widths
stack_std <- lapply(stack_grobs, function(x) {
  x$widths <- std_width
  return(x)})
# arrange ordinations
g <- ggarrange(plotlist=stack_std,
               nrow = 3,
               ncol = 2,
               labels = as.list( letters[1:length(eco_plots)]) )


ggsave("outputs/stacked_barplots/full_data_subset/facet_stacked_barplot.pdf",
       plot = g, 
       width = 10,
       height = 15)
