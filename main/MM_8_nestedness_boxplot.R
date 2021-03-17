###This script is used to generate a simple boxplot of nestedness indices for each site
# nestedness outputs were concatenated using cat *.out > All_Waimea_Sites.out

library(ggplot2)
library(gridExtra)
library(ggpubr)

source("../src/get_colors.R")

h2_out  <-read.table("../data/processed/nestedness_outputs/all_samples_empo_3_H2.out",
                     sep = "\t",
                     col.names = c("File","Index","Value"))

nodf_out <- read.table("../data/processed/nestedness_outputs/all_samples_empo_3_NODF.out",
                       sep = "\t",
                       col.names = c("File","Index","Value"))

nest_out <- rbind(h2_out, nodf_out)
nest_out$Group <- sub(".+\\/(.+).csv", "\\1",nest_out$File)

# pull out nestedness values for habi-sites (i.e. separated values per site and habitat)
# this means terrestrial and riverine are separated out
# contains both rarefied and urarefied data
habisites <- nest_out[grepl("habisite", nest_out$Group), ]

# clean up table contents
habisites$rarefied <- ifelse(grepl("rar",habisites$Group), yes = T, no = F)
habisites$HabiSite <- sub(".*habisite_(.+)_empo_3", "\\1",habisites$Group)
habisites$Value    <- as.numeric( gsub("\\[(.+)\\]","\\1", habisites$Value) )
habisites$Habitat  <- sub("(.+)-.+","\\1", habisites$HabiSite)
habisites$Site     <- sub(".+-(.+)","\\1", habisites$HabiSite)



# boxplots for each nestedness index

# set up comparisons for means
my_comparisons <- list( c(1, 2), c(2, 3), c(1, 3) )

# plot_nest_box() just subses the data to a specified index and makes a boxplot
plot_nest_box <- function(an_index, nest_data, description) {
  # subset by index
  nest_ind <- nest_data[nest_data$Index == an_index, ]
  # plot boxplot by Habitat
  p <-
    ggplot(nest_ind, aes(x = Habitat, y = Value, fill = Habitat)) +
    geom_boxplot()+
    geom_jitter(width = .2)+
    theme_pubr()+
    ylab(an_index)+ 
    labs(title = an_index, caption = "T-test")+
    
    stat_compare_means(comparisons = my_comparisons, method = "t.test") # Add pairwise comparisons p-value
  
  aov_out   <- aov(lm(Value ~ Habitat, nest_ind))
  tukey_out <- TukeyHSD(aov_out)
  
  # write out tukey test and individual plot
  capture.output(file = paste0("outputs/nestedness_boxplots/", description, "_tukey_test.txt"),
                 print( tukey_out)
                 )
  ggsave(filename = paste0("outputs/nestedness_boxplots/",description,"_empo3_",an_index,"_boxplot.pdf"),
         width = 4, height = 6)
  return(p)
}

rar_habisite_plots <- lapply(X =  unique(habisites$Index),
                            FUN = plot_nest_box,
                            nest_data = habisites[habisites$rarefied == T, ],
                            description = "rarefied")
               
raw_habisite_plots <- lapply(X =  unique(habisites$Index),
                            FUN = plot_nest_box,
                            nest_data = habisites[habisites$rarefied == F, ],
                            description = "raw")

ggsave(plot = rar_habisite_plots[[1]],
       filename = "outputs/nestedness_boxplots/rarified_")


g1 <-arrangeGrob(rar_habisite_plots[[1]], rar_habisite_plots[[2]],
                nrow = 1, ncol = 2)

ggsave("outputs/nestedness_boxplots/rarefied_empo_3_agg_habitat_nest_boxplot.pdf", plot = g1, width = 8, height = 6)

g2 <-arrangeGrob(raw_habisite_plots[[1]], raw_habisite_plots[[2]],
                 nrow = 1, ncol = 2)

ggsave("outputs/nestedness_boxplots/raw_count_empo_3_agg_habitat_nest_boxplot.pdf", plot = g2, width = 8, height = 6)

write.csv(habisites, "../data/processed/H2_NODF_empo3_agg_habitat_site.csv")

## Generate Vegan nestedness plot to provide visual alongside boxplots

# The first part is memory intensive and was run on a high performance computing cluster
# Plotting, which occurs after the stop point, can be done locally.

library(data.table)



# Generate data on HPC server --------------------------------------------
# read_hpc_abund() is a function for reading and cleaning C-MAIKI pipeline OTU table
read_hpc_abund <- function(abundance_file){
  
  # read in
  abundance_table <- fread(abundance_file,nrows = 1000)
  
  # drop extra columns
  abundance_table[, c('label','numOtus') :=NULL]i
  
  # transpose, keep OTU names
  abundance_table <- transpose(abundance_table, keep.names = "OTU_ID", make.names = "Group")
  
  setkey(abundance_table,"OTU_ID")
  
  return(abundance_table)
}


# identify data files

otu_file <-   "abundance_table_100.shared"
meta_file  <- "combined_prep_sample_metadata_16s.csv"

# read in all OTUs
otus <- read_hpc_abund(otu_file, header = T)

# read in metadata
sam_dat <- fread(meta_file)

# re-classify trophic level
sam_dat[trophic %in% c("Omnivore","Detritivore","Herbivore"), trophic := "Herbivore/Omnivore"]

# sum reads by trophic level
t_otus <- transpose(otus, keep.names = "sequencing_id", make.names = "OTU_ID")
t_otus[ , sequencing_id := as.integer(sequencing_id)]
troph_dat <- merge(sam_dat[ , .(sequencing_id, trophic)], t_otus, by = "sequencing_id")
troph_dat[ , sequencing_id := NULL]
troph_dat <- troph_dat[ , lapply(.SD, function(x) ifelse(any(x > 1), 1, 0)), by = trophic]


# reorder matrix by trophic levels
trophic_order <- c("Environmental",
                   "PrimaryProducer",
                   "Herbivore/Omnivore",
                   "Carnivore")
troph_otus <- as.matrix(troph_dat, rownames = "trophic")
troph_otus <- troph_otus[trophic_order, ]
troph_otus[1:3,1:5]

# write out data file
write.csv(troph_otus,"all_ASVs_trophic_matrix.csv", row.names = T)

stop()

# Plot on local computer ----------------------------------------------
library(data.table)
library(vegan)
troph_otus <-as.matrix(fread("../data/processed/all_ASVs_trophic_matrix.csv"))



# rownames of trophic presence absence matrix
trophic_order <- c("Environmental",
                   "PrimaryProducer",
                   "Herbivore/Omnivore",
                   "Carnivore")

row.names(troph_otus) <- trophic_order

# vegan nestedness plot

troph_nest <- nestedtemp(troph_otus)

png(file  = "outputs/nestedness_boxplots/vegan_nestedness_heatmap.png",
    width = 10, height = 4, units = "in", res = 700)
par(mar = c(4,10,4,4))
plot(troph_nest, kind = "incidence", col =c("white","darkgray"), names = c(T,F))
dev.off()

pdf(file  = "outputs/nestedness_boxplots/vegan_nestedness_heatmap.pdf",
    width = 40, height = 40)
par(mar = c(4,10,4,4))
plot(troph_nest, kind = "incidence", col =c("white","darkgray"), names = c(T,F))
dev.off()

