## set up environment
library(ape)
library(fields)
library(vegan)
library(philr)
library(Rcpp)
library(data.table)
library(parallel)

# source Jack Darcy's functions for collectors curves
source("../src/filtering_functions.r")


# Read in raw otu table and sample data

otu_file         <- "../data/processed/cleaned/raw_mm_16s_hiseqsabundance_table.csv"
sample_data_file <- "../data/processed/cleaned/mm_16s_hiseqs_metadata_table.csv"

sample_dat     <- fread(sample_data_file)

otu_table_raw  <- fread(otu_file)

# Make clean numeric abundance data.frame with samples as columns and ASVs as rows
sample_ids <- paste0("x",otu_table_raw$Group)
otu_table_raw[ , Group := NULL]
otu_mat <- t(setDF(otu_table_raw, rownames = sample_ids))

# Drop samples with less than 1000 reads
otu_mat <- otu_mat[ , colSums(otu_mat) > 10000]



figures_folder_fp <- "outputs/rarefaction"


# condense otu table by removing zeroes
# WARNING this speeds up collectors curves but rows no longer corespond to ASVs across samples
otu_mat_cond <- condense_otu_table(otu_mat)

# calculate collection curves for all samples
col_curve_list_raw <- multiple_collectors_curves(as.data.frame(otu_mat_cond),
                                                          met="sobs",
                                                          maxd=30000,
                                                          stepsize=500,
                                                          phy=NULL)


saveRDS(col_curve_list_raw, "../data/processed/binned_collectors_curve.rds")
# plot
pdf(paste(figures_folder_fp, "/all_samples_binned_collectors_curves.pdf", sep=""), useDingbats=F)
plot_multiple_binned_collectors_curves(col_curve_list_raw, ylab="Observed Richness of ASVs",
                                       n_bins=6,
                                       bin_d=20000,
                                       )
dev.off()

# calculate collection curves for each habitat
for(a_habitat in unique(sample_dat$habitat)){
  
  # pull out samples from a habitat
  
  habitat_samples <- sample_dat[habitat == a_habitat, sequencing_id]
  habitat_samples <- paste0("x",habitat_samples)
  
  habitat_otus <- otu_mat_cond[ , colnames(otu_mat_cond) %in% habitat_samples ]
  
  # generate collectors curves
  col_curve_list_habitat <- multiple_collectors_curves(as.data.frame(habitat_otus),
                                                   met="sobs",
                                                   maxd=30000,
                                                   stepsize=500,
                                                   phy=NULL)
  
  pdf(paste(figures_folder_fp, "/", a_habitat, "_binned_collectors_curves.pdf", sep=""), useDingbats=F)
  plot_multiple_binned_collectors_curves(col_curve_list_habitat, ylab="Observed Richness of ASVs",
                                         n_bins=6,
                                         bin_d= 20000,
  )
  dev.off()
  
  
}

# do a single collectors curve using vegan
all_samples <- rowSums(otu_mat)

# all_curve <- alpha_rare_curve_1sample(all_samples, stepsize  = 1e5, maxd = 1e8, seed = 2020)
# saveRDS(all_curve, "../data/processed/all_sample_curve.rds")

all_curve <-readRDS("../data/processed/all_sample_curve.rds")
# plot single curve

pdf(paste0(figures_folder_fp,"/all_samples_collectors_curve.pdf"))
plot(NULL, ylim=c(0, 4e5), xlim=c(0, 1e8), xlab="Number of sequences", ylab = "Number of ASVs")
points(all_curve$alpha ~ all_curve$step,  type = "l")
dev.off()
