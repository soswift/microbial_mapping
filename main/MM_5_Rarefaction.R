# Microbial Mapping - Rarefaction curves
# Document Summary
# This code is for running rarefaction curves on full the Waimea Dataset
# Written as an rscript instead of notebook, because analysis will be run on server
#
# For more information about iNext see:
# http://chao.stat.nthu.edu.tw/wordpress/wp-content/uploads/software/iNEXT_UserGuide.pdf



library(data.table)
library(iNEXT)
library(ggplot2)


# Define Functions ------------


## print.size()  prints size and dimensions of object
print.size <- function(object, description, units = "MB") {
  size <- object.size(object)
  message(paste(description, "size = ", sep = "\n"))
  print(size, units = units)
  
  if (is.vector(object)) {
    message("length = ")
    length(object)
  } else {
    message("dimensions =")
    print(paste(dim(object)))
  }
}


## read_hpc_abund(): reads in abundance file from 16S metaflowmics pipeline.
# Cleans up columns, transposes so OTUs are rows and samples are columns.
read_hpc_abund <- function(abundance_file){
  
  # read in
  abundance_table <- fread(abundance_file)
  
  # drop extra columns
  abundance_table[, c('label','numOtus') :=NULL]
  
  # transpose, keep OTU names
  abundance_table <- transpose(abundance_table, keep.names = "OTU_ID", make.names = "Group")
  
  setkey(abundance_table,"OTU_ID")
  
  return(abundance_table)
}

## set_standardID(): Identify the metadata column that contains sample names
# In all other functions, OTUs and metadata will always be linked using column standardID
set_standardID <- function(metadata, ID_column){
  ids <- metadata[[ID_column]]
  metadata[, standardID := ids]
}

## sub_by_samples(): subsets an OTU table to match a vector of sample names.
# Assumes that OTU table columns are samples, rows are otus.
# If a column of OTU names should kept, indicate the column name. 

sub_by_samples <- function(samples, otu_table, OTU_name_column) {
  # if OTU name column specified, add it to the vector of sample names
  if(!is.null(OTU_name_column)){
    samples <- c(OTU_name_column, samples)
  }
  
  # subset columns
  new_otu  <- otu_table[ , ..samples]
  
  # remove OTUs with rowSums that = 0
  new_otu[, Sums := rowSums(.SD), .SDcols = !OTU_name_column]
  new_otu <- new_otu[Sums != 0, .SD, .SDcols = !"Sums"]
  
  return(new_otu)
}

## make_pa(): converts a list to presence absence.
# Use lapply to convert a full OTU abundance table to presence absence data

make_pa <- function(Col) {
  Col[Col > 0] <- 1
  return(Col)
}

## abund_to_incidence(): converts OTU abundance data to incidence counts for use with iNext.
# Incidence counts are the number of samples taken followed by number of times each species was observed.
# Assumes OTUs are rows, samples are columns. Optionally exclude column containing OTU names.
# By default, all samples are included. To subset by certain samples, pass a list to samples.

abund_to_incidence <- function(samples = "All", otu_table, OTU_name_column = "OTU_ID"){
  
  # Keep all samples unless provided vector of sample names
  if(length(samples) == 1 && samples == "All"){
    otu_table <- otu_table
    
  } else{
    otu_table <- sub_by_samples(samples = samples,
                           otu_table = otu_table,
                           OTU_name_column = OTU_name_column)
  }
  
  # convert OTU table to presence absence
  otus <- otu_table[, lapply(.SD, make_pa ), .SDcols = !OTU_name_column]
  
  # count n observations for each OTU
  counts <- rowSums(otus)
  
  # count total samples (i.e. columns)
  total_samples <- ncol(otus)
  
  # put total samples in front of sample counts
  incidence_list <- list( c((total_samples),counts))
  
  return(incidence_list)
}


# Define Inputs ------------

# Toy data version
data_file <- "../../toy_data/toy_data/toy_hpc_data.csv"
metadata_file <- "../../toy_data/toy_data/toy_hpc_metadata.csv"

# # Full data version
# hpc_abund_tsv <- "inputs/waimea_16S_v3/abundance_table_100.shared"
# hpc_meta_csv  <- "inputs/waimea_16S_v3/combined_prep_sample_metadata_16s.csv"


# Read in ----------

# Read in toy data
otus <- fread(data_file, header = T)
meta <- fread(metadata_file)

# Read in full data
# otus <- read_hpc_abund(hpc_abund_tsv)
# meta <- fread(hpc_meta_csv)
# meta <- meta[sequencing_id %in% colnames(otus)]


# set key for linking metadata and OTUs
# in this case, sample ids are in the "sequencing_id" column of the metadata

set_standardID(meta, "sequencing_id")
meta$standardID[1:5]

# Rarefaction on full dataset ----------
# run iNext rarefaction on all samples, all OTUs

# calculate incidence counts
all_incidence <- abund_to_incidence(samples = "All",
                                    otu_table = otus,
                                    OTU_name_column = "OTU_ID")

# calculate diversity estimates
message("Running iNext on all samples")
start_time <- Sys.time()

all_div <- iNEXT(all_incidence[[1]],
                 q= 0,
                 datatype="incidence_freq")
end_time <- Sys.time()

message("Finished")
print(end_time - start_time)

# plot with free scale (each facet gets it's own scale)
p_all <- ggiNEXT(all_div, type=1) + 
  labs(title = paste("Sampling Curve-", "All Samples"),
                                      y= "ESV Richness")

ggsave( "all_rarefaction_curves.pdf", p_all)

# Rarefaction on each habitat type ----------

# calculate incidence counts by habitat
habitat_incidence <- meta[ , .(incidence = abund_to_incidence(samples = standardID,
                                                     otu_table = otus,
                                                     OTU_name_column = "OTU_ID")),
                          by = habitat]

habitat_list        <- habitat_incidence$incidence
names(habitat_list) <- habitat_incidence$habitat

message("OTUs counted per habitat")
lapply(habitat_list, length)

saveRDS(habitat_list, "habitat_incidence_counts.rds")

# calculate diversity estimates
message("starting rarefaction by habitat")


# loop through habitats
for (i in seq_along(habitat_list)) {
  
# start time
start_time <- Sys.time()

# get habitat name
habitat_name <- names(habitat_list[i])
print(habitat_name)

# run diversity calculations
habitat_div <- iNEXT(habitat_list[[i]],
                 q=0,
                 datatype="incidence_freq")

# save diversity calculations 
saveRDS( habitat_div, 
        file = paste(habitat_name, "_habitat_div.rds", sep = "") )

# plot with free scale
p_habitat <- ggiNEXT(habitat_div,
                     type=1) +
                   labs(title = paste("Sampling Curve-", habitat_name),
                        y= "ESV Richness")

ggsave( filename = paste( habitat_name, "_habitat_rarefaction_curves.pdf", sep = ""),
        p_habitat )

end_time <- Sys.time()
print(end_time - start_time)

}

# Rarefaction on non-marine samples -----------

# calculate incidence counts for non-marine samples 
no_marine_incidence <- meta[habitat != "Marine", .(incidence = abund_to_incidence(samples = standardID,
                                                          otu_table = otus,
                                                          OTU_name_column = "OTU_ID")) ]
no_marine_div <-  iNEXT(no_marine_incidence$incidence[[1]],
                        q=0,
                        datatype="incidence_freq")

p_no_marine <- ggiNEXT(no_marine_div,
                       type=1) +
                      labs(title = "Sampling Curve - Non-Marine Samples", y= "ESV Richness")
p_no_marine


ggsave( "no_marine_rarefaction_curves.pdf", p_no_marine)

