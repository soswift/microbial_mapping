### EMPO Data Exports for Nestedness Analysis
# Read in raw data set, export separate subsets based on EMP ontology.
# These tables will be used to run separate nestedness analyses for each group in the EMP ontology.

library(data.table)
# Read in raw otu table, drop unnecessary columns and transpose into long format
all_otus <- fread("../data/raw/sequence/waimea_16S_v1/abundance_table_100.shared",
                  drop = c("label","numOtus"))

all_otus <- transpose(all_otus,
                      keep.names = "OTU_ID",
                      make.names = "Group")

sample_dat <- read.csv("../data/processed/cleaned/mm_16s_hiseqs_metadata_table.csv")

# subset to relevant samples
sample_dat <- sample_dat[sample_dat$project_name %in% c("Keck", "JoshMS") ,]
sample_dat <- sample_dat[sample_dat$sequencing_id %in% colnames(all_otus) ,]

# sub_sum() sums read counts by samples, assuming sample ids are columns

sub_sum <- function(a_group, otu_table, sample_dat, group_column, key_column, otu_id_column){
  group_samples <- sample_dat[sample_dat[[group_column]] == a_group, key_column]
  group_sums    <- otu_table[ , rowSums(.SD), .SDcols = as.character(group_samples)]
  group_sums[is.na(group_sums)] <- 0
  names(group_sums) <- otu_table$OTU_ID
  return(group_sums)
}


# for each level of each category (trophic, habitat, site), write a table, agregated at the EMPO3 level

empo_agg_out <- list()

for(a_category in c("trophic", "habitat", "site")){
  
  # within a category, aggregate EMPO3 at each level
  for (a_level in unique(sample_dat[[a_category]])) {
    
    # subset sample data to a single level for a given column
    samples_subset <- sample_dat[sample_dat[[a_category]] == a_level , ]
    # subset the columns of the OTU table to match
    keep_samples <- c(samples_subset$sequencing_id, "OTU_ID")
    otus_subset    <- all_otus[ , ..keep_samples]
    
    # get a list of the unique EMPO 3 levels in this subset of samples
    empo_3_list <- as.character(unique(samples_subset$empo_3))
    
    # aggregate samples at each EMPO 3 level by summing read counts
    sums_otu_table <-
      sapply(empo_3_list,
             sub_sum,
             otu_table = otus_subset,
             sample_dat = samples_subset,
             group_column = "empo_3",
             key_column = "sequencing_id",
             otu_id_column = "OTU_ID")
    
    # Drop OTUS that have no reads and convert matrix to data.frame
    sums_otu_table <- sums_otu_table[rowSums(sums_otu_table[]) > 0 ,]
    # convert matrix to data.frame       
    sums_otu_table <- as.data.frame(sums_otu_table)
    # set column names to EMPO 3 names and add ASV names back in as "OTU_ID"
    colnames(sums_otu_table) <- empo_3_list
    sums_otu_table$OTU_ID    <- row.names(sums_otu_table)
    
    empo_agg_out[[paste0(a_level, "_", a_category, "_empo_3")]] <- sums_otu_table
  }
}

