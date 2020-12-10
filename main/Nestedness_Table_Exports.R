### EMPO Data Exports for Nestedness Analysis
# Read in raw data set, export separate subsets based on EMP ontology.
# These tables will be used to run separate nestedness analyses for each group in the EMP ontology.

library(data.table)
library(vegan)
# Read in raw otu table, drop unnecessary columns and transpose into long format
all_otus <- fread("../data/raw/sequence/waimea_16S_v1/abundance_table_100.shared",
                  drop = c("label","numOtus"), nrow = 100)

all_otus <- transpose(all_otus,
                      keep.names = "OTU_ID",
                      make.names = "Group")

sample_dat <- read.csv("../data/processed/cleaned/mm_16s_hiseqs_metadata_table.csv")

# subset to relevant samples
sample_dat <- sample_dat[sample_dat$project_name %in% c("Keck", "JoshMS") ,]
sample_dat <- sample_dat[sample_dat$sequencing_id %in% colnames(all_otus) ,]

# add metadata category for habitat and site combined (habisite) and all samples
# for habisite, consider beach water samples as riverine, since the question is about paired sites not habitat
sample_dat$habisite <- paste(sample_dat$habitat,
                             sample_dat$site_name,
                             sep = "-")
sample_dat[sample_dat$habisite == "Marine-Beach", "habisite"] <- "Riverine-Beach"

sample_dat$all <- "samples"
# sub_sum() sums read counts by samples, assuming sample ids are columns

sub_sum <- function(a_group, otu_table, sample_dat, group_column, key_column, otu_id_column){
  group_samples <- sample_dat[sample_dat[[group_column]] == a_group, key_column]
  group_sums    <- otu_table[ , rowSums(.SD), .SDcols = as.character(group_samples)]
  group_sums[is.na(group_sums)] <- 0
  names(group_sums) <- otu_table$OTU_ID
  return(group_sums)
}


# for each level of each category (trophic, habitat, habitat+site, site), write out two tables
# table 1: samples agregated at the EMPO3 level, raw counts
# table 2: samples aggregated at the EMPO3 level, rarefied to lowest EMPO3 sum

empo_agg_out <- list()

sums_out <- list()

for(a_category in c("all","habitat", "site_name", "habisite","trophic" )){
  
  ## within a category, aggregate EMPO3 at each level and write out table
  for (a_level in unique(sample_dat[[a_category]])) {
    
    # generate the name for this table
    table_name <- paste0( a_level, "_", a_category, "_empo_3")
    
    # subset sample data to a single level for a given category (i.e. data column)
    samples_subset <- sample_dat[sample_dat[[a_category]] == a_level , ]
    # subset the samples(i.e. columns) of the OTU table to match
    keep_samples <- c(samples_subset$sequencing_id, "OTU_ID")
    otus_subset    <- all_otus[ , ..keep_samples]
    
    # get a list of the unique EMPO 3 levels in this subset of samples
    empo_3_list <- as.character(unique(samples_subset$empo_3))
    
    # aggregate samples at each EMPO 3 level by summing read counts
    agg_otu_table <-
      sapply(empo_3_list,
             sub_sum,
             otu_table = otus_subset,
             sample_dat = samples_subset,
             group_column = "empo_3",
             key_column = "sequencing_id",
             otu_id_column = "OTU_ID")
    
    # Drop OTUS that have no reads and convert matrix to data.frame
    agg_otu_table <- agg_otu_table[rowSums(agg_otu_table[]) > 0 ,]
    # convert matrix to data.frame       
    agg_otu_table <- as.data.frame(agg_otu_table)
    # set column names to EMPO 3 names and add ASV names back in as "OTU_ID"
    colnames(agg_otu_table) <- empo_3_list
    
    ## Generate rarefied table
    # print out sums by EMPO 3 and save histogram to inform decisions about rarefaction depth
    empo3_sums <- colSums(agg_otu_table)
    
    print(table_name)
    print(empo3_sums)
    
    png(paste0(table_name, "_hist.png"))
    hist(empo3_sums,
         main = table_name)
    dev.off()
    
    # rarefy to lowest EMPO3 sum
    set.seed(2020)
    
    rar_agg_table <- rrarefy(t(as.matrix(agg_otu_table)),
                           sample = min(empo3_sums) )
    
    # add otu ids to rarefied aggregated table and store in output list
    rar_agg_table <- as.data.frame(t(rar_agg_table))
    rar_agg_table$OTU_ID <- row.names(agg_otu_table)
    empo_agg_out[[paste0("rar_", table_name)]]
    
    # add otu ids to raw aggregated table and store in output list
    agg_otu_table$OTU_ID <- row.names(agg_otu_table)
    empo_agg_out[[ table_name ]] <- agg_otu_table
    
    
  }
}
names(empo_agg_out)

for(a_table in names(empo_agg_out)){
  write.csv(empo_agg_out[[x]],
            file = paste0("nestedness_tables/",x,".csv"),
            row.names = F)
}
