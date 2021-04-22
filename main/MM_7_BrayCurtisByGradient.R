# Code by Anthony, modified by Sean
# Compare bray curtis distance between river and terrestrial sites across the elevation gradient
require(vegan)
require(geosphere)
library(data.table)

tm <- proc.time()

# source Jack's scripts
source("../src/filtering_functions.r")

# identify data files

otu_file         <- "../data/processed/cleaned/raw_mm_16s_hiseqsabundance_table.csv"
sample_data_file <- "../data/processed/cleaned/mm_16s_hiseqs_metadata_table.csv"

# read in data, samples as columns and ASVs as rows
# 'Group' in the OTU table matches 'sequencing_id' in the sample data
otu_table     <- transpose( fread(otu_file), make.names = "Group", keep.names = "ASV_name")
sample_dat    <- read.csv(sample_data_file, row.names = "sequencing_id")

# for this analysis, we will treat Beach samples as "Riverine" instead of "Marine"
sample_dat[sample_dat$site_name == "Beach" & sample_dat$habitat == "Marine", "habitat"] <- "Riverine"

# subset the data to only include the paired sites i.e. riverine and terrestrial sample
sample_dat <- sample_dat[sample_dat$habitat %in% c("Riverine","Terrestrial"), ]

# create single vector of unique habitat/site combinations
sample_dat$habisite <- paste(sample_dat$habitat, sample_dat$site_name, sep = "_")

# sub_sum() quickly calculates sums for ASVs (rows) for sample(columns) that occur in a habitat-site
sub_sum <- function(a_habisite, otu_table, sample_dat){
              habisite_samples <- row.names( sample_dat[sample_dat$habisite == a_habisite , ])
              habisite_sums    <- otu_table[ , rowSums(.SD), .SDcols = habisite_samples]
              habisite_sums[is.na(habisite_sums)] <- 0
              return(habisite_sums)
}

# calculate ASV sums for each habitat-site
sums_otu_table <-
  sapply(unique(sample_dat$habisite),
         sub_sum,
         otu_table = otu_table,
         sample_dat = sample_dat)

# reassign ASV_names in case we need them later
row.names(sums_otu_table) <- otu_table$ASV_name

# rarefy to the minimum sampling depth of a habitat-site for fair comparison
rar_otu_table <-rarefy_otu_table_cpp(sums_otu_table,
                                     depth = min(colSums(sums_otu_table)))

# finally, we can generate bray-curtis distance matrix for each habisite
bc_dist <-vegdist(t(rar_otu_table))

# now, we'll compare bray curtis distance between habitats (terrestrial/riverine) for each site
# extract bray curtis distances by site
unique_sites <- unique(sample_dat$site_name)


# for each site, pull out the maximum bray curtis distance
max_dists <- sapply(
  unique_sites,
  FUN = function(x) {
    # get names of habisites for each site
    site_names <- grep(x, labels(bc_dist), value = T)
    
    # pull out maximum bray curtis distances between the habitat-sites
    max_dist <- max(as.matrix(bc_dist)[site_names, site_names])
    names(max_dist) <- x
    return(max_dist)
  }
)

# organize bray-curtis distances with site information that can be used for plotting
site_summary <-setDT(sample_dat)[ , .N, by = .(site_name,
                                               rainfall,
                                               elevation,
                                               latitude_deg,
                                               longitude_deg)]

site_summary <- merge(site_summary,
                      data.table( site_name = names(max_dists),
                                  bc_dist = max_dists))

# identify point on Waimea shore from which to measure distance along the gradient
shore <- c(-158.063539,21.640629)

# calculate distance of sites from shore
site_summary$shore_dist <-  geosphere::distGeo(shore, as.matrix(site_summary[, .(longitude_deg,latitude_deg)] ) )

# check variation in bray-curtis distance across the sampling gradient
ggplot(site_summary, aes(x = site_summary[["shore_dist"]], y = site_summary[["bc_dist"]]))+
          geom_point()+
          labs(y = "Bray-Curtis Dissimilarity Between Habitats", x= "Distance from Shore")+
          theme_minimal()
ggsave("outputs/bray_curtis_gradient/bray_curtis_by_shore_distance.pdf")
proc.time() - tm

