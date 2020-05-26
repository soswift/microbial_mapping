

# read and clean 16S ------------------------------------------------------
# setwd("~/Desktop/R/CMAIKI_clean_and_query/bacterial_16S/")
# clean 16S pipeline outputs for R analyses

clean_16S_tables <- function(abundance_file = NULL,
                             taxonomy_file  = NULL,
                             metadata_file  = NULL,
                             description    = NULL,
                             output_dir = "./",
                             id_column = NULL,
                             cull = list(min.num = NULL, min.abund = NULL, min.single.abund = NULL)){
  
  # Requirements
  require(dplyr)
  require(tidyr)
  require(data.table)
  
  ########## START
  
  #### *** Test Case COMMENT OUT UNLESS TESTING ***

   abundance_file = "../data/raw/abundance_table_97.shared"
   taxonomy_file  = "../data/raw/annotations_97.taxonomy"
   metadata_file  = "../data/processed/interim/combined_prep_sample_metadata_16S.csv"

   description = "MM_16S_hiseqs"
   output_dir = "../data/processed/cleaned"
   id_column = "sequencing_id"
   
   cull = list(min.num = 5, min.abund = 500, min.single.abund = 10 )
  #### ** TEST END ***
   
  # check input files
  if(!file.exists(abundance_file)){
    warning("Abundance file doesn't exist")
  }
  if(!file.exists(taxonomy_file)){
    warning("Taxonomy file doesn't exist")
  }
  if(!file.exists(metadata_file)){
    warning("Metadata file doesn't exist")
  }
  
  
  ## Read in tables
  message("reading abundance")
  abund <- fread( abundance_file, header = T, sep = "\t")
  
  message("reading taxonomy")
  tax   <- fread( taxonomy_file, header = T, sep = "\t")

  message("reading metadata")
  fullmeta  <- fread( metadata_file, header = T,
                      sep = "," )
  ## Clean tables
  # Clean taxonomy table:
  #Split taxonomic annotations (col = Taxonomy, generates warning because of trailing ";")
  tax <- tax %>% 
            separate( Taxonomy,
                      c( "kingdom", "phylum", "class", "order","family","genus" ),
                      sep = ";")
  
  # Clean abundance table: Cut down id name and delete numOtus
  abund$Group <- sub(".*(1\\d{5}).*","\\1", abund$Group, perl = T)
  abund       <- abund[,!"numOtus"]
  
  

  
  ## Check if joins work:
  # Check abundance and taxonomy have the same OTUs
  
  if( !all( colnames(abund)[-c(1,2)] %in% tax$OTU)){ # test
          warning("*** WARNING Abundance and taxonomy have different OTUs")   # no
  } 
  else {
    message("OKAY Abundance and taxonomy have the same OTUs")                 # yes
    }         
  
  ## Culling abundance file to manageable size

    # define cull.otu
  
    cull.otu = function(abund.dt, min.num = 1, min.abund = 1, min.single.abund) {
    #Inputs:  relabund.df = dataframe containing ONLY relative abundance data, no metadata or other info. Samples in rows and OTUs in columns.
    
    #min.num = the minimum number of samples an OTU needs to be present in to not be culled.
    
    #min.abund = the minimum relative abundance an OTU needs to have in (the min.num) samples to not be culled.
    
    #min.single.abund = the minimum relative abundance an OTU needs to have in a SINGLE sample to not be culled.
      
      ## Test case
      
      
      
      min.abund = 2
      min.num = 1
      min.single.abund = 1
    
      
  abund.dt <- data.table(x = 1:4, y = 2:4, z = 1)
 
  names_vec <- names(abund.dt)  
  names_vec
  abund.dt
  #for (otu in names_vec) {
      otu = names_vec[1]
       if(abund.dt[get(otu) < min.abund, length(get(otu))] < min.num | abund.dt[get(otu) < min.single.abund, length(get(otu))] == 0 ) {
        abund.dt[ , get(otu):=NULL]
   #     }
  }
  abund.dt
  }
  
  
  #If cull is a list, cut down the abundance, and taxonomy files based on listed values
   if(is.list(cull) & !all(is.null(cull))){
   
    abund <- cull.otu(abund, min.num = cull[1], min.abund = cull[2], min.single.abund = cull[3])
    tax   <- tax[OTU %in% colnames(abund)]
   
   } 
    else{
    message("Skipping cull step")
  }

  
  # Clean metadata table: Subset to include only samples that appear in the abundance table
  meta        <- fullmeta[fullmeta[[id_column]] %in% abund[["Group"]]]
  
  
  not_in_abund <-  fullmeta[!(fullmeta[[id_column]] %in% abund[["Group"]])]
  
  if ( abund[!(Group %in% fullmeta[[id_column]]), Group] ){
    warning("Some abundance samples are missing from the metadata")
  }
  
  
  
  # Make all tables into a list
  result <- list( abund, tax, meta, fullmeta, not_in_abund)
  names(result) <- c("abundance", "taxonomy", "metadata","full_run_metadata","failed_samples")
  
  # Write out all tables as .csv and return a "result" list of tables
  lapply(1:length(result),
         function(x) fwrite(result[[x]],
                               file = paste(output_dir,"/",description, "_", names(result[x]),"_", "table.csv", sep = "")))
  # write message
  message(paste("tables written out as csv files starting with ", description, "...", sep = ""))
  return(result)
  ########## END

  
  
  }



# summarize 16S sequencing success --------------------------------------------------
pass_fail <- function(clean_tables_list = NULL, id_column = NULL){

          
            # shorten variable
            run   <- clean_tables_list 
          
            # pull out samples that are in the abundance table, assume these sequence well
            good_samples <- run$full_run_metadata[ run$full_run_metadata[, id_column] %in% run$abundance$Group,]
            
            # pull out samples that are not in the abundance table, assume these did not sequence well
            bad_samples  <- run$full_run_metadata[!(run$full_run_metadata[, id_column] %in% run$abundance$Group)]
            
            result<-list(good_samples,bad_samples)
            names(result) <- c("good_samples","bad_samples")
            message(paste("for ",
                    deparse(substitute(clean_tables_list)), # run name, ideally
                    ", ",
                    nrow(good_samples), # number of good samples
                    " samples sequenced well and ",
                    nrow(bad_samples), " did not", sep = "")) # number of bad samples
          return(result)
          ########### END
}


