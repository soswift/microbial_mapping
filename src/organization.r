## set up environment
	library(ape)
	library(fields)
	library(vegan)
	library(philr)
	library(Rcpp)

	source("filtering_functions.r") # takes a second because of compilation of c++ code

	# NOTE: for this to work, qiime1 and vsearch have to be working at the system level.
	# I used bioconda for this, so before running this script i run:
	# conda activate vsearch
	# conda activate qiime1

## define filepaths for inputs and outputs
	otutable_raw_fp <- "../../clustering/vsearch95_otutable.txt"
	metadata_raw_fp <- "../../metadata/Endophyte_metadata_cleaned_GIS.txt"
	taxonomy_raw_fp <- "../../ghosttree/tax_unite_fungi_gt/fungi_only_repset_tax_assignments.txt"
	ghosttree_raw_fp <- "../../ghosttree/ghost_tree_otuids.nwk"
	figures_folder_fp <- "./figures"
	workspace_folder_fp <- "./workspaces"
	plant_tree_fp <- "../../metadata/Supertree/supertree.nwk"

	# tar command to make tar.gz with above inputs for Anthony:
	# tarcmd <- paste("tar -cvzf inputs.tar.gz", otutable_raw_fp, metadata_raw_fp, taxonomy_raw_fp, ghosttree_raw_fp, plant_tree_fp)
	# system(tarcmd)

## make folders for figures and workspaces
	dir.create(figures_folder_fp)
	dir.create(workspace_folder_fp)

## Read in files
	# read in raw otutable
	otutable_raw <- t(read.table(otutable_raw_fp, header=T, 
		sep='\t', row.names=1))
	# read in metadata
	metadata_raw <- read.table(metadata_raw_fp, sep='\t', header=T, stringsAsFactors=F)
	# read in taxonomy, m
	taxonomy <- read.table(taxonomy_raw_fp, sep='\t', header=F, stringsAsFactors=F)

	# read in ghosttree
	ghosttree_raw <- read.tree(ghosttree_raw_fp)
	# read in plant supertree
	supertree <- read.tree(plant_tree_fp)
	# add header for taxonomy
	colnames(taxonomy) <- c("otu", "tax", "eval", "hit")

## Filter otutable taxonomy to only include OTUs that had matches in UNITE
	# function that subsets x to match the order and content of sampleIDs_to_keep
	# its name is a shortened "multiple which"
	# basiclly, mch(x, y) gives you the indices of x that correspond to their matches in y. example:
	#   x <- c(5,8,4,3,1,6,7,2)
	#   y <- c(1,2,3)
	#   # mch(x,y) should give us the indices of 1,2,3 in x, which are 5,8,4. x[c(5,8,4)] gives us 1,2,3.
	#   mch(x,y)
	#   x[mch(x,y)]
	# don't bother trying to understand it, it just works really well.
	mch <- function(unsorted, template){ sapply(X=template, FUN=function(x){which(x==unsorted)[1]}) }
	taxonomy <- taxonomy[mch(taxonomy$otu, rownames(otutable_raw)), ]
	# check:
	all(taxonomy$otu == rownames(otutable_raw))
	# logical vector, which rows in the otutable have a fungal match in UNITE?
	fungi_logical <- grepl("Fungi", taxonomy$tax)
	# drop non-fungi from otutable_raw and from taxonomy
	otutable_raw <- otutable_raw[fungi_logical, ]
	taxonomy <- taxonomy[fungi_logical, ]
	all(taxonomy$otu == rownames(otutable_raw))
	# still TRUE

## check if all unite otus made it into ghost tree
	all(rownames(otutable_raw) %in% ghosttree_raw$tip.label)
	# they did! If not, I'd have to drop them:
	# otutable_raw <- otutable_raw[ ! rownames(otutable_raw) %in% ghosttree_raw$tip.label, ]

## check that objects match up as expected
	# sanitize names of otutable_raw and metadata, so they can match up
	colnames(otutable_raw) <- make.names(colnames(otutable_raw))
	rownames(metadata_raw) <- make.names(gsub(metadata_raw$SampleID, pattern="\\.", replace=""))
	# check
	all(rownames(otutable_raw) %in% taxonomy[[1]])
	# drop old sampleid column from metadata to avoid confusion
	metadata_raw <- metadata_raw[colnames(metadata_raw) != "SampleID"]

## drop samples with no geographic information (i.e. no metadata)
	samps_with_geo <- rownames(metadata_raw)[ ! is.na(metadata_raw$Lon) ]
	otutable_raw <- otutable_raw[, colnames(otutable_raw) %in% samps_with_geo]

## drop samples with absurdly low depth (<= 150 seqs)
	otutable_raw <- otutable_raw[, colSums(otutable_raw) > 150 ]

## removing these samples will make some otus empty, so get rid of them
	otutable_raw <- otutable_raw[ rowSums(otutable_raw) > 0 ,]

## ensure otutable_raw and metadata_raw match each other
	metadata_raw <- metadata_raw[mch(rownames(metadata_raw), colnames(otutable_raw)), ]
	all(rownames(metadata_raw) == colnames(otutable_raw))

## make collectors curve figures
	col_curve_list_raw <- multiple_collectors_curves_parallel(otutable_raw, met="pd", maxd=3500, stepsize=50, phy=ghosttree_raw, ncores=10)
	pdf(paste(figures_folder_fp, "/binned_collectors_curves_pd.pdf", sep=""), useDingbats=F)
	plot_multiple_binned_collectors_curves(col_curve_list_raw, ylab="Phylogenetic Diversity of OTUs", n_bins=6, bin_d=1500)
	dev.off()

	col_curve_list_raw_sobs <- multiple_collectors_curves_parallel(otutable_raw, met="sobs", maxd=3500, stepsize=50, phy=NULL, ncores=10)
	pdf(paste(figures_folder_fp, "/binned_collectors_curves_sobs.pdf", sep=""), useDingbats=F)
	plot_multiple_binned_collectors_curves(col_curve_list_raw_sobs, ylab="Number of OTUs observed", n_bins=6, bin_d=1500)
	dev.off()
	
	# combined 2-up figure
	pdf(paste(figures_folder_fp, "/binned_collectors_curves_2up.pdf", sep=""), useDingbats=F)
	par(mfrow=c(2,1))
	plot_multiple_binned_collectors_curves(col_curve_list_raw, ylab="Phylogenetic Diversity of OTUs", n_bins=6, bin_d=1500)
	plot_multiple_binned_collectors_curves(col_curve_list_raw_sobs, ylab="Number of OTUs observed", n_bins=6, bin_d=1500)
	dev.off()


## Rarefy otutable to depth of 1500 sequences/sample
	# otutable_raw_r means otutable_raw_rarefied
	set.seed(12345)
	otutable_raw_r <- rarefy_otu_table_cpp(otutable_raw, depth=1500)


##
##
##    SAVE 1
##
##
## save workspace - enables script to be started from this point
	save(list=ls(), file=paste(workspace_folder_fp, "/workspace_organization_backup_1.rdata", sep="/"))

## plot of richness by frequency
	pdf(paste(figures_folder_fp, "species_observed_histograms.pdf", sep="/"), useDingbats=F)
	# otutable_raw_r
	sobs_otutable_raw_r <- colSums(otutable_raw_r > 0)
	hist(sobs_otutable_raw_r, breaks=max(sobs_otutable_raw_r), col="pink", xlab="Number of zOTUs observed", ylab="Frequency (number of samples)", main="otutable_raw_r")
	dev.off()

## make unifrac beta-diversity matrix - REQUIRES QIIME 1 TO BE INSTALLED AND ACTIVE
	# for UniFrac, first write rarefied otu table to biom
	write_biom_jld(otutable_raw_r, paste(workspace_folder_fp, "otutable_raw_r.biom", sep="/"))
	# do UniFrac with QIIME, read in matrix
	unifrac_command <- paste(
		"beta_diversity.py -i",
		paste(workspace_folder_fp, "otutable_raw_r.biom", sep="/"),
		"-t",
		ghosttree_raw_fp,
		"-o",
		paste(workspace_folder_fp, "otutable_raw_r_unifrac", sep="/")
	)
	system(unifrac_command)
	unifrac_file2get <- paste(workspace_folder_fp, "otutable_raw_r_unifrac", "weighted_unifrac_otutable_raw_r.txt", sep="/")
	unifrac_beta_distmat <- read.table(unifrac_file2get, sep='\t', row.names=1, header=T)
	
## ensure all data structures match each other
	# and also get rid of non-dicots (only a few left in by accident)
	# and also a sample where "Asteraceae" (not a genus) was recorded as the plant genus
	sampleIDs_to_keep <- colnames(otutable_raw_r)[
		colnames(otutable_raw_r) %in% rownames(metadata_raw)[! is.na(metadata_raw$PlantGenus) ] & 
		colnames(otutable_raw_r) %in% rownames(metadata_raw)[! metadata_raw$PlantGenus %in% c("Peperomia", "Smilax", "Asteraceae") ] & 
		colnames(otutable_raw_r) %in% rownames(unifrac_beta_distmat)
	]
	# function that subsets x to match the order and content of sampleIDs_to_keep:
	mch <- function(z, x=sampleIDs_to_keep){ sapply(X=x, FUN=function(x){which(x==z)[1]}) }
	# sort everything to match sampleIDs_to_keep
	otutable_raw_r <- otutable_raw_r[,mch(colnames(otutable_raw_r))]
	otutable_raw <- otutable_raw[,mch(colnames(otutable_raw))]
	metadata_raw <- metadata_raw[mch(rownames(metadata_raw)),]
	unifrac_beta_distmat <- unifrac_beta_distmat[ mch(rownames(unifrac_beta_distmat)), mch(colnames(unifrac_beta_distmat)) ]

	# two otutables - one rarefied, one not.
	# renaming to retain my sanity.
	otutable_rare <- otutable_raw_r
	otutable_raw  <- otutable_raw # no change here
	rm(otutable_raw_r)


	# drop empty species from both otutables
	otutable_rare <- otutable_rare[rowSums(otutable_rare) > 0,]
	otutable_raw <-  otutable_raw[rowSums(otutable_raw) > 0,]

	# separate taxonomies for rare/raw otutables
	taxonomy_rare <- taxonomy[mch(taxonomy$otu, rownames(otutable_rare)),]
	taxonomy_raw <- taxonomy[mch(taxonomy$otu, rownames(otutable_raw)),]
	rm(taxonomy)

	# check sorting
	all( colnames(otutable_rare) == rownames(metadata_raw))
	all( colnames(otutable_rare) == rownames(unifrac_beta_distmat))
	all( colnames(otutable_rare) == colnames(unifrac_beta_distmat))
	all( rownames(otutable_rare) == taxonomy_rare$otu)
	all( rownames(otutable_raw) == taxonomy_raw$otu)
	all( colnames(otutable_rare) == colnames(otutable_raw))

## make host genus cophenetic matrix
	supertree_genera <- as.vector(sapply(X=supertree$tip.label, FUN=function(x) strsplit(x, split="_")[[1]][1] ))
	supertree$tip.label <- supertree_genera

	# are all genera in metadata_raw are in the tree?
	all( unique(metadata_raw$PlantGenus) %in% supertree_genera )
		# TRUE! All good

	# make a distance matrix containing pairwise average phylogenetic distances between each pair of genera
	# first, make a tree that only has one tip per genus
	tip_inds_2keep <- sapply(X=unique(metadata_raw$PlantGenus), FUN=function(x){which(supertree$tip.label == x)[1]})
	tips_tf <- rep(F, length(supertree$tip.label)); tips_tf[tip_inds_2keep] <- TRUE
	tip_inds_2drop <- which(!tips_tf)
	supertree <- drop.tip(supertree, tip_inds_2drop)
	rm(tip_inds_2keep, tips_tf, tip_inds_2drop)

	mini_cophenetic <- cophenetic.phylo(supertree)

	# remove NA rows and columns from mini_cophenetic
	mini_cophenetic <- mini_cophenetic[!is.na(rownames(mini_cophenetic)), !is.na(colnames(mini_cophenetic))]

	# expand mini coph to a full sized version
	plant_distmat <- matrix(nrow(unifrac_beta_distmat), ncol(unifrac_beta_distmat), data=NA)
	rownames(plant_distmat) <- colnames(plant_distmat) <- rownames(metadata_raw)
	for(i in 1:ncol(mini_cophenetic)){for(j in 1:nrow(mini_cophenetic)){
		val_ij <- mini_cophenetic[i,j]
		# is i,j in the LT?
		if(!is.na(val_ij)){
			genus_i <- colnames(mini_cophenetic)[i]
			genus_j <- colnames(mini_cophenetic)[j]
			plant_distmat[ metadata_raw$PlantGenus %in% genus_i,	metadata_raw$PlantGenus %in% genus_j ] <- val_ij
			plant_distmat[ metadata_raw$PlantGenus %in% genus_j, metadata_raw$PlantGenus %in% genus_i ] <- val_ij
		}
	}}

	# check:
	all(colnames(plant_distmat) == colnames(unifrac_beta_distmat))
	all(rownames(plant_distmat) == rownames(unifrac_beta_distmat))

##
##
##    SAVE 2
##
##
## save workspace again
	save(list=ls(), file=paste(workspace_folder_fp, "/workspace_organization_backup_2.rdata", sep="/"))

## make bray-curtis distance matrices
	bc_beta_distmat <- as.matrix(vegdist(t(otutable_rare), method="bray"))

	# quick check to make sure it's sorted right:
	all(rownames(bc_beta_distmat) == rownames(unifrac_beta_distmat))
	# true = all good!

## make philr beta diversity matrix, because compositions
	# doesn't require rarefied table
	# renaming all nodes of ghosttree to make sense:
	ghosttree_raw$node.label <- paste0("n", 1:length(ghosttree_raw$node.label))
	# following philr vignette: browseVignettes("philr")

	# function to make philr distmat
		# requires no zeroes (add 1?)
		# requires otutable with species as columns
	make_philr_dm <- function(tab, tre=ghosttree_raw){
		contrasts <- philr(
			df=tab,
			tree=keep.tip(tre, colnames(tab)),
			part.weights='enorm.x.gm.counts', 
			ilr.weights='blw.sqrt'
		)
		dm <- as.matrix(dist(contrasts, method="euclidean"))
		return( dm/max(dm) )
	}

	philr_beta_distmat_raw  <- make_philr_dm(t(otutable_raw + 0.25))
	philr_beta_distmat_rare <- make_philr_dm(t(otutable_rare+ 0.25))


	plot(as.dist(philr_beta_distmat_raw) ~ as.dist(unifrac_beta_distmat), pch=".")
	plot(as.dist(philr_beta_distmat_rare) ~ as.dist(unifrac_beta_distmat), pch=".")
	plot(as.dist(philr_beta_distmat_rare) ~ as.dist(philr_beta_distmat_raw), pch=".")


## make geographic distance matrix out of lat/longs	
	geo_distmat <- distcalc(metadata_raw$Lat, metadata_raw$Lon, rownames(metadata_raw))
	# check:
	all(colnames(geo_distmat) == colnames(unifrac_beta_distmat))
	all(rownames(geo_distmat) == rownames(unifrac_beta_distmat))

## Aspect is in degrees, so euclidean distance (gdm)  is not appropriate. That's because 358 degrees is close to 0 degres.
	# to solve this, I'm using a custom distance metric (arc distance) and building a distance matrix, which
	# will be fed into GDM.
	degdist <- function(d1, d2){
		a <- d1 - d2
		a <- abs( ((a + 180) %% 360) - 180 )
		return(a)
	}
	aspect_distmat <- sapply(X=metadata_raw$Aspect, FUN=degdist, d2=metadata_raw$Aspect)

## add juliandate column to metadata_raw
	metadata_raw$JulianDate <- as.integer(as.Date(metadata_raw$CollectionDate, format="%m/%d/%y"))

## make "run" distance matrix. 1 if samples are from different runs, 0 if from same run.
	identity_dm <- function(x){
		output <- matrix(1, length(x), length(x))
		for(i in 1:length(x)){for(j in 1:length(x)){
			if(x[i] == x[j]){
				output[i,j] <- 0
			}
		}}
		return(output)
	}
	run_distmat <- identity_dm(metadata_raw$Run)

## save "final" objects in workspace
	metadata <- metadata_raw
	otutable <- otutable_rare
	taxonomy <- taxonomy_rare
	save(file=paste(workspace_folder_fp, "/organized_data.rdata", sep="/"), list=c(
		"metadata", 
		"otutable",
		"otutable_raw",
		"taxonomy",
		"taxonomy_raw",
		"plant_distmat", 
		"unifrac_beta_distmat",
		"bc_beta_distmat", 
		"philr_beta_distmat_raw",
		"philr_beta_distmat_rare",
		"geo_distmat",
		"aspect_distmat",
		"run_distmat",
		"supertree"
	))


## plot histograms to visualize distribution and range of each predictor variable
	pdf(paste(figures_folder_fp, "variable_distributions.pdf", sep="/"), useDingbats=F)
	
	df <- data.frame(
		NDVI=metadata$NDVI,
		RelativeHumidity=metadata$RelativeHumidity,
		Elevation=metadata$Elevation,
		WetCanopyEvap=metadata$WetCanopyEvap,
		Temperature=metadata$Temperature,
		SolarRadiation=metadata$SolarRadiation,
		LeafAreaIndex=metadata$LeafAreaIndex,
		Rainfall=metadata$Rainfall,
		Aspect=metadata$Aspect,
		Slope=metadata$Slope,
		Evapotranspiration=metadata$Evapotranspiration,
		CloudFrequency=metadata$CloudFrequency
	)
	df$Elevation <- df$Elevation - min(df$Elevation)

	kern <- "gaussian"
	bgcolor="gray"
	fillcolor="black"
	linecolor="white"
	labelcolor="black"
	labelbg="gray"
	overlap=0.3
	labelpad=0.40



	# reverse df order so entered columns order is top2bottom
	df <- df[rev(1:ncol(df))]

	# remove non-numeric cols from df
	df <- df[, apply(X=df, FUN=is.numeric, MAR=2)]
	# figure out plot dimensions
	n_plots <- ncol(df)
	heights <- (0:(n_plots-1)) * (1 - overlap)
	ymax <- max(heights) + 1

	# make plot
	par(mar=c(0, 0, 0, 0))
	plot(NULL, ylab="", xlab="", main="", bty="n", xaxt="n", yaxt="n", 
		xlim=c((-1 * labelpad),1), ylim=c(0, ymax+0.1), asp=0.09 )
	rect((-1 * labelpad), 0, 1, ymax+0.1, col=bgcolor, border=NA)

	rescale01 <- function(x){x <- x - min(x); return( x * (1/max(x)) )}
	mins <- maxs <- rep(0, ncol(df))
	# add dists for each var, top to bottom
	for(i in n_plots:1){
		dens_i <- density(df[,i][!is.na(df[,i])], kernel=kern)
		mins[i] <- min(df[,i]); maxs[i] <- max(df[,i])
		dens_i$x <- rescale01(dens_i$x)
		dens_i$y <- rescale01(dens_i$y)


		polygon(x=dens_i$x, y=dens_i$y + heights[i], col=fillcolor, border=NA)
		lines(x=dens_i$x, y=dens_i$y + heights[i], col=linecolor, lwd=2)
		lines(x=c((-1*labelpad), 1), y=rep(heights[i], 2), col=linecolor, lwd=2)
	}
	lines(x=c(0,0), y=c(0, ymax), col=linecolor, lwd=2)
	text(x=-1 * labelpad + 0.01, y=heights+(1-overlap-0.1), labels=colnames(df), col=labelcolor, adj=c(0,1) )
	
	minsp <- signif( mins, 2)
	maxsp <- signif( maxs, 2)
	text(x=0 + 0.01, y=heights+(1-overlap-0.1), labels=minsp, col=labelcolor, adj=c(0,1) )
	text(x=1 - 0.01, y=heights+(1-overlap-0.1), labels=maxsp, col=labelcolor, adj=c(1,1) )

	dev.off()


## write out mini metadata for supplemental table
	# just a mini metadata that only has sampleid, host genus, island, lat, lng.
	sup_metadata <- data.frame(
		SampleName=rownames(metadata), 
		HostGenus=metadata$PlantGenus,
		Island=metadata$Island,
		Latitude=metadata$Lat,
		Longitude=metadata$Lon
	)	
	write.table(file=paste(figures_folder_fp, "supplemental_metadata.txt", sep="/"), x=sup_metadata, sep='\t', row.names=F, quote=F)

## unifrac vs bray comparison stuff
	unifrac_beta_distmat <- unifrac_beta_distmat/max(unifrac_beta_distmat)
	pdf("figures/UnifVsBray.pdf", useDingbats = FALSE)
	plot(as.dist(unifrac_beta_distmat) ~ as.dist(bc_beta_distmat), pch=20, cex=0.10)
	hist(as.dist(unifrac_beta_distmat), main="UniFrac")
	hist(as.dist(bc_beta_distmat), main="Bray-Curtis")
	dev.off()



 


