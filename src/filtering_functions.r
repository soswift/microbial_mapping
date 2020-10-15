## function to rarefy one sample
	rarefy_sample <- function(abunds, names=NULL, depth){
		# this does rarefaction by indexing observations to species
		names_alias <- 1:length(abunds)
		obs_alias <- factor(unlist(apply(X=data.frame(names_alias, abunds), MAR=1, FUN=function(x){rep(x[1], x[2])} )), levels=names_alias)
		table_rare <- table( sample(obs_alias, size=depth, replace = FALSE) )
		# because obs_alias is a factor, table_rare will ALWAYS have colnames == names_alias
		return(as.vector(table_rare))
	}

## rarefy OTU table - wrapper for rarefy_sample
	rarefy_otu_table <- function(otu_table, depth, seed=12345){
		set.seed(seed)
		# remove columns of OTU table with insufficient sample depth
		otu_table <- otu_table[ , colSums(otu_table) >= depth]

		# rarefy remaining columns
		output_table <- apply(X=otu_table, MAR=2, FUN=rarefy_sample, names=row.names(otu_table), depth=depth)
		rownames(output_table) <- rownames(otu_table)
		return(output_table)
	}

# Above function re-written in c++ to go faster because why wait?
	# x: vector of species counts (i.e. column from site-by-species matrix)
	# d: depth to which you'd like to rarefy
	# note that this will just take d down to sum(x) if d > x.
	library(Rcpp)
	cppFunction('IntegerVector rarefy_sample_cpp(const IntegerVector x, const int d) {
		// x is the "sample", i.e. column from otutable
		// need to expand x to be an array of size sum(x)
		int xsum = sum(x);
		IntegerVector xexp(xsum);
		int nspec = x.size();
		int pos = 0;
		int ind = 0;
		// this next bit "unpacks" x. Its EXACTLY what qiime1 does to rarefy.
		// for each position in x
		for(int i = 0; i < nspec; i++){
			int si = x[i];
			// for each observation at x[i]
			for(int j = 0; j < si; j++){
				xexp[pos] = ind;
				pos++;
			}
			ind++;
		}
		// shuffle xexp
		std::random_shuffle(xexp.begin(), xexp.end());
		// build output
		IntegerVector output(x.size());
		for(int i = 0; i < d; i++){
			// prevent errors by making sure d is reasonable!
			if(i < xsum){
				output[xexp[i]]++;
			}else{
				break;
			}
		}
		return output;
	}')



## rarefy OTU table C++ wrapper (in R, but runs cpp function)
	rarefy_otu_table_cpp <- function(otu_table, depth, seed=12345, rarefun=rarefy_sample_cpp){
		set.seed(seed)
		# make sure otu_table is a df, otherwise this won't work!
		if(is.data.frame(otu_table)){
			# ok woops I wrote this backwards lol
		}else if(is.matrix(otu_table)){
			otu_table <- as.data.frame(otu_table)
		}else{
			stop("otu_table needs to be a data.frame object.")
		}
		# todo: add a check that everything is positive integers.
		# remove columns of OTU table with insufficient sample depth
		otu_table <- otu_table[ , colSums(otu_table) >= depth]
		# rarefy each column using cpp function
		# this works because data.frames in R are secretly lists of columns
		# parallelizing this resulted in ignoring seed for some reason, so whatever.
		# would need to write a wrapper that takes a daughter seed and I'm too lazy.
		output <- simplify2array(lapply(X=as.list(otu_table), FUN=rarefun, d=depth))
		rownames(output) <- rownames(otu_table)
		return(output)
	}


## phylodiversity function
	one_samp_pd <- function(counts, names, tree){ 
		require(ape)
		# which tree tips are in names where counts > 0
		ids_overzero <- names[counts > 0]
		tips2keep <- tree$tip.label[ tree$tip.label %in% ids_overzero ]
		tips2drop <- tree$tip.label[ ! tree$tip.label %in% tips2keep ]
		# remove tips not in ids_overzero
		tree <- drop.tip(tree, tips2drop)
		# return pd of pruned tree
		return(sum(tree$edge.length))
	}

## rarefaction curve for one sample
	# abunds - column vector of species abundances
	# stepsize - how often to report diversity (number of seqs sampled)
	# met - alpha div metric. either "sobs" or "pd"
	# phy - needed if met="pd". a phylogenetic tree 
	# ids - column vector of species names/ids. only needed if met="pd"
	# maxd - max depth to construct curve (will stop making curve at maxd seqs sampled)
	alpha_rare_curve_1sample <- function(x, ids=NULL, stepsize=50, met="sobs", phy=NULL, maxd=NULL, seed=NULL){
		if(!is.null(seed)){ set.seed(seed) }
		abunds <- x
		output_adiv <- NULL # just going to concatenate to this
		output_step <- NULL 
		resampled_abunds <- rep(0, length(abunds))
		# take care of null value for maxd
		if(is.null(maxd)){
			maxd <- sum(abunds)
		}else if(maxd > sum(abunds)){
			maxd <- sum(abunds)
		}

		# calculate number of steps
		nsteps <- floor(maxd / stepsize)

		# initialize output data frame
		output <- data.frame(alpha=rep(0, nsteps), step=rep(0, nsteps))

		# x is going to be re-created as x2 via resampling, so first initialize x2 as all zeroes
		x2 <- x-x

		# for each step, take a sample of depth stepsize, then subtract that sample out from x
		# so that in the next step, those obs are still gone
		for(s in 1:nsteps){
			# deltax is just the difference between x[s-1] and x[s]
			deltax <- rarefy_sample_cpp(x, d=stepsize)
			x2 <- x2 + deltax
			x <- x - deltax
			output$step[s] <- s * stepsize

			if(met=="sobs"){
				output$alpha[s] <- sum(x2 > 0)
			}else if(met=="pd"){
				output$alpha[s] <- one_samp_pd(x2, ids, phy)
			}
		}
		return(output)
	}


##  make distance matrix out of lat/longs	
	library(fields)
	distcalc <- function(lat, lng, sampIDs){
		require(fields)
		longlats <- data.frame(lng, lat)
		rownames(longlats) <- sampIDs
		input <- as.matrix(longlats)
		distmat=rdist.earth(input, miles=F, R=NULL) 
		return(distmat)
	}

## faith's phylogenetic diversity
	pd <- function(tree, otutable){
		tips_in_table <- rownames(otutable)[ rowSums(otutable) > 0 ]
		unused_tips <- tree$tip.label[ ! tree$tip.label %in% tips_in_table ]
		tree <- drop.tip(tree, unused_tips)
		return(sum(tree$edge.length))

	}

## convert lat long (ddg) coordinates to utm
	LongLatToUTM<-function(lon_ddg, lat_ddg, zone=4){
		x <- lon_ddg; y <- lat_ddg
		xy <- data.frame(ID = 1:length(x), X = x, Y = y)
		coordinates(xy) <- c("X", "Y")
		proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
		res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
		return(as.data.frame(res))
}


## write data frame to biom file - requires BIOM installed at system lvl
	write_biom_jld <- function(df, outname){
		df <- data.frame("#SampleID"=rownames(df), df)
		somenumbers <- paste(sample(0:9, 12, replace=T), collapse="")
		tempfilename <- paste("temp_tsv_", somenumbers, "_", Sys.getpid(), ".txt.tmp", sep="")
		write.table(file=tempfilename, x=df, sep='\t', quote=F, row.names=F)
		command <- paste("biom convert -i", tempfilename, "--to-json -o", outname)
		system(command)
		system(paste("rm", tempfilename))
	}

## read data frame from biom file - requires BIOM installed at system lvl
	read_biom_jld <- function(biom_fp){
		somenumbers <- paste(sample(0:9, 12, replace=T), collapse="")
		tempfilename <- paste("temp_tsv_", somenumbers, "_", Sys.getpid(), ".txt.tmp", sep="")
		command <- paste("biom convert -i", biom_fp, "--to-tsv -o", tempfilename)
		system(command)
		df <- read.table(tempfilename, skip=1, comment.char="", header=T, sep="\t", stringsAsFactors=FALSE)
		rownames(df) <- df[,1]
		df <- df[, 2:ncol(df)]
		system(paste("rm", tempfilename))
		return(df)
	}


## col2pastel
	# a simple function that pastelizes a color
	col2pastel <- function(cols, pwr=0.3){
		apply(X=rgb2hsv(col2rgb(cols)), MAR=2, FUN=function(x){hsv(x[1], x[2]*pwr, x[3])})
	}
	# example:
	# par(mfrow=c(3,1))
	# raw_cols <- c("red", "blue","green")
	# barplot(1:3, col=raw_cols)
	# barplot(1:3, col=col2pastel(raw_cols, pwr=0.6))
	# barplot(1:3, col=col2pastel(raw_cols, pwr=0.3))

## apply2jagged
	# this function applies another function ACROSS values in a jagged array, omitting NAs
	# for example, consider the jagged list
	# foo <- list(a=1:4, b=1:2, c=1:6, d=1:5)
	# a | 1 2 3 4
	# b | 1 2
	# c | 1 2 3 4 5 6
	# d | 1 2 3 4 5
	# now you want to take a sum of each of those COLUMNS
	# apply2jagged(X=foo, FUN=sum)
	apply2jagged <- function(X, FUN, ...){
		# transform X into a matrix by padding vectors with NAs
		# done with this little interior pad function that takes vector x and pads it to length l
		padNA <- function(x, l){  c(x, rep(NA, l-length(x) ))  }
		maxlength <- max(sapply(X=X, FUN=length))
		x_matrix <- sapply(X=X, FUN=padNA, l=maxlength)
		return(apply(X=x_matrix, MAR=1, FUN=function(x){ FUN(as.vector(na.omit(x)), ...)}))
	}


##  condense_otu_table removes all zeroes prior to generating a collectors curve
	condense_otu_table <- function(x){
	condense_samp <- function(s, l){
	  saz <- (s > 0)
	  nnz <- sum(saz)
	  c(s[saz], rep(0, l - nnz))
	}
	ocmax <- max(colSums(x > 0))
	apply(X=x, MAR=2, FUN=condense_samp, l=ocmax)
	}
	
## mclapply2 function
	# this is a version of mclapply that has a PROGRESS INDICATOR
	# THANKS SO MUCH to wannymahoots on stackoverflow, what a HERO.
	# this is over twice as fast as pblapply from the pbapply package.
	# https://stackoverflow.com/questions/10984556/is-there-way-to-track-progress-on-a-mclapply
	mclapply2 <- function(X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L), mc.cleanup = TRUE, mc.allow.recursive = TRUE, mc.progress=TRUE, mc.style=3) {
		if(version$major==3 && version$minor==5.1){
			message("Ignore \"cannot wait for child\" warnings. Only in R 3.5.1, they are not a problem.")
		}
		if (!is.vector(X) || is.object(X)) X <- as.list(X)

		if (mc.progress) {
			f <- fifo(tempfile(), open="w+b", blocking=T)
			p <- parallel:::mcfork()
			pb <- txtProgressBar(0, length(X), style=mc.style)
			setTxtProgressBar(pb, 0) 
			progress <- 0
			if (inherits(p, "masterProcess")) {
				while (progress < length(X)) {
					readBin(f, "double")
					progress <- progress + 1
					setTxtProgressBar(pb, progress) 
				}
				cat("\n")
				parallel:::mcexit()
			}
		}
		tryCatch({
			result <- mclapply(X, ..., function(...) {
					res <- FUN(...)
					if (mc.progress) writeBin(1, f)
					res
				}, 
				mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,
				mc.silent = mc.silent, mc.cores = mc.cores,
				mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive
			)

		}, finally = {
			if (mc.progress) close(f)
		})
		result
	}

## multiple_collectors_curves()
	# this function generates collectors curves for a whole siteXspecies table
	# the output is a list of data.frame objects, each with columns alpha and step.
	# otu_table - siteXspecies table. if met="pd", all rownames MUST be ids present in phy's tip labels
	# stepsize - how often to report diversity (number of seqs sampled)
	# met - alpha div metric. either "sobs" or "pd". 
		# Maybe more if more options are added to alpha_rare_curve_1sample()
		# note that met="pd" takes a LONG time
	# max_d - max depth to construct curve (will stop making curve at maxd seqs sampled)
	# phy - needed if met="pd". a phylogenetic tree.
	# seed - just a random seed so this thing is repeatable
	multiple_collectors_curves <- function(siteXspecies, met="sobs", maxd=3500, stepsize=50, phy=NULL, seed=12345){
		set.seed(seed)
		output <- apply( X=siteXspecies, MAR=2, FUN=alpha_rare_curve_1sample, ids=rownames(otu_table), stepsize=stepsize, met=met, phy=phy, maxd=maxd)
		return(output)
	}

	# same function as above but parallelized
	# siteXspecies MUST be a data.frame for this to work, with COLUMNS AS SAMPLES AND ROWS AS SPECIES
	multiple_collectors_curves_parallel <- function(siteXspecies, met="sobs", maxd=3500, stepsize=50, phy=NULL, seed=12345, ncores=8){
		require(parallel)
		set.seed(seed)
		# not all processes will finish at the same time, so separate seeds are generated PER SAMPLE to make this repeatable
		seeds <- replicate(as.integer(paste(sample(0:9, 5), collapse="")), n=ncol(siteXspecies))
		# reconfigure siteXspecies to be a list, where each item in the list is the vector of occurrences, and also the seed (as a lit)
		siteXspecies_list <- list()
		for(i in 1:ncol(siteXspecies)){ siteXspecies_list[[i]] <- list(seed=seeds[i], obs=siteXspecies[,i]) }
		# do parallel operations
		#output <- mclapply(
		output <- mclapply2(
			X=siteXspecies_list, 
			FUN=function(x){alpha_rare_curve_1sample(
				x$obs, ids=rownames(siteXspecies), stepsize=stepsize, met=met, phy=phy, maxd=maxd, seed=x$seed
			)}, 
			mc.cores=ncores
		)
		return(output)
	}


## multiple_binned_collectors_curves
	# this function plots lots of collectors curves, as generate by multiple_collectors_curves()
	# these curves are all plotted in the background, and are binned by richness and averaged
	# ylabel - label for the plot's Y-axis. Should describe met.
	# n_bins - how many bins do you want?
	# bin_d - the value of d (x-axis) where curves are binned and aggregated by met.
	# ids - column vector of species names/ids. only needed if met="pd"
	# col_1 - color for odd bins
	# col_2 - color for even bins
	plot_multiple_binned_collectors_curves <- function(curves_list, n_bins=6, bin_d=1500, col_1="blue", col_2="red", col_0="gray", ...){
		# initialize plot
		max_val <- max( sapply(X=curves_list, FUN=function(x){max(x$alpha) } ))
		max_steps <- max( sapply(X=curves_list, FUN=function(x){max(x$step) } ))
		plot(NULL, ylim=c(0, max_val), xlim=c(0, max_steps), xlab="Number of sequences", ...)
		# i don't remember what this is for, so it's commented out
		# plot_fake_var <- sapply(X=curves_list, FUN=function(x){ points(x$alpha ~ x$step, col="gray", type="l")})

		# bin samples based on value at bin_d
		vals_at_bin_n  <- unlist(sapply(
			X=curves_list, 
			FUN=function(x, d){
				if(d>=min(x$step) & d <=max(x$step)){
					return(x$alpha[ which.min(abs(x$step - d)) ])
				}else{
					return(NA)
				}
			}, 
			d=bin_d 
		))

		# bounds for binning
		bin_bound_lo <- seq(from=min(na.omit(vals_at_bin_n)), to=max(na.omit(vals_at_bin_n)), length.out=n_bins+1)
		# there's an "extra" lower bound at the end of bin_bound_lo so that x<bin_bound_lo[i+1] makes sense later on
		# add a little bit to this value just so that x< makes sense for the top bin
		bin_bound_lo[n_bins+1] <- bin_bound_lo[n_bins+1] + 0.10 * bin_bound_lo[n_bins+1]

		# assign a bin to each sample. bin 0 is for samples that don't fit in any bin (maybe they never reached bin_d)
		samp_bins <- rep(0, length(vals_at_bin_n))
		for(b in 1:n_bins){ samp_bins[ vals_at_bin_n >= bin_bound_lo[b] & vals_at_bin_n < bin_bound_lo[b+1] ] <- b }

		# plot curve for each sample
		for(i in 1:length(samp_bins)){
			if(samp_bins[i]==0){ 
				col_i <- col2pastel(col_0, 0.3)
			}else if(samp_bins[i]%%2==1){
				col_i <- col2pastel(col_1, 0.3)
			}else{
				col_i <- col2pastel(col_2, 0.3)
			}
			points(curves_list[[i]]$alpha~curves_list[[i]]$step, col=col_i, type="l")
		}

		# plot averages for each bin
		for(b in 1:n_bins){
			# pull out curves that are in bin b
			rare_curves_i <- curves_list[ vals_at_bin_n >= bin_bound_lo[b] & vals_at_bin_n < bin_bound_lo[b+1] ]
			# plot them
			#lapply(X=rare_curves_i, FUN=function(x){ points(x$alpha ~ x$step, col=bin_cols_lite[b], type="l")})
			# calculate averages, ONLY if averaging 3+ curves
			mean3plus <- function(x){if(length(x) >= 3){return(mean(x))}else{return(NA)}}
			av_alpha_b <- apply2jagged(X=lapply(X=rare_curves_i, FUN=function(x){x$alpha}), FUN=mean3plus)
			av_step_b <- apply2jagged(X=lapply(X=rare_curves_i, FUN=function(x){x$step}), FUN=mean)

			# bin color
			col_b <- col_1
			if(b%%2==0){col_b <- col_2}

			# plot bin average
			points(av_alpha_b ~ av_step_b, lwd=2, col=col_b, type="l")
		}
	}
