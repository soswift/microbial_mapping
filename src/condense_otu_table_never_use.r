# an ungodly function that's only good for making
# an otutable smaller for species-observed collector's
# curves and NOTHING ELSE
# x - an otutable full of ints where cols are samps
# 	and rows are species
condense_otu_table <- function(x){
	condense_samp <- function(s, l){
		saz <- (s > 0)
		nnz <- sum(saz)
		c(s[saz], rep(0, l - nnz))
	}
	ocmax <- max(colSums(x > 0))
	apply(X=x, MAR=2, FUN=condense_samp, l=ocmax)
}



# shitty test code
# shity otutable with 250000 OTUs and 1600 samples 
shit <- matrix(sample(c(rep(0,11),1), 4E8, replace=T), ncol=1600)
shit2 <- condense_otu_table(shit)
dim(shit2)
