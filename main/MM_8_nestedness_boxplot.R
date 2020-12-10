###This script is used to generate a simple boxplot of nestedness indices for each site
# nestedness outputs were concatenated using cat *.out > All_Waimea_Sites.out

library(ggplot2)
library(gridExtra)
library(ggpubr)

h2_out  <-read.table("../data/processed/nestedness_outputs/all_samples_empo_3_H2.out", sep = "\t", col.names = c("File","Index","Value"))
nodf_out <- read.table("../data/processed/nestedness_outputs/all_samples_empo_3_NODF.out", sep = "\t", col.names = c("File","Index","Value"))

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
plot_nest_box <- function(an_index, nest_data) {
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
  print( tukey_out)
  return(p)
}

rar_habisite_plots <- lapply(X =  unique(habisites$Index),
                            FUN = plot_nest_box,
                            nest_data = habisites[habisites$rarefied == T, ])

raw_habisite_plots <- lapply(X =  unique(habisites$Index),
                            FUN = plot_nest_box,
                            nest_data = habisites[habisites$rarefied == F, ])


g1 <-arrangeGrob(rar_habisite_plots[[1]], rar_habisite_plots[[2]],
                nrow = 1, ncol = 2)

ggsave("outputs/nestedness_boxplots/rarefied_empo_3_agg_habitat_nest_boxplot.pdf", plot = g1, width = 8, height = 6)

g2 <-arrangeGrob(raw_habisite_plots[[1]], raw_habisite_plots[[2]],
                 nrow = 1, ncol = 2)

ggsave("outputs/nestedness_boxplots/raw_count_empo_3_agg_habitat_nest_boxplot.pdf", plot = g2, width = 8, height = 6)

write.csv(habisites, "../data/processed/H2_NODF_empo3_agg_habitat_site.csv")


