###This script is used to generate a simple boxplot of nestedness indices for each site
# nestedness outputs were concatenated using cat *.out > All_Waimea_Sites.out

library(ggplot2)
library(gridExtra)
library(ggpubr)

nest_out  <-read.table("../data/processed/nestedness_outputs/All_Waimea_Sites.out", sep = "\t", col.names = c("HabiSite","Index","Value"))

# clean up table contents
nest_out$HabiSite <- sub(".+\\/(.+).csv", "\\1",nest_out$HabiSite)
nest_out$Value    <- as.numeric( gsub("\\[(.+)\\]","\\1", nest_out$Value) )
nest_out$Habitat  <- sub("(.+)-.+","\\1", nest_out$HabiSite)
nest_out$Site     <- sub(".+-(.+)","\\1", nest_out$HabiSite)
nest_out$Index[nest_out$Index == "X"] <- "web asymmetry"


# boxplots for each nestedness index

# set up comparisons for means
my_comparisons <- list( c(1, 2), c(2, 3), c(1, 3) )

index_plots <- lapply(
  X =  unique(nest_out$Index),
  FUN = function(an_index, data) {
    # subset by index
    nest_ind <- nest_out[nest_out$Index == an_index, ]
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
  },
  data = nest_out
)

g <-arrangeGrob(index_plots[[1]], index_plots[[2]], index_plots[[3]], index_plots[[4]],index_plots[[5]],
                nrow = 3, ncol = 2)

ggsave("outputs/nestedness_boxplots/nested_habitat_boxplot.pdf", plot = g, width = 10, height = 15)
