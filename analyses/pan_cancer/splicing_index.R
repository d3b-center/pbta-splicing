################################################################################
# splicing_index.R
# script that takes in "results/splicing_index.wdPSI_per_sample.txt" data file 
# and computes relative proportion of aberrant splicing changes in samples
# written by Ammar Naqvi
#
# usage: Rscript splicing_index.R
################################################################################

library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)

dataDir = "~/Desktop/AS-DMG/analyses/pan_cancer/results/"

## more accurate version -- dPSI version
file <- "splicing_index.wdPSI_per_sample.txt"
splice_index  = read.delim(paste0(dataDir, file), sep = "\t", header=TRUE,row.names=1)


## violin plot version
e <- ggplot(splice_index, aes(x = hist, y = splice_index))


e + geom_violin(trim = FALSE) + 
  stat_summary(
    fun.data = "mean_sdl",  fun.args = list(mult = 1), 
    geom = "pointrange", color = "black"
  )

e + geom_violin(aes(fill = hist), trim = FALSE) + 
  geom_boxplot(width = 0.2) +
  #scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07")) +
  geom_point(color="black", size=1, position = position_jitter(w=0.02)) +
  theme(legend.position = "none") + theme_Publication()

plots_dir <- "/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/plots/"

# Save plot as PNG
#png(file.path(plots_dir, "splice_index_initial_v2.png"))

dev.copy(png,'/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/plots/splice_index_initial_v2.png', width=1550, height=1494, res=120)
dev.off()

## Make CDF plot for splicing index
splice_index <- splice_index %>%
  as.data.frame(stringsAsFactors = FALSE)


# Set up the data.frame for plotting
si_cdf <- splice_index %>%
  
  # We only really need these two variables from data.frame
  dplyr::transmute(
    group = hist,
    number = as.numeric(splice_index)
  ) %>%
  
  # Group by specified column
  dplyr::group_by(group) %>%
  
  # Only keep groups with the specified minimum number of samples
  dplyr::filter(dplyr::n() > 1) %>%
  
  # Calculate group median
  dplyr::mutate(
    group_median = median(number, na.rm = TRUE),
    group_rank = rank(number, ties.method = "first") / dplyr::n(),
    sample_size = paste0("n = ", dplyr::n())
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(group = reorder(group, group_median))

# Set y_lim as min and max if it was not set
#if (is.na(y_lim)) {
# y_lim = c(min(df$number),
#           max(df$number))
#}

si_cdf %>%
  # Now we will plot these as cumulative distribution plots
  ggplot2::ggplot(ggplot2::aes(
    x = group_rank,
    y = number
  )) +
  
  ggplot2::geom_point(color = "black") +
  
  # Add summary line for median
  ggplot2::geom_segment(
    x = 0, xend = 1, color = "blue",
    ggplot2::aes(y = group_median, yend = group_median)
  ) +
  
  # Separate by histology
  ggplot2::facet_wrap(~ group + sample_size, nrow = 1, strip.position = "bottom") +
  ggplot2::theme_classic() +
  ggplot2::xlab("Histology") +
  ggplot2::ylab("Splicing Index") +
  
  # Making it pretty
  ggplot2::theme(legend.position = "none") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    strip.placement = "outside",
    strip.text = ggplot2::element_text(size = 10, angle = 90, hjust = 1),
    strip.background = ggplot2::element_rect(fill = NA, color = NA)
  ) +
  ggplot2::ggtitle("Splicing Index")
dev.copy(png,'/Users/naqvia/Desktop/AS-DMG/analyses/pan_cancer/plots/splice_index_CDF_v2.png', width=1550, height=1494, res=120)
dev.off()


##theme for all plots
theme_Publication <- function(base_size=15, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            # legend.margin = unit(0.5, "cm"),
            legend.margin = margin(5,5,5,5),
            legend.title = element_text(face="bold"),
            #plot.margin=unit(c(10,5,5,5),"mm"),
            plot.margin=unit(c(10,5,5,10),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}


