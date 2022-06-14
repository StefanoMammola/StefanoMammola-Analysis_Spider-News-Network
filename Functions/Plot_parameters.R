###############################################################

## Mammola, S. et al. (2022) The global spread of misinformation on spiders. Current Biology.

###############################################################

## ------------------------------------------------------------------------
# 'Custom plot parameters for ggplot2'
## ------------------------------------------------------------------------

## Authors: Stefano Mammola
## Software: R (v. R 4.1.0) and R studio (v. 1.4.1103)

################
# axis labels ##
################
axis_labels_plot_model1 <- c("Year of publication",
                          "Type [Magazine]",
                          "Type [Online only]",
                          "Circulation [Internat.]",
                          "Circulation [National]",
                          "Event [Bite]",
                          "Event [Deadly Bite]",
                          "Species photo [Present]",
                          "Bite photo [Present]",
                          "Errors [Present]",
                          "Expert [Doctor]",
                          "Expert [Arachnologist]",
                          "Expert [Others]")

axis_labels_plot_model2 <- c("Prop. of sensationalistic news", 
                          "Prop. of news with errors",
                          "N° of deadly spiders",
                          "Internet users",
                          "Language [Arabic]",
                          "Language [Others]",
                          "Language [Russian]",
                          "Language [Spanish]",
                          "Node match: Language")

title_plot_model1 <- "Drivers of sensationalism"
title_plot_model2 <- "Drivers of probability of connection\nin the network"

xlab_plot_model1 <- expression(paste("Odds ratio" %+-% "Standard Error"))
xlab_model2      <- expression(paste("Mean effect size (95% confidence interval)"))

###########
# Themes ##
###########
theme_custom <- function(){
  theme_bw() +
    theme(
      axis.text = element_text(size = 12), 
      axis.title = element_text(size = 14),
      axis.line.x = element_line(color="black"), 
      axis.line.y = element_line(color="black"),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),                                          
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),  
      plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
      plot.title = element_text(size = 15, vjust = 1, hjust = 0),
      legend.text = element_text(size = 12),          
      legend.title = element_blank(),                              
      legend.position = c(0.95, 0.15), 
      legend.key = element_blank(),
      legend.background = element_rect(color = "black", 
                                       fill = "transparent", 
                                       size = 2, linetype = "blank"))}

theme_map_custom <- theme( 
  axis.line=element_blank(),axis.text.x=element_blank(),
  axis.text.y=element_blank(),axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),legend.position="none",
  panel.background=element_rect(fill ="white", colour="white"),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_rect(fill ="white", colour="white"))

# -- End of the script -- #