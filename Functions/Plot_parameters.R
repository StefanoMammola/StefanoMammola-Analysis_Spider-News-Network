###############################################################

## ------------------------------------------------------------------------
# 'Custom plot parameters for ggplot2'
## ------------------------------------------------------------------------

## Authors: Stefano Mammola
## Software: R (v. R 4.1.0) and R studio (v. 1.4.1103)

###########
# Colors ##
###########
color_event   <- c("turquoise3", "orangered", "grey10")

################
# axis labels ##
################
axis_labels_sjPlot_m1 = c("Year of publication",
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

axis_labels_sjPlot_m2 = c("Year of publication",
                          "Type [Magazine]",
                          "Type [Online only]",
                          "Circulation [Internat.]",
                          "Circulation [National]",
                          "Event [Bite]",
                          "Event [Deadly Bite]",
                          "Sensationalism [Yes]",
                          "Expert [Doctor]",
                          "Expert [Arachnologist]",
                          "Expert [Others]")

axis_labels_sjPlot_m3 = c("N° of deadly spiders",
                          "Press freedom",
                          "N° of newspaper",
                          "Internet users",
                          "Sensationalistic news [Proportion]",
                          "News with errors [Proportion]",
                          "Lenguage [Spanish]",
                          "Lenguage [Russian]",
                          "Lenguage [Arabic]",
                          "Lenguage [English]")

axis_labels_ergm1 <- c("Edges", 
                       "nodeMatch [Language]", 
                       "Language [English]",
                       "Language [Others]",
                       "Language [Russian]",
                       "Lenguage [Spanish]",
                       "Sensationalism",
                       "Errors",
                       "Number of News Articles")

title_sjPlot_m1    = "Drivers of sensationalism"
title_sjPlot_m2    = "Drivers of errors [any error type]"
title_sjPlot_m3    = "Drivers of country centrality\nin the network"
title_ergm1        = "Drivers of country connectivity\nin the network"

xlab_sjPlot_m1_m2 = expression(paste("Odds ratio" %+-% "Standard Error"))
xlab_sjPlot_m3    = expression(paste("Incidence odds ratio" %+-% "Standard Error"))
xlab_ergm1        = expression(paste("Estimated coefficient" %+-% "Standard Error"))

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
                                       size = 2, linetype = "blank"))

}

theme_map_custom <- theme( 
  axis.line=element_blank(),axis.text.x=element_blank(),
  axis.text.y=element_blank(),axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),legend.position="none",
  panel.background=element_rect(fill ="gray40", colour="gray40"),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_rect(fill ="gray40", colour="gray40")
)
