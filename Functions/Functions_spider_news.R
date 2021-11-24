###############################################################

## ------------------------------------------------------------------------
# 'Custom functions for analyses'
## ------------------------------------------------------------------------

## Authors: Stefano Mammola
## Software: R (v. R 4.1.0) and R studio (v. 1.4.1103)

##########################################################
# Function to extract characteristics of individual news #
##########################################################

#id:: a subset of the original database
extractID_attributes <- function(id){
  
  require("BAT") ; require("tidyverse")
  
  ID              <- as.character(id$ID_Event)[1] 
  n               <- nrow(id) 
  year            <- min(id$yr, na.rm = TRUE) 
  lon             <- id$lon[1]
  lat             <- id$lat[1]
  
  Temporal_span   <- ifelse(is.na(min(id$Year_news)), 
                            NA,  
                            length(seq(min(id$Year_news,na.rm = TRUE), max(id$Year_news,na.rm = TRUE), "days"))) 
  
  Country_event   <- ifelse(length(table(id$Country_event)) == 0,
                            NA,
                            names(sort(table(id$Country_event), decreasing = TRUE))[1])
  
  Genus           <- names(sort(table(id$Genus), decreasing = TRUE))[1] 
  Sensationalism  <- sum(id$Sensationalism, na.rm = TRUE)/nrow(id) 
  Taxonomic_error <- sum(id$Taxonomic_error, na.rm = TRUE)/nrow(id) 
  Venom_error     <- sum(id$Venom_error, na.rm = TRUE)/nrow(id) 
  Anatomy_error   <- sum(id$Anatomy_error, na.rm = TRUE)/nrow(id) 
  Photo_error     <- sum(id$Photo_error, na.rm = TRUE)/nrow(id) 
  Bite            <- sum(id$Bite, na.rm = TRUE)/nrow(id) 
  Death           <- sum(id$Death, na.rm = TRUE)/nrow(id)
  
  Distance_news <- ifelse(nrow(id) == 1,
                          NA,
                          id %>% dplyr::select(Species,Bite,Death,Figure_species,
                                        Figure_bite,Taxonomic_error,Venom_error,
                                        Anatomy_error,Photo_error) %>% BAT::gower(convert=c(1:9)) 
                          
  ) #get news distance
  
  Distance_information <- ifelse(is.na(Distance_news), 0, mean(as.matrix(Distance_news)[,1]))
  
  return( data.frame(ID, n, year, Temporal_span, lon, lat,
                     Country_event, Genus, Bite, Death,
                     Sensationalism,
                     Distance_information,
                     Taxonomic_error, Venom_error, Anatomy_error, Photo_error))
  
}

##############################################
# Function to extract node traits in a graph #
##############################################

# For a discussion:
# https://medium.com/@615162020004/social-network-analysis-with-r-centrality-measure-86d7fa273574

# Graph should be a tidygraph dataset

NodeTraitGet <- function(Graph, mode = "in", dir = TRUE){
  
  list(
    ID          = Graph %>% activate(nodes) %>% as.data.frame %>% pull(1),
    Degree      = degree(Graph, mode = mode),
    Strength    = strength(Graph, mode = mode),
    Eigenvector = Graph %>% eigen_centrality(directed = dir) %>% extract2("vector"),
    Betweenness = Graph %>% betweenness(directed = dir),
    Closeness   = Graph %>% closeness(mode = mode)
    
    # Eigenvector_Weighted = Graph %>% 
    #   eigen_centrality(weights = Graph %>% activate(edges) %>% pull(weight)) %>% 
    #   extract2("vector"),
    #<< 
    
    #Clustering = Graph %>% transitivity(type = "local"),
    
  )
}

#########################################
# Function to extract graph-level stats #
#########################################

# Graph should be an iGraph object

NetworkTraitGet <- function(Graph){
  
  data.frame(
    
    Size = vcount(Graph),
    
    Diameter = diameter(Graph),
    
    MeanDegree = mean(igraph::degree(Graph)),
    
    DegreeVariance = sd(igraph::degree(Graph)),
    
    Components = igraph::components(Graph)$no,
    
    Transitivity = transitivity(Graph),
    
    Density = ggregplot::Prev(get.adjacency(Graph, sparse = F) > 0),
    
    LouvainModularity = Graph %>% cluster_louvain %>% membership %>% modularity(Graph, .)
    
  ) %>% return
  
}

#################################
# Other miscellaneous functions #
#################################

range01 <- function(x){(x-min(x))/(max(x)-min(x))} #Range
