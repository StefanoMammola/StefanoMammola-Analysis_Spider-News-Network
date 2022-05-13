###############################################################

## Mammola, S. et al. (2022) The global spread of (mis-)information on spiders. Current Biology.

###############################################################

## ------------------------------------------------------------------------
# 'Custom functions for analyses'
## ------------------------------------------------------------------------

## Authors: Stefano Mammola
## Software: R (v. R 4.1.0) and R studio (v. 1.4.1103)

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
    
    LouvainModularity = Graph %>% cluster_louvain %>% membership %>% modularity(Graph, .),
    
    RealizedConnectance = ecount(Graph) / ((vcount(Graph)*(vcount(Graph)-1))/2)
    
  ) %>% return
  
}

#################################
# Other miscellaneous functions #
#################################

range01 <- function(x){(x-min(x))/(max(x)-min(x))} #Range
