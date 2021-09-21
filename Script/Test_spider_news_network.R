###############################################################

## Flow of spider related information
## Mammola, S. et al. 2021

###############################################################

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

## Authors: Stefano Mammola
## Last update: 29 Jul 2021, Helsinki, Finland
## Software: R (v. R 4.1.0) and R studio (v. 1.4.1103)

###############################################################

# clean the workspace -----------------------------------------------------

rm(list=ls())

# Loading R package -------------------------------------------------------

library("Amelia")        # for exploring missing data
library("BAT")
library("bipartite")     # for network
library("plyr")         # for data wrangling and pipe function
library("geosphere")     # for spherical operations
library("GGally")        # for network
library("ggalt")         # for plotting
library("ggplot2")       # for plotting
library("ggthemes")      # for plotting
library("gridExtra")     # for plates
library("igraph")        # for network (measures of centrality)
library("lme4")          # for glmm
library("maps")          # for the world map
library("network")       # for network analyses
library("performance")   # model validation
library("PupillometryR") # for plotting
library("sna")           # for network plotting
library("tidyverse")     # for tidy operations
library("tidygraph")
library("magrittr")
library("ggraph")

# library(tidyverse); library(igraph); library(ggraph); library(tidygraph); library(ggregplot)
# library(fs); library(magrittr); library(cowplot)

# Creating useful functions ------------------------------------------------

db_i <- db %>% filter(ID_Event == "R70") %>%
  arrange(Year_news) %>% droplevels()

a <- db_i %>% select(Species,Bite,Death,Figure_species,
                Figure_bite,Taxonomic_error,Venom_error,
                Anatomy_error,Photo_error) %>% BAT::gower(convert=c(1:9))

#cleaver way to do???
DeerLocations %<>% group_by(Code) %>% 
  summarise_at(c("Easting", "Northing"), ~mean(.x, na.rm = T)) %>% 
  ungroup

# mean(as.matrix(a)[,1])

##########################################################
# Function to extract characteristics of individual news #
##########################################################

#Id should be a subset of the original database
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
                          id %>% select(Species,Bite,Death,Figure_species,
                                        Figure_bite,Taxonomic_error,Venom_error,
                                        Anatomy_error,Photo_error) %>% BAT::gower(convert=c(1:9)) 
                          
  ) #get news distance
  
  Distance_information <- ifelse(is.na(Distance_news), NA, mean(as.matrix(Distance_news)[,1]))
  
  # lm?
  
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

NodeTraitGet <- function(Graph, mode = "in", dir = TRUE){
  
  list(
    ID = Graph %>% activate(nodes) %>% as.data.frame %>% pull(1),
    Degree = degree(Graph, mode = mode),
    Strength = strength(Graph, mode = mode),
    Eigenvector = Graph %>% eigen_centrality(directed = dir) %>% extract2("vector"),
    Betweenness = Graph %>% betweenness(directed = dir),
    Closeness = Graph %>% closeness(mode = mode)
    
    # Eigenvector_Weighted = Graph %>% 
    #   eigen_centrality(weights = Graph %>% activate(edges) %>% pull(weight)) %>% 
    #   extract2("vector"),
    #<< 
    
    #Clustering = Graph %>% transitivity(type = "local"),
  
  )
}

##############################################
# Function to extract graph-level stats ######
##############################################

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


# Plot style. Modified from:
# https://ourcodingclub.github.io/tutorials/dataviz-beautification-synthesis/

theme_custom <- function(){
  theme_bw() +
    theme(#text = element_text(family = "Arial"),
      axis.text = element_text(size = 14), 
      axis.title = element_text(size = 18),
      axis.line.x = element_line(color="black"), 
      axis.line.y = element_line(color="black"),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),                                          
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),  
      plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
      plot.title = element_text(size = 18, vjust = 1, hjust = 0),
      legend.text = element_text(size = 12),          
      legend.title = element_blank(),                              
      legend.position = c(0.95, 0.15), 
      legend.key = element_blank(),
      legend.background = element_rect(color = "black", 
                                       fill = "transparent", 
                                       size = 2, linetype = "blank"))
}

#Range
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

###############################################################

## Data preparation:

###############################################################

# Loading the Database ----------------------------------------------------

db <- read.csv(file = "/Users/stefanomammola/Desktop/course on NETWORKS/Me playing around/Data_spider_news_global.csv",sep='\t', dec='.',header=T,as.is=F)

str(db)
dim(db)

#Converting factors to chr
db$Title <- as.character(db$Title)
db$Notes <- as.character(db$Notes)

#Sort
db <- db %>% arrange(Country_search) 

# Creating new variables --------------------------------------------------

# Rename the ID_event
levels(db$ID_Event) <- c(paste(rep("R",nlevels(db$ID_Event)),1:nlevels(db$ID_Event),sep=""))  

#Type of Event
db <- db %>% mutate(TypeEvent = Bite + Death)

db$TypeEvent <- as.factor(db$TypeEvent) ; levels(db$TypeEvent) <- c("Encounter","Bite","Deadly bite")

#Total number of Errors
db <- db %>% mutate(TotalError = tidyr::replace_na(as.numeric(Taxonomic_error),0) + 
                      tidyr::replace_na(as.numeric(Venom_error),0)     + 
                      tidyr::replace_na(as.numeric(Anatomy_error),0)  + 
                      tidyr::replace_na(as.numeric(Photo_error),0))

#Expert
db <- db %>% mutate(TotalExpert = tidyr::replace_na(as.numeric(Expert_arachnologist),0) + 
                      tidyr::replace_na(as.numeric(Expert_doctor),0)        + 
                      tidyr::replace_na(as.numeric(Expert_others),0))

db$TotalExpert <- ifelse(db$TotalExpert > 0 , 1 , 0) #converting to binary

#ID for the event
db$ID_event <- as.factor(paste(db$Location_event,db$Year_event,sep="_"))

#Date
db$Year_news   <-  as.Date(paste(db$d,db$m,db$yr,sep="/"),format='%d/%m/%Y')

#Database only with distinct event
db_unique_event <- distinct(db, ID_event, .keep_all = TRUE) 

#Database only with distinct news
db_unique_news  <- distinct(db, ID, .keep_all = TRUE) 

###########################################################################
###########################################################################
# Network analyses --------------------------------------------------------
###########################################################################
###########################################################################

# TYPE OF DATA
# I have collected media news on spiders published between 2010 - 2020.
# Each news reports about a specific spider bite event. There can be multiple articles published in different countries talking about the same event.

# QUESTIONS / Are they edge-level, node-level, network-level? 

# First question: (Node-level) Why are some countries more influential than others in 
# generating spider-related content taken up in  other places?

# Second question: (Edge-level) Why determines the strength of specific connection (news flow) among countries?

# Third question: (Network level) 
# Is the network diffuse or country are more isolated based? [not sure how to answer, maybe some measure edge density?]
# If I slice the network by years, does the density (connectiveness?) change over time?

## What is the unit of analysis in each case?

# Nodes: each Country

# Edged: flow of news between countries (n° of news). directional.

# To answer the node-level question, I can get variable for each country (e.g. number of venomous species,
# )

# TYPE OF NETWORK: here I'm not fully sure.

# Network: Adjacency matrix, simmetrical (=unipartite). 
# Or should it be a edge list cause I want to associate variables to each country?


# What processes structure your network?

#-- Country level meta-data
  
# Node-level/edge-level/neither?

#-- Missing nodes (namely, unsampled news and unsample countries). Higher probability of missing nodes (unsampled news) in the past.

# Which of these are known, and which are hypotheses?

#-- The existence of missing node is known. Country level difference is an hypothesis.

# What sampling biases are you vulnerable to? How would you account for them? 
  
#-- Hard to say... probably I would randomly remove or add nodes in the observed network, and see if the resulting properties of the network deviate significantly from the observed one

##### Attribute for each event

sort(table(db$ID_Event), decreasing = TRUE)

#Extract ID attributes for each ID event
ID_attributes <- db %>% filter(ID_Event == levels(ID_Event)[1]) %>%
                 arrange(Year_news) %>% droplevels() %>% extractID_attributes #first event

for(i in 2:nlevels(db$ID_Event)) #takes about 1 minutes!
  ID_attributes <- db %>% filter(ID_Event == levels(ID_Event)[i]) %>% arrange(Year_news) %>% droplevels() %>%
                   extractID_attributes %>% bind_rows(ID_attributes) #all others

ID_attributes %<>% mutate_all(function(x) ifelse(is.infinite(x), 0, x)) 

ID_attributes %>% head(20)

write.table(ID_attributes,"ID_attributes.txt")
ID_attributes <- read.table("ID_attributes.txt")

ggplot(data = ID_attributes, aes(x= Distance_information, y= Temporal_span)) + geom_point()
ggplot(data = ID_attributes, aes(x= lon, y= Temporal_span)) + geom_point()
ggplot(data = ID_attributes, aes(x= lat, y= Temporal_span)) + geom_point()

##### NEWS-LEVEL NETWORK ####

# Generate a bi-partite network -------------------------------------------

# Get subset database
db_graph <- db %>% dplyr::select(Country_search,ID_Event) %>% table() 
  
# Network
Graph_bipartite <- igraph::graph_from_incidence_matrix(db_graph, directed = TRUE)

V(Graph_bipartite)$type
print(Graph_bipartite, e=TRUE, v=TRUE)
vcount(Graph_bipartite) ; ecount(Graph_bipartite)

#Graph_bipartite <- igraph::simplify(Graph_bipartite, edge.attr.comb = "sum") 

# Get attribute table
Graph_tbl_bip <- tidygraph::as_tbl_graph(Graph_bipartite, directed = TRUE)  

# Collapse it into an unipartite 
Graph_unipartite_full <- igraph::bipartite_projection(Graph_bipartite)

# Country -----------------------------------------------------------------

Graph_unipartite <- Graph_unipartite_full$proj1

#summary stats
NetworkTraitGet(Graph_unipartite)

E(Graph_unipartite)$weight

# Get the adjacency matrix
Graph_adj_matrix <- Graph_unipartite %>% get.adjacency(attr = "weight", sparse = FALSE)

vcount(Graph_unipartite) ; ecount(Graph_unipartite)

# Get attribute table
Graph_tbl_uni <- Graph_unipartite %>% as_tbl_graph(directed = TRUE) %>% 
  activate(edges) %>% #%>% mutate(weight = 1) 
  igraph::simplify(edge.attr.comb = "sum") %>% 
  as_tbl_graph

# Calculating node level traits 

Node_attributes <- Graph_tbl_uni %>% NodeTraitGet(mode = "in", dir = TRUE) %>% bind_cols()

# Storing & cleaning
Sens <- db %>% dplyr::select(Country_search,Sensationalism,TotalError,Lenguage,lon3,lat3) %>% group_by(Country_search) %>% 
                              summarise(Sens = sum(Sensationalism, na.rm = TRUE)/length(Sensationalism),
                   TotalError = sum(TotalError, na.rm = TRUE)/length(TotalError),
                   Lenguage = names(sort(table(Lenguage), decreasing = TRUE))[1],
                   lon = mean(lon3),
                   lat = mean(lat3))

Sens$Lenguage <- ifelse(Sens$Lenguage %in% names(which(table(Sens$Lenguage)>9)), Sens$Lenguage , "Others")
#Sens$Lenguage <- ifelse(Sens$Lenguage == "English","English", "Non-English")

N <- db %>% group_by(Country_search) %>% count()

Node_attributes <- data.frame(Node_attributes,
                              Sens[,-1],
                              N = N[,2])

Graph_tbl_uni <- Graph_tbl_uni %>% tidygraph::activate(nodes) %>% 
  left_join(Node_attributes, by = c("name" = "ID"))

Graph_tbl_uni %>% activate(nodes) %>% as.data.frame

Graph_tbl_uni %>% activate(nodes) %>% as.data.frame %>% 
  ggplot(aes(y=Strength,x=Sens)) + geom_point()

Graph_tbl_uni %>% activate(nodes) %>% as.data.frame %>% 
  ggplot(aes(y=Strength,x=n)) + geom_point()

Graph_tbl_uni %>% activate(nodes) %>% as.data.frame %>% 
  ggplot(aes(y=Strength,x=Lenguage)) + geom_boxplot()

SpatialLayout <- Node_attributes %>% dplyr::select(lon,lat) %>% as.matrix #geo-coordinates
Layout2 <- layout_with_kk(Graph_tbl_uni) # Kamada-Kawai
Layout3 <- layout_with_mds(Graph_tbl_uni) # Multi-dimensional scaling

#Comparing 3 visualisations
Graph_tbl_uni %>% igraph::simplify(edge.attr.comb = "sum") %>% ggraph::ggraph(SpatialLayout) +
  geom_edge_density(fill="orange", alpha=1) +
  geom_edge_fan(aes(width=weight),color="gray60", alpha=0.1) +
  geom_node_point(col="grey30", alpha = .8, 
                  aes(size=n,fill=Lenguage), shape = 21) + 
  geom_node_text(aes(label = name), size=2, color="gray10", repel=TRUE) +
  scale_fill_manual(values = c("blue", "orange", "turquoise","purple", "grey15"))+
  theme_void() + theme(legend.position = "bottom",legend.direction = "vertical")+ coord_fixed()# add edges to the plot geom_node_point()

Graph_tbl_uni %>% igraph::simplify(edge.attr.comb = "sum") %>% ggraph::ggraph(Layout2) +
  geom_edge_density(fill="orange", alpha=0.9) +
  geom_edge_fan(aes(width=weight),color="gray70", alpha=0.1) +
  geom_node_point(col="grey30", alpha = .8, 
                  aes(size=n,fill=Lenguage), shape = 21) + 
  scale_fill_manual(values = c("blue", "orange", "turquoise","purple", "grey15"))+
  geom_node_text(aes(label = name), size=2, color="gray10", repel=TRUE) +
  theme_void() + theme(legend.position = "bottom",legend.direction = "vertical")+ coord_fixed()# add edges to the plot geom_node_point()

Graph_tbl_uni %>% igraph::simplify(edge.attr.comb = "sum") %>% ggraph::ggraph(Layout3) +
  geom_edge_density(fill="orange", alpha=0.9) +
  geom_edge_fan(aes(width=weight),color="gray70", alpha=0.1) +
  geom_node_point(col="grey30", alpha = .8, 
                  aes(size=n,fill=Lenguage), shape = 21) + 
  scale_fill_manual(values = c("blue", "orange", "turquoise","purple", "grey15"))+
  geom_node_text(aes(label = name), size=2, color="gray10", repel=TRUE) +
  theme_void() + theme(legend.position = "bottom")+ coord_fixed()# add edges to the plot geom_node_point()

# plot removing disconnected nodes
Graph_tbl_uni %>% activate(nodes) %>% 
  mutate(Degree = degree(.)) %>% 
  filter(Degree > 0)  %>% 
  #actual plot
  ggraph::ggraph(layout="kk") +
  #geom_edge_density(fill="orange", alpha=0.1) +
  #geom_edge_fan(color="gray30", width=0.1, alpha=0.7) +
  geom_edge_arc(aes(width = weight),col = "grey90", alpha=0.05) +
  geom_node_point(col="grey5", alpha = .9, 
                  aes(size=Closeness, fill = Lenguage), shape = 21) + 
  geom_node_text(aes(label = name), size=3, col="gray10", repel=TRUE) +
  theme_void() + theme(legend.position = "bottom") + coord_fixed()# add edges to the plot geom_node_point()

# Try to fit a model ------------------------------------------------------

# Exponential Random Graph Models

# Getting an adj matrix, weighted by the number of events, for the country level projected network
AdjMatrix <- Graph_tbl_uni %>% get.adjacency(attr = "weight", sparse = FALSE) 

#Convert to incidence
AdjMatrix[AdjMatrix>0] <- 1 ; AdjMatrix %>% dim #check it's a square

colnames(AdjMatrix) <- rownames(AdjMatrix) <- Graph_tbl_uni %>% activate(nodes) %>% pull(name) #assign row and col names

#Response variables
ResponseNetwork <- AdjMatrix %>% as.matrix %>% network::network(directed = FALSE)

#Adding node-level attributes
Graph_tbl_uni %>% as.data.frame %>% colnames

ResponseNetwork %v% "n" <- Graph_tbl_uni %>% as.data.frame %>% pull(n)

ResponseNetwork %v% "Lenguage" <- Graph_tbl_uni %>% as.data.frame %>% pull(Lenguage) # should be as.character

ResponseNetwork %v% "Sens" <- Graph_tbl_uni %>% as.data.frame %>% pull(Sens)

ResponseNetwork %v% "TotalError" <- Graph_tbl_uni %>% as.data.frame %>% pull(TotalError)

Spatial <- SpatialLayout %>% dist %>% as.matrix

##temporal network is a list of network

#Model fit
ergm0 <- ergm::ergm(ResponseNetwork ~ edges, estimate = "MLE")
ergm1 <- ergm::ergm(ResponseNetwork ~ edges + nodefactor("Lenguage"), estimate = "MLE")
ergm2 <- ergm::ergm(ResponseNetwork ~ edges + nodefactor("Lenguage") + nodecov("Sens"), estimate = "MLE")
ergm3 <- ergm::ergm(ResponseNetwork ~ edges + nodefactor("Lenguage") + nodecov("Sens") + nodecov("TotalError"), estimate = "MLE")
ergm4 <- ergm::ergm(ResponseNetwork ~ edges + nodefactor("Lenguage") + nodecov("Sens") + nodecov("TotalError") + edgecov(Spatial), estimate = "MLE")
ergm5 <- ergm::ergm(ResponseNetwork ~ edges + edgecov(Spatial), estimate = "MLE")

#nodematch("Lenguage", diff = TRUE) # give individual estimate

list(ergm0, ergm1, ergm2, ergm3, ergm4) %>% map_dbl(BIC) %>% plot

ERGM1 <- ergm::ergm(ResponseNetwork ~
                edges + nodecov("Sens") + nodecov("TotalError") +
                nodefactor("Lenguage") + nodematch("Lenguage") + dyadcov(Spatial),
                estimate = "MLE")

#Model validation

ergm3 %>% ergm::gof() %>% plot

# Simulating ####

BlankNetwork <- ResponseNetwork
BlankNetwork[] <- 0
# BlankNetwork %>% as.network(directed = T)

NSims <- 1000

Sims <- simulate(ERGM1, 
                 nsim = NSims, 
                 basis = BlankNetwork,
                 # reference = ~Poisson("random"),
                 # response = "value"
)

Sims[[1]] %>% as.matrix %>% #c %>% 
  Prev

Sims %>% 
  map_dbl(~.x %>% as.matrix %>% c %>% Prev) %>% 
  qplot() + 
  geom_vline(xintercept = Prev(as.matrix(ResponseNetwork)),
             colour = "red")

qplot(1:NSims, 
      Sims %>% 
        map_dbl(~.x %>% as.matrix %>% c %>% Prev)) + 
  geom_hline(yintercept = Prev(as.matrix(ResponseNetwork)),
             colour = "red")

qplot(1:NSims, 
      Sims %>% 
        map_dbl(~.x %>% as.matrix %>% c %>% Prev)) + 
  ylab("prevalence")+
  geom_hline(yintercept = Prev(as.matrix(ResponseNetwork)),
             colour = "red") + 
  lims(x = c(250, 750))

#Plotting
EstimateDF <- 
  ERGM1 %>% summary %>% extract2("coefficients") %>% 
  as.data.frame %>% 
  rownames_to_column("Variable")

ERGM1 %>% 
  confint %>% # Grabbing the 95% confidence intervals
  as.data.frame %>% 
  rename(Lower = 1, Upper = 2)

EstimateDF %<>% # Bind them together
  bind_cols(ERGM1 %>% confint %>% as.data.frame %>% 
              rename(Lower = 1, Upper = 2))

EstimateDF %>% head

# EstimateDF %>% ggplot(aes(Variable, Estimate)) +
#   geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.3) +
#   geom_point()

EstimateDF %>% ggplot2::ggplot(aes(Variable, Estimate)) +
  geom_hline(lty = 2, yintercept = 0) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  geom_point() + theme_bw() +
  coord_flip()

#### Count data

http://statnet.org/Workshops/valued.html#2_network_and_edge_attributes
library("ergm.count")

AdjMatrix <- Graph_tbl_uni %>% get.adjacency(attr = "weight", sparse = FALSE) 
AdjMatrix %>% as.matrix %>% dim #square

ResponseNetwork <- as.network(AdjMatrix %>% as.matrix, 
                              directed = TRUE, 
                              matrix.type = "a", 
                              ignore.eval = FALSE, 
                              names.eval = "weight")  # Important! )

#Adding node-level attributes
Graph_tbl_uni %>% as.data.frame %>% colnames

ResponseNetwork %v% "n" <- Graph_tbl_uni %>% as.data.frame %>% pull(n)

ResponseNetwork %v% "Lenguage" <- Graph_tbl_uni %>% as.data.frame %>% pull(Lenguage) # should be as.character

ResponseNetwork %v% "Sens" <- Graph_tbl_uni %>% as.data.frame %>% pull(Sens)

ResponseNetwork %v% "lon" <- Graph_tbl_uni %>% as.data.frame %>% pull(lon)

ResponseNetwork %v% "lat" <- Graph_tbl_uni %>% as.data.frame %>% pull(lat)

ResponseNetwork %v% "TotalError" <- Graph_tbl_uni %>% as.data.frame %>% pull(TotalError)

as.matrix(ResponseNetwork, attrname = "weight")

#
ergm::summary_formula(ResponseNetwork ~ sum, response = "weight")

set.edge.value(ResponseNetwork,
               "weight",
               c(AdjMatrix))

library("ergm.count")

NMCMC <- 10000

test <- ergm(ResponseNetwork ~ sum + #nonzero +  
                 nodecov("Sens") + 
               # nodecov("TotalError") +
                 nodefactor("Lenguage"), control = control.ergm(
                 parallel =10, #parallel.type="PSOCK",
                 MCMC.samplesize = NMCMC,
                 MCMLE.maxit = 50),
                      response = "weight", reference = ~ Poisson)


summary(test)


mcmc.diagnostics(test)

summary(test)

EstimateDF <- 
  test %>% summary %>% extract2("coefficients") %>% 
  as.data.frame %>% 
  rownames_to_column("Variable")

test %>% 
  confint %>% # Grabbing the 95% confidence intervals
  as.data.frame %>% 
  rename(Lower = 1, Upper = 2)

EstimateDF %<>% # Bind them together
  bind_cols(test %>% confint %>% as.data.frame %>% 
              rename(Lower = 1, Upper = 2))

EstimateDF %>% head

# EstimateDF %>% ggplot(aes(Variable, Estimate)) +
#   geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.3) +
#   geom_point()

EstimateDF %>% ggplot(aes(Variable, Estimate)) +
  geom_hline(lty = 2, yintercept = 0) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  geom_point() + theme_bw() +
  coord_flip()



# Event -----------------------------------------------------------------

Graph_unipartite2 <- igraph::simplify(Graph_unipartite_full$proj2, edge.attr.comb = "sum") 

#summary stats
NetworkTraitGet(Graph_unipartite2)

#plot(Graph_unipartite$proj1, edge.arrow.size = 0.5, edge.arrow.mode = "-")
#plot(Graph_unipartite$proj2, edge.arrow.size = 0.5, edge.arrow.mode = "-")

E(Graph_unipartite2)$weight

# Get the adjacency matrix
Graph_adj_matrix2 <- Graph_unipartite2 %>% get.adjacency(sparse = FALSE)

#Graph_unipartite$proj1
#Graph_unipartite$proj2
vcount(Graph_unipartite2)
ecount(Graph_unipartite2)

# Get attribute table
Graph_tbl_uni2 <- Graph_unipartite2 %>% as_tbl_graph(directed = TRUE) %>% 
  activate(edges) %>% #%>% mutate(weight = 1) 
  igraph::simplify(edge.attr.comb = "sum") %>% 
  as_tbl_graph

# Calculating node level traits

Node_attributes2 <- Graph_tbl_uni2 %>% NodeTraitGet(mode = "in", dir = TRUE) %>% bind_cols()

Node_attributes2 %>% head()

#Add node attributes
Graph_tbl_uni2 <-  Graph_tbl_uni2 %>% tidygraph::activate(nodes) %>% 
  left_join(Node_attributes2, by = c("name" = "ID"))

#Add news attributes
Graph_tbl_uni2 <-  Graph_tbl_uni2 %>% tidygraph::activate(nodes) %>% 
  left_join(ID_attributes, by = c("name" = "ID"))

# Check
Graph_tbl_uni2 %>% activate(nodes) %>% as.data.frame %>% 
  ggplot(aes(y=Strength,x=lat)) + geom_point()

# plot
Graph_tbl_uni2 %>% activate(nodes) %>% 
  mutate(Degree = degree(.)) %>% 
  filter(Degree > 0) %>% 
  ggraph::ggraph(layout="kk") +
  #geom_edge_fan0(col = "grey90", alpha= .05) +
  geom_edge_density(fill="orange", alpha=0.9) +
  geom_node_point(col="grey5", alpha = .9, fill="blue", shape = 21) + 
  theme_void() + coord_fixed()# add edges to the plot geom_node_point()


# Network by year ---------------------------------------------------------

# Temporal ERGM
# https://www.rdocumentation.org/packages/tergm/versions/3.7.0

db$yr_f <- as.factor(db$yr)

db$Year_news_f   <-  as.factor(paste(db$m,db$yr,sep="/"))

list_db_yr <- list()

for(i in levels(db$Year_news_f))
  list_db_yr[[i]] <- db %>% filter(Year_news_f == i)

NetworkTraitGet2 <- function(Graph){
  
  data.frame(
    
    Size = vcount(Graph),
    
    Diameter = diameter(Graph),
    
    MeanDegree = mean(igraph::degree(Graph)),
    
    DegreeVariance = sd(igraph::degree(Graph)),
    
    Components = igraph::components(Graph)$no,
    
    Transitivity = transitivity(Graph),
    
    Density = Prev(get.adjacency(Graph, sparse = F) > 0),
    
    #LouvainModularity = Graph %>% cluster_louvain %>% membership %>% modularity(Graph, .)
    
  ) %>% return
  
}

NetworkTraitYEAR <- list_db_yr %>% map(~.x %>% dplyr::select(Country_search,Country_event) %>% table() %>% 
                       igraph::graph_from_incidence_matrix(directed = TRUE) %>% bipartite_projection %>% 
                       extract2(1) %>% NetworkTraitGet) %>% bind_rows(.id = "ID") %>% mutate_at("ID", as.numeric)

NetworkTraitYEAR$ID <- levels(db$Year_news_f)

ggplot(data=NetworkTraitYEAR, aes(y = MeanDegree, x = ID)) + geom_point()

# Country search ----------------------------------------------------------

# Countries of the Search
Country <- data.frame(table(db$Country_search))
colnames(Country) <- c("Country", "N_news")

# Adding the coordinate of each country
Country <- unique(dplyr::left_join(x  = Country, 
                                   y  = data.frame(Country = db$Country_search, long = db$lon3, lat = db$lat3), 
                                   by ="Country", copy = FALSE)) 

##Calculating country connections
all_pairs <- data.frame(long1 = NA, 
                        long2 = NA, 
                        lat1  = NA, 
                        lat2  = NA, 
                        lwd   = NA, 
                        countrySearch = NA, 
                        country = NA) 

for (i in 1:nlevels(db$Country_search)) {
  
  #select first country
  country_i <- db[db$Country_search == as.character(unique(db$Country_search)[i]),]
  
  #remove potential NAs
  country_i <- country_i[1:table(db$Country_search)[i],]
  
  country_i$Country_event <- droplevels(country_i$Country_event)
  
  country_i <- country_i[country_i$Country_event != as.character(unique(db$Country_search)[i]), ]
  
  country_i <- subset(country_i, !is.na(lon2) & !is.na(lat2))
  
  
  if(nrow(country_i) < 1){
    NULL
  }
  
  else {
    len <- length(table(droplevels(country_i$Country_event)))
    
    all_pairs2 <- data.frame( long1 = rep(country_i$lon3[1], len),
                              long2 = c(unique(country_i$lon2)),
                              lat1 = rep(country_i$lat3[1], len),
                              lat2 = c(unique(country_i$lat2)),
                              lwd = as.numeric(table(droplevels(country_i$Country_event))),
                              countrySearch = rep(as.character(unique(db$Country_search)[i]),len),
                              country= names(table(droplevels(country_i$Country_event)))
    ) 
    
    all_pairs  <- rbind(all_pairs,all_pairs2)
    
  }
}

all_pairs <- na.omit(all_pairs)

# Countries with no news --------------------------------------------------

Botwsana   <- c(23.90,  -21.60)
Iceland    <- c(-19.22,  64.82)
#Montenegro <- c(19.17,   42.83)

NoNews <- rbind(Botwsana,
                Iceland) %>% 
  as.data.frame()

colnames(NoNews) <- c("long","lat")

# Making the map ----------------------------------------------------------

# https://datascience.blog.wzb.eu/2018/05/31/three-ways-of-visualizing-a-graph-on-a-map/

# https://www.r-graph-gallery.com/how-to-draw-connecting-routes-on-map-with-r-and-great-circles.html

world<-map_data("world")

(map2 <- ggplot() +
   geom_map(map = world, data = world,
            aes(long, lat, map_id = region), 
            color = "gray50", fill = "grey70", size = 0.3) +
   
   geom_curve(aes(x = jitter(long1,0.0001), 
                  y = jitter(lat1,0.0001), 
                  xend = jitter(long2, 0.0001), 
                  yend = jitter(lat2, 0.0001),  # draw edges as arcs
                  size = lwd),
              data = all_pairs, curvature = 0.22,
              alpha = 0.2,  color = "orange") +
   
   geom_point(data = Country, 
              aes(x = long, y = lat),
              alpha = 0.7, colour = "black",fill="blue",
              size = range01(sqrt(Country$N_news))*13,
              shape = 21,stroke = 0.8)+
   
   geom_point(data = NoNews, 
              aes(x = long, y = lat),
              alpha = 0.7, colour = "black",fill="red",
              size = 2,
              shape = 21,stroke = 0.8)+
   
   scale_size_continuous("Number of connections:", breaks=c(1,5,10,15))+
   
   theme_map()+
   # theme(legend.position = "bottom",
   #       legend.text = element_text(size = 12),
   #       legend.title = element_text(size = 12),
   theme(
     axis.line=element_blank(),axis.text.x=element_blank(),
     axis.text.y=element_blank(),axis.ticks=element_blank(),
     axis.title.x=element_blank(),
     axis.title.y=element_blank(),legend.position="none",
     panel.background=element_rect(fill = "black", colour = "black"),
     panel.border=element_blank(),panel.grid.major=element_blank(),
     panel.grid.minor=element_blank(),
     plot.background= element_rect(fill = "black", colour = "black"))
 
)

pdf("Figure_conn.pdf", width = 14, height = 8)
map2
dev.off()



