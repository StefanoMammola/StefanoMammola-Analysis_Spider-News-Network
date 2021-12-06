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
library("Amelia")       
library("bipartite")     
library("dplyr")         
library("geosphere")     
library("GGally")        
library("ggalt")         
library("ggplot2")
library("ggraph")
library("ggthemes")      
library("gridExtra")
library("ggregplot")
library("igraph")        
library("lme4")
library("MASS")
library("magrittr")
library("network")       
library("performance")   
library("PupillometryR") 
library("rworldmap")
library("scatterpie")
library("sna")           
library("tidygraph")
library("tidyverse")     

# Loading useful functions ------------------------------------------------

source("Functions/Functions_spider_news.R")

# Loading plot parameters -------------------------------------------------

source("Functions/Plot_parameters.R")

###############################################################

## Data preparation:

###############################################################

# Loading the Database ----------------------------------------------------

db <- read.csv(file = "Data/Data_spider_news_global.csv", sep = '\t', dec = '.', header = TRUE, as.is = FALSE)

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

#Date
db$Year_news   <-  as.Date(paste(db$d,db$m,db$yr,sep="/"),format='%d/%m/%Y')

# Subset databases --------------------------------------------------------

#Database only with distinct event
db_unique_event <- distinct(db, ID_Event, .keep_all = TRUE) 

#Database only with distinct news
db_unique_news  <- distinct(db, ID, .keep_all = TRUE) 

###############################################################

## Summary statistics:

###############################################################

#Number of unique News
nrow(db_unique_news)

#Number of unique EVents
nrow(db_unique_event)

# % of news with error
sum(ifelse(db_unique_news$TotalError > 0 , 1 , 0),na.rm = TRUE)/nrow(db_unique_news)*100
sum(db_unique_news$Taxonomic_error,na.rm = TRUE)/nrow(db_unique_news)*100
sum(db_unique_news$Venom_error,na.rm = TRUE)/nrow(db_unique_news)*100
sum(db_unique_news$Anatomy_error,na.rm = TRUE)/nrow(db_unique_news)*100
sum(as.numeric(db_unique_news$Photo_error),na.rm = TRUE)/nrow(db_unique_news)*100

# % of sensationalistic news
sum(na.omit(db_unique_news$Sensationalism)) / nrow(db) * 100

# % of news with Expert
sum(db_unique_news$TotalExpert, na.rm = TRUE) / nrow(db) * 100
sum(db_unique_news$Expert_arachnologist, na.rm = TRUE) / nrow(db) * 100
sum(db_unique_news$Expert_doctor, na.rm = TRUE) / nrow(db) * 100
sum(db_unique_news$Expert_others, na.rm = TRUE) / nrow(db) * 100

#Number of unique species involved
nlevels(droplevels(db$Species))
as.character(sort(unique(droplevels(db$Species)))) #id of species
sort(table(db$Species))

#Number of country and news country
nlevels(droplevels(db_unique_news$Country_search))  # +2 with no news (Botswana / Iceland)
sort(table(db_unique_news$Country_search))

#Number of lenguages
nlevels(droplevels(db$Lenguage)) 

#Distinct event
table(db_unique_news$TypeEvent)

# Create attribute table for Network analysis -----------------------------

#############################
# At the Country level ######
#############################

# Loading country-level attributes

CountryAttributes <- read.csv(file = "Data/CountryAttributes_Demography.csv", sep = '\t', dec = ',', header = TRUE, as.is = FALSE)

# Summarizing country attributes from the Spider news database
country_attr <- db %>% dplyr::select(Country_search,
                                     Sensationalism,
                                     TotalError,
                                     Lenguage,
                                     lon3,
                                     lat3) %>% group_by(Country_search) %>% 
  summarise(N = length(Country_search),
            Sensationalism = sum(Sensationalism, na.rm = TRUE)/length(Sensationalism), #proportion of sensationalistic news/country
            TotalError = sum(TotalError, na.rm = TRUE)/length(TotalError), #proportion of news with error/country
            Lenguage = names(sort(table(Lenguage), decreasing = TRUE))[1], #main country lenguage
            lon = mean(lon3), #Longitude
            lat = mean(lat3)) #latitude

# All language with less than 10 country become "Others":
country_attr$Lenguage <- ifelse(country_attr$Lenguage %in% names(which(table(country_attr$Lenguage)>9)), country_attr$Lenguage , "Others")

# Add Countries with no News
NoNews <- data.frame(Country_search = c("Botswana", "Iceland"),
                        N = c(0,0),
                        Sensationalism = c(0,0),
                        TotalError = c(0,0),
                        Lenguage = c("English","English"),
                        lon = c(23.90,-19.22),
                        lat = c(-21.60,64.82))

country_attr <- rbind(country_attr,NoNews)

# Check that the two databases match
unique(CountryAttributes$Country %in% country_attr$Country_search) #should be TRUE

# Merge and clean
country_attr <- country_attr %>% dplyr::left_join(CountryAttributes, by = c("Country_search" = "Country")) %>% 
    dplyr::select(Country_search,
                  N,
                  Sensationalism,
                  TotalError,
                  Lenguage,
                  lon, 
                  lat, 
                  ISA, #number of arachnologist in the country.
                  Education_index = Education_Index_2019,  
                  Internet_users  = Internet_users_2018,
                  N_newspapers    = Number_newspapers_Mean1996.2005,
                  Press_Freedom = Press_Freedom_Index,
                  N_Spiders = N_Spiders_WSC,
                  N_Deadly_Spiders = N_Deadly)
                  
rm(CountryAttributes) #clean

#########################
# At the ID_Event level #
#########################

sort(table(db$ID_Event), decreasing = TRUE)

#Extract ID attributes for each ID event
ID_attr <- db %>% filter(ID_Event == levels(ID_Event)[1]) %>%
  arrange(Year_news) %>% droplevels() %>% extractID_attributes #first event

#takes about 1 minutes::
for(i in 2:nlevels(db$ID_Event)) 
  ID_attr <- db %>% filter(ID_Event == levels(ID_Event)[i]) %>% arrange(Year_news) %>% droplevels() %>%
  extractID_attributes %>% bind_rows(ID_attr) #all others

ID_attr %<>% mutate_all(function(x) ifelse(is.infinite(x), NA, x)) #convert infinite to NA

# All genus with less than 50 occurrence become "Others":
ID_attr$Genus <- ifelse(ID_attr$Genus %in% names(which(table(ID_attr$Genus)>50)), ID_attr$Genus, "Others")

#Replace missing data in temporal span (assuming it is 1)
ID_attr$Temporal_span <- ifelse(is.na(ID_attr$Temporal_span) == TRUE, 1, ID_attr$Temporal_span)

rm(i) #clean

###############################################################

## Figure map:

###############################################################

pie_1  <- data.frame(table(db$Country_search,db$Sensationalism)) ; levels(pie_1$Var2) <- c("Non_Sensationalist", "Sensationalist")
radius <- data.frame(log(table(db$Country_search)+1)) ; rownames(radius) <- NULL
n      <- data.frame(table(db$Country_search)) ; rownames(radius) <- NULL

pie <- data.frame(Country = levels(pie_1$Var1),
                  n =  n$Freq,
                  Non_Sensationalist = pie_1$Freq[1:nlevels(pie_1$Var1)],
                  Sensationalist  = pie_1$Freq[c(nlevels(pie_1$Var1)+1):(2*c(nlevels(pie_1$Var1)))],
                  radius      = (radius$Freq)*0.8)

pie <- unique(dplyr::left_join(x  = pie, 
                               y  = data.frame(Country = db$Country_search, lon = db$lon3, lat = db$lat3), 
                               by = "Country", copy = FALSE)) 

net_map <- getMap(resolution="high") %>% 
        ggplot() +
        geom_polygon(aes(long,lat, group=group), colour="grey5", fill="grey5", size = 0.25)+
  ylim(-56.8,90)+ xlim(-180,195)+
  geom_curve(aes(x = jitter(lon3,0.0001), 
                 y = jitter(lat3,0.0001), 
                 xend = jitter(lon, 0.0001), 
                 yend = jitter(lat, 0.0001)),
             data = db, curvature = 0.24,lwd=0.4,
             alpha = 0.05,  color = "aquamarine")+#"#FFEA46FF")+
 
  geom_point(data = db_unique_event, 
               aes(x = lon, y = lat), size = 0.3,
               alpha = 0.9, color="aquamarine",
               stroke = 0.8) +
  ggthemes::theme_map() + theme_map_custom

figure_1 <- net_map + scatterpie::geom_scatterpie(data = pie, aes(x=lon, y= lat, group = Country, r = radius),
                                            cols = c("Non_Sensationalist","Sensationalist"), color = "white", lwd=0.05, alpha=.8) +
    theme(legend.position = c(0.03, 0.2),
          legend.text = element_text(size = 10),
          legend.background = element_rect(fill = "gray40", colour="transparent"))+
    geom_scatterpie_legend(pie$radius,
                           x= -165,
                           y= -45, n = 2,
                           labeller = function (x) x=c(min(pie$n), max(pie$n)))+
    scale_fill_manual("",labels = c("Not Sensationalist","Sensationalist"), values = c("blue","darkmagenta"))

ggplot2::ggsave("Figures/Figure_1.pdf", 
                figure_1, 
                device = cairo_pdf,
                units = "cm",
                width = 20,
                height = 9)

rm(pie, pie_1, n, radius, net_map) #clean

###############################################################

## Analysis #1: What factors explain sensationalism and error?

###############################################################

# Data preparation --------------------------------------------

# Converting factors as factor
db$Sensationalism       <- as.factor(db$Sensationalism)
db$TotalError.01        <- as.factor(ifelse(db$TotalError > 0 , 1 , 0))
db$Expert_arachnologist <- as.factor(db$Expert_arachnologist)
db$Expert_doctor        <- as.factor(db$Expert_doctor)
db$Expert_others        <- as.factor(db$Expert_others)
db$TotalExpert          <- as.factor(db$TotalExpert)
db$Figure_species       <- as.factor(db$Figure_species)
db$Figure_bite          <- as.factor(db$Figure_bite)

# Reneaming low-frequency factors
db$Family2 <- db$Family ; table(db$Family2)
rename <- names(which(table(db$Family2) < 20))
levels(db$Family2)[levels(db$Family2) %in% rename]  <- "Other"

db$Country_search2 <- db$Country_search ; table(db$Country_search2)
rename <- names(which(table(db$Country_search2) < 50))
levels(db$Country_search2)[levels(db$Country_search2) %in% rename]  <- "Other"

db$Lenguage2 <- db$Lenguage ; table(db$Lenguage2)
rename <- names(which(table(db$Lenguage2) < 50))
levels(db$Lenguage2)[levels(db$Lenguage2) %in% rename]  <- "Other"

# Setting baseline for factors
db <- within(db, Family2 <- relevel(Family2, ref = "Other"))
db <- within(db, Country_search2 <- relevel(Country_search2, ref = "Other"))
db <- within(db, Circulation <- relevel(Circulation, ref = "Regional"))
db <- within(db, Type_of_newspaper <- relevel(Type_of_newspaper, ref = "Traditional newspaper"))
db <- within(db, Lenguage2<- relevel(Lenguage2, ref = "Other"))

rm(rename)

# Data exploration --------------------------------------------------------

# Creating a database for models
db_m <- db %>% dplyr::select(
  Sensationalism,
  TotalError.01, 
  Type_of_newspaper, 
  Continent,
  Circulation,
  yr,
  m,
  TypeEvent,
  Figure_species, 
  Figure_bite,
  Expert_doctor, 
  Expert_arachnologist, 
  Expert_others, 
  ID_Event, 
  ID,
  Family2,
  Genus,
  Country_search2,
  Lenguage2,
  lon2,
  lat2) %>% mutate_if(is.numeric, scale)

#Balanced levels of factors: dependent variables
table(db_m$Sensationalism) #ok
table(db_m$TotalError.01)  #ok

# Missing data
Amelia::missmap(db_m)
db_m <- na.omit(db_m) #Omitting missing data

# Final sample size
nrow(db_m)

# Fitting the model: Sensationalism ---------------------------------------

# Fit the model
m1 <- lme4::glmer(Sensationalism ~ yr + Type_of_newspaper + Circulation + TypeEvent + Figure_species + Figure_bite +
                    TotalError.01 + Expert_doctor + Expert_arachnologist + Expert_others +
                    (1|ID_Event) + (1|Genus) + (1|Country_search2) + (1|Lenguage2), 
                  data    = db_m, 
                  family  = binomial(link = "logit"),
                  control = glmerControl(optimizer="bobyqa"))

# Check model
performance::check_model(m1, check = c("vif", "reqq"))

# Model residuals vs temporal and spatial factors
R_m1 <- residuals(m1)

par(mfrow = c(2,2),mar = c(2,2,2,2))
plot(db_m$lon2, R_m1, xlab = "Longitude", ylab = "Residuals",  main = "Residuals vs Longitude")
plot(db_m$lat2, R_m1, xlab = "Latitude",  ylab = "Residuals",  main = "Residuals vs Latitude")
boxplot(R_m1 ~ as.factor(db_m$m), xlab = NULL, ylab = "Residuals",  main = "Residuals vs month")
boxplot(R_m1 ~ as.factor(db_m$yr), xlab = NULL,  ylab = "Residuals",  main = "Residuals vs year")

# Interpret the model
(fit <- parameters::model_parameters(m1))
performance::r2(m1)
fit$Parameter

# color for plot 
color_p_value <- ifelse(round(fit$p[2:14],3) > 0.05, 
                        "grey5", 
                        ifelse(fit$Coefficient[2:14] > 0, "darkmagenta", "blue"))
                        
# Plot
(plot_model1 <- sjPlot::plot_model(m1, 
                                   title = title_sjPlot_m1,
                                   sort.est = FALSE,  
                                   vline.color = "grey80",
                                   colors = "grey5",#rev(color_p_value),
                                   axis.title = xlab_sjPlot_m1_m2, #"Odds ratios ",
                                   axis.labels = rev(axis_labels_sjPlot_m1), #rev cause it traspose x and y axes
                                   show.values = FALSE, 
                                   value.offset = .3, 
                                   se = TRUE, 
                                   show.p = FALSE) + ylim(-.8,1.4) + theme_custom())

rm(fit,R_m1) #clean

# Fitting the model: Errors ---------------------------------------
m2 <- lme4::glmer(TotalError.01 ~ yr + Type_of_newspaper + Circulation + TypeEvent +
                      Sensationalism + Expert_doctor + Expert_arachnologist + Expert_others +
                    (1|ID_Event) + (1|Genus) + (1|Lenguage2) + (1|Country_search2), 
                  data = db_m, 
                  family = binomial(link = "logit"),
                  control = glmerControl(optimizer="bobyqa"))

# Check model
performance::check_model(m2, check = c("vif", "reqq"))

# Model residuals vs temporal and spatial factors
R_m2 <- residuals(m2)

par(mfrow = c(2,2),mar = c(2,2,2,2))
plot(db_m$lon2, R_m2, xlab = "Longitude", ylab = "Residuals",  main = "Residuals vs Longitude")
plot(db_m$lat2, R_m2, xlab = "Latitude",  ylab = "Residuals",  main = "Residuals vs Latitude")
boxplot(R_m2 ~ as.factor(db_m$m), xlab = NULL, ylab = "Residuals",  main = "Residuals vs month")
boxplot(R_m2 ~ as.factor(db_m$yr), xlab = NULL,  ylab = "Residuals",  main = "Residuals vs year")

# Interpret the model
(fit <- parameters::model_parameters(m2))
performance::r2(m2)

# Plot
(plot_model2 <- sjPlot::plot_model(m2, 
                                   title = title_sjPlot_m2,
                                   sort.est = FALSE,  
                                   vline.color = "grey80",
                                   color = "grey5",
                                   axis.title = xlab_sjPlot_m1_m2, 
                                   axis.labels = rev(axis_labels_sjPlot_m2), 
                                   show.values = FALSE, 
                                   value.offset = .3, 
                                   se = TRUE, 
                                   show.p = FALSE) + ylim(-.8,1.4) +theme_custom())

rm(fit,R_m1) #clean

# Figure 2 ----------------------------------------------------------------

pdf(file = "Figures/Figure_2.pdf", width = 12, height =5)
ggpubr::ggarrange(plot_model1,plot_model2,
                  common.legend = FALSE,
                  hjust = -4,
                  align = "v",
                  labels = c("A","B"),
                  ncol=2)
dev.off()

###############################################################

## Analysis #2: Network analysis

###############################################################

# Constructing a network -----------------------------------------
Graph_bipartite <- db %>% 
            dplyr::select(Country_search,ID_Event) %>% 
            table() %>% 
            igraph::graph_from_incidence_matrix(directed = TRUE) %>% 
            tidygraph::as_tbl_graph(directed = TRUE) 
 
# Collapse it into an unipartite 
Graph_unipartite_full <- igraph::bipartite_projection(Graph_bipartite)

# Analysing country-level network properties ------------------------------

# Question 1a: (Node-level)
# Why are some countries more influential than others in generating 
# spider-related content taken up in  other places?

# Question 1b: (Edge-level)
#Why determines the existence of specific connection (news flow) 
# among countries?

#--------------------------------------------------------------------------

# Takes the unipartite project graph
Graph_unipartite <- Graph_unipartite_full$proj1  %>% as_tbl_graph(directed = TRUE) %>% 
  activate(edges) %>% #%>% mutate(weight = 1) 
  igraph::simplify(edge.attr.comb = "sum") %>% 
  as_tbl_graph

# Summary stats
NetworkTraitGet(Graph_unipartite)

# Calculating node level traits 
node_trait <- Graph_unipartite %>% 
                   NodeTraitGet(mode = "in", dir = FALSE) %>% bind_cols()
node_trait

# Combine them with the Country attributes
country_attr <- country_attr %>% dplyr::left_join(node_trait, by = c("Country_search" = "ID"))
country_attr[80:81,15:19] <- 0 #give values of zero for Network properties for Iceland and Botswana

# Attach them to the Graph_tbl_uni
Graph_unipartite <- Graph_unipartite %>% tidygraph::activate(nodes) %>% 
  left_join(country_attr, by = c("name" = "Country_search"))

########################
# Plotting the network #
########################

# Plot the network [comparing three layouts]
SpatialLayout <- country_attr[1:79,] %>% dplyr::select(lon,lat) %>% as.matrix #geo-coordinates

#SpatialLayout
(plot_network <- Graph_unipartite %>% 
  igraph::simplify(edge.attr.comb = "sum") %>% 
  ggraph::ggraph(SpatialLayout) +
  #geom_edge_density(fill="orange", alpha=1) +
  geom_edge_fan(aes(width=weight),color="gray60", alpha=0.1) +
  geom_node_point(col="grey30", alpha = .8, 
                  aes(size=N,fill=Lenguage), shape = 21) + 
  geom_node_text(aes(label = name), size=2, color="gray10", repel=TRUE) +
  scale_fill_manual(values = c("blue", "orange", "turquoise","purple", "grey15")) +
  theme_void() + theme(legend.position = "top",legend.direction = "vertical") + coord_fixed())

###################################
# Regression Models #
###################################

db_m3_cat <- country_attr %>% dplyr::select(Degree, Country_search, Lenguage) 
db_m3_con <- country_attr %>% dplyr::select(N, ISA, Sensationalism, TotalError, Education_index,
                                            Internet_users,N_newspapers, N_Spiders, N_Deadly_Spiders) 

db_m3_con <- data.frame(db_m3_con)

# Check distributon
par(mfrow = c(3,3),mar = c(2,2,2,2))
for(i in 1:ncol(db_m3_con))
  dotchart(as.numeric(db_m3_con[,i]), main = colnames(db_m3_con)[i])

#log transform
db_m3_con$N <- log(db_m3_con$N+1)
db_m3_con$ISA <- log(db_m3_con$ISA+1)
db_m3_con$N_newspapers <- log(db_m3_con$N_newspapers+1)

# Check collinearity
psych::pairs.panels(db_m3_con) 

#Education index and internewt users 
#Number of newspaper with N 
#ISA and N

#Check balance of levels
db_m3_cat$Lenguage %>% table

# scale all var
db_m3_con <- db_m3_con %>% apply(MARGIN = 2, scale) #scale cont variable

#merge
db_m3 <- data.frame(db_m3_cat, db_m3_con) ; rm(db_m3_cat, db_m3_con)

# Remove missing values
db_m3_noNA <- na.omit(db_m3)

# Set reference level
db_m3_noNA$Lenguage <- as.factor(db_m3_noNA$Lenguage)

# Fit the model -----------------------------------------------------------

# # Set formula
# formula_m3 <- as.formula("Degree ~ TotalError + Sensationalism + Internet_users + N_Spiders + (1|Lenguage)")
# 
# # Fit the model
# m3 <- lme4::glmer(formula_m3, 
#             data = db_m3_noNA, 
#             family = poisson(link = "log"),
#             control = glmerControl(optimizer="bobyqa"))
# 
# # Check model
# performance::check_overdispersion(m3) #Overdispersed
# 
# # Switch to negative binomial
# m3 <- lme4::glmer.nb(formula_m3, data = db_m3_noNA)
# 
# summary(m3)
# performance::check_model(m3)
# performance::r2(m3)
# 
# (plot_model3 <- sjPlot::plot_model(m3))

# Set formula
formula_m3 <- as.formula("Degree ~ Lenguage + Sensationalism + TotalError + Internet_users + N_Deadly_Spiders")

# Fit the model
m3 <- glm(formula_m3, family = "poisson", data = db_m3_noNA)

# Check model
performance::check_overdispersion(m3) #Overdispersed

# Switch to negative binomial
m3 <- MASS::glm.nb(formula_m3, data = db_m3_noNA)

# Check model
performance::check_model(m3)

# Interpret the model
parameters::model_parameters(m3)
performance::r2(m3)

# Plot
(plot_model3 <- sjPlot::plot_model(m3, 
                                   title = title_sjPlot_m3,
                                   sort.est = FALSE,  
                                   vline.color = "grey80",
                                   color = "grey5",
                                   axis.title = xlab_sjPlot_m3, 
                                   #axis.labels = rev(axis_labels_sjPlot_m3), 
                                   show.values = TRUE, 
                                   value.offset = .3, 
                                   se = TRUE, 
                                   show.p = FALSE) + theme_custom())

###################################
# Exponential Random Graph Models #
###################################

# Get the adjacency matrix
AdjMatrix <- Graph_unipartite %>% get.adjacency(attr = "weight", sparse = FALSE)

AdjMatrix[AdjMatrix>0] <- 1  

dim(AdjMatrix) #check if it's a square

#Assign row and col names
colnames(AdjMatrix) <- rownames(AdjMatrix) <- Graph_unipartite %>% activate(nodes) %>% pull(name)

#Response variables
ResponseNetwork <- AdjMatrix %>% as.matrix %>% network::network(directed = FALSE)

#Adding node-level attributes
ResponseNetwork %v% "N"              <- Graph_unipartite %>% as.data.frame %>% dplyr::mutate(N = log(N+1)) %>% pull(N)
ResponseNetwork %v% "Lenguage"       <- Graph_unipartite %>% as.data.frame %>% pull(Lenguage) # should be as.character
ResponseNetwork %v% "Sensationalism" <- Graph_unipartite %>% as.data.frame %>% pull(Sensationalism)
ResponseNetwork %v% "TotalError"     <- Graph_unipartite %>% as.data.frame %>% pull(TotalError)

#Data exploration

# Collinearity 
psych::pairs.panels(Graph_unipartite %>% as.data.frame %>% select(N, Sensationalism, TotalError))

# Check association between Lenguage and continous var
Graph_unipartite %>% as.data.frame %>% ggplot(aes(x=Lenguage,y=lon)) + geom_boxplot() + theme_custom()
Graph_unipartite %>% as.data.frame %>% ggplot(aes(x=Lenguage,y=lat)) + geom_boxplot() + theme_custom()
Graph_unipartite %>% as.data.frame %>% ggplot(aes(x=Lenguage,y=Sensationalism)) + geom_boxplot() + theme_custom()
Graph_unipartite %>% as.data.frame %>% ggplot(aes(x=Lenguage,y=TotalError)) + geom_boxplot() + theme_custom()
Graph_unipartite %>% as.data.frame %>% ggplot(aes(x=Lenguage,y=N)) + geom_boxplot() + theme_custom()

# Model fit
ergm0 <- ergm::ergm(ResponseNetwork ~ edges, estimate = "MLE")
ergm1 <- ergm::ergm(ResponseNetwork ~ edges + edgecov(SpatialLayout), estimate = "MLE")
ergm2 <- ergm::ergm(ResponseNetwork ~ edges + edgecov(SpatialLayout) + nodecov("Sensationalism"), estimate = "MLE")
ergm3 <- ergm::ergm(ResponseNetwork ~ edges + edgecov(SpatialLayout) + nodecov("Sensationalism") + nodecov("TotalError"), estimate = "MLE")
ergm4 <- ergm::ergm(ResponseNetwork ~ edges + edgecov(SpatialLayout) + nodecov("Sensationalism") + nodecov("TotalError") + nodecov("N"), estimate = "MLE")
ergm5 <- ergm::ergm(ResponseNetwork ~ edges + nodematch("Lenguage", diff = FALSE) + nodecov("Sensationalism") + nodecov("TotalError") + nodecov("N"), estimate = "MLE")
#ergm6 <- ergm::ergm(ResponseNetwork ~ edges + nodematch("Lenguage", diff = F) + nodefactor("Lenguage") + nodecov("Sensationalism") + nodecov("TotalError") + nodecov("N"), estimate = "MLE")

list(ergm0, ergm1, ergm2, ergm3, ergm4, ergm5) %>% map_dbl(BIC) %>% plot

ergm_BIC <- ergm5 ; rm(ergm0, ergm1, ergm2, ergm3, ergm4, ergm5)

# Model validation

# Simulating from the model

BlankNetwork   <- ResponseNetwork ; BlankNetwork[] <- 0 #generate an empty network with the same dim as ResponseNetwork 
NSims <- 1000 #N simulation

Sims <- stats::simulate(ergm_BIC, nsim = NSims, basis = BlankNetwork)

ResponseNetwork %>% as.matrix %>% Prev # average probability of having an edge connecting two nodes

Sims %>% map_dbl(~.x %>% as.matrix %>% c %>% Prev) %>% 
  qplot() + geom_vline(xintercept = Prev(as.matrix(ResponseNetwork)), colour = "red") + theme_custom()

qplot(1:NSims, Sims %>% map_dbl(~.x %>% as.matrix %>% c %>% Prev)) + ylab("prevalence") + xlab("simulation ID") +
  geom_hline(yintercept = Prev(as.matrix(ResponseNetwork)), colour = "red") + theme_custom()

rm(BlankNetwork, NSims)#clean

# Interpret the model

EstimateDF <- 
  ergm_BIC %>% summary %>% extract2("coefficients") %>% 
  as.data.frame %>% 
  rownames_to_column("Variable")

ergm_BIC %>% 
  confint %>% # Takes the 95% confidence intervals
  as.data.frame %>% 
  rename(Lower = 1, Upper = 2)

EstimateDF %<>% # Bind them together
  bind_cols(ergm_BIC %>% confint %>% as.data.frame %>% 
              rename(Lower = 1, Upper = 2))

EstimateDF %>% head

EstimateDF$Variable <- axis_labels_ergm1

(plot_ergm1 <- EstimateDF %>% ggplot2::ggplot(aes(Variable, Estimate)) +
  geom_hline(lty = 1, col = "grey60", yintercept = 0) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0) +
  labs(y     = xlab_ergm1_ergm2, 
       title = title_ergm1) +
      geom_point() + theme_custom() + theme(axis.title.y=element_blank()) + coord_flip())

# Analysing ID_Event-level network properties -----------------------------

# Question 1a: (Node-level)

# Question 1b: (Edge-level)

#--------------------------------------------------------------------------

# Takes the unipartite project graph
Graph_unipartite2 <- Graph_unipartite_full$proj2 %>% as_tbl_graph(directed = FALSE) %>% 
  activate(edges) %>% #%>% mutate(weight = 1) 
  igraph::simplify(edge.attr.comb = "sum") %>% 
  as_tbl_graph#note that it loses the directionality

# Summary stats
NetworkTraitGet(Graph_unipartite2)
# Size Diameter MeanDegree DegreeVariance Components Transitivity   Density LouvainModularity
# 2644       16   271.1135       250.9714         10    0.7854129 0.1025391         0.5366639

# Calculating node level traits (warning: takes about 1 minute)
node_trait2 <- Graph_unipartite2 %>% 
  NodeTraitGet(mode = "in", dir = FALSE) %>% bind_cols()

# Combine them with the ID_attr
ID_attr <- ID_attr %>% left_join(node_trait2, by = c("ID" = "ID"))

# Attach them to the Graph_tbl_uni
Graph_unipartite2 <- Graph_unipartite2 %>% tidygraph::activate(nodes) %>% left_join(ID_attr, by = c("name" = "ID"))

###################################
# Regression Models #
###################################
db_m3_cat <- ID_attr %>% dplyr::select(Degree, Country_event, Genus) 
db_m3_con <- ID_attr %>% dplyr::select(n, Temporal_span, Sensationalism, Distance_information, Bite, Death) %>% apply(MARGIN = 2, scale)

db_m3 <- data.frame(db_m3_cat, db_m3_con) ; rm(db_m3_cat, db_m3_con)

m4 <- lme4::glmer(Degree ~ n + Bite + Death + Temporal_span + Sensationalism +
                        Genus + (1|Country_event),
                  family = "poisson",
                  data = db_m3,
                  control = glmerControl(optimizer="bobyqa")
)

performance::check_overdispersion(m4)

m4 <- lme4::glmer.nb(Degree ~ n + Bite + Death + Temporal_span + Sensationalism + Genus + (1|Country_event),
                  data = db_m3,
                  control = glmerControl(optimizer="bobyqa")
)

summary(m4)

sjPlot::plot_model(m4, 
                   #title = title_sjPlot_m2,
                   sort.est = FALSE,  
                   vline.color = "grey80",
                   color = "grey5",
                   #axis.title = xlab_sjPlot_m1_m2, 
                   #axis.labels = rev(axis_labels_sjPlot), 
                   show.values = TRUE, 
                   value.offset = .3, 
                   se = TRUE, 
                   show.p = TRUE) + theme_custom()
)


###################################
# Exponential Random Graph Models #
###################################

#Convert the Adjacency matrix to incidence
AdjMatrix2 <- Graph_unipartite2 %>% get.adjacency(sparse = FALSE)

AdjMatrix2[AdjMatrix2>0] <- 1  

dim(AdjMatrix2) #check it's a square

#Assign row and col names
colnames(AdjMatrix2) <- rownames(AdjMatrix2) <- Graph_unipartite2 %>% activate(nodes) %>% pull(name)

#Response variables
ResponseNetwork2 <- AdjMatrix2 %>% as.matrix %>% network::network(directed = FALSE)

#Adding node-level attributes
Graph_unipartite2 %>% as.data.frame %>% colnames

ResponseNetwork2 %v% "year" <- Graph_unipartite2 %>% as.data.frame %>% pull(year)
ResponseNetwork2 %v% "Bite" <- Graph_unipartite2 %>% as.data.frame %>% pull(Bite) 
ResponseNetwork2 %v% "Death" <- Graph_unipartite2 %>% as.data.frame %>% pull(Death)
ResponseNetwork2 %v% "Genus" <- Graph_unipartite2 %>% as.data.frame %>% pull(Genus)
ResponseNetwork2 %v% "Sensationalism" <- Graph_unipartite2 %>% as.data.frame %>% pull(Sensationalism)
ResponseNetwork2 %v% "Temporal_span" <- Graph_unipartite2 %>% as.data.frame %>% pull(Temporal_span)
ResponseNetwork2 %v% "Distance_information" <- Graph_unipartite2 %>% as.data.frame %>% pull(Distance_information)

Spatial2 <- ID_attr %>% dplyr::select(lon,lat) %>% as.matrix

#Data exploration

# Collinearity 
psych::pairs.panels(Graph_unipartite2 %>% as.data.frame %>% select(year, Bite, Death, Sensationalism, Temporal_span, Distance_information))

# Model fit
ergm0 <- ergm::ergm(ResponseNetwork2 ~ edges, estimate = "MLE")
ergm1 <- ergm::ergm(ResponseNetwork2 ~ edges + nodefactor("Genus"), estimate = "MLE")
ergm2 <- ergm::ergm(ResponseNetwork2 ~ edges + nodefactor("Genus") + nodecov("Sensationalism"), estimate = "MLE")
ergm3 <- ergm::ergm(ResponseNetwork2 ~ edges + nodefactor("Genus") + nodecov("Sensationalism") + nodecov("Temporal_span"), estimate = "MLE")
ergm4 <- ergm::ergm(ResponseNetwork2 ~ edges + nodefactor("Genus") + nodecov("Sensationalism") + nodecov("Temporal_span") + nodecov("Distance_information"), estimate = "MLE")
ergm5 <- ergm::ergm(ResponseNetwork2 ~ edges + nodefactor("Genus") + nodematch("Genus") + nodecov("Sensationalism") + nodecov("Temporal_span") + nodecov("Distance_information"), estimate = "MLE")

list(ergm0, ergm1, ergm2, ergm3, ergm4, ergm5) %>% map_dbl(BIC) %>% plot

ergm_BIC2 <- ergm1 # best model

# Model validation

# Simulating from the model

BlankNetwork   <- ResponseNetwork2 ; BlankNetwork[] <- 0 #generate an empty network with the same dim as ResponseNetwork 
NSims <- 1000 #N simulation

Sims <- stats::simulate(ergm_BIC2, nsim = NSims, basis = BlankNetwork)

ResponseNetwork2 %>% as.matrix %>% Prev # average probability of having an edge connecting two nodes

Sims %>% map_dbl(~.x %>% as.matrix %>% c %>% Prev) %>% 
  qplot() + geom_vline(xintercept = Prev(as.matrix(ResponseNetwork2)), colour = "red") + theme_custom()

qplot(1:NSims, Sims %>% map_dbl(~.x %>% as.matrix %>% c %>% Prev)) + ylab("prevalence") + xlab("simulation ID") +
  geom_hline(yintercept = Prev(as.matrix(ResponseNetwork2)), colour = "red") + theme_custom()

# Interpret the model
summary(ergm_BIC2)

EstimateDF2 <- 
  ergm_BIC2 %>% summary %>% extract2("coefficients") %>% 
  as.data.frame %>% 
  rownames_to_column("Variable")

ergm_BIC2 %>% 
  confint %>% # Takes the 95% confidence intervals
  as.data.frame %>% 
  rename(Lower = 1, Upper = 2)

EstimateDF2 %<>% # Bind them together
  bind_cols(ergm_BIC2 %>% confint %>% as.data.frame %>% 
              rename(Lower = 1, Upper = 2))

EstimateDF2 %>% head

EstimateDF2 %>% ggplot2::ggplot(aes(Variable, Estimate)) +
  geom_hline(lty = 1, col = "grey60", yintercept = 0) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  geom_point() + theme_custom() + coord_flip()

#########################################

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
                        lat1 = NA, 
                        lat2 = NA, 
                        lwd = NA, 
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

# Calculating measures of Centrality --------------------------------------

# https://medium.com/@615162020004/social-network-analysis-with-r-centrality-measure-86d7fa273574

Net_2 <- igraph::graph.data.frame(Network_matrix_igraph, directed=TRUE)
plot(Net_2,edge.arrow.size=0.5,edge.arrow.mode = "-")

# Degree Centrality

# Definition: Degree centrality assigns an importance score based purely on the number of links held by each node.
# What it tells us: How many direct, ‘one hop’ connections each node has to other nodes within the network.
# When to use it: For finding very connected individuals, popular individuals, individuals who are likely to hold most information or individuals who can quickly connect with the wider network.

Centrality_Degree <- data.frame(igraph::degree(Net_2, mode = "out"))
Centrality_Degree <- data.frame(country = as.character(rownames(Centrality_Degree)),
                                n = Centrality_Degree[,1]) %>% arrange(country) 

# Betweenness Centrality

# Definition: Betweenness centrality measures the number of times a node lies on the shortest path between other nodes.
# What it tells us: This measure shows which nodes act as ‘bridges’ between nodes in a network. It does this by identifying all the shortest paths and then counting how many times each node falls on one.
# When to use it: For finding the individuals who influence the flow around a system.

Centrality_Betweeness <- data.frame(igraph::betweenness(Net_2))
Centrality_Betweeness <- data.frame(country = as.character(rownames(Centrality_Betweeness)),
                                    n = Centrality_Betweeness[,1]) %>% arrange(country) 

# Closeness Centrality

# Definition: This measure scores each node based on their ‘closeness’ to all other nodes within the network.
# What it tells us: This measure calculates the shortest paths between all nodes, then assigns each node a score based on its sum of shortest paths.
# When to use it: For finding the individuals who are best placed to influence the entire network most quickly.

Centrality_Closeness <- data.frame(igraph::closeness(Net_2))
Centrality_Closeness <- data.frame(country = as.character(rownames(Centrality_Closeness)),
                                   n = Centrality_Closeness[,1]) %>% arrange(country) 

# Storing & cleaning

Centrality <- data.frame(country     = Centrality_Degree$country,
                         Cdegree     = Centrality_Degree$n,
                         Cbetween    = Centrality_Betweeness$n,
                         Cclosenness = Centrality_Closeness$n)

rm(Centrality_Degree,Centrality_Betweeness,Centrality_Closeness,Net_2)

# Creating a connection plot ----------------------------------------------

# Countries of the Search
Country <- data.frame(table(db$Country_search),Centrality[,2:4])
colnames(Country) <- c("Country", "N_news","Cdegree","Cbetween","Cclosenness")

# Adding the coordinate of each country
Country <- unique(dplyr::left_join(x = Country, 
                                   y = data.frame(Country = db$Country_search, long = db$lon3, lat = db$lat3), 
                                   by = "Country", copy = FALSE)) 


##Calculating country connections
all_pairs <- data.frame(long1 = NA, 
                        long2 = NA, 
                        lat1 = NA, 
                        lat2 = NA, 
                        lwd = NA, 
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
Montenegro <- c(19.17,   42.83)

NoNews <- rbind(Botwsana,
                Iceland,
                Montenegro) %>% 
  as.data.frame()

colnames(NoNews) <- c("long","lat")

# Making the map ----------------------------------------------------------

# https://datascience.blog.wzb.eu/2018/05/31/three-ways-of-visualizing-a-graph-on-a-map/

# https://www.r-graph-gallery.com/how-to-draw-connecting-routes-on-map-with-r-and-great-circles.html

(map2 <- ggplot() +
   geom_map(map = world, data = world,
            aes(long, lat, map_id = region), 
            color = "gray50", fill = "grey70", size = 0.3) +
   
   geom_curve(aes(x = jitter(long1,0.0001), 
                  y = jitter(lat1,0.0001), 
                  xend = jitter(long2, 0.0001), 
                  yend = jitter(lat2, 0.0001),    # draw edges as arcs
                  size = lwd),
              data = all_pairs, curvature = 0.34,
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



# Wordcloud ---------------------------------------------------------------

library(tidytext)
library(R.utils)
library(wordcloud)

db_english <- db[db$Lenguage == "English",]

# A list of boring and non-useful words, bundled with `tidytext`
data(stop_words)
Title_sensationalistic <- db_english[db_english$Sensationalism == 1,] %>%
  mutate(title = as.character(Title)) %>%
  unnest_tokens(output = title_word,
                input = Title) %>%
  anti_join(stop_words, by = c("title_word" = "word")) %>%
  count(title_word, sort = TRUE) 

Title_sensationalistic[2:60,] %>% with(wordcloud(words = title_word, 
                                                 freq = n, 
                                                 max.words = 200,
                                                 scale=c(4,.2),
                                                 random.color=TRUE, color = c("orange","orange","darkorange","black","black")))


Title_non <- db_english[db_english$Sensationalism == 0,] %>%
  mutate(title = as.character(Title)) %>%
  unnest_tokens(output = title_word,
                input = Title) %>%
  anti_join(stop_words, by = c("title_word" = "word")) %>%
  count(title_word, sort = TRUE) 

Title_non[2:60,] %>% with(wordcloud(words = title_word, 
                                    freq = n, 
                                    max.words = 200,
                                    scale=c(4,.2),
                                    random.color=TRUE, color = c("aquamarine3","aquamarine4","darkblue","black")))






## https://cedricscherer.netlify.app/2019/05/17/the-evolution-of-a-ggplot-ep.-1/

## load fonts
font_add_google("Poppins", "Poppins")
font_add_google("Roboto Mono", "Roboto Mono")
showtext_auto()




ggplot(db, aes(x =TotalError , y =Country_search)) +
  stat_summary(fun = mean, geom = "point", size = 5) +
  geom_jitter(size = 2, alpha = 0.25, width = 0.2) +
  theme_light(base_size = 15) +
  labs(caption = "Data: UNESCO Institute for Statistics")+
  
  
  geom_curve(
    # data = arrows, aes(x = x1, xend = x2,
    #                    y = y1, yend = y2),
    arrow = arrow(length = unit(0.08, "inch")), size = 0.5,
    color = "gray20", curvature = -0.3#
  ) +
  
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    plot.caption = element_text(size = 11, color = "gray50"),
    panel.grid = element_blank()
  )





scale_fill_manual(labels=c("  Latrodectus tredecimguttatus", "  Loxosceles rufescens"),values=c("black", "#E69F00"))+
  scale_x_continuous(breaks = c(2010:2019), labels = as.character(2010:2019))+ 
  labs(title="Annual distribution of media reports",x=NULL, y = "Density")+
  
  #annotate("text",label="Fred Vargas' book\n'Quand sort la recluse'",hjust = 0,x=2015,y=0.4,size=4)+
  #annotate("segment", x = 2017, xend = 2017, y = 0.25, yend = 0.35, colour = "grey20", size=0.8)+
  theme_classic()+
  theme(legend.position = c(0.2, 0.7),
        legend.text = element_text(size=10,face = "italic"),
        legend.title = element_text(size=0))








table(db_unique_ID$Species) # 5 Others; 22 Latrodectus tredecimguttatus;  71 Loxosceles rufescens 
table(db_unique_ID$Bite)  # 69 bites
table(db_unique_ID$Death) # 3 deaths

Death_db=db_unique_ID[db_unique_ID$Death== "Yes",]

# merging poorly represented species in the category "others"

db$Species2 = db$Species
levels(db$Species2) = c("Others","Others","Latrodectus tredecimguttatus","Loxosceles rufescens","Others")

## Temporal trends

db_ll <- db[db$Species == "Loxosceles rufescens" | db$Species == "Latrodectus tredecimguttatus",]
db_ll$Species= droplevels(db_ll$Species)

p1 = ggplot(db_ll[db_ll$yr<2020,], aes(x=yr,fill=Species))+
  geom_density(alpha=0.6,adjust=1.5)+
  scale_fill_manual(labels=c("  Latrodectus tredecimguttatus", "  Loxosceles rufescens"),values=c("black", "#E69F00"))+
  scale_x_continuous(breaks = c(2010:2019), labels = as.character(2010:2019))+ 
  labs(title="Annual distribution of media reports",x=NULL, y = "Density")+
  
  #annotate("text",label="Fred Vargas' book\n'Quand sort la recluse'",hjust = 0,x=2015,y=0.4,size=4)+
  #annotate("segment", x = 2017, xend = 2017, y = 0.25, yend = 0.35, colour = "grey20", size=0.8)+
  theme_classic()+
  theme(legend.position = c(0.2, 0.7),
        legend.text = element_text(size=10,face = "italic"),
        legend.title = element_text(size=0))

p2 = ggplot(db_ll, aes(x=m, fill=Species))+
  geom_density(alpha=0.6)+
  scale_fill_manual(labels=c("  Latrodectus tredecimguttatus", "  Loxosceles rufescens"),values=c("black", "#E69F00"))+
  scale_x_continuous(breaks = 1:12, labels = array(month.abb))+ 
  labs(title="Seasonal distribution of media reports",x=NULL, y = "Density")+
  theme_classic()+
  theme(legend.position = "none")

pdf("/Users/stefanomammola/Desktop/Figure_temporal.pdf", width = 14, height = 5)

grid.arrange(p1,p2,ncol=2)

dev.off()

###### CONTENT OF REPORTS

CONT1= data.frame(table(db_unique_ID$Species,db_unique_ID$Bite))
DEATH1=data.frame(table(db_unique_ID$Species,db_unique_ID$Death))

CONT1[6,3] <- 59
levels(CONT1$Var2) = c("Encounter","Bite")

DEATH1=data.frame(Var1=DEATH1[1:3,1],Var2=rep("Deadly bite",3),Freq=c(0,0,3))

CONT1 = rbind(CONT1,DEATH1)

c1 = ggplot(CONT1, aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity",color="black",position=position_dodge(), alpha=1)+
  geom_text(aes(label=Freq), vjust=-1, color="black",
            position = position_dodge(0.9), size=3.5)+
  
  scale_fill_manual("",labels=c("Encounter", "Bite", "Deadly bite"),values=c("grey20", "deepskyblue4","orange"))+
  
  labs(title="Type of event",subtitle =  "[Only unique event are counted]",x=NULL, y = "Frequency")+
  theme_classic()+
  theme(legend.position = c(0.15, 0.5),
        legend.text = element_text(size=10,face = "plain"),
        legend.title = element_text(size=12,face="bold"), 
        axis.text.x = element_text(size=12,face = c("plain","italic","italic")))


c1


CONT2= data.frame(table(db$Species,db$Figure_species))

c2 = ggplot(CONT2, aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity",color="black",position=position_dodge(), alpha=1)+
  geom_text(aes(label=Freq), vjust=-1, color="black",
            position = position_dodge(0.9), size=3.5)+
  
  scale_fill_manual("",labels=c("Without", "With"),values=c("grey20", "deepskyblue3"))+
  
  labs(title="Photography of the species",x=NULL, y = "Frequency")+
  theme_classic()+
  theme(legend.position = c(0.15, 0.5),
        legend.text = element_text(size=10,face = "plain"),
        legend.title = element_text(size=12,face="bold"), 
        axis.text.x = element_text(size=12,face = c("plain","italic","italic")))

c2

CONT3= data.frame(table(db$Species2,db$Figure_bite))

c3 = ggplot(CONT3, aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity",color="black",position=position_dodge(), alpha=1)+
  geom_text(aes(label=Freq), vjust=-1, color="black",
            position = position_dodge(0.9), size=3.5)+
  
  scale_fill_manual("",labels=c("Without", "With"),values=c("grey20", "deepskyblue3"))+
  
  labs(title="Photographs of the bite",x=NULL, y = "Frequency")+
  theme_classic()+
  theme(legend.position = c(0.15, 0.5),
        legend.text = element_text(size=10,face = "plain"),
        legend.title = element_text(size=12,face="bold"), 
        axis.text.x = element_text(size=12,face = c("plain","italic","italic")))

c3

CONT4= data.frame(table(db$Species2,db$Expert))

c4 = ggplot(CONT4, aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity",color="black",position=position_dodge(), alpha=1)+
  geom_text(aes(label=Freq), vjust=-1, color="black",
            position = position_dodge(0.9), size=3.5)+
  
  scale_fill_manual("",labels=c("No", "Yes"),values=c("grey20", "deepskyblue3"))+
  
  labs(title="Expert consulted",x=NULL, y = "Frequency")+
  theme_classic()+
  theme(legend.position = c(0.15, 0.5),
        legend.text = element_text(size=10,face = "plain"),
        legend.title = element_text(size=12,face="bold"), 
        axis.text.x = element_text(size=12,face = c("plain","italic","italic")))

c4

CONT5= data.frame(table(db$Species2,db$Sensationalism_consensus))

c5 = ggplot(CONT5, aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity",color="black",position=position_dodge(), alpha=1)+
  geom_text(aes(label=Freq), vjust=-1, color="black",
            position = position_dodge(0.9), size=3.5)+
  
  scale_fill_manual("",labels=c("No", "Yes"),values=c("grey20", "deepskyblue3"))+
  
  labs(title="Sensationalistic media report",x=NULL, y = "Frequency")+
  theme_classic()+
  theme(legend.position = c(0.15, 0.5),
        legend.text = element_text(size=10,face = "plain"),
        legend.title = element_text(size=12,face="bold"), 
        axis.text.x = element_text(size=12,face = c("plain","italic","italic")))

c5
str(db)

## Facebook share
c6 = ggplot(data=db, aes(y=log(Share_FB+1), x=Species2))+
  labs(title="Share on social media", subtitle = "[Expressed as the number of shares on Facebook; Nanni et al., 2020]",x=NULL, y = "N° of shares [logarithm]")+
  geom_boxplot(fill=c("darkseagreen2", "grey60", "#E69F00"),alpha=0.8)+ 
  geom_jitter(alpha=0.6,col="grey20",shape=16, position=position_jitter(0.2))+
  theme_classic()+
  theme(axis.text.x = element_text(size=12,face = c("plain","italic","italic")))

c6

pdf("/Users/stefanomammola/Desktop/Fiugre_content.pdf", width = 21, height = 15.5)

grid.arrange(c1,c6,c2,c3,c4,c5, ncol=2)

dev.off()

######## QUALITY OF PAPERS

db$Species3 = db$Species2
levels(db$Species3) = c("Others" , "L. tredecimguttatus", "L. rufescens")    

BAR4= data.frame(table(db$Species3,db$Error_figures))

e2 = ggplot(BAR4, aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity",color="black",position=position_dodge(), alpha=1)+
  geom_text(aes(label=Freq), vjust=-1, color="black",
            position = position_dodge(0.9), size=3.5)+
  
  scale_fill_manual("",labels=c("No", "Yes"),values=c("grey20", "deepskyblue3"))+
  labs(title="Error in photographs", subtitle =  "[i.e., wrong species depicted]",x=NULL, y = "Frequency")+
  theme_classic()+
  theme(legend.position = c(0.15, 0.8),
        legend.text = element_text(size=10,face = "italic"),
        legend.title = element_text(size=12,face="bold"), 
        axis.text.x = element_text(size=12,face = c("plain","italic","italic")))

e2

BAR5= data.frame(table(db$Species3,db$Taxonomic_error))

e3 = ggplot(BAR5, aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity",color="black",position=position_dodge(), alpha=1)+
  geom_text(aes(label=Freq), vjust=-1, color="black",
            position = position_dodge(0.9), size=3.5)+
  scale_fill_manual(labels=c("No", "Yes"),values=c("grey20", "deepskyblue3"))+
  labs(title="Error related to taxonomy", subtitle =  "[e.g., '... a spider is an insect...' ]",x=NULL, y = "Frequency")+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12,face = c("plain","italic","italic")))

e3

BAR6= data.frame(table(db$Species3,db$Venom_error))

e4 = ggplot(BAR6, aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity",color="black",position=position_dodge(), alpha=1)+
  geom_text(aes(label=Freq), vjust=-1, color="black",
            position = position_dodge(0.9), size=3.5)+
  
  scale_fill_manual(labels=c("No", "Yes"),values=c("grey20", "deepskyblue3"))+
  labs(title="Error related to venom", subtitle =  "[e.g., '... the spiders injected a venom sac...' ]",x=NULL, y = "Frequency")+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12,face = c("plain","italic","italic")))

e4


BAR7= data.frame(table(db$Species3,db$Sting_error))

e5 = ggplot(BAR7, aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity",color="black",position=position_dodge(), alpha=1)+
  geom_text(aes(label=Freq), vjust=-1, color="black",
            position = position_dodge(0.9), size=3.5)+
  ylim(0,140)+
  
  scale_fill_manual(labels=c("No", "Yes"),values=c("grey20", "deepskyblue3"))+
  labs(title="Error related to spider anatomy", subtitle = "[e.g., '... spiders sting ...' ]",x=NULL, y = "Frequency")+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12,face = c("plain","italic","italic")))

e5

pdf("/Users/stefanomammola/Desktop/Fiugre_errors.pdf", width = 19, height = 10)

lay <- rbind(c(1,1,2,3),
             c(1,1,4,5))
grid.arrange(e1,e2,e3,e4,e5, layout_matrix = lay)

dev.off()

############## SHARE ON SOCIAL MEDIA

dotchart(db$Share_FB)

#Zero inflatio
sum(db_ll$Share_FB==0)/nrow(db)*100

db_glm=db_ll[db_ll$Share_FB<15000,]
db_glm=db_glm[db_glm$yr<2020,]

db_glm$yr = as.numeric(db_glm$yr)
db_glm$ID_event = as.factor(db_glm$ID_event)

ggplot(data=db_glm, aes(y=Share_FB, x=Sensationalism_consensus))+
  labs(x="Sensationalistic new", y = "Shares on Facebook")+
  geom_boxplot(fill=c("grey20", "deepskyblue3"))+theme_classic()

ggplot(data=db_glm, aes(y=Share_FB, x=Species2))+
  labs(x="Species considered in the new", y = "Shares on Facebook")+
  geom_boxplot(fill=c("grey50", "#E69F00"))+theme_classic()

ggplot(data=db_glm, aes(y=Share_FB, x=Figure_species))+
  labs(x="News with a figure of the species", y = "Shares on Facebook")+
  geom_boxplot(fill=c("grey20", "deepskyblue3"))+theme_classic()

ggplot(data=db_glm, aes(y=Share_FB, x=Bite))+
  labs(x="News with a figure of the bite", y = "Shares on Facebook")+
  geom_boxplot(fill=c("grey20", "deepskyblue3"))+theme_classic()

ggplot(data=db_glm, aes(y=Share_FB, x=Expert))+
  labs(x="News when an expert was contacted", y = "Shares on Facebook")+
  geom_boxplot(fill=c("grey20", "deepskyblue3"))+theme_classic()

ggplot(data=db_glm, aes(y=Share_FB, x=Circulation))+
  labs(x="Journal' circulation", y = "Shares on Facebook")+
  geom_boxplot()+theme_classic()

ggplot(data=db_glm[db_glm$yr<2020,], aes(y=Share_FB, x=as.factor(as.character(yr))))+
  labs(x=NULL, y = "Shares on Facebook")+
  geom_point()+theme_classic()


## Regression model Poisson

M0 <- glmmadmb(Share_FB ~ Circulation + Bite + Sensationalism_consensus + Species2 + Figure_species + Figure_bite + Expert + (1|Newspaper), data=db_glm,family="poisson")

summary(M0)

overdisp_fun(M0)

## Regression model Negative binomial

M0 <- glmmadmb(Share_FB ~ Bite + Circulation + yr + m + I(m^2) + Sensationalism_consensus + Species2 + Figure_species + Figure_bite + Expert + (1|Newspaper) + (1|ID_event), data=db_glm, family="nbinom")

summary(M0)

M1 <- glmmadmb(Share_FB ~  Bite +  Circulation + yr +Sensationalism_consensus + Species2 + Figure_species + Figure_bite + Expert + (1|Newspaper)+ (1|ID_event), data=db_glm, family="nbinom")

summary(M1)

M2 <- glmmadmb(Share_FB ~  Circulation + yr + Sensationalism_consensus + Species2 + Figure_species + Figure_bite + Expert + (1|Newspaper)+ (1|ID_event), data=db_glm[db_glm$yr<2020,], family="nbinom")

summary(M2)

M3 <- glmmadmb(Share_FB ~ Circulation + yr + Sensationalism_consensus + Figure_species + Figure_bite + Expert + (1|Newspaper)+ (1|ID_event), data=db_glm[db_glm$yr<2020,], family="nbinom")

summary(M3)

M4 <- glmmadmb(Share_FB ~ Circulation + yr + Sensationalism_consensus + Figure_bite + Expert + (1|Newspaper) + (1|ID_event), data=db_glm[db_glm$yr<2020,], family="nbinom")

summary(M4)

M5 <- glmmadmb(Share_FB ~ Circulation + yr + Sensationalism_consensus + Expert + (1|Newspaper)+ (1|ID_event), data=db_glm, family="nbinom")

summary(M5)

M6 <- glmmadmb(Share_FB ~ Circulation + yr + Sensationalism_consensus + (1|Newspaper)+ (1|ID_event), data=db_glm[db_glm$yr<2020,], family="nbinom")

summary(M6)

M7 <- glmmadmb(Share_FB ~ yr + Sensationalism_consensus + (1|Newspaper)+ (1|ID_event), data=db_glm[db_glm$yr<2020,], family="nbinom")

summary(M7)


plot(M7)

E7 = resid(M7, type="response")
boxplot(E7 ~ db_glm[db_glm$yr<2020,]$ID_event) 
boxplot(E7 ~ db_glm[db_glm$yr<2020,]$Newspaper) 
plot(E7 ~ db_glm[db_glm$yr<2020,]$yr) 
boxplot(E7 ~ db_glm[db_glm$yr<2020,]$Sensationalism_consensus) 


VISREG=visreg(M7,"Sensationalistic", gg=TRUE,
              line=list(col="deepskyblue3",size=1.5),
              fill=list(fill="grey70",alpha=0.6),
              points=list(size=2,col="grey25", pch=16,alpha=0.6)) +  
  labs(x="Sensationalistic new", y = "Residuals (Shares on Facebook)")+
  theme_bw()


### model selection

AICs <- AIC(M0,M1,M2,M3,M4,M5,M6,M7)
MyDf <- AICs[,1]
AICsNum <- AICs[,2]
minAW <- min(AICsNum)
Delta <- AICsNum-minAW
RL <- exp(-0.5 * Delta)
wi <- round(RL / sum(RL),3)
Z <- data.frame(Model=rownames(AICs),AICs, AICsNum, Delta, wi)
Z <- Z[order(Z$AIC),]
Z


# Model df     AIC AICsNum Delta    wi
# M7    M7  6 3057.76 3057.76  0.00 0.414
# M5    M5  8 3059.16 3059.16  1.40 0.206
# M6    M6  7 3059.72 3059.72  1.96 0.156
# M4    M4  9 3060.40 3060.40  2.64 0.111
# M3    M3 10 3061.48 3061.48  3.72 0.065
# M2    M2 11 3063.16 3063.16  5.40 0.028
# M1    M1 12 3064.68 3064.68  6.92 0.013
# M0    M0 14 3065.66 3065.66  7.90 0.008


# Call:
#   glmmadmb(formula = Share_FB ~ yr + Sensationalism_consensus + 
#              (1 | Newspaper) + (1 | ID_event), data = db_glm[db_glm$yr < 
#                                                                2020, ], family = "nbinom")
# 
# AIC: 3057.8 
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                 -2688.362    434.440   -6.19  6.1e-10 ***
#   yr                              1.334      0.215    6.20  5.6e-10 ***
#   Sensationalism_consensusYes     1.150      0.503    2.29    0.022 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Number of observations: total=312, Newspaper=166, ID_event=92 
# Random effect variance(s):
#   Group=Newspaper
# Variance StdDev
# (Intercept)    6.554   2.56
# Group=ID_event
# Variance   StdDev
# (Intercept) 3.247e-05 0.005698



M7 <- glm.nb(Share_FB ~ yr + Sensationalism_consensus, data=db_glm)

MyData1 <- data.frame(yr=seq(from =range(db_glm$yr)[1],
                             to = range(db_glm$yr)[2],
                             length = 100),
                      Sensationalism_consensus="No")

MyData2 <- data.frame(yr=seq(from =2010,
                             to = 2019,
                             length = 100),
                      Sensationalism_consensus="Yes")

P1 <- predict(M7, newdata = MyData1,se=TRUE,type="response") 
P2 <- predict(M7, newdata = MyData1,se=TRUE,type="response",interval = "confidence")
MyData1$mu    <- P1$fit  #Fitted values
MyData1$selow <- P1$fit - 2 * P1$se.fit  #lower bound
MyData1$seup  <- P1$fit + 2 * P1$se.fit   #lower bound

MyData2$mu    <- P1$fit  #Fitted values
MyData2$selow <- P1$fit - 2 * P1$se.fit  #lower bound
MyData2$seup  <- P1$fit + 2 * P1$se.fit   #lower bound

MyData =rbind(MyData1,MyData2)

ggplot() + 
  xlab("Year") +
  ylab("Shares on Facebook")+
  #facet_grid(~ Sensationalistic)+
  
  scale_x_continuous(breaks = c(2010:2019), labels = as.character(c(2010:2019)))+ 
  geom_point(data = db_glm, aes(y = Share_FB, x = yr),shape = 16, size = 3, col="black",alpha=0.6) + 
  #geom_point(data = db_glm, aes(y = Share_FB, x = yr, col=Species2),shape = 16, size = 2) + 
  scale_color_manual(values=c("black", "#E69F00"))+
  geom_ribbon(data = MyData2, aes(x = yr, ymax = seup,ymin = selow ),fill="deepskyblue3",alpha = 0.6)+
  geom_line(data = MyData2,aes(x = yr, y = mu)) +
  
  theme_classic()+theme(legend.position = c(0.25, 0.8),
                        legend.text = element_text(size=12,face = "italic"),
                        legend.title = element_text(size=0,face="bold"))


ggplot(data=db_glm, aes(y=Share_FB, x=Sensationalism_consensus))+
  labs(x="Sensationalistic content", y = "Shares on Facebook")+
  geom_boxplot(fill=c("grey40", "deepskyblue3"),alpha=0.8)+ 
  geom_jitter(alpha=0.6,col="orange",shape=16, position=position_jitter(0.4))+
  
  theme_classic()




##############

library(maptools)
library(ggplot2)
library(ggalt)
library(ggthemes)
library(tibble)
library(viridis)

# get italy region map
italy_map <- map_data("italy")

# your data will need to have these region names
print(unique(italy_map$region))

# we'll simulate some data for this
set.seed(1492)
choro_dat <- data_frame(region=unique(italy_map$region),
                        value=sample(100, length(region)))

# we'll use this in a bit
italy_proj <- "+proj=aea +lat_1=38.15040684902542
+lat_2=44.925490198742295 +lon_0=12.7880859375"

gg <- ggplot()

# lay down the base layer
gg <- gg + geom_map(data=italy_map, map=italy_map,
                    aes(long, lat, map_id=region),
                    color="#b2b2b2", size=0.1, fill=NA)



# fill in the regions with the data
gg <- gg + geom_map(data=choro_dat, map=italy_map,
                    #aes(fill=value, map_id=region)
                    aes(map_id=region),
                    fill="grey70",
                    color="#b2b2b2", size=0.1)

# great color palette (use a better legend title)
gg <- gg + scale_fill_viridis(name="Scale title")

# decent map projection for italy choropleth
gg <- gg + coord_proj(italy_proj)

gg <- gg + geom_point(data=db_unique_ID, aes(x=x, y=y),color = "black",  size=2.4)
gg <- gg + geom_point(data=db_unique_ID, aes(x=x, y=y,col=Species), size=2)
gg <- gg +   scale_color_manual(values=c("blue","grey20", "#E69F00"))

# good base theme for most maps
gg <- gg + theme_classic()

# move the legend
gg <- gg + theme(legend.position = c(0.2, 0.2),
                 legend.text = element_text(size=12,face = "italic"),
                 legend.title = element_text(size=0,face="bold"))


gg



##############


library(tidyverse)
library(viridis)
library(hrbrthemes)
library(mapdata)
library(dplyr)
library(ggplot2)


library(rgdal)
library(raster)
boundaries<-readOGR("/Users/stefanomammola/Desktop/BIOCLIM/SHAPE","GSHHS_i_L1") 
boundaries=crop(boundaries,extent(0,40,34,50))


library(rworldmap)
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

Italy <- world %>% filter(sovereignt == "Italy")


# plot
db %>%
  #filter(homecontinent=='Europe') %>%
  ggplot( aes(x=x, y=y)) + 
  #geom_sf(data = Italy)+
  geom_hex(bins=39,alpha=1) +
  #geom_polygon(data=boundaries,aes(x=long, y=lat, group=group),fill=NA, size=.4,color="grey20",lty=1, alpha = 1)+
  ggplot2::annotate("text", x = -27, y = 72, label="Where spider-related events took place", colour = "black", size=5, alpha=1, hjust=0) +
  ggplot2::annotate("segment", x = -27, xend = 10, y = 70, yend = 70, colour = "black", size=0.2, alpha=1) +
  theme_void() +
  xlim(0, 40) +
  ylim(34, 50) +
  scale_fill_viridis(
    option="B",
    breaks=c(10,20,30,40),
    name="N° of events", 
    guide = guide_legend( keyheight = unit(2.5, units = "mm"), keywidth=unit(10, units = "mm"), label.position = "bottom", title.position = 'top', nrow=1) 
  )  +
  ggtitle( "" ) +
  theme(
    legend.position = c(0.8, 0.09),
    legend.title=element_text(color="black", size=8),
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(size= 13, hjust=0.1, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
  ) 




library(ggplot2)
library(maptools)

# Load the dataset of points



italy_map <- map_data("italy")
colnames(italy_map)[1] <- "x"
colnames(italy_map)[2] <- "y"

p  <- ggplot(data=db, aes(x=x,y=y)) + stat_binhex(bins=80, binwidth=c(1,1)) +
  geom_map(data=italy_map, map=italy_map,
           aes(x, y, map_id=region),
           color="#b2b2b2", size=1, fill=NA)+
  scale_fill_gradientn(colours=c('light gray','blue'),name='Frequency',na.value=NA) +
  coord_equal() +
  labs(x=NULL, y=NULL) +
  theme_bw() +
  theme(legend.position=c(0.075,0.28))



hist(rnorm(10))

