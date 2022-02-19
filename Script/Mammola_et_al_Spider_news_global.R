###############################################################

## Mammola, S. et al. (2022) The global spread of misinformation on spiders

###############################################################

## ------------------------------------------------------------------------
# 'R script to reproduce the analyses'
## ------------------------------------------------------------------------

## Authors: Stefano Mammola
## Software: R (v. R 4.1.0) and R studio (v. 1.4.1103)

###############################################################

# clean the workspace -----------------------------------------------------

rm(list=ls())

# Loading R package -------------------------------------------------------
library("Amelia")
library("BAT")
library("bipartite")     
library("dplyr")  
library("ergm")
library("flextable")
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
library("magrittr")
library("MuMIn")
library("network")       
library("performance")   
library("PupillometryR") 
library("rworldmap")
library("scatterpie")
library("sna")           
library("tidygraph")
library("tidyverse")   
library("tidylog")

# Loading useful functions ------------------------------------------------

source("Functions/Functions_spider_news.R")

# Loading plot parameters -------------------------------------------------

source("Functions/Plot_parameters.R")

# Set seed ----------------------------------------------------------------

set.seed(123)

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

#Number of Languages
nlevels(droplevels(db$Language)) 

#Distinct event
table(db_unique_news$TypeEvent)

# Create attribute table for Network analysis -----------------------------

#############################
# At the Country level ######
#############################

# Loading country-level attributes
CountryAttributes <- read.csv(file = "Data/CountryAttributes_Demography.csv", sep = '\t', dec = ',', header = TRUE, as.is = FALSE)

#Proportion of missing data?
for (i in 2:(ncol(CountryAttributes)-1))
  message(paste("-- Column ", 
                colnames(CountryAttributes)[i], 
                " contains ", 
                round((sum(is.na(CountryAttributes[, i]) == TRUE)/nrow(CountryAttributes)) * 100, 2), 
                "% of missing data.", sep = ""))

#Collinearity?
psych::pairs.panels(CountryAttributes[,3:13])

# Summarizing country attributes from the Spider news database
country_attr <- db %>% dplyr::select(Country_search,
                                     Sensationalism,
                                     TotalError,
                                     Language,
                                     lon3,
                                     lat3) %>% group_by(Country_search) %>% 
  summarise(N = length(Country_search),
            Sensationalism = sum(Sensationalism, na.rm = TRUE)/length(Sensationalism), #proportion of sensationalistic news/country
            TotalError = sum(TotalError, na.rm = TRUE)/length(TotalError), #proportion of news with error/country
            Language = names(sort(table(Language), decreasing = TRUE))[1], #main country Language
            lon = mean(lon3), #Longitude
            lat = mean(lat3)) #latitude

# Add Countries with no News
NoNews <- data.frame(Country_search = c("Botswana", "Iceland"),
                        N = c(0,0),
                        Sensationalism = c(0,0),
                        TotalError = c(0,0),
                        Language = c("English","Icelandic"),
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
                  Language,
                  lon, 
                  lat, 
                  ISA, #number of arachnologist in the country.
                  HDI = Human_Development_Index_.HDI._2019,
                  Education_index  = Education_Index_2019,  
                  Internet_users   = Internet_users_2018,
                  N_newspapers     = Number_newspapers_Mean1996.2005,
                  Press_Freedom    = Press_Freedom_Index,
                  N_Spiders        = N_Spiders,
                  N_Deadly_Spiders = N_Deadly)

rm(CountryAttributes) #clean

# All language with less than 10 country become "Others":
country_attr$Language <- ifelse(country_attr$Language %in% names(which(table(country_attr$Language)>9)), country_attr$Language , "Others")

# Predict missing data

country_attr[,9:13] <- BAT::fill(data.frame(country_attr[,9:13]))

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

# net_map <- getMap(resolution="high") %>% 
#         ggplot() +
#         geom_polygon(aes(long,lat, group=group), colour="grey5", fill="grey5", size = 0.25)+
#   ylim(-56.8,90)+ xlim(-180,195)+
#   geom_curve(aes(x = jitter(lon3,0.0001), 
#                  y = jitter(lat3,0.0001), 
#                  xend = jitter(lon, 0.0001), 
#                  yend = jitter(lat, 0.0001)),
#              data = db, curvature = 0.24,lwd=0.4,
#              alpha = 0.05,  color = "aquamarine")+#"#FFEA46FF")+
#  
#   geom_point(data = db_unique_event, 
#                aes(x = lon, y = lat), size = 0.3,
#                alpha = 0.9, color="aquamarine",
#                stroke = 0.8) +
#   ggthemes::theme_map() + theme_map_custom
# 
# figure_1 <- net_map + scatterpie::geom_scatterpie(data = pie, aes(x=lon, y= lat, group = Country, r = radius),
#                                             cols = c("Non_Sensationalist","Sensationalist"), color = "white", lwd=0.05, alpha=.8) +
#     theme(legend.position = c(0.03, 0.2),
#           legend.text = element_text(size = 10),
#           legend.background = element_rect(fill = "gray40", colour="transparent"))+
#     geom_scatterpie_legend(pie$radius,
#                            x= -165,
#                            y= -45, n = 2,
#                            labeller = function (x) x=c(min(pie$n), max(pie$n)))+
#     scale_fill_manual("",labels = c("Not Sensationalist","Sensationalist"), values = c("blue","darkmagenta"))
# 
# ggplot2::ggsave("Figures/Figure_1.pdf", 
#                 figure_1, 
#                 device = cairo_pdf,
#                 units = "cm",
#                 width = 20,
#                 height = 9)
# 
# rm(pie, pie_1, n, radius, net_map) #clean

net_map <- getMap(resolution="high") %>% 
  ggplot() +
  geom_polygon(aes(long,lat, group=group), colour="grey15", fill="grey15", size = 0.25)+
  ylim(-56.8,90)+ xlim(-180,195)+
  geom_curve(aes(x = jitter(lon3,0.0001), 
                 y = jitter(lat3,0.0001), 
                 xend = jitter(lon, 0.0001), 
                 yend = jitter(lat, 0.0001)),
             data = db, curvature = 0.24,lwd=0.4,
             alpha = 0.05,  color = "darkcyan")+#"#FFEA46FF")+
  
  geom_point(data = db_unique_event, 
             aes(x = lon, y = lat), size = 0.3,
             alpha = 0.9, color="darkcyan",
             stroke = 0.8) +
  ggthemes::theme_map() + theme_map_custom

figure_1 <- net_map + scatterpie::geom_scatterpie(data = pie, aes(x=lon, y= lat, group = Country, r = radius),
                                                  cols = c("Non_Sensationalist","Sensationalist"), color = "white", lwd=0.05, alpha=.8) +
  theme(legend.position = c(0.03, 0.2),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = "white", colour="transparent"))+
  geom_scatterpie_legend(pie$radius,
                         x= -165,
                         y= -45, n = 2,
                         labeller = function (x) x=c(min(pie$n), max(pie$n)))+
  scale_fill_manual("",labels = c("Not Sensationalist","Sensationalist"), values = c("blue","darkmagenta"))

ggplot2::ggsave("Figures/Figure_1.pdf", 
                figure_1, 
                device = cairo_pdf,
                units = "cm",
                width = 25,
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

db$Language2 <- db$Language ; table(db$Language2)
rename <- names(which(table(db$Language2) < 50))
levels(db$Language2)[levels(db$Language2) %in% rename]  <- "Other"

# Setting baseline for factors
db <- within(db, Family2 <- relevel(Family2, ref = "Other"))
db <- within(db, Country_search2 <- relevel(Country_search2, ref = "Other"))
db <- within(db, Circulation <- relevel(Circulation, ref = "Regional"))
db <- within(db, Type_of_newspaper <- relevel(Type_of_newspaper, ref = "Traditional newspaper"))
db <- within(db, Language2<- relevel(Language2, ref = "Other"))

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
  Language2,
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

formula_m1 <- as.formula("Sensationalism ~ yr + Type_of_newspaper + Circulation + TypeEvent + Figure_species + Figure_bite +
                    TotalError.01 + Expert_doctor + Expert_arachnologist + Expert_others +
                    (1|ID_Event) + (1|Genus) + (1|Country_search2) + (1|Language2)")
# Fit the model
m1 <- lme4::glmer(formula_m1, 
                  data    = db_m, 
                  family  = binomial(link = "logit"),
                  control = glmerControl(optimizer="bobyqa"))

# Check model
performance::check_model(m1, check = c("vif"))
performance::check_model(m1, check = c("reqq"))

R_m1 <- residuals(m1)

par(mfrow = c(2,2),mar = c(2,2,2,2))
plot(db_m$lon2, R_m1, xlab = "Longitude", ylab = "Residuals",  main = "Residuals vs Longitude")
plot(db_m$lat2, R_m1, xlab = "Latitude",  ylab = "Residuals",  main = "Residuals vs Latitude")
boxplot(R_m1 ~ as.factor(db_m$m),  xlab = NULL,  ylab = "Residuals",  main = "Residuals vs month")
boxplot(R_m1 ~ as.factor(db_m$yr), xlab = NULL,  ylab = "Residuals",  main = "Residuals vs year")

# Interpret the model
(fit_m1 <- parameters::model_parameters(m1))
performance::r2(m1)

# Save the table
sjPlot::tab_model(m1, file = "Tables/Table S1.docx")

# Plot
Estimates_m1 <- 
  m1 %>% 
  summary %>% 
  magrittr::extract2("coefficients") %>% # extract estimates
  as.data.frame %>% rownames_to_column("Variable") %>% 
  dplyr::filter(!row_number() %in% 1) %>%  #remove intercept
  dplyr::rename(SE = 3, z = 4, p = 5) #rename

Estimates_m1$Variable <- axis_labels_sjPlot_m1

(plot_model1 <- ggplot2::ggplot(data = Estimates_m1, aes(Variable, Estimate)) +
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    geom_errorbar(aes(ymin = Estimate-SE, ymax = Estimate+SE), width = 0, col = "grey10") +
    
    geom_text(aes(Variable, Estimate), 
              label = round(Estimates_m1$Estimate,2), 
              vjust = -1, size = 3) +
    geom_point(size = 2, pch = 21, col = "grey10", fill = "grey20") +
    labs(title = title_sjPlot_m1,
         y = xlab_sjPlot_m1_m2,
         x = NULL)+
    theme_custom() + coord_flip())

rm(R_m1) #clean

# Fitting the model: Errors ---------------------------------------

formula_m2 <- as.formula("TotalError.01 ~ yr + Type_of_newspaper + Circulation + TypeEvent +
                         Sensationalism + Expert_doctor + Expert_arachnologist + Expert_others +
                         (1|ID_Event) + (1|Genus) + (1|Language2) + (1|Country_search2)")
# Fit the model
m2 <- lme4::glmer(formula_m2, 
                  data = db_m, 
                  family = binomial(link = "logit"),
                  control = glmerControl(optimizer="bobyqa"))

# Check model
performance::check_model(m2, check = c("vif"))
performance::check_model(m2, check = c("reqq"))

R_m2 <- residuals(m2)

par(mfrow = c(2,2),mar = c(2,2,2,2))
plot(db_m$lon2, R_m2, xlab = "Longitude", ylab = "Residuals",  main = "Residuals vs Longitude")
plot(db_m$lat2, R_m2, xlab = "Latitude",  ylab = "Residuals",  main = "Residuals vs Latitude")
boxplot(R_m2 ~ as.factor(db_m$m), xlab = NULL, ylab = "Residuals",  main = "Residuals vs month")
boxplot(R_m2 ~ as.factor(db_m$yr), xlab = NULL,  ylab = "Residuals",  main = "Residuals vs year")

# Interpret the model
(fit_m2 <- parameters::model_parameters(m2))
performance::r2(m2)

#save the table
sjPlot::tab_model(m2, file = "Tables/Table S2.docx")

#Plot
Estimates_m2 <- 
  m2 %>% 
  summary %>% 
  magrittr::extract2("coefficients") %>% # extract estimates
  as.data.frame %>% rownames_to_column("Variable") %>% 
  dplyr::filter(!row_number() %in% 1) %>%  #remove intercept
  dplyr::rename(SE = 3, z = 4, p = 5) #rename

Estimates_m2$Variable <- axis_labels_sjPlot_m2

(plot_model2 <- ggplot2::ggplot(data = Estimates_m2, aes(Variable, Estimate)) +
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    geom_errorbar(aes(ymin = Estimate-SE, ymax = Estimate+SE), width = 0, col = "grey10") +
    
    geom_text(aes(Variable, Estimate), 
              label = round(Estimates_m2$Estimate,2), 
              vjust = -1, size = 3) +
    geom_point(size = 2, pch = 21, col = "grey10", fill = "grey20") +
    labs(title = title_sjPlot_m2,
         y = xlab_sjPlot_m1_m2,
         x = NULL)+
    theme_custom() + coord_flip())

rm(R_m2) #clean

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
round(NetworkTraitGet(Graph_unipartite),2)

# Calculating node level traits 
node_trait <- Graph_unipartite %>% 
                   NodeTraitGet(mode = "in", dir = FALSE) %>% bind_cols()

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

#layout <- create_layout(Graph_unipartite, layout = 'igraph', algorithm = 'kk')

#SpatialLayout
(plot_network <- Graph_unipartite %>% 
  igraph::simplify(edge.attr.comb = "sum") %>% 
  ggraph::ggraph(SpatialLayout) +
  #geom_edge_density(fill="orange", alpha=1) +
  geom_edge_fan(aes(width=weight),color="darkcyan", alpha=0.05) +
  geom_node_point(col="grey30", alpha = .8, 
                  aes(size=N,fill=Language), shape = 21) + 
  geom_node_text(aes(label = name), size=2, color="gray10", repel=TRUE) +
  scale_fill_manual(values = c("blue", "orange", "pink","purple", "grey15")) +
  theme_void() + theme(legend.position = "bottom",legend.direction = "vertical") + coord_fixed())

# Figure S1 ---------------------------------------------------------------

pdf(file = "Figures/Figure_S1.pdf", width = 10, height = 7)
plot_network
dev.off()

###################################
# Regression Models #
###################################

db_m3_cat <- country_attr %>% dplyr::select(Degree, Country_search, Language) #Dependent variable + categorical
db_m3_con <- country_attr %>% dplyr::select(N, ISA, Sensationalism, TotalError, Education_index, HDI,
                                            Internet_users, N_newspapers, Press_Freedom, N_Spiders, N_Deadly_Spiders) #Continuous variables

db_m3_con <- data.frame(db_m3_con)

# Check distributon
par(mfrow = c(3,3),mar = c(2,2,2,2))
for(i in 1:ncol(db_m3_con))
  dotchart(as.numeric(db_m3_con[,i]), main = colnames(db_m3_con)[i])

#log transform
db_m3_con$N            <- log(db_m3_con$N+1)
db_m3_con$ISA          <- log(db_m3_con$ISA+1)
db_m3_con$N_newspapers <- log(db_m3_con$N_newspapers+1)

# Check collinearity
psych::pairs.panels(db_m3_con) 

#Education index, with HDI and internet users  
#Number of newspaper with N 
#ISA and N
#Dealy spider and N spiders


# Check collinearity vs categorical
(collinearity <- country_attr %>% dplyr::select(Language, N, ISA, Sensationalism, TotalError, Education_index, HDI,
                               Internet_users, N_newspapers, Press_Freedom, N_Spiders, 
                               N_Deadly_Spiders) %>% 
       GGally::ggpairs(aes(colour = Language, alpha = 0.4)))


# Figure S2 ---------------------------------------------------------------

pdf(file = "Figures/Figure_S2.pdf", width = 18, height =14)
psych::pairs.panels(db_m3_con) 
dev.off()

# #

#Check balance of levels
db_m3_cat$Language %>% table

# Scale all var
db_m3_con <- db_m3_con %>% apply(MARGIN = 2, scale) #scale cont variable

# Merge
db_m3 <- data.frame(db_m3_cat, db_m3_con) ; rm(db_m3_cat, db_m3_con)

# Remove missing values (if any)
db_m3_noNA <- na.omit(db_m3)

# Set reference level
db_m3_noNA$Language <- as.factor(db_m3_noNA$Language)

db_m3_noNA$Language <- factor(db_m3_noNA$Language, levels = c("Others", "English", "Arabic", "Russian", "Spanish"))

# Fit the model -----------------------------------------------------------

# Set formula
formula_m3 <- as.formula("Degree ~ Sensationalism + TotalError + Internet_users +
                         Press_Freedom + N_Spiders + (1|Language)")

# Fit the model
m3 <- lme4::glmer(formula_m3, 
                  data    = db_m3_noNA, 
                  family  = poisson,
                  control = glmerControl(optimizer="bobyqa"))

# Check model
performance::check_overdispersion(m3)

# Fit the model (nd)
m3 <- lme4::glmer.nb(formula_m3, 
                  data    = db_m3_noNA, 
                  control = glmerControl(optimizer="bobyqa"))

# Check model
performance::check_model(m3, check = c("vif"))
performance::check_model(m3, check = c("reqq"))

# Interpret the model
(fit_me <- parameters::model_parameters(m3))
performance::r2(m3)

# Save the table
sjPlot::tab_model(m3, file = "Tables/Table S3.docx")

# Plot
Estimates_m3 <- 
  m3 %>% 
  summary %>% 
  magrittr::extract2("coefficients") %>% # extract estimates
  as.data.frame %>% rownames_to_column("Variable") %>% 
  dplyr::filter(!row_number() %in% 1) %>%  #remove intercept
  dplyr::rename(SE = 3, z = 4, p = 5) #rename

Estimates_m3$Variable <- axis_labels_sjPlot_m3

(plot_model3 <- ggplot2::ggplot(data = Estimates_m3, aes(Variable, Estimate)) +
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    geom_errorbar(aes(ymin = Estimate-SE, ymax = Estimate+SE), width = 0, col = "grey10") +
    geom_text(aes(Variable, Estimate), 
              label = round(Estimates_m3$Estimate,2), 
              vjust = -1, size = 3) +
    geom_point(size = 2, pch = 21, col = "grey10", fill = "grey20") +
    labs(title = title_sjPlot_m3,
         y = xlab_sjPlot_m3,
         x = NULL)+
    theme_custom() + coord_flip())

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

levels(db_m3_noNA$Language)[2] <- "AAA_English"

#Adding node-level attributes
ResponseNetwork %v% "Language"         <- as.character(db_m3_noNA[1:79,]$Language)
ResponseNetwork %v% "Sensationalism"   <- db_m3_noNA[1:79,]$Sensationalism
ResponseNetwork %v% "TotalError"       <- db_m3_noNA[1:79,]$TotalError
ResponseNetwork %v% "Internet"         <- db_m3_noNA[1:79,]$Internet_users
ResponseNetwork %v% "Freedom"          <- db_m3_noNA[1:79,]$TotalError
ResponseNetwork %v% "N_Spiders"        <- db_m3_noNA[1:79,]$N_Spiders

# Model fit
# ergm0 <- ergm::ergm(ResponseNetwork ~ edges, estimate = "MLE")
# ergm1 <- ergm::ergm(ResponseNetwork ~ edges + edgecov(SpatialLayout), estimate = "MLE")
# ergm2 <- ergm::ergm(ResponseNetwork ~ edges + edgecov(SpatialLayout) + nodecov("Sensationalism"), estimate = "MLE")
# ergm3 <- ergm::ergm(ResponseNetwork ~ edges + edgecov(SpatialLayout) + nodecov("Sensationalism") + nodecov("TotalError"), estimate = "MLE")
# ergm4 <- ergm::ergm(ResponseNetwork ~ edges + edgecov(SpatialLayout) + nodecov("Sensationalism") + nodecov("TotalError") + nodecov("N_Spiders"), estimate = "MLE")
# ergm5 <- ergm::ergm(ResponseNetwork ~ edges + edgecov(SpatialLayout) + nodecov("Sensationalism") + nodecov("TotalError") + nodecov("N_Spiders") + nodefactor("Language"), estimate = "MLE")
# ergm6 <- ergm::ergm(ResponseNetwork ~ edges + edgecov(SpatialLayout) + nodecov("Sensationalism") + nodecov("TotalError") + nodecov("N_Spiders") + nodefactor("Language") + nodecov("Internet"), estimate = "MLE")
# ergm7 <- ergm::ergm(ResponseNetwork ~ edges + edgecov(SpatialLayout) + nodecov("Sensationalism") + nodecov("TotalError") + nodecov("N_Spiders") + nodefactor("Language") + nodecov("Internet") + nodecov("Freedom"), estimate = "MLE")
# ergm8 <- ergm::ergm(ResponseNetwork ~ edges + edgecov(SpatialLayout) + nodecov("Sensationalism") + nodecov("TotalError") + nodecov("N_Spiders") + nodefactor("Language") + nodematch("Language") + nodecov("Internet"), estimate = "MLE")
ergm <- ergm::ergm(ResponseNetwork ~ edges + nodecov("Sensationalism") + nodecov("TotalError") + nodecov("N_Spiders") + nodecov("Internet") + nodefactor("Language") + nodematch("Language"), estimate = "MLE")

# list(ergm0, ergm1, ergm2, ergm3, ergm4, ergm5, ergm6, ergm7, ergm8, ergm9) %>% map_dbl(BIC) %>% plot
# ergm_BIC <- ergm

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

summary(ergm)

# Interpret the model
EstimateDF <- 
  ergm %>% summary %>% extract2("coefficients") %>% 
  as.data.frame %>% 
  rownames_to_column("Variable")

ergm %>% 
  confint %>% # Takes the 95% confidence intervals
  as.data.frame %>% 
  rename(Lower = 1, Upper = 2)

EstimateDF %<>% # Bind them together
  bind_cols(ergm %>% confint %>% as.data.frame %>% 
              rename(Lower = 1, Upper = 2))

# Clean the table
EstimateDF_tab <- EstimateDF

CI <- paste0(round(EstimateDF_tab$Lower,2), 
             rep(" — ", nrow(EstimateDF_tab)), 
             round(EstimateDF_tab$Upper,2))

EstimateDF_tab <- EstimateDF_tab[ , -c(4,7,8)]
EstimateDF_tab <- data.frame(EstimateDF_tab[,1:3], CI, EstimateDF_tab[,4:5])
rownames(EstimateDF_tab) <- NULL 

for(s in c(2,3,5))
  EstimateDF_tab[,s] <- round(EstimateDF_tab[,s],2)

colnames(EstimateDF_tab)[c(3,5,6)] <- c("SE", "z", "p")

EstimateDF_tab$p <- ifelse(EstimateDF_tab$p < 0.001, "<0.001", as.character(round(EstimateDF_tab$p,3)))  
  
#Save the table
write.table(EstimateDF_tab,"Tables/Table S4.csv") ; rm(EstimateDF_tab, CI)

#Remove edges for plot
Estimate_m4 <- EstimateDF[2:nrow(EstimateDF),]

Estimate_m4$Variable <- c("Prop. of sensationalistic news", 
                         "Prop. of news with errors",
                         "N° of spiders",
                         "Internet users",
                         "Language [Arabic]",
                         "Language [Others]",
                         "Language [Russian]",
                         "Language [Spanish]",
                         "Node match: Language")

# Plot
(plot_ergm1 <- Estimate_m4 %>% ggplot2::ggplot(aes(Variable, Estimate)) +
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, col = "grey10") +
    geom_text(aes(Variable, Estimate), 
              label = round(Estimate_m4$Estimate,2), 
              vjust = -1, size = 3) +
    geom_point(size = 2, pch = 21, col = "grey10", fill = "grey20") +  
    labs(y = xlab_ergm1, 
         title = title_ergm1) +
    geom_point() + theme_custom() + theme(axis.title.y=element_blank()) + coord_flip())

# Figure 2 ----------------------------------------------------------------

pdf(file = "Figures/Figure_2.pdf", width = 16, height =12)
ggpubr::ggarrange(plot_model1,plot_model2,plot_model3,plot_ergm1,
                  common.legend = FALSE,
                  hjust = -5,
                  align = "hv",
                  labels = c("A", "B", "C", "D"),
                  ncol=2, nrow=2)
dev.off()

# With count data ---------------------------------------------------------

#Generalized exponential random graph model
library("ergm.count")

AdjMatrix2 <- Graph_unipartite %>% get.adjacency(attr = "weight", sparse = FALSE) 
AdjMatrix2 %>% as.matrix %>% dim #square

ResponseNetwork2 <- as.network(AdjMatrix2 %>% as.matrix, 
                              directed = TRUE, 
                              matrix.type = "a", 
                              ignore.eval = FALSE, 
                              names.eval = "weight")  # Important! )

levels(db_m3_noNA$Language)[2] <- "AAA_English"

#Adding node-level attributes
ResponseNetwork2 %v% "Language"         <- as.character(db_m3_noNA[1:79,]$Language)
ResponseNetwork2 %v% "Sensationalism"   <- db_m3_noNA[1:79,]$Sensationalism
ResponseNetwork2 %v% "TotalError"       <- db_m3_noNA[1:79,]$TotalError
ResponseNetwork2 %v% "Internet"         <- db_m3_noNA[1:79,]$Internet_users
ResponseNetwork2 %v% "Freedom"          <- db_m3_noNA[1:79,]$TotalError
ResponseNetwork2 %v% "N_Spiders"        <- db_m3_noNA[1:79,]$N_Spiders

set.edge.value(ResponseNetwork,
               "weight",
               c(AdjMatrix))

ergm_count  <- ergm(ResponseNetwork2 ~ 
                       sum + 
                       nodecov("Sensationalism", form = "sum") + 
                       nodecov("TotalError", form = "sum") +
                       nodecov("Internet", form = "sum") +
                       nodecov("N_Spiders", form = "sum") +
                       nodematch("Language", diff = FALSE, form = "sum"),
                     control = control.ergm(
                       parallel = 10, #change with your parallel
                       MCMC.samplesize = 10000,
                       MCMLE.maxit = 100),
                     response = "weight", 
                     reference = ~ Poisson)

# Check model
mcmc.diagnostics(ergm_count) #good mixing of chains

# Summary stats
summary(ergm_count)

# Plotting
EstimateDF2 <- 
  ergm_count %>% summary %>% extract2("coefficients") %>% 
  as.data.frame %>% 
  rownames_to_column("Variable")

ergm_count %>% 
  confint %>% # Grabbing the 95% confidence intervals
  as.data.frame %>% 
  rename(Lower = 1, Upper = 2)

EstimateDF2 %<>% # Bind them together
  bind_cols(ergm_count %>% confint %>% as.data.frame %>% 
              rename(Lower = 1, Upper = 2))

# Clean the table
EstimateDF_tab2 <- EstimateDF2

CI2 <- paste0(round(EstimateDF_tab2$Lower,2), 
             rep(" — ", nrow(EstimateDF_tab2)), 
             round(EstimateDF_tab2$Upper,2))

EstimateDF_tab2 <- EstimateDF_tab2[ , -c(4,7,8)]
EstimateDF_tab2 <- data.frame(EstimateDF_tab2[,1:3], CI2, EstimateDF_tab2[,4:5])
rownames(EstimateDF_tab2) <- NULL 

for(s in c(2,3,5))
  EstimateDF_tab2[,s] <- round(EstimateDF_tab2[,s],2)

colnames(EstimateDF_tab2)[c(3:5,6)] <- c("SE", "CI", "z", "p")

EstimateDF_tab2$p <- ifelse(EstimateDF_tab2$p < 0.001, "<0.001", as.character(round(EstimateDF_tab2$p,3))) 

#Save the table
write.table(EstimateDF_tab2,"Tables/Table S5.csv") ; rm(EstimateDF_tab2, CI2)

#Remove edges for plot
Estimate_m5 <- EstimateDF2[2:nrow(EstimateDF2),]

Estimate_m5$Variable <- c("Prop. of sensationalistic news",
                         "Prop. of news with errors",
                         "Internet users",
                         "N° of spiders",
                         "Node match: Language")

# Plot
(plot_ergm2 <- Estimate_m5 %>% ggplot2::ggplot(aes(Variable, Estimate)) +
    geom_hline(lty = 3, size = 0.7, col = "grey50", yintercept = 0) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0, col = "grey10") +
    geom_text(aes(Variable, Estimate), 
              label = round(Estimate_m5$Estimate,2), 
              vjust = -1, size = 3) +
    geom_point(size = 2, pch = 21, col = "grey10", fill = "grey20") +
    labs(y     = xlab_ergm1, 
         title = title_ergm2) +
    geom_point() + theme_custom() + theme(axis.title.y=element_blank()) + coord_flip())

  
# Figure S3 ---------------------------------------------------------------

pdf(file = "Figures/Figure_S3.pdf", width = 7, height = 5)
plot_ergm2
dev.off()

# end of the analysis

###############################################################

### Supplementary analysis on network of news ###

###############################################################

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


# Analysing ID_Event-level network properties -----------------------------

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

