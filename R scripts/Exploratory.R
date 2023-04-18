#load packages
library(sf)
library(rgdal)
library(ggplot2)
library(here)
library(readxl)
library(dplyr)
library(lubridate)
library(tidyr)
library(CARBayes)
library(CARBayesST)
library(GGally)
library(spdep)
library (readr)
library(coda)

##################################
### READ DATA AND MANIPULATION ###
##################################

#shapefile to get the map of municipalities of Cuba
cuba_shape <- here("C:/Users/u0149006/OneDrive - KU Leuven/KU Leuven/ASSISTENTSCHAP/Doctoraat/WP CUBA/Data/Map Cuba/Map Cuba/Munic.shp") %>%
  st_read()
cuba_shape

#plot municipalities
ggplot(data = cuba_shape) +
  geom_sf()

#population by municipality
Population_by_municipality <- read_excel("C:/Users/u0149006/OneDrive - KU Leuven/KU Leuven/ASSISTENTSCHAP/Doctoraat/WP CUBA/Data/Population by municipality.xlsx")

#data 
data_Cuba_COVID <- read_delim("C:/Users/u0149006/OneDrive - KU Leuven/KU Leuven/ASSISTENTSCHAP/Doctoraat/WP CUBA/Data/data Cuba COVID until 311021 by municipality ENG.txt", 
                              delim = ";", escape_double = FALSE, trim_ws = TRUE)

#determining spatial neighborhood matrix
centroids <- st_centroid(st_geometry(cuba_shape), of_largest_polygon = T)

nb_q <- poly2nb(st_geometry(cuba_shape), queen = TRUE)

plot(st_geometry(cuba_shape), border='grey', main='Queen contiguity')
plot(nb_q, centroids, col='red', add=TRUE)

#Isla de la Juventud ID=168 -> problem with continuity
no_neighs <- which(card(nb_q)==0)
#  Artificialy create a continuity 
k1nb <- knn2nb(knearneigh(centroids, k=1))
k1nb[no_neighs] <- as.integer(47)
k1nb[no_neighs]
nb_q[no_neighs] <-  k1nb[no_neighs]
nb_q[no_neighs]
attr(nb_q, "sym") <- is.symmetric.nb(nb_q, force=TRUE)
nb_q <- make.sym.nb(nb_q)
#Build new neighborhood matrix
W <- nb2mat(nb_q, style = "B")
W_list <- nb2listw(nb_q, style = "B")

plot(st_geometry(cuba_shape), border='grey', main='Neighborhoud structure')
plot(nb_q, centroids, col='red', add=TRUE)

#convert daily data to weekly data to take into account week/weekend effect
Sys.getlocale("LC_TIME")
Sys.setlocale("LC_TIME", "English")
data_Cuba_COVID$date_daily <- as.POSIXct(data_Cuba_COVID$Date, format = "%d-%b-%y")
data_Cuba_COVID$year <- format(data_Cuba_COVID$date_daily, format = "%Y")

#select only data of 2021
data_Cuba_COVID_2021 <- data_Cuba_COVID[data_Cuba_COVID$year == "2021",]

data_Cuba_COVID_2021$date_weekly <- paste(year(data_Cuba_COVID_2021$date_daily), week(data_Cuba_COVID_2021$date_daily), sep = "-")
data_Cuba_COVID_2021$week <- week(data_Cuba_COVID_2021$date_daily)

weekly_data_Cuba_COVID_2021 <- data_Cuba_COVID_2021 %>%
  group_by(ID, Municipality, date_weekly, week) %>%
  summarise(Imported = sum(Imported), Deaths = sum(Death), Cases = sum(confirmed_Cases),
            Cases_under18 = sum(age_under_18), Cases_19_39 = sum(age_19_39), 
            Cases_40_59 = sum(age_40_59), Cases_older60 = sum(age_60),
            Native = sum(native)) %>%
  ungroup()

#merge data set with population 
weekly_data_Cuba_COVID_2021_population <- left_join(weekly_data_Cuba_COVID_2021, Population_by_municipality, 
                           by = c("ID" = "ID"))
weekly_data_Cuba_COVID_2021_population <- weekly_data_Cuba_COVID_2021_population[order(weekly_data_Cuba_COVID_2021_population$ID, weekly_data_Cuba_COVID_2021_population$week),]

#COVID Data per 10 000 inhabitants
weekly_data_Cuba_COVID_2021_population$Imported_std <- (weekly_data_Cuba_COVID_2021_population$Imported / weekly_data_Cuba_COVID_2021_population$Total)*10000
weekly_data_Cuba_COVID_2021_population$Deaths_std <- (weekly_data_Cuba_COVID_2021_population$Deaths / weekly_data_Cuba_COVID_2021_population$Total)*10000

#summary statistics
mean(weekly_data_Cuba_COVID_2021_population$Imported_std)
min(weekly_data_Cuba_COVID_2021_population$Imported_std)
max(weekly_data_Cuba_COVID_2021_population$Imported_std)
median(weekly_data_Cuba_COVID_2021_population$Imported_std)

mean(weekly_data_Cuba_COVID_2021_population$Deaths_std)
min(weekly_data_Cuba_COVID_2021_population$Deaths_std)
max(weekly_data_Cuba_COVID_2021_population$Deaths_std)
median(weekly_data_Cuba_COVID_2021_population$Deaths_std)

###############################################
## EXPLORING THE SPATIAL AND TEMPORAL TRENDS ##
###############################################

#exploring spatial trend by looking at average for each municipality
Weekly_summaries_avg <- summarise(group_by(weekly_data_Cuba_COVID_2021_population,ID), 
                              Imported_std_avg = mean(Imported_std), 
                              Death_std_avg = mean(Deaths_std))

Weekly_summaries_avg <- summarise(group_by(weekly_data_Cuba_COVID_2021_population,ID), 
                                  Imported_std_avg = median(Imported_std), 
                                  Death_std_avg = median(Deaths_std))

cuba_df <- cuba_shape %>%
  left_join(Weekly_summaries_avg, by = c("ID" = "ID"))

ggplot(data = cuba_df) +
  geom_sf(aes(fill = Imported_std_avg)) +
  scale_fill_gradient(name="Imported cases", low = "white", high = "gray0")+ 
  ggtitle("Average number of imported cases per 10 000 inhabitants") + labs(fill = "Imported cases")

ggplot(data = cuba_df) +
  geom_sf(aes(fill = Death_std_avg)) +
  scale_fill_gradient(name="Deaths", low = "white", high = "gray0")+ 
  ggtitle("Average number of deaths per 10 000 inhabitants") + labs(fill = "Deaths")


#exploring temporal trend
ggplot(weekly_data_Cuba_COVID_2021_population, aes(x=week, y=Imported_std, color=as.factor(ID))) + 
  geom_point() +
  geom_smooth(se=FALSE) + 
  theme(legend.position="none") + ylab("Number of imported cases per 10 000 inhabitants")

ggplot(weekly_data_Cuba_COVID_2021_population, aes(x=week, y=Deaths_std, color=as.factor(ID))) + 
  geom_point() +
  geom_smooth(se=FALSE) + 
  theme(legend.position="none") + ylab("Number of deaths per 10 000 inhabitants")

#create a lag for analyses
weekly_data_Cuba_COVID_2021_population_lag <- weekly_data_Cuba_COVID_2021_population %>% 
  group_by(ID) %>%
  mutate(deaths_lag1 = lead(Deaths, order_by = ID, n=1),
         deaths_lag2 = lead(Deaths, order_by = ID, n=2),
         deaths_lag3 = lead(Deaths, order_by = ID, n=3),
         deaths_lag4 = lead(Deaths, order_by = ID, n=4),
         deaths_lag5 = lead(Deaths, order_by = ID, n=5),
         deaths_lag6 = lead(Deaths, order_by = ID, n=6),
         deaths_lag7 = lead(Deaths, order_by = ID, n=7),
         deaths_lag8 = lead(Deaths, order_by = ID, n=8),
         deaths_lag9 = lead(Deaths, order_by = ID, n=9),
         deaths_lag10 = lead(Deaths, order_by = ID, n=10),
         deaths_lag11 = lead(Deaths, order_by = ID, n=11),
         deaths_lag12 = lead(Deaths, order_by = ID, n=12))
