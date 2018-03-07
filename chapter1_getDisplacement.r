## Load libraries

library(tidyverse)
library(motus)

library(geosphere)
library(circular)

source("loadMotusData.r")
source("latLonDist.r")

########
### Import datasets into list
########

# Set viewport to 2 x 2 grid
par(mfrow=c(2,2))

# Set timezone to GMT
Sys.setenv(TZ = 'GMT')

# Define Project ID
projectID <- 109

# Local address of database
subfolder <- "data/"

# Create a new SQLITE database?
newDatabase <- F

# Load tag hits
tagHits <- loadMotusData(projectID, database, newDatabase)

# Create a vectors of BP and Seal sites
bpsites <- c('BPHill','BPhill','BPLH','BPwestomni', 'BPLab', 'BPNorth', 'BPnorthomni', 'BPnorthyagi')
sealsites <- c('SealSouth','Seal South','SealWest','SealNorth')
shitsites <- read_csv(paste0(database, '109-shitSites.csv', collapse = ''), col_names = 'site')

fixedHits <- tagHits %>%
  rename(recvDeployName = recvDepName)  %>% 
  filter(runLen > 2, 
         !recvDeployName %in% shitsites$site, 
         !is.na(tagDeployID) & motusTagID != 21273,
         !is.na(lat),
         !is.na(lon),
         lat > 41, lon >-71) %>% 
  mutate(recvDeployName = ifelse(recvDeployName %in% bpsites, 'Bon Portage', ifelse(recvDeployName %in% sealsites, 'Seal Island', recvDeployName)))

# Get list of ages
tagAges <- fixedHits %>% group_by(markerNumber) %>% summarise(age = paste(age[1], sep = ''))

# Load list of flight runIDs
tagMeta <- read_csv(paste0(subfolder, "Project109-sigplots-Flights.csv", collapse = '')) %>%
  select(markerNumber, bandsite) %>%
  left_join(tagAges, by = 'markerNumber')

# Only select receivers from BP and Seal
allSites <- siteTrans(fixedHits, latCoord = 'lat', lonCoord = 'lon')

# Plot of latitude over julian date
allSites %>%
  rename(markerNumber = tagDeployID) %>%
  left_join(tagMeta, by = 'markerNumber') %>%
  ggplot(aes(as.integer(format(ts.y, "%j")), lat.y, group = markerNumber, color = age)) +
  geom_line() +
  scale_x_continuous()+
  xlab('Julian Date')+
  ggtitle('Latitude by julian date')


# Plot julian date over longitude
allSites %>%
  rename(markerNumber = tagDeployID) %>%
  left_join(tagMeta, by = 'markerNumber') %>%
  ggplot(aes(lon.y, as.integer(format(ts.y, "%j")), group = markerNumber, color = age)) +
  geom_line() +
  scale_y_continuous()+
  ylab('Julian Date')+
  ggtitle('Julian date by longitude')

dCumumlative <- allSites %>%
  group_by(tagDeployID) %>%
  summarise(distCumulative = sum(dist))

dNet <- allSites %>%
  group_by(tagDeployID) %>%
  summarise(lat.x = lat.x[which.min(ts.x)], 
    lon.x = lon.x[which.min(ts.x)], 
    lat.y = lat.y[which.max(ts.y)], 
    lon.y = lon.y[which.max(ts.y)], 
    recvDeployName.x = recvDeployName.x[which.min(ts.x)], 
    recvDeployName.y = recvDeployName.y[which.max(ts.y)]) %>%
  rowwise() %>%
  mutate(dist = distHaversine(c(lon.x, lat.x), c(lon.y, lat.y)))

displacement <- allSites %>%
  group_by(tagDeployID) %>%
  summarise(distCumulative = sum(dist),
            lat.x = lat.x[which.min(ts.x)], 
            lon.x = lon.x[which.min(ts.x)], 
            lat.y = lat.y[which.max(ts.y)], 
            lon.y = lon.y[which.max(ts.y)], 
            recvDeployName.x = recvDeployName.x[which.min(ts.x)], 
            recvDeployName.y = recvDeployName.y[which.max(ts.y)],
            ts.y = max(ts.y),
            tot_ts = difftime(max(ts.y), min(ts.x), units = 'hours')) %>%
  rowwise() %>%
  mutate(distNet = latLonDist(lat.x, lon.x, lat.y, lon.y),
         netVCum = distNet/distCumulative,
         rate = distNet/as.integer(tot_ts)) %>%
  rename(markerNumber = tagDeployID) %>%
  filter(!recvDeployName.y %in% c('Bon Portage', 'Seal Island', 'Shag Harbour 2')) %>%
  left_join(tagMeta, by = 'markerNumber')


# Net v. cumulative displacement

displacement %>%
  ggplot(aes(netVCum)) +
  geom_histogram()

displacement %>%
  filter(!is.na(age)) %>%
  ggplot(aes(netVCum)) +
  geom_histogram() +
  facet_grid(age~bandsite)

displacement %>%
  filter(!is.na(age)) %>%
  ggplot(aes(bandsite, netVCum)) +
  geom_boxplot() +
  facet_grid(age~.)


##
## Models
##


# Test: does net v. cumulative displacement differ among ages and/or band site?

glm1 <- displacement %>% glm(formula = netVCum ~ age * bandsite, family = gaussian)

plot(glm1)

summary(glm1)

anova(glm1, test = "F")

################## Departures

# Load list of flight runIDs
flights.raw <- read_csv(paste0(subfolder, "Project109-sigplots-Flights.csv", collapse = ''))

flights.raw <- flights.raw %>%
  mutate(depRunID = ifelse(is.na(depRunID), max(flights.raw[flights.raw$markerNumber == markerNumber, 5:16]), depRunID))

# Get flight meta data
flight.meta <- flights.raw %>% select(bandsite, markerNumber, depRunID)

# Select all hits which correspond to each flight runID
flights.df <- tagHits[tagHits$runID %in% as.numeric(unlist(flights.raw[5:16])),]

# Group flights by runID and select hitID and ts that corresponds to max signals strength during that run
flights <- flights.df %>%
  group_by(runID) %>%
  summarise(ts = ts[which.max(sig)], markerNumber = markerNumber[1])

# Create a new tibble of just departure flights and join it with flight meta data
depFlights <- flights %>%
  left_join(flight.meta, by = 'markerNumber') %>%
  group_by(markerNumber) %>%
  summarise(depRunID = ifelse(is.na(depRunID), runID[which.max(ts)], max(depRunID))[1], 
            ts = ts[runID==depRunID]) %>%
  select(markerNumber, ts.dep = ts)

# Make new tibble with rate of displacement using departure time as start time
displacementRate <- displacement %>% 
  left_join(depFlights, by = 'markerNumber') %>%
  mutate(rate = distNet/as.integer(difftime(ts.y, ts.dep, unit = 'hours'))) %>%
  filter(year(ts.dep) == year(ts.y))

# Rate (log)

displacementRate %>%
  ggplot(aes(log(rate))) +
  geom_histogram()

displacementRate %>%
  filter(!is.na(age)) %>%
  ggplot(aes(log(rate))) +
  geom_histogram() +
  facet_grid(age~bandsite)

displacementRate %>%
  filter(!is.na(age)) %>%
  ggplot(aes(bandsite, log(rate))) +
  geom_boxplot() +
  facet_grid(age~.)



# Test: does rate of movement differ among ages and/or band site

glm2 <- displacementRate %>% glm(formula = log(rate) ~ age * bandsite, family = Gamma)

plot(glm2)

summary(glm2)

anova(glm2, test = "F")
