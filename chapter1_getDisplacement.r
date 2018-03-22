## Load libraries

library(tidyverse)
library(motus)

library(geosphere)
library(circular)

library(mgcv)
library(gamm4)
library(lmtest)

library(RNCEP) # weather data
library(survival) # survival (proportional hazards)

library(itsadug)

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
tagHits <- loadMotusData(projectID, subfolder)

# Create a vectors of BP and Seal sites
bpsites <- c('BPHill','BPhill','BPLH','BPwestomni', 'BPLab', 'BPNorth', 'BPnorthomni', 'BPnorthyagi')
sealsites <- c('SealSouth','Seal South','SealWest','SealNorth')
shitsites <- read_csv(paste0(subfolder, '109-shitSites.csv', collapse = ''), col_names = 'site')

fixedHits <- tagHits %>%
#  rename(recvDeployName = recvDepName)  %>% 
  filter(runLen > 2, 
         !recvDeployName %in% shitsites$site, 
         !is.na(tagDeployID) & motusTagID != 21273,
         !is.na(lat),
         !is.na(lon)
         ) %>% 
  mutate(recvDeployName = ifelse(recvDeployName %in% bpsites, 'Bon Portage', ifelse(recvDeployName %in% sealsites, 'Seal Island', recvDeployName)))

# Remove the unnecessary and HUGE dataframe
rm(tagHits)
# Clear out the memory of unnecessary garbage
gc()

# Get list of ages
tagMeta <- fixedHits %>%
  mutate(bandsite = ifelse(depLat > 43.45, 'BP', 'SEAL')) %>% 
  group_by(markerNumber, tagDeployID) %>% 
  summarise(age = as.factor(age[1]), bandsite = as.factor(bandsite[1]), tagDeployStart = tagDeployStart[1])

#write_rds(tagMeta, 'tagMeta.rds')

################## Departures

# Load list of flight runIDs
flights.raw <- read_csv(paste0(subfolder, "Project109-sigplots-Flights.csv", collapse = ''))

flights.raw <- flights.raw %>%
  mutate(depRunID = ifelse(is.na(depRunID), max(flights.raw[flights.raw$markerNumber == markerNumber, 5:16]), depRunID))

# Get flight meta data
flight.meta <- flights.raw %>% select(bandsite, markerNumber, depRunID)

# Select all hits which correspond to each flight runID
flights.df <- fixedHits[fixedHits$runID %in% as.numeric(unlist(flights.raw[5:16])),]

# Group flights by runID and select hitID and ts that corresponds to max signals strength during that run
flights <- flights.df %>%
  group_by(runID) %>%
  summarise(ts = ts[which.max(sig)], markerNumber = markerNumber[1])

# Create a new tibble of just departure flights and join it with flight meta data
depFlights <- flights %>%
  left_join(flight.meta, by = 'markerNumber') %>%
  group_by(markerNumber) %>%
  summarise(depRunID = ifelse(is.na(depRunID), runID[which.max(ts)], max(depRunID))[1], 
            ts.dep = ts[runID==depRunID])# %>%
  select(markerNumber, ts.dep = ts)

################## END Departures

# Create a list of sites each bird located outside of their natal island
sitesVisited <- fixedHits %>%
  filter(
    lat > 42, lon > -71 & lon < 0,
    (lat <  43.39 | lat > 43.48) &
      (lon < -66.03 | lon > -65.74) &
      !recvDeployName %in% c('Bon Portage', 'Seal Island', 'Shag Harbour 2')
  ) %>%
  filter(spEN == "Swainson's Thrush") %>%
  left_join(depFlights, by = "markerNumber")  %>%
  filter(year(ts) == year(ts.dep)) %>%
  group_by(markerNumber) %>%
  summarise(nSites = length(unique(recvDeployName)),
            siteNames = paste0(unique(recvDeployName), collapse = ', ')) %>%
  left_join(tagMeta, by = 'markerNumber')

# Get some stats on the sites visited
sitesVisited %>%
  group_by(bandsite) %>%
  #group_by(age) %>%
  #group_by(age, bandsite) %>%
  summarise(mean = mean(nSites), sd = sd(nSites), n = n()) %>% View()

sitesVisited %>%
  ggplot(aes(nSites)) +
  geom_histogram()

sitesVisited %>%
  mutate(bandsite = ifelse(bandsite == 'SEAL', 'SI', 'BP')) %>%
  ggplot(aes(bandsite, nSites)) +
  geom_boxplot() +
  ylab('Number of sites visited')+
  xlab('Banding site')+
  facet_grid(age~.)

# Test: does the number of sites visited differ among ages and/or band site?

glm0 <- sitesVisited %>% glm(formula = (nSites) ~ age * bandsite, family = poisson)

plot(glm0)

summary(glm0)

anova(glm0, test = "Chisq")

glm0 <- sitesVisited %>% filter(bandsite == 'SEAL') %>% glm(formula = (nSites) ~ age, family = poisson)

plot(glm0)

summary(glm0)

anova(glm0, test = "Chisq")


glm0 <- sitesVisited %>% filter(bandsite == 'BP') %>% glm(formula = (nSites) ~ age, family = poisson)

plot(glm0)

summary(glm0)

anova(glm0, test = "Chisq")

# Create a new tibble of site transitions
allSites <- fixedHits %>%
  filter(#lat > 25, lon < 0
         lat > 42, lon > -71 & lon < 0
         ) %>%
  left_join(depFlights, by = "markerNumber")  %>%
  filter(recvDeployName != 'Shag Harbour 2', ts >= ts.dep, year(ts) == year(ts.dep)) %>%
  siteTrans(latCoord = 'lat', lonCoord = 'lon') %>%
  arrange(tagDeployID, ts.x) %>%
  left_join(tagMeta, by = "tagDeployID")
  #mutate(bandsite = ifelse(depLat > 43.45, 'BP', 'SEAL'))

gc()
  
#%>%
#  mutate(
#    isSimultaneous = (
#      recvDeployName.y == lag(recvDeployName.x) &
#      recvDeployName.x == lag(recvDeployName.y) & 
#      recvDeployName.y == lead(recvDeployName.x) &
#      recvDeployName.x == lead(recvDeployName.y) &
#      tagDeployID == lag(tagDeployID),
#    )
#  ) 

# Plot of latitude over julian date
allSites %>%
  ggplot(aes(as.integer(format(ts.y, "%j")), lat.y, group = markerNumber, color = age)) +
  geom_line() +
  scale_x_continuous()+
  xlab('Julian Date')+
  ggtitle('Latitude by julian date')


# Plot julian date over longitude
allSites %>%
  ggplot(aes(lon.y, as.integer(format(ts.y, "%j")), group = markerNumber, color = age)) +
  geom_line() +
  scale_y_continuous()+
  ylab('Julian Date')+
  ggtitle('Julian date by longitude')


# Plot latitude date over longitude
allSites %>%
  ggplot(aes(lon.y, lat.y, group = markerNumber, color = age)) +
  geom_line() +
  ggtitle('Latitude by longitude')+
  coord_cartesian()


# Create new dataframe and calculate cumulative and net displacement
displacement <- allSites %>%
 # filter(is.na(isSimultaneous) | !isSimultaneous) %>% #View()
  group_by(tagDeployID) %>%
  summarise(distCumulative = sum(dist),
            lat.x = lat.x[which.min(ts.x)], 
            lon.x = lon.x[which.min(ts.x)], 
            lat.y = lat.y[which.max(ts.y)], 
            lon.y = lon.y[which.max(ts.y)], 
            recvDeployName.x = recvDeployName.x[which.min(ts.x)], 
            recvDeployName.y = recvDeployName.y[which.max(ts.y)],
            ts.x = min(ts.x),
            ts.y = max(ts.y),
            ts.diff = difftime(max(ts.y), min(ts.x), units = 'hours')) %>%
  rowwise() %>%
  mutate(distNet = latLonDist(lat.x, lon.x, lat.y, lon.y),
         netVCum = distNet/distCumulative,
         rate = distNet/as.integer(ts.diff)) %>%
  filter(
    !is.na(rate) & rate != Inf, 
    recvDeployName.x != recvDeployName.y,
    netVCum <= 1
  ) %>%
  left_join(tagMeta, by = 'tagDeployID')


# Net v. cumulative displacement

displacement %>%
  ggplot(aes(netVCum)) +
  geom_histogram()

displacement %>%
#  filter(!is.na(age)) %>%
  ggplot(aes(netVCum, group = age, color = age)) +
  geom_density() +
  facet_grid(.~bandsite)

displacement %>%
  filter(!is.na(age)) %>%
  ggplot(aes(bandsite, netVCum)) +
  geom_boxplot() +
  facet_grid(age~.)

displacement %>%
  filter(!is.na(age)) %>%
  ggplot(aes(bandsite, distNet)) +
  geom_boxplot() +
  facet_grid(age~.)

displacement %>%
  filter(!is.na(age)) %>%
  ggplot(aes(bandsite, distCumulative)) +
  geom_boxplot() +
  facet_grid(age~.)

displacement %>%
  group_by(age) %>%
  summarise(mean.n = mean(distNet), sd.n = sd(distNet),
            mean.c = mean(distCumulative), sd.c = sd(distCumulative),
            mean.nvc = mean(netVCum), sd.nvc = sd(netVCum)) %>%
  View()

displacement %>%
  group_by(bandsite) %>%
  summarise(mean.n = mean(distNet), sd.n = sd(distNet),
            mean.c = mean(distCumulative), sd.c = sd(distCumulative),
            mean.nvc = mean(netVCum), sd.nvc = sd(netVCum)) %>%
  View()

##
## Models
##

# Test: does net v. cumulative displacement differ among ages and/or band site?

glm1 <- displacement %>% glm(formula = log(netVCum) ~ age * bandsite, family = gaussian)

plot(glm1)

summary(glm1)

anova(glm1, test = "F")


# Test: does net displacement differ among ages and/or band site?

glm1.2 <- displacement %>% glm(formula = log(distNet) ~ age * bandsite, family = gaussian)

plot(glm1.2)

summary(glm1.2)

anova(glm1.2, test = "F")


# Test: does cumulative displacement differ among ages and/or band site?

glm1.3 <- displacement %>% glm(formula = log(distCumulative) ~ age * bandsite, family = gaussian)

plot(glm1.3)

summary(glm1.3)

anova(glm1.3, test = "F")

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
  ggplot(aes(log(rate), group = age, color = age)) +
  geom_density() +
  facet_grid(.~bandsite)

displacementRate %>%
  filter(!is.na(age)) %>%
  ggplot(aes(bandsite, log(rate))) +
  geom_boxplot() +
  facet_grid(age~.)


displacementRate %>%
  group_by(age) %>%
  summarise(mean = mean(rate), sd = sd(rate)) %>%
  View()

displacementRate %>%
  group_by(bandsite) %>%
  summarise(mean = mean(rate), sd = sd(rate)) %>%
  View()

# Test: does rate of movement differ among ages and/or band site

glm4 <- displacementRate %>% glm(formula = log(rate) ~ age * bandsite, family = Gamma)

plot(glm4)

summary(glm4)

anova(glm4, test = "F")


# Make a new tibble of displacement rates for each movement made between sites
rateFromDeparture <- allSites %>%
  left_join(depFlights, by = "markerNumber")  %>%
#  left_join(select(displacementRate, markerNumber, age, bandsite, ts.dep, tagDeployID), by = "tagDeployID") %>%
  ungroup() %>%
  mutate(ts.dep.diff = as.integer(difftime(ts.x, ts.dep, units = 'hours')), tagDeployID = as.factor(tagDeployID)) %>%
#         age = ifelse(age == 'Adult', 1, 0),
#         bandsite = ifelse(bandsite == 'BP', 1, 0)) %>% 
  filter(rate < 25, !is.na(age), !is.na(bandsite))

# Visualise smooth functions
  #geom_point(aes(color = age))+
  #    geom_line(aes(color = age, group = markerNumber))+
  #geom_smooth(aes(color = age, group = age))+
  #geom_smooth(aes(color = bandsite, group = bandsite))

rateFromDeparture %>% select(ts.dep.diff, age, bandsite, rate) %>% filter(ts.dep.diff<=671)

# Make generalized additive mixed models for each variable

gamm1 <- gamm(rate ~ s(ts.dep.diff), data = rateFromDeparture)

gamm2 <- gamm(rate ~ s(ts.dep.diff, age, bs = 'fs'), data = rateFromDeparture)

gamm3 <- gamm(rate ~ s(ts.dep.diff, age, bs = 'fs') + bandsite, data = rateFromDeparture)

gamm4 <- gamm(rate ~ s(ts.dep.diff, bandsite, bs = 'fs'), data = rateFromDeparture)

gamm5 <- gamm(rate ~ s(ts.dep.diff, bandsite, bs = 'fs') + age, data = rateFromDeparture)

# Set viewport to 2 x 2 grid
par(mfrow=c(2,3))
par(mfrow=c(1,1))


p.gamm <- rateFromDeparture %>%
#  ggplot(aes(ts.dep.diff, rate))
  ggplot(aes(ts.dep.diff, rate, color = age))

p.gamm +
  geom_point() +
  geom_smooth(method = 'gam', formula = y ~ s(x))

p.gamm +
  stat_smooth(method = 'gam', formula = y ~ s(x, age, bs = 'fs'))

summary(gamm2$gam)

plot(gamm1$gam, xlab = 'Time (hours since departure)', ylab = 's(Time)')
#mtext(text = "b)", side = 3, adj = 0, at = c(-150))

title(main = 'gamm1')

plot(gamm2$gam)
title(main = 'gamm2')

plot(gamm3$gam)
title(main = 'gamm3')

plot(gamm4$gam)
title(main = 'gamm4')

plot(gamm5$gam)
title(main = 'gamm5')

# Compare linear mixed-effects of models with log-ratio test
lrtest(gamm1$lme, gamm2$lme)
lrtest(gamm2$lme, gamm3$lme)
lrtest(gamm1$lme, gamm4$lme)
lrtest(gamm4$lme, gamm5$lme)
lrtest(gamm2$lme, gamm5$lme)

gamm4.1 <- gamm4(rate ~ s(ts.dep.diff), data = rateFromDeparture)

gamm4.2 <- gamm4(rate ~ s(ts.dep.diff, age, bs = 'fs'), data = rateFromDeparture)

gamm4.3 <- gamm4(rate ~ s(ts.dep.diff, age, bs = 'fs') + s(ts.dep.diff, bandsite, bs = 'fs'), data = rateFromDeparture)

# Compare linear mixed-effects of models with log-ratio test
lrtest(gamm4.1$mer, gamm4.2$mer)
lrtest(gamm4.2$mer, gamm4.3$mer)
lrtest(gamm4.1$mer, gamm4.3$mer)