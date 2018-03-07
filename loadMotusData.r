library(tidyverse) 
library(motus)

loadMotusData <- function (projectID, subfolder, newDatabase) {

hits.sql <- tagme(projectID, new = newDatabase, update = !newDatabase, forceMeta = TRUE, dir = subfolder)

hits.tbl <- tbl(hits.sql, "alltags")

tagMeta <- metadata(hits.sql, projectIDs = projectID)

hits.df <- select(hits.tbl, 
                  motusTagID, id = mfgID, hitID, runID, batchID, 
                  ts, sig, runLen, freqsd, sigsd, slop, burstSlop,
                  antType, antBearing, lat = gpsLat, lon = gpsLon, recv,
                  depLat = tagDeployLat, depLon = tagDeployLon, site = recvDeployName, markerNumber, spEN = speciesEN) %>% 
  distinct() %>% collect() %>% as.data.frame() %>%
  mutate(ts = as.POSIXct(ts, origin="1970-01-01"), 
         year = year(ts))

### Load RECEIVER METADATA ###
receiverData <- read_csv(paste0(subfolder, "receiver-deployments.csv", collapse = '')) %>%
  mutate(
    tsStart = as.POSIXct(tsStart,origin = '1970-01-01'),
    tsEnd = as.POSIXct(tsEnd,origin = '1970-01-01'),
    deploymentName = gsub("<[^\\]]*>", "", deploymentName, perl = TRUE)
  ) %>%
  select(deploymentName, recv = receiverID, latitude, longitude, recvDepStart = tsStart, recvDepEnd = tsEnd) 

### Load Tag deployment METADATA ###
tagDeploymentData <- read_csv(paste0(subfolder, "tag-deployments.csv", collapse = '')) %>%
  mutate(
    id = as.factor(mfgID),
    tsStart = as.POSIXct(tsStart,origin = '1970-01-01'),
    tsEnd = as.POSIXct(tsEnd,origin = '1970-01-01'),
    age = ifelse(age == 2, 'Hatch-year', 'Adult')
  ) %>%
  select(markerNumber, tsStart, tsEnd, age, sex, tagDeployID, weight, wing) 

# Fix raw data dates/times and ID
rawDataComb <- hits.df %>%
  mutate(
    ts = as.POSIXct(ts,origin = '1970-01-01'),
    date = format(ts, "%Y-%m-%d"),
    # Fix any missing lat/lon values with appropriate deployment from receiverData file
    lat = ifelse(lat == 0 | is.na(lat), receiverData[receiverData$recv == recv & receiverData$recvDepStart <= ts & (is.na(receiverData$recvDepEnd) | receiverData$recvDepEnd >= ts), ]$latitude, lat),
    lon = ifelse(lat == 0 | is.na(lon), receiverData[receiverData$recv == recv & receiverData$recvDepStart <= ts & (is.na(receiverData$recvDepEnd) | receiverData$recvDepEnd >= ts), ]$longitude, lon)
  ) %>%
  mutate(
    # If there are still some missing lat/lon values, use latest deployment of respective receiver (based on receiver) to get lat/lon
    lat = ifelse(lat == 0 | is.na(lat), receiverData[receiverData$deploymentName == site, ]$latitude[which.max(receiverData[receiverData$deploymentName == site, ]$recvDepStart)], lat),
    lon = ifelse(lat == 0 | is.na(lon), receiverData[receiverData$deploymentName == site, ]$longitude[which.max(receiverData[receiverData$deploymentName == site, ]$recvDepStart)], lon)
  ) %>%
  left_join(tagDeploymentData, by = 'markerNumber')

# Fix data Again
rawDataComb %>%
  select(id, hitID, runID, motusTagID,
         sig, markerNumber, batchID, 
         freqsd, sigsd, slop, burstSlop,
         antType, antBearing, 
         recv, tagDeployID,
         ts, date, site,
         age, sex, wing, weight, 
         depLon, depLat, lon, lat, 
         runLen, spEN) %>%
  mutate(
    date = format(ts,"%Y-%m-%d"),
    hour = format(ts,"%Y-%m-%d %H"),
    age = ifelse(is.na(age), 'Unknown', age),
    year = format(ts, '%Y'),
    recvDeployName = site,
    tagDeployIDmarkerNumber
  )
}
