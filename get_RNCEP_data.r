library(survival)
library(RNCEP)

bp <- c(43.4617827,-65.7567228)
seal <- c(43.4083358,-66.0250048)

site.ll <- tribble(~bandsite,  ~lat,        ~lon,
                  'BP',    43.4617827,  -65.7567228,
                  'SEAL',  43.4083358,  -66.0250048)

w.variable <- "tcdc.eatm"
w.level <- 'gaussian'
w.months.minmax <- c(8, 9)
w.years.minmax <- c(2016, 2017)
w.lat.southnorth <- c(43.39, 43.49)
w.lon.westeast <- c(-66.02, -65.74)

w.data <- NCEP.gather(w.variable, w.level, w.months.minmax, w.years.minmax, w.lat.southnorth, w.lon.westeast)


w.aggregate.data <- NCEP.aggregate(w.data)

NCEP.vis.area(w.aggregate.data)

depFlights$ts.dep

site.ll[site.ll$site == flightBySunset$bandsite,]$lat

flightByWeather <- flightBySunset %>% left_join(site.ll, by = 'bandsite')

get.interp.data <- function(w.variable, w.level, lats, lons, times) {
  NCEP.interp(w.variable, w.level, lats, lons, times)
}

flightBySunset.cc <- get.interp.data('tcdc.eatm', 'gaussian', flightBySunset$lat, flightBySunset$lon, flightBySunset$ts) # cloud cover

flightBySunset.vwind <- get.interp.data('vwnd.sig995', 'surface', flightBySunset$lat, flightBySunset$lon, flightBySunset$ts) # lat wind

flightBySunset.uwind <- get.interp.data('uwnd.sig995', 'surface', flightBySunset$lat, flightBySunset$lon, flightBySunset$ts) # lon wind


flightBySunset.wind.abs <- sqrt((flightBySunset.uwind^2) + (flightBySunset.vwind^2))

flightBySunset.wind.dir <- atan2((flightBySunset.uwind/flightBySunset.wind.abs), (flightBySunset.vwind/flightBySunset.wind.abs))*(180/pi) + 180



flightByWeather$vwind <- flightBySunset.vwind
flightByWeather$uwind <- flightBySunset.uwind
flightByWeather$wind.abs <- flightBySunset.wind.abs * (3.6) # Convert m/s to km/h
flightByWeather$wind.dir <- flightBySunset.wind.dir


flightByWeather$cc <- flightBySunset.cc



#### Cloud cover

flightByWeather %>%
  filter(isDep) %>% 
  ggplot(aes(bandsite, (cc)))+
  geom_boxplot()+
  facet_grid(age~.)+
  ylab('Cloud Cover (%)')+
  xlab('Banding site')

flightByWeather %>%
  filter(isDep) %>% 
  ggplot(aes((cc))) +
  geom_histogram()+
  facet_grid(age~bandsite)

glm7 <- flightByWeather %>% filter(isDep) %>% glm(formula = (cc) ~ age * bandsite, family = gaussian)

plot(glm7)

summary(glm7)

anova(glm7, test = 'F')


#### Wind Velocity

flightByWeather %>%
  filter(isDep) %>% 
  ggplot(aes(bandsite, log(wind.abs)))+
  geom_boxplot()+
  facet_grid(age~.)+
#  ylab('Cloud Cover (%)')+
  xlab('Banding site')

flightByWeather %>%
  filter(isDep) %>% 
  ggplot(aes(log(wind.abs))) +
  geom_histogram()+
  facet_grid(age~bandsite)

glm8 <- flightByWeather %>% filter(isDep) %>% glm(formula = log(wind.abs) ~ age * bandsite, family = gaussian)

plot(glm8)

summary(glm8)

anova(glm8, test = 'F')


#### Wind direction

flightByWeather %>%
  filter(isDep) %>% 
  ggplot(aes(bandsite, (wind.dir)))+
  geom_boxplot()+
  facet_grid(age~.)+
#  ylab('Cloud Cover (%)')+
  xlab('Banding site')

flightByWeather %>%
  filter(isDep) %>% 
  ggplot(aes((wind.dir))) +
  geom_histogram()+
  facet_grid(age~bandsite)

glm9.1 <- flightByWeather %>% filter(isDep) %>% glm(formula = (wind.dir) ~ age * bandsite, family = gaussian)
glm9.2 <- flightByWeather %>% filter(isDep) %>% glm(formula = (wind.dir) ~ age + bandsite, family = gaussian)
glm9.3 <- flightByWeather %>% filter(isDep) %>% glm(formula = (wind.dir) ~ 1, family = gaussian)

plot(glm9)

summary(glm9.1)
summary(glm9.2)
summary(glm9.3)

anova(glm9.1, glm9.3, test = 'F')

anova(glm9, update(glm9, .~ -age:bandsite), test = 'F')
anova(glm9, update(glm9, .~ -age:bandsite), test = 'F')
