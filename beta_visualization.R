
library(dplyr)
library(lubridate)
library(readr)
library(deSolve)
library(ggplot2)
library(scales)
library(tidyverse)
library(rstudioapi)


directory <- dirname(rstudioapi::getSourceEditorContext()$path) # doesn't work if ran in Terminal

#reading in dataset

cases <- read.csv(file.path(directory, "data/cases_by_status_and_phu.csv"))

#grouping all PHU's together

cases_ont <- cases %>% 
  group_by(FILE_DATE) %>% 
  summarise(active = sum(ACTIVE_CASES), resolved = sum(RESOLVED_CASES), deaths = sum(DEATHS))

#reading in dataset with new daily cases

new_cases <- read.csv(file.path(directory, "data/new_cases.csv"))
new_cases$date <- as.Date(new_cases$category, format = "%a %b %d %Y")
new_cases <- new_cases %>% filter(date >= as.Date('2020-04-01'))

#combining two datasets

df <- data.frame(cases_ont, new = new_cases$New.cases)

#making new variables to calculate proportions

N <- 14.57*1000000 # population of ontario

df <- df %>% 
  filter(FILE_DATE >= as.Date('2020-04-01') & FILE_DATE <= as.Date('2021-05-22')) %>%
  mutate(
    act.p = active / N,
    res.p = resolved/N,
    death.p = deaths/N,
    new.p = new/N
  )

#function to calculate beta

sir.calc <- function(date, t=1, N=14570000){
  Inew <- df[df$FILE_DATE == date, c("new")]  # number of new infections since last sample
  I <- df[df$FILE_DATE == date, c("active")] # number of infected individuals
  R<- df[df$FILE_DATE == date, c("resolved")] # number of recovered 
  D<- df[df$FILE_DATE == date, c("deaths")] # number of deaths
  S <- N - I - R - D
  beta <- (-1/t) * log(1-Inew*((1/I)+(1/S)))
}

#make new column for beta by date

df$beta <- sir.calc(df$FILE_DATE)

## Add columns
# Moving average
library(runner); library(psych)
df <- df %>% 
  mutate(Date = as.Date(FILE_DATE),
         roll_mean  = zoo::rollmean(beta, k = 7, fill = NA), # seven day moving average
         roll_geo = runner(beta, k = 7, idx = Date, f = geometric.mean))  # seven day geometric moving average

matplot(df$Date, df[c('roll_mean','roll_geo')], type = 'l', lty = 1) # not much difference between arithmetic and geometric moving average 

# Temperature data
temperature <- read.csv2('data/en_climate_daily_ON_6158355_2020_2021_P1D.csv', sep=',')
temperature$Date <- as.Date(temperature$Date.Time, format="%m/%d/%Y")
temperature$DailyMean <- as.numeric(temperature$Mean.Temp..Â.C.)

# merge temperature data to df
df <- df %>% left_join(temperature[c('Date','DailyMean')], by = 'Date')

#graphing proportion active cases and beta

coeff = 10

springIdx <- grep('2020-03-20|2020-06-20',df$FILE_DATE)
springIdx <- c(0, springIdx)
summerIdx <- grep('2020-06-20|2020-09-21',df$FILE_DATE)
fallIdx <- grep('2020-09-21|2020-12-20',df$FILE_DATE)
winterIdx <- grep('2020-12-20|2021-03-20',df$FILE_DATE)
springIdxNew <- grep('2021-03-20|2021-06-20',df$FILE_DATE)
springIdxNew <- c(springIdxNew, nrow(df))

ggplot(df, aes(x = FILE_DATE))+
  geom_point(aes(y = act.p)) +
  annotate("rect", xmin = springIdx[1], xmax = springIdx[2], ymin = -Inf, ymax = Inf,
           alpha = 0.2, fill = "green1") +
  annotate("rect", xmin = summerIdx[1], xmax = summerIdx[2], ymin = -Inf, ymax = Inf,
           alpha = 0.1,  fill = "salmon") +
  annotate("rect", xmin = fallIdx[1], xmax = fallIdx[2], ymin = -Inf, ymax = Inf,
           alpha = 0.3,  fill = "pink") +
  annotate("rect", xmin = winterIdx[1], xmax = winterIdx[2], ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "blue1") +
  annotate("rect", xmin = springIdxNew[1], xmax = springIdxNew[2], ymin = -Inf, 
           ymax = Inf, alpha = 0.2, fill = "green1") +
  geom_point(aes(y = act.p)) +
  geom_point(aes(y = roll_mean/coeff), color='violetred') +
  scale_y_continuous(name = 'proportion active cases', 
                     sec.axis = sec_axis(~.*coeff, name = 'beta'))+
  theme(axis.text.x = element_text(angle = 90))+
  geom_vline(xintercept = c('2020-09-08'), color = 'blue', size = .25) + #schools open
  geom_vline(xintercept = c('2020-09-17'), color = 'blue', size = .25) + #gathering limits reduced
  geom_vline(xintercept = c('2020-09-25'), color = 'blue', size = .25) + #restrictions on opening hours for bars/restaurants
  geom_vline(xintercept = c('2020-09-28'), color = 'blue', size = .25) + #indoor dining banned
  geom_vline(xintercept = c('2020-10-02'), color = 'blue', size = .25) + #capacity limits on bars/ gyms
  geom_vline(xintercept = c('2020-10-10'), color = 'blue', size = .25) + #Peel/ ottawa/ toronto back to stage 2
  geom_vline(xintercept = c('2020-10-16'), color = 'blue', size = .25) + #york moved to stage 2
  geom_vline(xintercept = c('2020-11-06'), color = 'blue', size = .25) + #Peel moved to red tier
  geom_vline(xintercept = c('2020-11-07'), color = 'blue', size = .25) + #Shift to 5 tier ssystem to determine restrictions
  geom_vline(xintercept = c('2020-11-14'), color = 'blue', size = .25) + #toronto moved to red tier
  geom_vline(xintercept = c('2020-11-16'), color = 'blue', size = .25) + # other regions moved to red tier
  geom_vline(xintercept = c('2020-11-23'), color = 'blue', size = .25) + #all non-essential businessess closed
  geom_vline(xintercept = c('2020-12-14'), color = 'blue', size = .25) + #york, windsor moved to lockdown
  geom_vline(xintercept = c('2020-12-15'), color = 'blue', size = .25) + #pfizer vaccine administration begins
  geom_vline(xintercept = c('2020-12-21'), color = 'blue', size = .25) + #hamilton moves to lockdown
  geom_vline(xintercept = c('2020-12-26'), color = 'red', size = 1.0) + #provincewide shutdown
  geom_vline(xintercept = c('2021-01-14'), color = 'blue', size = .25) + #additional stay at home orders 
  geom_vline(xintercept = c('2021-01-29'), color = 'blue', size = .25) + #new travel restrictions - individuals on foreign flights required to take PCR test
  geom_vline(xintercept = c('2021-02-10'), color = 'blue', size = .25) + #stayat home orders lifted in some areas
  geom_vline(xintercept = c('2021-03-03'), color = 'blue', size = .25) + #first AZ vaccines distributed in torotno windsor kingston
  geom_vline(xintercept = c('2021-03-05'), color = 'blue', size = .25) + #ontario exits stay at home orders.  toronto/ peel enter lockdown
  geom_vline(xintercept = c('2021-04-08'), color = 'red', size = 1.0) + #second stay at home order, vaccinations by postal codes
  scale_x_discrete(breaks = df$FILE_DATE[grepl('01$',df$FILE_DATE)]) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab("Date")



