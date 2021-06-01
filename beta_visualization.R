
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

#graphing proportion active cases and beta

coeff = 10

ggplot(df, aes(x = FILE_DATE))+
  geom_point(aes(y = act.p)) +
  geom_point(aes(y = beta/coeff)) +
  scale_y_continuous(name = 'proportion active cases', sec.axis = sec_axis(~.*coeff, name = 'beta'))+
  theme(axis.text.x = element_text(angle = 90))+
  geom_vline(xintercept = c('2020-09-08'), color = 'blue', size = .5) + #schools open
  geom_vline(xintercept = c('2020-09-17'), color = 'blue', size = .5) + #gathering limits reduced
  geom_vline(xintercept = c('2020-09-25'), color = 'blue', size = .5) + #restrictions on opening hours for bars/restaurants
  geom_vline(xintercept = c('2020-09-28'), color = 'blue', size = .5) + #indoor dining banned
  geom_vline(xintercept = c('2020-10-02'), color = 'blue', size = .5) + #capacity limits on bars/ gyms
  geom_vline(xintercept = c('2020-10-10'), color = 'blue', size = .5) + #Peel/ ottawa/ toronto back to stage 2
  geom_vline(xintercept = c('2020-10-16'), color = 'blue', size = .5) + #york moved to stage 2
  geom_vline(xintercept = c('2020-11-06'), color = 'blue', size = .5) + #Peel moved to red tier
  geom_vline(xintercept = c('2020-11-07'), color = 'blue', size = .5) + #Shift to 5 tier ssystem to determine restrictions
  geom_vline(xintercept = c('2020-11-14'), color = 'blue', size = .5) + #toronto moved to red tier
  geom_vline(xintercept = c('2020-11-16'), color = 'blue', size = .5) + # other regions moved to red tier
  geom_vline(xintercept = c('2020-11-23'), color = 'blue', size = .5) + #all non-essential businessess closed
  geom_vline(xintercept = c('2020-12-14'), color = 'blue', size = .5) + #york, windsor moved to lockdown
  geom_vline(xintercept = c('2020-12-15'), color = 'blue', size = .5) + #pfizer vaccine administration begins
  geom_vline(xintercept = c('2020-12-21'), color = 'blue', size = .5) + #hamilton moves to lockdown
  geom_vline(xintercept = c('2020-12-26'), color = 'blue', size = .5) + #provincewide shutdown
  geom_vline(xintercept = c('2021-01-14'), color = 'blue', size = .5) + #additional stay at home orders 
  geom_vline(xintercept = c('2021-01-29'), color = 'blue', size = .5) + #new travel restrictions - individuals on foreign flights required to take PCR test
  geom_vline(xintercept = c('2021-02-10'), color = 'blue', size = .5) + #stayat home orders lifted in some areas
  geom_vline(xintercept = c('2021-03-03'), color = 'blue', size = .5) + #first AZ vaccines distributed in torotno windsor kingston
  geom_vline(xintercept = c('2021-03-05'), color = 'blue', size = .5) + #ontario exits stay at home orders.  toronto/ peel enter lockdown
  geom_vline(xintercept = c('2021-04-08'), color = 'blue', size = .5) #second stay at home order, vaccinations by postal codes
