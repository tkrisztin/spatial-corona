rm(list=ls())

library(spdep)
library(MASS)
library(expm)
library(dplyr)
library(readxl)
library(spdep)
library(countrycode)
library(sp)
require(sf)
library(sna)
library(tidyr)
library(readr)
library(zoo)
library(ggplot2)
library(stringr)
library(scales)
require(dplyr)


#------------------------------------------------------------------------------------------
# Data aggregation

cases <- read_excel("./Data/Covid19_28Mar.xlsx", sheet = "confirmed") %>%
  rename("Date" = ...1) %>%
  pivot_longer(names_to = "countries",cols = -1) %>%
  mutate(
    iso3c = countrycode(countries,origin = "country.name",destination = "iso3c")
  ) %>%
  filter(!is.na(iso3c)) %>%
  pivot_wider(names_from = Date,values_fill = 0)


#############################################################################
##### other W data

Wdata <- read_excel("./Data/otherWdata.xlsx")[,c(1,2,4,5,7)]

Wdata$importer <- str_replace(Wdata$importer,pattern = "China HK SAR", replacement = "Hongkong")
Wdata$exporter <- str_replace(Wdata$exporter,pattern = "China HK SAR", replacement = "Hongkong")
Wdata$importer <- str_replace(Wdata$importer,pattern = "Untd Arab Em", replacement = "United Arab Emirate")
Wdata$exporter <- str_replace(Wdata$exporter,pattern = "Untd Arab Em", replacement = "United Arab Emirate")
Wdata$importer <- str_replace(Wdata$importer,pattern = "Cent.Afr.Rep", replacement = "Central African Republic")
Wdata$exporter <- str_replace(Wdata$exporter,pattern = "Cent.Afr.Rep", replacement = "Central African Republic")
Wdata$importer <- str_replace(Wdata$importer,pattern = "Dominican Rp", replacement = "Dominican Republic")
Wdata$exporter <- str_replace(Wdata$exporter,pattern = "Dominican Rp", replacement = "Dominican Republic")

Wdata$importer <- countrycode(Wdata$importer,origin = "country.name",destination = "iso3c")
Wdata$exporter <- countrycode(Wdata$exporter,origin = "country.name",destination = "iso3c")

W.language <- Wdata[,c(1,2,5)]
W.language <- na.omit(W.language)
W.language <- data.frame(spread(W.language, importer, LanguageTies)[,-1])
W.language[is.na(W.language)] <- 0
W.language <- as.matrix(W.language)
subs1 <- colnames(W.language)
subs1 <- subset(subs1, subs1!="SYR")

W.FTA <- Wdata[,c(1,2,4)]
W.FTA <- na.omit(W.FTA)
W.FTA <- data.frame(spread(W.FTA, importer, FTA)[,-1])
W.FTA[is.na(W.FTA)] <- 0
W.FTA <- as.matrix(W.FTA)

W.currency <- Wdata[,c(1,2,3)]
W.currency <- na.omit(W.currency)
W.currency <- data.frame(spread(W.currency, importer, SameCurrency)[,-1])
W.currency[is.na(W.currency)] <- 0
W.currency <- as.matrix(W.currency)

#############################################################################
##### flight data

flight <- read_csv("./Data/openflights_data.csv")[,3:5]

flight$source <- countrycode(flight$source,origin = "country.name",destination = "iso3c")
flight$destination <- countrycode(flight$destination,origin = "country.name",destination = "iso3c")
flight <- aggregate(flight, by=list(flight$source, flight$destination), FUN=mean)[,c(1,2,5)]


Wflight <- data.frame(spread(flight, Group.1, nr_flights)[,-1])
Wflight$NIU <- NULL
Wflight[is.na(Wflight)] <- 0
Wflight <- as.matrix(Wflight)

#############################################################################
##### flight data robust expo


Wflight.rob <- Wflight
Wflight.rob[Wflight.rob!=0] <- exp(Wflight.rob[Wflight.rob!=0])

#############################################################################
##### trade data

export <- read.csv("./Data/Trade_Export.csv")[,c(3,6,10)]
export <- aggregate(export$TradeValue,by=list(export$ReporterISO3,export$PartnerISO3),FUN=sum)
export <- export[as.character(export$Group.1)!=as.character(export$Group.2),]

import <- read.csv("./Data/Trade_Import.csv")[,c(3,6,10)]
import <- aggregate(import$TradeValue,by=list(import$ReporterISO3,import$PartnerISO3),FUN=sum)
import <- import[as.character(import$Group.1)!=as.character(import$Group.2),]

import <- import[,c(2,1,3)]
colnames(import) <- colnames(export)

trade <- rbind(export,import)
trade[,c(1,2)] <- apply(trade[,c(1,2)],2,as.character)
trade <- aggregate(trade$x, by=list(trade$Group.1, trade$Group.2), FUN = mean)
trade <- trade[trade$Group.1%in%cases$iso3c & trade$Group.2%in%cases$iso3c,]
trade.rob <- trade

trade <- trade %>%
  group_by(Group.1) %>%
  top_n(7, x) 


Wtrade <- data.frame(spread(trade, Group.1, x), stringsAsFactors = FALSE)
Wtrade <- Wtrade[Wtrade$Group.2 %in% colnames(Wtrade),]
Wtrade <- Wtrade[,colnames(Wtrade) %in% Wtrade$Group.2]
Wtrade[is.na(Wtrade)] <- 0
Wtrade <- as.matrix(Wtrade)


#############################################################################
##### trade data robust

Wtrade.rob <- data.frame(spread(trade.rob, Group.1, x), stringsAsFactors = FALSE)
Wtrade.rob <- Wtrade.rob[Wtrade.rob$Group.2 %in% colnames(Wtrade.rob),]
Wtrade.rob <- Wtrade.rob[,colnames(Wtrade.rob) %in% Wtrade.rob$Group.2]
Wtrade.rob[is.na(Wtrade.rob)] <- 0
Wtrade.rob <- as.matrix(Wtrade.rob)



Wtrade.rob.ex <- Wtrade.rob/1000000
Wtrade.rob.ex[Wtrade.rob.ex!=0] <- exp(Wtrade.rob.ex[Wtrade.rob.ex!=0])




subs2 <- colnames(Wtrade)
subs2 <- intersect(subs2, subs1)
subs2 <- intersect(subs2, cases$iso3c)

##### get population data from worldbank




###### Shape file

shp <- st_read("./Data/world_shape/CNTR_RG_60M_2016_3035.shp", stringsAsFactors=FALSE)

shp = shp %>% rename(iso3c = ISO3_CODE) %>%
  left_join(cases,by = join_by(iso3c)) %>%
  filter(iso3c %in% subs2)

coords <-st_coordinates(shp)
subs <- shp$iso3c

######################################
#### sorting matrices

rownames(Wflight) <- colnames(Wflight)
rownames(Wtrade) <- colnames(Wtrade)
rownames(Wflight.rob) <- colnames(Wflight.rob)
rownames(Wtrade.rob) <- colnames(Wtrade.rob)
rownames(Wtrade.rob.ex) <- colnames(Wtrade.rob.ex)
rownames(W.language) <- colnames(W.language)
rownames(W.FTA) <- colnames(W.FTA)
rownames(W.currency) <- colnames(W.currency)

Wflight <- Wflight[rownames(Wflight)%in%subs, colnames(Wflight)%in%subs]
Wtrade <- Wtrade[rownames(Wtrade)%in%subs, colnames(Wtrade)%in%subs]
W.language <- W.language[rownames(W.language)%in%subs, colnames(W.language)%in%subs]
W.FTA <- W.FTA[rownames(W.FTA)%in%subs, colnames(W.FTA)%in%subs]
W.currency <- W.currency[rownames(W.currency)%in%subs, colnames(W.currency)%in%subs]
Wflight.rob <- Wflight.rob[rownames(Wflight.rob)%in%subs, colnames(Wflight.rob)%in%subs]
Wtrade.rob <- Wtrade.rob[rownames(Wtrade.rob)%in%subs, colnames(Wtrade.rob)%in%subs]
Wtrade.rob.ex <- Wtrade.rob.ex[rownames(Wtrade.rob.ex)%in%subs, colnames(Wtrade.rob.ex)%in%subs]

Wtrade <- Wtrade[match(subs,rownames(Wtrade)),match(subs,colnames(Wtrade))]
Wflight <- Wflight[match(subs,rownames(Wflight)),match(subs,colnames(Wflight))]
W.language <- W.language[match(subs,rownames(W.language)),match(subs,colnames(W.language))]
W.FTA <- W.FTA[match(subs,rownames(W.FTA)),match(subs,colnames(W.FTA))]
W.currency <- W.currency[match(subs,rownames(W.currency)),match(subs,colnames(W.currency))]
Wtrade.rob <- Wtrade.rob[match(subs,rownames(Wtrade.rob)),match(subs,colnames(Wtrade.rob))]
Wflight.rob <- Wflight.rob[match(subs,rownames(Wflight.rob)),match(subs,colnames(Wflight.rob))]
Wtrade.rob.ex <- Wtrade.rob.ex[match(subs,rownames(Wtrade.rob.ex)),match(subs,colnames(Wtrade.rob.ex))]
##### design of W matrices

##Queen contiguity matrix
queen_nb <- poly2nb(shp, row.names = shp$data$Id, queen = T, snap = 0.0001)  #creates a neighborhoodlist
W.list.queen <- nb2listw(queen_nb, style = "W", zero.policy = TRUE) #creates a weights-list
W.queen <- listw2mat(W.list.queen) #creates a weigths matrix
# W.queen.sq <- W.queen %*% W.queen
# diag(W.queen.sq) <- 0
# W.queen.sq <- make.stochastic(W.queen.sq,'row')



wmats <- c("Queen","language","currency" ,"FTA" , "flight", "trade", "flight.rob", "trade.rob","trade.rob.ex")

W.language <- make.stochastic(W.language,'row')
Wflight <- make.stochastic(Wflight,'row')
Wtrade <- make.stochastic(Wtrade,'row')
W.FTA <- make.stochastic(W.FTA,'row')
W.currency <- make.stochastic(W.currency,'row')
Wflight.rob <- make.stochastic(Wflight.rob,'row')
Wtrade.rob <- make.stochastic(Wtrade.rob,'row')
Wtrade.rob.ex <- make.stochastic(Wtrade.rob.ex,'row')


W.list <- list()
W.list[[1]] <- W.queen
W.list[[2]] <- W.language
W.list[[3]] <- W.currency
W.list[[4]] <- W.FTA
W.list[[5]] <- Wflight
W.list[[6]] <- Wtrade
W.list[[7]] <- Wflight.rob
W.list[[8]] <- Wtrade.rob
W.list[[9]] <- Wtrade.rob.ex
names(W.list) = wmats

save(W.list,shp,file="./Estimation/xy.RData")
