chr.dir.r.files <- "m:/models/bacteria/hspf/bigelkhydrocal201601/r-files"
source(file = paste0(chr.dir.r.files, "/functions.R"))

chr.dir.hspf <- "m:/models/bacteria/hspf/bigelkhydrocal201601/hspf-files"
chr.file.hspf.out <- "behydcal.out"

## drainage area in sqr mi for Big Elk Creek at outlet
da.be <- 88.8

df.mod <- rpltgen(chr.dir = chr.dir.hspf, chr.file = chr.file.hspf.out)

## calculate model groups
## mlog log base 10 of stream flow
mlog <- log10(df.mod$Rch18.flow)

## mflow stream flow
mflow <- df.mod$Rch18.flow

## mbaseind, base flow index 
## baseflow seperation using USGS-HySep R version in DVstats
require(DVstats)
df.hysep88.8 <- hysep(Flow = df.mod$Rch18.flow , 
                      Dates = as.Date(df.mod$tmp.date), da = da.be)

## caclulate baseflow in the same way as calculated from the observed data
mbaseind <- sum(df.hysep88.8$BaseQ)/sum(df.hysep88.8$Flow)

## clean up
rm(df.hysep88.8)

## convert stream flow from cu ft / sec to ac-ft / day for use in volumes
## (1 cu ft / sec) * (86400 sec / day) * (1 ac-ft / 43559.9 cu ft)
df.mod <- cbind(df.mod, 
                flow.ac.ft =  86400 * (1 / 43559.9) * df.mod$Rch18.flow)

## need doBy package to sums for annual, summer and winter
require(doBy)

## mvol_ann - annual volumes in ac-ft
## create factor for year
df.mod <- cbind(df.mod, 
                fac.ann  = as.factor(
                  strftime(df.mod$tmp.date, format = "%Y")))
mvol_ann <- summaryBy(flow.ac.ft ~ fac.ann, data = df.mod, FUN = sum)

## create factor for month used in mvol_smr and mvol_wtr calculations
df.mod <- cbind(df.mod, 
                fac.mon  = as.factor(
                  strftime(df.mod$tmp.date, format = "%b")))

## create factor for month used in mvol_smr and mvol_wtr calculations

df.tmp <- data.frame(mon=as.character(df.mod$fac.mon), season = "none", 
                     stringsAsFactors = FALSE)

## summer season, summer is Jun, Jul and Aug
lng.smr <- grep("Jun|Jul|Aug", df.tmp$mon)

## winter season
lng.wtr <- grep("Dec|Jan|Feb", df.tmp$mon)

## assign summer and winter values to season. leave spring and fall as none
df.tmp$season[lng.smr] <- "summer"
df.tmp$season[lng.wtr] <- "winter"

## add season as factor to df.mod 
df.mod <- data.frame(df.mod, fac.season = as.factor(df.tmp$season))

## clean up
rm(df.tmp, lng.smr, lng.wtr)

df.vol.seasons <- summaryBy(flow.ac.ft ~ fac.ann + fac.season , data = df.mod, FUN = sum)


## mvol_smr - summer volumes in ac-ft
mvol_smr <- df.vol.seasons[as.character(df.vol.seasons$fac.season) == "summer", 
                           c("fac.ann", "flow.ac.ft.sum")]

## mvol_wtr
mvol_wtr <- df.vol.seasons[as.character(df.vol.seasons$fac.season) == "winter", 
                           c("fac.ann", "flow.ac.ft.sum")]
## storm information
## get storm dates from text file Information in this file from 
## Select_Storm_HydCal repo
## column 2 is the begin date of storm and column 8 is the end date of storm


df.strm.dates <- read.delim(file = paste0(chr.dir.r.files, "/dates_stm.dat"),
                            header = FALSE, sep = " ", 
                            stringsAsFactors = FALSE)[ , c(2, 8)]



## mpeak

## mvol_stm

## mtime

