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

## mvol_ann
## mvol_smr
## mvol_wtr

## mpeak
## mvol_stm

## mtime

