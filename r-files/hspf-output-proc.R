## load packages
library(DVstats, quietly = TRUE) # USGS-HySep R version in DVstats
library(doBy, quietly = TRUE) # need doBy package to sums for annual, summer and winter

chr.dir.r.files <- "m:/models/bacteria/hspf/bigelkhydrocal201601/r-files"
source(file = paste0(chr.dir.r.files, "/functions.R"))

chr.dir.hspf <- "m:/models/bacteria/hspf/bigelkhydrocal201601/hspf-files"
chr.file.hspf.out <- "behydcal.out"

## path to storm dates file
chr.dir.stm.dates <- "m:/models/bacteria/hspf/HydroCal201506/R_projs/Select_Storm_HydCal"

## drainage area in sqr mi for Big Elk Creek at outlet
da.be <- 88.8

df.mod <- rpltgen(chr.dir = chr.dir.hspf, chr.file = chr.file.hspf.out)

## calculate model groups
## mlog log base 10 of stream flow
mlog <- log10(df.mod$Rch18.flow)
mlog[mlog == -Inf] <- -6

## mflow stream flow
mflow <- df.mod$Rch18.flow

## mbaseind, base flow index 
## baseflow seperation using USGS-HySep R version in DVstats
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

## mvol_ann - annual volumes in ac-ft
## create factor for year
df.mod <- cbind(df.mod, 
                fac.ann  = as.factor(
                  strftime(df.mod$tmp.date, format = "%Y")))
mvol_ann <- as.numeric(
  summaryBy(flow.ac.ft ~ fac.ann, data = df.mod, FUN = sum)[ ,2])

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
mvol_smr <- as.numeric(df.vol.seasons[as.character(df.vol.seasons$fac.season) == "summer", 
                           "flow.ac.ft.sum"])

## mvol_wtr
mvol_wtr <- as.numeric(df.vol.seasons[as.character(df.vol.seasons$fac.season) == "winter", 
                           "flow.ac.ft.sum"])
## storm information
## get storm dates from text file Information in this file from 
## Select_Storm_HydCal repo
## column 2 is the begin date of storm and column 8 is the end date of storm
df.strm.dates.raw <- read.delim(file = paste0(chr.dir.stm.dates, "/dates_stm.dat"),
                            header = FALSE, sep = " ", 
                            stringsAsFactors = FALSE)[ , c(2, 8)]
## convert to POSIXct dates
df.strm.dates <- data.frame(apply(df.strm.dates.raw, MARGIN = 2, strptime, 
                                  format = "%m/%d/%Y"))
## set names
names(df.strm.dates) <- c("begin", "end")

## clean up
rm(df.strm.dates.raw)

## storm durations in days
df.strm.dur <- as.numeric(df.strm.dates$end - df.strm.dates$begin)

## mpeak
mpeak <- rep(-1, length(df.strm.dates$begin))

for(ii in 1:length(mpeak)) {
  mpeak[ii] <- max(df.mod[df.mod$tmp.date >= df.strm.dates$begin[ii] & 
               df.mod$tmp.date <= df.strm.dates$end[ii], ]$Rch18.flow)
}
rm(ii)

## mvol_stm in cu-ft for storm convert cu-ft/sec to cu-ft/day 
## using 1 day = 86400 s
mvol_stm <- rep(-1, length(df.strm.dates$begin))

for(ii in 1:length(mvol_stm)) {
  mvol_stm[ii] <- sum(df.mod[df.mod$tmp.date >= df.strm.dates$begin[ii] & 
                            df.mod$tmp.date <= 
                              df.strm.dates$end[ii], ]$flow) * 
    (df.strm.dur[ii] * 86400)
}
## mtime - % exceedance for flow, using 1%, 5%, 10%, 25%, 50%, 75%, 95%, 99%
## this is different than what Cadmus using in tsproc which is the fraction
## of time the flow is above some value. I am not going to use tsproc when
## doinmg the calculations. I will use R script

## percents used
tmp.per <- c(0.0001, 0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99)

mtime <- as.numeric(quantile(x = df.mod$Rch18.flow, probs = tmp.per))

## clean up
rm(tmp.per)


## get observational groups names

# get rows where the block names are
chr.dir.pst <- "m:/models/bacteria/hspf/bigelkhydrocal201601/pest-files"
str.control <- scan(paste0(chr.dir.pst,"/control.pst"), sep = "\n", 
                    what = "character", quiet = TRUE)
tmp.blk.hd <- grep("\\*", str.control)
str.obs.grp.names <- 
  str.control[(tmp.blk.hd[grep("[Oo]bs.*[Gg]roups", 
                               str.control[tmp.blk.hd])] + 1):
                (tmp.blk.hd[grep("[Oo]bs.*[Gg]roups", 
                                 str.control[tmp.blk.hd]) + 1] - 1)]

lng.0.pd <- nchar(as.character(max(sapply(mget(str.obs.grp.names), length))))

lng.0.pd <- nchar(as.character(sapply(mget(str.obs.grp.names), length)))

to.df.cur.data <- function(x) {
  n.pad <- nchar(as.character(length(get(x))))
  y <- data.frame(
    name = paste0(x, "_", 
                  sprintf(fmt = paste0("%0", n.pad, "d"), 
                          1:length(get(x)))), 
    val = get(x), stringsAsFactors = FALSE)
  return(y)
} 

tmp.blk.data <- do.call(rbind, lapply(str.obs.grp.names, FUN = to.df.cur.data))


##
## write output to filed format text file

## get length of longest variable name and add 5
#lng.name <- max(nchar(attr(tmp.blk.data, "names"))) + 5
lng.name <- max(nchar(tmp.blk.data$name)) + 5



## write output to chracater vector. the format is 2s, variable name left 
## justified and width 5 plus length of longest variable name and value as 
## 1.5E+00. 
## total width of variable name and value is length of longest variable name + 5
## + 11 = length of longest variable name + 16
chr.mod.output <- paste0(
  sprintf(paste0("  %-", lng.name, "s"), tmp.blk.data$name), 
  sprintf("%.5E", tmp.blk.data$val))

write.table(data.frame(out=chr.mod.output), 
            file = paste0(chr.dir.pst, "/model.out"), 
            quote = FALSE, col.names = FALSE, row.names = FALSE,
            sep = "\n")


##
## write model.ins file
lng.var.val <- nchar(chr.mod.output[1])

chr.mod.ins <- c("pif $", paste0(
  sprintf("l1  [%s]",tmp.blk.data$name), lng.name, ":",lng.var.val + 1))

write.table(data.frame(ins=chr.mod.ins), 
            file = paste0(chr.dir.pst, "/model.ins"), 
            quote = FALSE, col.names = FALSE, row.names = FALSE,
            sep = "\n")


