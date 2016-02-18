## analyze pest-hspf results

## load packages
library(ggplot2)
library(reshape)
require(DVstats)

## working path 
chr.dir <- "M:/Models/Bacteria/HSPF/bigelkHydroCal201601"
## path to storm dates file
chr.dir.stm <- "M:/Models/Bacteria/HSPF/HydroCal201506/R_projs/Select_Storm_HydCal"

## load functions
source("m:/Models/Bacteria/LDC/Calculations/Rscripts/LDC Functions.R")
source(file=paste0(chr.dir.stm,"/devel/functions.R"))


## read residuals file
chr.res <- scan(file = paste0(chr.dir,"/pest-files/control.res"),
                what = "character", sep = "\n", quiet = TRUE)
## replace * with x in field names
chr.res[1] <- gsub("\\*","x",chr.res[1])

## replace spaces among columns with comma
chr.res <- gsub("( ){1,}",",",chr.res)

## convert chracter vector to data.frame
df.res <- data.frame(do.call(rbind,strsplit(chr.res,split = ",")), 
                     stringsAsFactors = FALSE)

## remove first column becuase it is empty
df.res <- df.res[ ,-1]

## first row is names for columns
names(df.res) <- df.res[1, ]

## discard first row
df.res <- df.res[-1, ]

## get simulation dates from UCI
chr.sim.dates <- gsub("([aA-zZ ])|(00\\:00)|(24\\:00)","",
     grep("START", scan(file = paste0(chr.dir, "/hspf-files/bigelk.uci"), 
                   sep = "\n", what = "character", quiet = TRUE), value = TRUE))
dte.str <- as.Date(substr(chr.sim.dates, start =  1, stop = 10), fmt = "%Y/%m/%d")
dte.end <- as.Date(substr(chr.sim.dates, start = 11, stop = 20), fmt = "%Y/%m/%d")

## create date sequence
dte.flows <- seq(from = as.Date(dte.str), to = as.Date(dte.end), by = "day")

## get mlog time-series
df.mlog <- data.frame(dates = dte.flows, 
                      df.res[grep("mlog", as.character(df.res$Group)), ])
df.mlog[, 4:12] <- sapply(df.mlog[ , 4:12], as.numeric)
## add year as factor
df.mlog <- cbind(df.mlog, year = factor(format(df.mlog$dates, "%Y")))
### add month
df.mlog <- cbind(df.mlog, 
                 month = 
                   factor(
                     format(df.mlog$dates, "%b"),
                     levels = c("Oct", "Nov", "Dec", "Jan","Feb",
                                "Mar", "Apr", "May","Jun", "Jul", 
                                "Aug", "Sep")))
## add season
chr_season <- function(chr.month) {
  if(chr.month %in% c("Dec", "Jan", "Feb")) season  <- "winter"
  if(chr.month %in% c("Mar", "Apr", "May")) season   <- "spring"
  if(chr.month %in% c("Jul", "Jun", "Aug")) season   <- "summer"
  if(chr.month %in% c("Sep", "Oct", "Nov")) season <- "fall"
  return(season)
}
df.mlog <- cbind(df.mlog, 
                 season = sapply(as.character(df.mlog$month), chr_season))

## assign flow zones used in LCD
## add ecdf to mlog
mlog.ecdf <- ecdf(df.mlog$Modelled)
df.mlog <- cbind(df.mlog, exceed = 100 * ( 1 - round(mlog.ecdf(df.mlog$Modelled),3)) )
## function taken from the LDC functions
get.flow.zone <- function(flow.exceed) {
  flow.zone <- NA
  if(is.na(flow.exceed) == TRUE) {return(flow.zone)}
  if(flow.exceed >= 0 & flow.exceed <= 10) {
    flow.zone <- "high"
  } else if(flow.exceed > 10 & flow.exceed <= 40) {
    flow.zone <- "transitional"
  } else if(flow.exceed > 40 & flow.exceed <= 60) {
    flow.zone <- "typical"
  } else if(flow.exceed > 60 & flow.exceed <= 90) {
    flow.zone <- "dry"
  } else if(flow.exceed > 90 & flow.exceed <= 100) {
    flow.zone <- "low"
  }
  return(flow.zone)
}
df.mlog <- cbind(df.mlog, 
                 flw.zn = factor(
                   sapply(df.mlog$exceed,get.flow.zone), 
                   levels = c("high", "transitional", "typical", 
                              "dry", "low")))
## boxplot of weight x residuals
p.mlog.bar.wt.rs.all <- ggplot(data = df.mlog, 
                               aes(x=factor(0), y = WeightxResidual)) + 
  geom_boxplot()
plot(p.mlog.bar.wt.rs.all)

## boxplot of weight x residuals by year
p.mlog.bar.wt.rs.yr <- ggplot(data = df.mlog, 
                               aes(x=year, y = WeightxResidual)) + 
  geom_boxplot()
plot(p.mlog.bar.wt.rs.yr)

## boxplot of weight x residuals by season
p.mlog.bar.wt.rs.sn <- ggplot(data = df.mlog, 
                              aes(x=season, y = WeightxResidual)) + 
  geom_boxplot()
plot(p.mlog.bar.wt.rs.sn)

## boxplot of weight x residuals by ldc flow zone
p.mlog.bar.wt.rs.fz <- ggplot(data = df.mlog, 
                              aes(x=flw.zn, y = WeightxResidual)) + 
  geom_boxplot()
plot(p.mlog.bar.wt.rs.fz)


## get mflow time-series
df.mflow <- data.frame(dates = dte.flows, 
                      df.res[grep("mflow", as.character(df.res$Group)), ])
df.mflow[, 4:12] <- sapply(df.mflow[ , 4:12], as.numeric)
## add year as factor
df.mflow <- cbind(df.mflow, year = factor(format(df.mflow$dates, "%Y")))
### add month
df.mflow <- cbind(df.mflow, 
                 month = 
                   factor(
                     format(df.mflow$dates, "%b"),
                     levels = c("Oct", "Nov", "Dec", "Jan","Feb",
                                "Mar", "Apr", "May","Jun", "Jul", 
                                "Aug", "Sep")))
## add season
df.mflow <- cbind(df.mflow, 
                 season = sapply(df.mflow$month, chr_season))

## assign flow zones used in LCD
## add ecdf to mlog
mflow.ecdf <- ecdf(df.mflow$Modelled)
df.mflow <- cbind(df.mflow, exceed = 100 * ( 1 - round(mflow.ecdf(df.mflow$Modelled),3)) )
df.mflow <- cbind(df.mflow, 
                 flw.zn = factor(
                   sapply(df.mflow$exceed,get.flow.zone), 
                   levels = c("high", "transitional", "typical", 
                              "dry", "low")))
## boxplot of weight x residuals
p.mflow.bar.wt.rs.all <- ggplot(data = df.mflow, 
                               aes(x=factor(0), y = WeightxResidual)) + 
  geom_boxplot()
plot(p.mflow.bar.wt.rs.all)

## boxplot of weight x residuals by year
p.mflow.bar.wt.rs.yr <- ggplot(data = df.mflow, 
                              aes(x=year, y = WeightxResidual)) + 
  geom_boxplot()
plot(p.mflow.bar.wt.rs.yr)

## boxplot of weight x residuals by season
p.mflow.bar.wt.rs.sn <- ggplot(data = df.mflow, 
                              aes(x=season, y = WeightxResidual)) + 
  geom_boxplot()
plot(p.mflow.bar.wt.rs.sn)

## boxplot of weight x residuals by ldc flow zone
p.mflow.bar.wt.rs.fz <- ggplot(data = df.mflow, 
                              aes(x=flw.zn, y = WeightxResidual)) + 
  geom_boxplot()
plot(p.mflow.bar.wt.rs.fz)

## get mtime
df.mtime <- data.frame(
  probs = eval(
    parse(text = 
            gsub("(tmp.per)|(<-)|(( ){1, })", "", 
                 grep("tmp.per <-", 
                      scan(
                        file = paste0(chr.dir, "/r-files/hspf-output-proc.R"), 
                        sep = "\n", 
                        what = "character", 
                        quiet = TRUE), 
                      value = TRUE)))), 
  df.res[grep("mtime", as.character(df.res$Group)), ])
df.mtime[ , 4:12] <- sapply(df.mtime[, 4:12], as.numeric)
## create a long format table for mtime with factor indicating if from obs or model
df.mtime.lg <- data.frame(rbind(cbind(src = "obs", value = df.mtime[ , 5], df.mtime[ ,c(1,6:12)]),
                          cbind(src = "model", value = df.mtime[ , 4], df.mtime[ ,c(1,6:12)])),
                          stringsAsFactors = FALSE)
df.mtime.lg$src <- as.factor(df.mtime.lg$src)
levels(df.mtime.lg$src) <- c("Obs", "Model")
df.mtime.lg[ , -1] <- sapply(df.mtime.lg[ , -1], as.numeric)

## estimate fdc
## stations at outlet is LASAR 34453
tmp.one.station <- 34453
tmp.ss.est.fn <- paste0("st",tmp.one.station,".xml")
tmp.fdc.ss.est <- fdc.ss.estimate(ss.fn=tmp.ss.est.fn, ss.path=get.path("StreamStatsBacteria"))
tmp.fdc.ss.name <- paste0("fdc.ss.st",tmp.one.station)
eval(parse(text=paste0(tmp.fdc.ss.name," <- tmp.fdc.ss.est")))
rm(list=ls(pattern="^tmp\\.*")) ## clean up

## create three difference data.frames for plots
df.mtime.obs <- df.mtime.lg[df.mtime.lg$src == "Obs", c("probs", "value")]
df.mtime.mod <- df.mtime.lg[df.mtime.lg$src == "Model", c("probs", "value", "Residual","WeightxResidual")]
df.mtime.eq <- fdc.ss.st34453

df.mtime.obs$probs <- 100 * (1 - df.mtime.obs$probs)
df.mtime.mod$probs <- 100 * (1 - df.mtime.mod$probs)
df.mtime.eq$FDPercent <- 100 * df.mtime.eq$FDPercent

names(df.mtime.obs) <- c("x", "y")
names(df.mtime.mod) <- c("x", "y", "res", "wxres")
names(df.mtime.eq)  <- c("x", "y", "ymin", "ymax")

df.mtime.obs <- df.mtime.obs[order(df.mtime.obs$x), ]
df.mtime.mod <- df.mtime.mod[order(df.mtime.mod$x), ]

df.mtime.all <- melt(list(obs = df.mtime.obs, 
                          mod = df.mtime.mod, 
                          eq = df.mtime.eq),
                     id.vars = "x")

df.mtime.all$L1 <- factor(df.mtime.all$L1, levels = c("obs", "mod", "eq"))


## plot mtime and related data
p.mtime00 <- ggplot(data = df.mtime.all) + 
  scale_y_log10() + 
  scale_colour_manual(name = "", breaks = c("obs", "mod", "eq"), 
                        labels = c("Obs", "Model", "USGS Eq"), 
                      values = c("blue", "green", "black")) +
  scale_shape_manual(name = "", breaks = c("obs", "mod", "eq"), 
                     labels = c("Obs", "Model", "USGS Eq"),
                     values = c(21, 17, 15)) +
    xlab("Percent Time Greater") + ylab("Mean Daily Flow (cfs)")

## add obs and mod flow data along with reg eq
p.mtime00 <- p.mtime00 + 
  geom_line(data = df.mtime.all[as.character(df.mtime.all$variable) == "y", ],
            aes(x = x, y = value, colour = L1))

p.mtime00 <- p.mtime00 + 
  geom_point(data = df.mtime.all[as.character(df.mtime.all$variable) == "y", ],
             aes(x = x, y = value, colour = L1, shape = L1), size = 4)

## add error bars from reg eq
p.mtime00 <- p.mtime00 + 
  geom_errorbar(data = df.mtime.all[df.mtime.all$L1 == "eq", ], 
                aes(x = df.mtime.all[df.mtime.all$L1 == "eq" & 
                                       as.character(df.mtime.all$variable) == "ymin", 'x'],
    ymin = df.mtime.all[df.mtime.all$L1 == "eq" & 
                          as.character(df.mtime.all$variable) == "ymin", 'value'],
    ymax = df.mtime.all[df.mtime.all$L1 == "eq" & 
                          as.character(df.mtime.all$variable) == "ymax", "value"]
    ))
plot(p.mtime00)

## get mvol_ann
chr.yrs <- unique(format(dte.flows, "%Y"))
df.mvol_ann <- data.frame(year = factor(chr.yrs), 
                       df.res[grep("mvol_ann", as.character(df.res$Group)), ])
df.mvol_ann[, 4:12] <- sapply(df.mvol_ann[ , 4:12], as.numeric)

## boxplot of weight x residuals
p.mvol_ann.bar.wt.rs.all <- ggplot(data = df.mvol_ann, 
                                aes(x=factor(0), y = WeightxResidual)) + 
  geom_boxplot()
plot(p.mvol_ann.bar.wt.rs.all)

## scatter plot of weight x residuals by year
p.mvol_ann.pnt.wt.rs.yr <- ggplot(data = df.mvol_ann, 
                               aes(x=as.numeric(as.character(df.mvol_ann$year)), 
                                   y = WeightxResidual)) +
  geom_point(shape = 1, size = 4)
plot(p.mvol_ann.pnt.wt.rs.yr)

## bar plot of weight x residuals by year
p.mvol_ann.bar.wt.rs.yr <- ggplot(data = df.mvol_ann, 
                                  aes(x=as.numeric(as.character(df.mvol_ann$year)), 
                                      y = WeightxResidual)) + xlab("year") +
  geom_bar(stat = "identity", fill = "blue", position=position_dodge())
plot(p.mvol_ann.bar.wt.rs.yr)

## get mvol_smr
df.mvol_smr <- data.frame(
  year = factor(
    unique(
      format(df.mlog$dates[
        grep("summer", as.character(df.mlog$season))], "%Y"))), 
  df.res[grep("mvol_smr", as.character(df.res$Group)), ])
df.mvol_smr[, 4:12] <- sapply(df.mvol_smr[ , 4:12], as.numeric)

## boxplot of weight x residuals
p.mvol_smr.bar.wt.rs.all <- ggplot(data = df.mvol_smr, 
                                   aes(x=factor(0), y = Residual)) + 
  geom_boxplot()
plot(p.mvol_smr.bar.wt.rs.all)

## scatter plot of weight x residuals by year
p.mvol_smr.pnt.wt.rs.yr <- 
  ggplot(data = df.mvol_smr,
         aes(
           x=as.numeric(as.character(df.mvol_smr$year)),
           y = Residual)) + 
  xlab("year") + ylab("residual (ac-ft)") + geom_point(shape = 1, size = 4)
plot(p.mvol_smr.pnt.wt.rs.yr)

## bar plot of weight x residuals by year
p.mvol_smr.bar.wt.rs.yr <- ggplot(data = df.mvol_smr,
                                  aes(x=as.numeric(
                                    as.character(df.mvol_smr$year)),
                                    y = Residual)) + 
  xlab("year") + ylab("residual (ac-ft)") +
  geom_bar(stat = "identity", fill = "blue", position=position_dodge())
plot(p.mvol_smr.bar.wt.rs.yr)

## get mvol_wtr
df.mvol_wtr <- data.frame(
  year = factor(
    unique(
      format(df.mlog$dates[
        grep("winter", as.character(df.mlog$season))], "%Y"))), 
  df.res[grep("mvol_wtr", as.character(df.res$Group)), ])
df.mvol_wtr[, 4:12] <- sapply(df.mvol_wtr[ , 4:12], as.numeric)

## boxplot of weight x residuals
p.mvol_wtr.bar.wt.rs.all <- ggplot(data = df.mvol_wtr, 
                                   aes(x=factor(0), y = Residual)) + 
  xlab("year") + ylab("residual (ac-ft)") +
  geom_boxplot()
plot(p.mvol_wtr.bar.wt.rs.all)

## scatter plot of weight x residuals by year
p.mvol_wtr.pnt.wt.rs.yr <- 
  ggplot(data = df.mvol_wtr,
         aes(
           x=as.numeric(as.character(df.mvol_wtr$year)),
           y = Residual)) + 
  xlab("year") + ylab("residual (ac-ft)") + geom_point(shape = 1, size = 4)
plot(p.mvol_wtr.pnt.wt.rs.yr)

## bar plot of weight x residuals by year
p.mvol_wtr.bar.wt.rs.yr <- ggplot(data = df.mvol_wtr,
                                  aes(x=as.numeric(
                                    as.character(df.mvol_wtr$year)),
                                    y = Residual)) + 
  xlab("year") + ylab("residual (ac-ft)") +
  geom_bar(stat = "identity", fill = "blue", position=position_dodge())
plot(p.mvol_wtr.bar.wt.rs.yr)

## get storms
## storm information
## get storm dates from text file Information in this file from 
## Select_Storm_HydCal repo
## column 2 is the begin date of storm and column 8 is the end date of storm
df.strm.dates.raw <- read.delim(file = paste0(chr.dir.stm, "/dates_stm.dat"),
                                header = FALSE, sep = " ", 
                                stringsAsFactors = FALSE)[ , c(2, 8)]
## convert to POSIXct dates
df.strm.dates <- data.frame(apply(df.strm.dates.raw, MARGIN = 2, strptime, 
                                  format = "%m/%d/%Y"))
## set names
names(df.strm.dates) <- c("begin", "end")

## peaks
df.storms.peak <- df.res[ df.res$Group == "mpeak", ]
df.storms.peak[ ,3:11] <- sapply(df.storms.peak[ ,3:11],as.numeric)

## combine dates with storms
df.storms.peak <- cbind(df.strm.dates, df.storms.peak)
## add year from year storm starts
df.storms.peak <- cbind(df.storms.peak, 
                        year = factor(format(df.storms.peak$begin, "%Y")))
## add month from month storm starts
df.storms.peak <- cbind(df.storms.peak, 
                   month = 
                     factor(
                       format(df.storms.peak$begin, "%b"),
                       levels = c("Oct", "Nov", "Dec", "Jan",
                                  "Feb", "Mar","Apr", "May","Jun", 
                                  "Jul", "Aug", "Sep")))
## add season from month storm starts
df.storms.peak <- cbind(df.storms.peak, 
                   season = factor(
                     sapply(as.character(df.storms.peak$month), chr_season)))
## add flow regime for peak flow
mflow.ecdf.mod <- ecdf(df.mflow$Modelled)
mflow.ecdf.obs <- ecdf(df.mflow$Measured)
df.storms.peak <- cbind(df.storms.peak, 
                        mod.exceed = 
                          100 * ( 1 - round(
                            mflow.ecdf.mod(
                              df.storms.peak$Modelled),3)),
                        obs.exceed = 
                          100 * ( 1 - round(
                            mflow.ecdf.obs(
                              df.storms.peak$Measured),3))
                        )
df.storms.peak <- cbind(df.storms.peak, 
                        mod.flw.zn = factor(
                          sapply(df.storms.peak$mod.exceed,get.flow.zone),
                          levels = c("dry", "low", "typical", 
                                     "transitional", "high")),
                        obs.flw.zn = factor(
                          sapply(df.storms.peak$obs.exceed,get.flow.zone),
                          levels = c("dry", "low", "typical", 
                                     "transitional", "high"))
                        
                        )
## volumes
df.storms.vol <- df.res[ df.res$Group == "mvol_stm", ]
df.storms.vol[ ,3:11] <- sapply(df.storms.vol[ ,3:11],as.numeric)
## combine dates with storms
df.storms.vol <- cbind(df.strm.dates, df.storms.vol)
## add info from peaks data.frame
df.storms.vol <- cbind(df.storms.vol, df.storms.peak[ , 14:20])

## get precip, obs and model flow for each storms and make plots as in 
## the Storm_Select script

# get precip data
source(file=paste0(chr.dir.stm,"/devel/get-precip-data.R"),
       local = TRUE)
# use max precip between the two gages (src) for daily precip
tmp.daily.precip.max.stations <- 
  summaryBy(prec.sum ~ date_org, df.daily.precip, FUN = max)

# get precip for dates in flow ts
df.daily.precip <- tmp.daily.precip.max.stations[do.call(rbind,lapply(strftime(df.mflow$dates, format = "%Y%m%d"), grep, strftime(tmp.daily.precip.max.stations$date_org, format = "%Y%m%d"))), ]
names(df.daily.precip) <- c("dates", "daily.precip")
# ## last day of precip not included in simulation
# df.daily.precip.max.stations <- 
#   df.daily.precip.max.stations[-1 * 
#                                  length(df.daily.precip.max.stations$date), ]
# names(df.daily.precip.max.stations) <- c("date_org", "daily.precip")
# ## use mflow dates and drop date_org
# df.daily.precip <- data.frame(dates = df.mflow$dates, 
#                          daily.precip = 
#                            df.daily.precip.max.stations$prec.sum.max)
# baseflow seperation using USGS-HySep R version
df.hysep88.8.obs <- hysep(Flow = df.mflow$Measured, 
                      Dates = as.Date(df.mflow$dates), da = 88.8)
# baseflow seperation using USGS-HySep R version
df.hysep88.8.mod <- hysep(Flow = df.mflow$Modelled, 
                          Dates = as.Date(df.mflow$dates), da = 88.8)


