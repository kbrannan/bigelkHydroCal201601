## analyze pest-hspf results

## load packages
library(ggplot2)

## working path 
chr.dir <- "M:/Models/Bacteria/HSPF/bigelkHydroCal201601"

## read residuals file
chr.res <- scan(file = paste0(chr.dir,"/pest-files/control.res"),
                what = "character", sep = "\n")
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
                   sep = "\n", what = "character"), value = TRUE))
dte.str <- as.POSIXct(substr(chr.sim.dates, start =  1, stop = 10), fmt = "%Y/%m/%d")
dte.end <- as.POSIXct(substr(chr.sim.dates, start = 11, stop = 20), fmt = "%Y/%m/%d")

## create date sequence
dte.flows <- seq(from = dte.str, to = dte.end, by = "day")

## get mlog time-series
df.mlog <- data.frame(dates = dte.flows, 
                      df.res[grep("mlog", as.character(df.res$Group)), ])
df.mlog[, 4:12] <- sapply(df.mlog[ , 4:12], as.numeric)
## add year as factor
df.mlog <- cbind(df.mlog, year = factor(format(df.mlog$dates, "%Y")))
### add month
df.mlog <- cbind(df.mlog, 
                 month = as.factor(format(df.mlog$dates, "%b")))

levels(df.mlog$month) <- c("Oct", "Nov", "Dec", "Jan","Feb", 
                           "Mar", "Apr", "May","Jun", "Jul", "Aug", "Sep")


## add season
chr_season <- function(chr.month) {
  if(chr.month %in% c("Dec", "Jan", "Feb")) season  <- "winter"
  if(chr.month %in% c("Mar", "Apr", "May")) season   <- "spring"
  if(chr.month %in% c("Jul", "Jun", "Aug")) season   <- "summer"
  if(chr.month %in% c("Sep", "Oct", "Nov")) season <- "fall"
  return(season)
}
df.mlog <- cbind(df.mlog, 
                 season = sapply(df.mlog$month, chr_season))

## bar plot of weight residuals
p.mlog.bar.wt.rs.all <- ggplot(data = df.mlog, 
                               aes(x=factor(0), y = WeightxResidual)) + 
  geom_boxplot()
plot(p.mlog.bar.wt.rs.all)

summary(df.mlog$WeightxResidual)


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
## plot mtime ADD information from USGS regression equations
p.mtime00 <- ggplot(data = df.mtime.lg, 
                    aes(x = 100 * (1 - probs), colour = src)) + 
  scale_y_log10() + 
  scale_colour_discrete(name = "", breaks = c("obs", "model"), 
                        labels = c("Obs", "Model")) +
  xlab("Percent Time Greater") + ylab("Mean Daily Flow (cfs)") 
## add flow data
p.mtime00 <- p.mtime00 + geom_line(aes(y = value))
## add weighted residuals
p.mtime00 <- p.mtime00 + geom_point(data = df.mtime.lg[df.mtime.lg$src == "obs", ], 
                                 aes(x = 100 * (1 - probs), y = value, 
                                     size = WeightxResidual)) + 
  scale_size_continuous(name = "Weighted Resudual", range = c(4,8))
plot(p.mtime00)

## 

