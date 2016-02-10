## analyze pest-hspf results

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


