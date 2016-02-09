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

head(df.res)
