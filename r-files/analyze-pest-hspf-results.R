## analyze pest-hspf results

## working path 
chr.dir <- "M:/Models/Bacteria/HSPF/bigelkHydroCal201601"

## read residuals file
chr.res <- scan(file = paste0(chr.dir,"/pest-files/control.res"),
                what = "character", sep = "\n")
head(chr.res)

chr.res[1] <- gsub("\\*","x",chr.res[1])

chr.res <- gsub("( ){1,}",",",chr.res)

df.res <- data.frame(do.call(rbind,strsplit(chr.res,split = ",")))[ , -1]

names(df.res) <- df.res[1, ]

df.res <- df.res[ ,-1]
