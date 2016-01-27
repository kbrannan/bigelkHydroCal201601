## path and file name for uci file
chr.dir <- "M:/Models/Bacteria/HSPF/bigelkHydroCal201601/hspf-files"
chr.uci.file <- "bigelk.uci"

## read uci file as character vector
chr.uci <- scan(file = paste0(chr.dir,"/",chr.uci.file), sep = "\n", 
                what = "character")

## get network block
lng.str <- grep("^SCHEMATIC$", chr.uci)
lng.end <- grep("^END SCHEMATIC$", chr.uci)
chr.network.blk <- chr.uci[(lng.str + 1):(lng.end - 1)]

## get the IMPLND areas
lng.implnd <- grep("IMPLND", chr.network.blk)
chr.implnd <- do.call(rbind, strsplit(chr.network.blk[lng.implnd], 
                                      split = " {1,}"))

paste0("total area of IMPLND is ", sum(as.numeric(chr.implnd[ , 3])), " ac")

