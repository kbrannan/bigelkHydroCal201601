chr.dir.r.files <- "m:/models/bacteria/hspf/bigelkhydrocal201601/r-files"
source(file = paste0(chr.dir.r.files, "/functions.R"))

chr.dir.hspf <- "m:/models/bacteria/hspf/bigelkhydrocal201601/hspf-files"
chr.file.hspf.out <- "behydcal.out"

df.mod <- rpltgen(chr.dir = chr.dir.hspf, chr.file = chr.file.hspf.out)

## calculate model groups
## mlog log base 10 of stream flow
mlog <- log10(df.mod$Rch18.flow)

## mflow stream flow
mflow <- df.mod$Rch18.flow

## mbaseind


## mpeak
## mvol_ann
## mvol_smr
## mvol_wtr
## mvol_stm
## mtime

