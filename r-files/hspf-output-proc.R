chr.dir.r.files <- "m:/models/bacteria/hspf/bigelkhydrocal201601/r-files"
source(file = paste0(chr.dir.r.files, "/functions.R"))

chr.dir.hspf <- "m:/models/bacteria/hspf/bigelkhydrocal201601/hspf-files"
chr.file.hspf.out <- "behydcal.out"

junk <- rpltgen(chr.dir = chr.dir.hspf, chr.file = chr.file.hspf.out)
