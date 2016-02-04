## get the calibrated parameter values from the Cadmus work and isert them into
## the current PEST control.file

# zipfile path and name
tmp.zip.dir <- "//deqhq1/tmdl/TMDL_WR/MidCoast/Models/Bacteria/HSPF/Hydro Calibration/Files from Cadmus"
tmp.zip <- "Final_Deliverables_EPA_July2012.zip"

# get list of files in the zipfile
tmp.fns <- unzip(zipfile = paste0(tmp.zip.dir,"/", tmp.zip), list = TRUE)

# see how many occurances there are of the control file in the zipfile
tmp.par <- grep("*\\.par", tmp.fns$Name, value = TRUE)
tmp.only.par <- tmp.par[-1 * grep("parset", tmp.par)]

# Use this one
tmp.file <- "Final_Deliverables_EPA_July2012/PEST_end/calib.par"

# read the PEST par file as a character vector
str.par <- scan(unz(paste0(tmp.dir,"/", tmp.zip),
                    tmp.file)
                , sep = "\n", what = character())
str.par

# get current control file
tmp.pst.dir <- "//deqhq1/tmdl/TMDL_WR/MidCoast/Models/Bacteria/HSPF/bigelkHydroCal201601/pest-files"
str.cur.pst <- scan(paste0(tmp.pst.dir,"/control.pst"), sep = "\n", what = character())

grep("\\*", str.cur.pst, value = TRUE)

## insert block of calib pars into the control
lng.cur.st <- grep("\\* parameter data" ,str.cur.pst)
lng.cur.ed <- lng.cur.st + min(grep("\\* " , 
                                    str.cur.pst[(lng.cur.st + 1):length(str.cur.pst)]))

str.cur.pst[lng.cur.st:lng.cur.ed]

str.new.pst <- c(str.cur.pst[1:lng.cur.st], 
                     str.par[2:length(str.par)],
                     str.cur.pst[lng.cur.ed:length(str.cur.pst)])

str.new.pst[lng.cur.st:lng.cur.ed]

write.table(str.new.pst, file = paste0(tmp.pst.dir,"/control-cadmus.pst"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)



