## get the calibrated parameter values from the Cadmus work and isert them into
## the current PEST control.file

# zipfile path and name
tmp.dir <- "//deqhq1/tmdl/TMDL_WR/MidCoast/Models/Bacteria/HSPF/Hydro Calibration/Files from Cadmus"
tmp.zip <- "Final_Deliverables_EPA_July2012.zip"

# get list of files in the zipfile
tmp.fns <- unzip(zipfile = paste0(tmp.dir,"/", tmp.zip), list = TRUE)

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
