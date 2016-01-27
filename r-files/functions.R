
rpltgen <- function(chr.dir = "m:/models/bacteria/hspf/",
                        chr.file = "beflhyd.out") {
  ## funtion that reads PLTGEN output file from HSPF simulations
  ## created 2016-01-27 by Kevin Brannan
  ## based on function read.pltgen in the GitHub repo
  ## https://github.com/kbrannan/construct_control
  ## input:
  ## chr.dir is the path to the pltgen file
  ## chr.file is the name of the pltgen file to read
  
  ## set the working directory
  ##setwd("m:/models/bacteria/hspf/HydroCal201506/R_Projs/construct_control")
  
  
  ## set path for the PLTGEN file
  ##chr.dir <- paste0(getwd(), "/ODEQ_hspf/extended period/")
  
  ## set name for the PLTGEN file
  ## chr.file <- "beflhyd.out"
  
  ## read in the PLTGEN file
  chr.pltgen <- scan(file = paste0(chr.dir,"/", chr.file), sep = "\n", 
                     what = "character")
  
  ## get first line of data. the "-1.0000000E+30" is a flag for no data and 
  ## should only occur on the first day for a daily tiime step aggregation of
  ## hourly data. take min just in case
  lng.str <- min(min(grep("-1\\.0000000E\\+30{1,}", chr.pltgen)) + 1)
  
  
  
  ## get data only from orginal PLTGEN file
  str.data <- chr.pltgen[lng.str:length(chr.pltgen)]
  
  
  ## get variable names from PLTGEN file
  ## get line where the variable name list starts
  lng.name.str <- grep("^.*Label( ){1,}LINTYP", 
                       chr.pltgen[1:lng.str]) + 1
  
  ## get line where the variable name list ends
  lng.name.end <- min(
    grep(
      paste0(
        gsub(" Label.*$","",chr.pltgen[lng.name.str -1]),"$"), 
      chr.pltgen[lng.name.str:lng.str])) + lng.name.str - 2
  
  gsub(paste0(
    chr.pltgen[lng.name.end + 1], " "), "", chr.pltgen[lng.name.str:lng.name.end])
  
  
  ## get variable names
  str.var.names <- do.call(rbind,
                           lapply(gsub(paste0(
                             chr.pltgen[lng.name.end + 1], " "), "", 
                             chr.pltgen[lng.name.str:lng.name.end]),
                                  function(x){
                                    y <- gsub("( ){3,}.*","", 
                                              gsub("^( ){1,}To( ){1,}" ,
                                                   "", x))
                                    z <- gsub(" ", ".", y)
                                    return(z)
                                  }
                           )
  )
  
  ## get data in dtat frame
  ## fisrt get data as characters
  df.data.chr <- data.frame(do.call(rbind,
                                    strsplit(x = gsub("( ){1,}24( ){1,}0", " ",
                                                      gsub("^( ){1,}To( ){1,}", "",
                                                           str.data))
                                             , split = "( ){1,}")),
                            stringsAsFactors = FALSE)
  
  ## rename variables in the data frame to names in PLTGEN file
  names(df.data.chr) <- c("year", "month", "day", str.var.names)
  
  ## create data frame with date and numeric values
  tmp.flows <- df.data.chr[, -1 * c(1, 2, 3)]
  for(ii in 1:length(tmp.flows[1, ])) {
    tmp.flows[ , ii] <- as.numeric(tmp.flows[ , ii])
  }
  tmp.date <- as.POSIXct(apply(df.data.chr[, 1:3], 
                               1, function(x) paste0(x, collapse = "-")))
  
  ## create data.frame and return
  df.data <- cbind(tmp.date, tmp.flows)
  return(df.data)
  
  
}