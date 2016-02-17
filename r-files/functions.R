
rpltgen <- function(chr.dir = "m:/models/bacteria/hspf/",
                        chr.file = "beflhyd.out") {
  ## funtion that reads PLTGEN output file from HSPF simulations
  ## created 2016-01-27 by Kevin Brannan
  ## based on function read.pltgen in the GitHub repo
  ## https://github.com/kbrannan/construct_control
  ## input:
  ## chr.dir is the path to the pltgen file
  ## chr.file is the name of the pltgen file to read
  
  ## read in the PLTGEN file
  chr.pltgen <- scan(file = paste0(chr.dir,"/", chr.file), sep = "\n", 
                     what = "character", quiet = TRUE)
  
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
  
  ## get data in data frame

  ## get rid of leading characters
  str.data <- gsub(paste0(chr.pltgen[lng.name.end + 1], 
                          "( ){1, }"), "", str.data)

  ## get rid of trailing spaces
  str.data <- gsub("( ){1, }$", "", str.data)

  ## extract data from character vector
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

storms_plot_to_file <- function(dte.stms = dte.stms, 
                                dte.flows  = dte.flows,
                                obs.flow   = obs.flow,
                                mod.flow   = mod.flow,
                                obs.bflow  = obs.bflow,
                                mod.bflow  = mod.bflow,
                                obs.pflow  = obs.pflow,
                                mod.pflow  = mod.pflow,
                                precip     = precip,
                                out.file  = "strmInvdPlots.pdf") {
  
  # plot flow and precip times series with storm hydrographs highlighted for 
  # each individula storm and send figure to pdf file
  # input:
  #
  # dte.stms - data.frame of POSIXct dates that correspond to first and last 
  #            days of storms
  # dte.flows - vector of POSIXct dates that correspond to flow values
  # obs.flow - vector of observed daily flow values in cfs
  # mod.flow - vector of modeled daily flow values in cfs
  # obs.bflow - vector of observed daily base flow values in cfs
  # mod.pflow - vector of modeled daily base flow values in cfs
  # obs.pflow - vector of observed peak flow values in cfs for storms
  # mod.pflow - vector of modeled peak flow values in cfs for storms
  #   precip - vector of numeric precip values in inches
  # out.file - string for pdf filename (and path) where figure is sent
  #
  # output:
  # No output in R environment. Figure is sent to pdf file at location sepcified
  # in out.file vairable
  
  # creating temporary data sets for flow and precip  
  df.f <- data.frame(dates = dates, flow = flow)
  df.p <- data.frame(date = dates, p = precip)
  
  # creating temporary data sets for storm list data set
  tmp.peaks <- pot.strm$convex
  tmp.rises <- pot.strm$concave
  tmp.rises.sel <- pot.strm$rises.sel
  tmp.pot.strms <- pot.strm$pot.strm
  strm.nums <- as.numeric(unique(as.character(tmp.pot.strms$strm.num)))
  
  # open pdf file for output
  pdf(file = out.file, width = 11, height = 8.5, onefile = TRUE)
  
  # loop to print a figure for each storm to a signle page in "out.file"
  for(ii in 1:(length(strm.nums))) {
    # get info for storm ii
    x <- tmp.pot.strms[tmp.pot.strms$strm.num == strm.nums[ii], ]
    
    # set y-limits for current storm
    tmp.ylims <- c(10 ^ (floor(log10(min(x$flow))   - 1)), 
                   10 ^ (ceiling(log10(max(x$flow)) + 1)))
    
    # set x-limits for current storm
    tmp.xlims <- c(min(x$date) - 1, max(x$date) + 1)
    
    # subset precip and flow for current storm
    tmp.p <- df.p[df.p$date >= tmp.xlims[1] & df.p$date <= tmp.xlims[2], ]
    tmp.f <- df.f[df.f$date >= tmp.xlims[1] & 
                    df.f$date <= tmp.xlims[2], ]
    
    # set plot area matrix for 2 rows and one column along with other pars
    par(mfrow = c(2, 1), tck = 0.01,  mar = c(0, 1, 0, 0.5), 
        oma = c(7, 5, 7, 2))
    
    # vary height of plots. Make precip hight smaller than flow
    layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE), heights = c(1, 3), 
           widths = c(1, 1))
    
    # precip plot set up, don't plot data
    plot(x = tmp.p$date, y = tmp.p$p, xlab = "",pch = "", xlim = tmp.xlims,
         ylim = c(0, max(c(tmp.p$p, 0.1))), xaxt = "n")
    
    # title for plot is current sorm
    title(xlab = "", ylab = "", main = paste0("storm num ", strm.nums[ii]), 
          outer = TRUE, line = 3)
    
    # plot vertical lines for each precip obs
    lines(x = tmp.p$date, y = tmp.p$p, type = "h")
    
    # add grid lines in plot for dates
    grid(nx = 30, ny = NULL)
    
    # set up plot for flow, don't plot data. Plot flow data after the storm
    # polygons are plotted so flow lines are borders of polygons  
    plot(x = tmp.f$date, y = tmp.f$flow, type = "l", log = "y", lty = "blank",
         xlim = tmp.xlims, ylim = tmp.ylims, xaxt = "n")
    
    # storm flow as a filled in polygon
    polygon(x = x$date, y = x$flow, col = "yellow", lty = "blank")
    
    # flow data for year plotted over storm polygon
    lines(x = tmp.f$date, y = tmp.f$flow, type = "l", col = "blue")
    
    # points for rises ploted over storm polygon and flow data
    points(x = tmp.rises$date, y = tmp.rises$flow)
    
    # points for peaks ploted over storm polygon and flow data
    points(x = tmp.peaks$date, y = tmp.peaks$flow)
    
    # add grid lines in plot for dates
    grid(nx = 30, ny = NULL)
    
    # format x-axis on flow plot to include day, month and year
    axis.Date(side = 1, x = tmp.f$date, format = "%m-%d-%Y")
  }
  
  # close file done
  dev.off()
}