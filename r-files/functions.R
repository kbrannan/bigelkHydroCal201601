
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
                                storms.peak  = storms.peak,
                                storms.vol = storms.vol, 
                                precip     = precip,
                                out.file  = "strmInvdPlots.pdf") {
  
  # plot flow and precip information foreach individula storm and send
  # figure to pdf file
  #
  # input:
  #
  #    dte.stms - data.frame of POSIXct dates that correspond to first and last 
  #               days of storms
  #   dte.flows - vector of POSIXct dates that correspond to flow values
  #    obs.flow - vector of observed daily flow values in cfs
  #    mod.flow - vector of modeled daily flow values in cfs
  #   obs.bflow - vector of observed daily base flow values in cfs
  # storms.peak - data.frame for peak flow data for each storm
  #  storms.vol - data.frame for flow volume data for each storm
  #      precip - vector of numeric precip values in inches
  #    out.file - string for pdf filename (and path) where figure is sent
  #
  # output:
  # No output in R environment. Figure is sent to pdf file at location sepcified
  # in out.file vairable

  # temerpory assignments for function development
#   dte.stms <- df.storms.peak[ , c("begin", "end")]
#   dte.flows  <- df.mflow$dates
#   obs.flow   <- df.mflow$Measured
#   mod.flow   <- df.mflow$Modelled
#   obs.bflow  <- df.hysep88.8.obs$BaseQ
#   mod.bflow  <- df.hysep88.8.mod$BaseQ
#   storms.peak <- df.storms.peak
#   storms.vol <- df.storms.vol
#   obs.pflow  <- storms.peak$Measured
#   mod.pflow  <- storms.peak$Modelled
#   obs.vflow  <- storms.vol$Measured
#   mod.vflow  <- storms.vol$Modelled
#   precip     <- df.daily.precip$daily.precip
#   out.file <- "M:/Models/Bacteria/HSPF/bigelkHydroCal201601/indvstrms.pdf"
  
  
    
  # creating temporary data sets for storm list data set
  strm.nums <- 1:length(dte.stms$begin)
  
  # open pdf file for output
  pdf(file = out.file, width = 11, height = 8.5, onefile = TRUE)
  
  # loop to print a figure for each storm to a signle page in "out.file"
    for(ii in 1:(length(strm.nums))) {
      # get info for storm ii
      ##ii <- 25
      ## begin row for storm in flow time-series
      lng.begin <- grep(strftime(dte.stms$begin[ii], format = "%Y%m%d"), 
           strftime(dte.flows, format = "%Y%m%d"))
      lng.end   <- grep(strftime(dte.stms$end[ii], format = "%Y%m%d"), 
                        strftime(dte.flows, format = "%Y%m%d"))
      if(length(lng.begin) > 1 | length(lng.end) > 1) print(paste0("Too many matches for dates for storm ", ii))
      ## expand range for plotting
      lng.ex <- 2 # number of days before start and after end
      lng.begin.ex <- max(lng.begin - lng.ex, 1)
      lng.end.ex   <- min(lng.end + lng.ex, length(dte.flows))
      ## dates for storm time series
      tmp.dte.stm.ts <- dte.flows[lng.begin:lng.end]
      # get data for current storm
      tmp.dte.flows  <- dte.flows[lng.begin.ex:lng.end.ex]
      tmp.obs.flow   <- obs.flow[lng.begin.ex:lng.end.ex]
      tmp.mod.flow   <- mod.flow[lng.begin.ex:lng.end.ex]
      tmp.obs.bflow  <- obs.bflow[lng.begin.ex:lng.end.ex]
      tmp.mod.bflow  <- mod.bflow[lng.begin.ex:lng.end.ex]
      tmp.dte.obs.pflow <- min(tmp.dte.stm.ts[df.mflow$Measured[lng.begin:lng.end] == storms.peak$Measured[ii]])
      tmp.dte.mod.pflow <- min(tmp.dte.stm.ts[df.mflow$Modelled[lng.begin:lng.end] == storms.peak$Modelled[ii]])
      tmp.precip     <- df.daily.precip$daily.precip[lng.begin.ex:lng.end.ex]
      
  
      
      
      
      # set y-limits for current storm
      tmp.ylims <- 
        c(10 ^ (floor(log10(min(tmp.obs.bflow, tmp.mod.bflow)))), 
          10 ^ (ceiling(log10(max(tmp.obs.flow, tmp.mod.flow,
                                  tmp.obs.pflow, tmp.mod.pflow)) )))
      
      # set x-limits for current storm
      tmp.xlims <- c(dte.flows[lng.begin.ex], dte.flows[lng.end.ex])
      
      # set plot area matrix for 2 rows and one column along with other pars
      par(mfrow = c(2, 1), tck = 0.01,  mar = c(0, 1, 0, 0.5), 
          oma = c(7, 5, 7, 2))
      
      # vary height of plots. Make precip hight smaller than flow
      layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE), heights = c(1, 3), 
             widths = c(1, 1))
      
      # precip plot set up, don't plot data
      plot(x = tmp.dte.flows, y = tmp.precip, xlab = "",pch = "", xlim = tmp.xlims,
           ylim = c(0, max(c(tmp.precip, 0.1))), xaxt = "n")
      
      # title for plot is current storm
      title(xlab = "", ylab = "", 
            main = paste0("storm num ", strm.nums[ii],
            "  Season = ", storms.peak$season[ii],
            "\nObs Flow-zone = ", storms.peak$obs.flw.zn[ii],
            "  Mod Flow-zone = ", storms.peak$mod.flw.zn[ii],
            "\nObs peak flow = ", sprintf(fmt = "%-.2f cfs", storms.peak$Measured[ii]),
            " Mod peak flow = ", sprintf(fmt = "%-.2f cfs", storms.peak$Modelled[ii]),
            "  Obs storm vol = ", sprintf(fmt = "%-.2E ac-ft", storms.vol$Measured[ii]),
            " Mod storm vol = ", sprintf(fmt = "%-.2E ac-ft", storms.vol$Modelled[ii])),
            outer = TRUE, line = 3)
      
      # plot vertical lines for each precip obs
      lines(x = tmp.dte.flows, y = tmp.precip, type = "h")
      
      # add grid lines in plot for dates
      grid(nx = 30, ny = NULL)
      
      # set up plot for flow, don't plot data. Plot flow data after the storm
      # polygons are plotted so flow lines are borders of polygons  
      plot(x = tmp.dte.flows, y = tmp.obs.flow, type = "l", log = "y", lty = "blank",
           xlim = tmp.xlims, ylim = tmp.ylims, xaxt = "n")
      
      # storm flow as a filled in polygon
      #polygon(x = tmp.dte.flows, y = tmp.obs.flow, col = "yellow", lty = "blank")
      
      # flow data storm polygon
      lines(x = tmp.dte.flows, y = tmp.obs.flow, type = "l", col = "blue")
      lines(x = tmp.dte.flows, y = tmp.mod.flow, type = "l", col = "red")
      
      # base flow data storm polygon
      lines(x = tmp.dte.flows, y = tmp.obs.bflow, type = "l", lty = "dashed", col = "blue")
      lines(x = tmp.dte.flows, y = tmp.mod.bflow, type = "l", lty = "dashed", col = "red")
      
  
      # points for peaks ploted over storm polygon and flow data
        points(x = tmp.dte.obs.pflow, y = storms.peak$Measured[ii], col = "blue", 
               pch = 2, cex = 1.1)
        points(x = tmp.dte.mod.pflow, y = storms.peak$Modelled[ii], col = "red", 
               pch = 2, cex = 1.1)
  
      # add grid lines in plot for dates
      grid(nx = 30, ny = NULL)
      
      # format x-axis on flow plot to include day, month and year
      axis.Date(side = 1, x = tmp.dte.flows, format = "%m-%d-%Y")
    }
    
    # close file done
    dev.off()
}

storms_plot_to_list <- function(dte.stms = dte.stms, 
                                dte.flows  = dte.flows,
                                obs.flow   = obs.flow,
                                mod.flow   = mod.flow,
                                obs.bflow  = obs.bflow,
                                mod.bflow  = mod.bflow,
                                storms.peak  = storms.peak,
                                storms.vol = storms.vol, 
                                precip     = precip) {
  
  # plot flow and precip information for each individul storm and store
  # in a list
  # input:
  #
  #    dte.stms - data.frame of POSIXct dates that correspond to first and last 
  #               days of storms
  #   dte.flows - vector of POSIXct dates that correspond to flow values
  #    obs.flow - vector of observed daily flow values in cfs
  #    mod.flow - vector of modeled daily flow values in cfs
  #   obs.bflow - vector of observed daily base flow values in cfs
  # storms.peak - data.frame for peak flow data for each storm
  #  storms.vol - data.frame for flow volume data for each storm
  #      precip - vector of numeric precip values in inches
  #    out.file - string for pdf filename (and path) where figure is sent
  #
  # output:
  # my.plots <- list of plots capterted using recordPlot function
  
  # temerpory assignments for function development
#       dte.stms <- df.storms.peak[ , c("begin", "end")]
#       dte.flows  <- df.mflow$dates
#       obs.flow   <- df.mflow$Measured
#       mod.flow   <- df.mflow$Modelled
#       obs.bflow  <- df.hysep88.8.obs$BaseQ
#       mod.bflow  <- df.hysep88.8.mod$BaseQ
#       storms.peak <- df.storms.peak
#       storms.vol <- df.storms.vol
#       obs.pflow  <- storms.peak$Measured
#       mod.pflow  <- storms.peak$Modelled
#       obs.vflow  <- storms.vol$Measured
#       mod.vflow  <- storms.vol$Modelled
#       precip     <- df.daily.precip$daily.precip

  
  ##dev.new(noRStudioGD = FALSE)

  ## create list to put plots in
  my.plots <- vector(length(dte.stms$begin), mode='list')
  
  # creating temporary data sets for storm list data set
  strm.nums <- 1:length(dte.stms$begin)

  # loop to print a figure for each storm to a signle page in "out.file"
  for(ii in 1:(length(strm.nums))) {
    # get info for storm ii
    ##ii <- 25
    ## begin row for storm in flow time-series
    lng.begin <- grep(strftime(dte.stms$begin[ii], format = "%Y%m%d"), 
                      strftime(dte.flows, format = "%Y%m%d"))
    lng.end   <- grep(strftime(dte.stms$end[ii], format = "%Y%m%d"), 
                      strftime(dte.flows, format = "%Y%m%d"))
    if(length(lng.begin) > 1 | length(lng.end) > 1) print(paste0("Too many matches for dates for storm ", ii))
    ## expand range for plotting
    lng.ex <- 2 # number of days before start and after end
    lng.begin.ex <- max(lng.begin - lng.ex, 1)
    lng.end.ex   <- min(lng.end + lng.ex, length(dte.flows))
    ## dates for storm time series
    tmp.dte.stm.ts <- dte.flows[lng.begin:lng.end]
    # get data for current storm
    tmp.dte.flows  <- dte.flows[lng.begin.ex:lng.end.ex]
    tmp.obs.flow   <- obs.flow[lng.begin.ex:lng.end.ex]
    tmp.mod.flow   <- mod.flow[lng.begin.ex:lng.end.ex]
    tmp.obs.bflow  <- obs.bflow[lng.begin.ex:lng.end.ex]
    tmp.mod.bflow  <- mod.bflow[lng.begin.ex:lng.end.ex]
    tmp.dte.obs.pflow <- min(tmp.dte.stm.ts[df.mflow$Measured[lng.begin:lng.end] == storms.peak$Measured[ii]])
    tmp.dte.mod.pflow <- min(tmp.dte.stm.ts[df.mflow$Modelled[lng.begin:lng.end] == storms.peak$Modelled[ii]])
    tmp.precip     <- df.daily.precip$daily.precip[lng.begin.ex:lng.end.ex]
    
    
    
    
    
    # set y-limits for current storm
    tmp.ylims <- 
      c(10 ^ (floor(log10(min(tmp.obs.bflow, tmp.mod.bflow)))), 
        10 ^ (ceiling(log10(max(tmp.obs.flow, tmp.mod.flow)) )))
    
    # set x-limits for current storm
    tmp.xlims <- c(dte.flows[lng.begin.ex], dte.flows[lng.end.ex])
    
    # set plot area matrix for 2 rows and one column along with other pars
    par(mfrow = c(2, 1), tck = 0.01,  mar = c(0, 1, 0, 0.5), 
        oma = c(7, 5, 7, 2))
    
    # vary height of plots. Make precip hight smaller than flow
    layout(matrix(c(1, 1, 2, 2), 2, 2, byrow = TRUE), heights = c(1, 3), 
           widths = c(1, 1))
    
    # precip plot set up, don't plot data
    plot(x = tmp.dte.flows, y = tmp.precip, xlab = "",pch = "", xlim = tmp.xlims,
         ylim = c(0, max(c(tmp.precip, 0.1))), xaxt = "n")
    
    # title for plot is current storm
    title(xlab = "", ylab = "", 
          main = paste0("storm num ", strm.nums[ii],
                        "  Season = ", storms.peak$season[ii],
                        "\nObs Flow-zone = ", storms.peak$obs.flw.zn[ii],
                        "  Mod Flow-zone = ", storms.peak$mod.flw.zn[ii],
                        "\nObs peak flow = ", sprintf(fmt = "%-.2f cfs", storms.peak$Measured[ii]),
                        " Mod peak flow = ", sprintf(fmt = "%-.2f cfs", storms.peak$Modelled[ii]),
                        "  Obs storm vol = ", sprintf(fmt = "%-.2E ac-ft", storms.vol$Measured[ii]),
                        " Mod storm vol = ", sprintf(fmt = "%-.2E ac-ft", storms.vol$Modelled[ii])),
          outer = TRUE, line = 3)
    
    # plot vertical lines for each precip obs
    lines(x = tmp.dte.flows, y = tmp.precip, type = "h")
    
    # add grid lines in plot for dates
    grid(nx = 30, ny = NULL)
    
    # set up plot for flow, don't plot data. Plot flow data after the storm
    # polygons are plotted so flow lines are borders of polygons  
    plot(x = tmp.dte.flows, y = tmp.obs.flow, type = "l", log = "y", lty = "blank",
         xlim = tmp.xlims, ylim = tmp.ylims, xaxt = "n")
    
    # storm flow as a filled in polygon
    #polygon(x = tmp.dte.flows, y = tmp.obs.flow, col = "yellow", lty = "blank")
    
    # flow data storm polygon
    lines(x = tmp.dte.flows, y = tmp.obs.flow, type = "l", col = "blue")
    lines(x = tmp.dte.flows, y = tmp.mod.flow, type = "l", col = "red")
    
    # base flow data storm polygon
    lines(x = tmp.dte.flows, y = tmp.obs.bflow, type = "l", lty = "dashed", col = "blue")
    lines(x = tmp.dte.flows, y = tmp.mod.bflow, type = "l", lty = "dashed", col = "red")
    
    
    # points for peaks ploted over storm polygon and flow data
    points(x = tmp.dte.obs.pflow, y = storms.peak$Measured[ii], col = "blue", 
           pch = 2, cex = 1.1)
    points(x = tmp.dte.mod.pflow, y = storms.peak$Modelled[ii], col = "red", 
           pch = 2, cex = 1.1)
    
    # add grid lines in plot for dates
    grid(nx = 30, ny = NULL)
    
    # format x-axis on flow plot to include day, month and year
    axis.Date(side = 1, x = tmp.dte.flows, format = "%m-%d-%Y")
    
    ## store current plot in list
    my.plots[[ii]] <- recordPlot()
    
    ## reset graphics
    ##graphics.off()
  }
##  dev.off()
return(my.plots)
}

storm_plot <- function(lng.stm,
                        dte.stms = dte.stms,
                        dte.flows  = dte.flows,
                        obs.flow   = obs.flow,
                        mod.flow   = mod.flow,
                        obs.bflow  = obs.bflow,
                        mod.bflow  = mod.bflow,
                        storms.peak  = storms.peak,
                        storms.vol = storms.vol,
                        precip     = precip) {
  
  # plot flow and precip information for each individul storm and store
  # in a list
  # input:
  #
  #     lng.stm - storm number (row number in storm info data.frames)
  #    dte.stms - data.frame of POSIXct dates that correspond to first and last 
  #               days of storms
  #   dte.flows - vector of POSIXct dates that correspond to flow values
  #    obs.flow - vector of observed daily flow values in cfs
  #    mod.flow - vector of modeled daily flow values in cfs
  #   obs.bflow - vector of observed daily base flow values in cfs
  # storms.peak - data.frame for peak flow data for each storm
  #  storms.vol - data.frame for flow volume data for each storm
  #      precip - vector of numeric precip values in inches
  #
  # output:
  # p.storm <- list of precip and flow ggplot objects
  
    
  # temporary assignments for function development
#   dte.stms <- df.storms.peak[ , c("begin", "end")]
#   dte.flows  <- df.mflow$dates
#   obs.flow   <- df.mflow$Measured
#   mod.flow   <- df.mflow$Modelled
#   obs.bflow  <- df.hysep88.8.obs$BaseQ
#   mod.bflow  <- df.hysep88.8.mod$BaseQ
#   storms.peak <- df.storms.peak
#   storms.vol <- df.storms.vol
#   obs.pflow  <- storms.peak$Measured
#   mod.pflow  <- storms.peak$Modelled
#   obs.vflow  <- storms.vol$Measured
#   mod.vflow  <- storms.vol$Modelled
#   precip     <- df.daily.precip$daily.precip
# 
#   lng.stm <- 4

  
  ## load pakage
  require(ggplot2, quietly = TRUE)
  require(gridExtra, quietly = TRUE)

  # title for plot is current storm
  tmp.title <- paste0("storm num ", lng.stm,
                     "  Season = ", storms.peak$season[lng.stm],
                     "\nObs Flow-zone = ", storms.peak$obs.flw.zn[lng.stm],
                     "  Mod Flow-zone = ", storms.peak$mod.flw.zn[lng.stm],
                     "\nObs peak flow = ", sprintf(fmt = "%-.2f cfs", storms.peak$Measured[lng.stm]),
                     " Mod peak flow = ", sprintf(fmt = "%-.2f cfs", storms.peak$Modelled[lng.stm]),
                     "  Obs storm vol = ", sprintf(fmt = "%-.2E ac-ft", storms.vol$Measured[lng.stm]),
                     " Mod storm vol = ", sprintf(fmt = "%-.2E ac-ft", storms.vol$Modelled[lng.stm]))

  ## get data for storm plot  
  ## begin row for storm in flow time-series
  lng.begin <- grep(strftime(dte.stms$begin[lng.stm], format = "%Y%m%d"), 
                    strftime(dte.flows, format = "%Y%m%d"))
  lng.end   <- grep(strftime(dte.stms$end[lng.stm], format = "%Y%m%d"), 
                    strftime(dte.flows, format = "%Y%m%d"))
  if(length(lng.begin) > 1 | length(lng.end) > 1) print(paste0("Too many matches for dates for storm ", lng.stm))
  ## expand range for plotting
  lng.ex <- 2 # number of days before start and after end
  lng.begin.ex <- max(lng.begin - lng.ex, 1)
  lng.end.ex   <- min(lng.end + lng.ex, length(dte.flows))
  ## dates for storm time series
  tmp.dte.stm.ts <- dte.flows[lng.begin:lng.end]
  # get data for current storm
  tmp.dates  <- dte.flows[lng.begin.ex:lng.end.ex]
  tmp.data <- rbind(
    data.frame(date = tmp.dates,
               value = df.daily.precip$daily.precip[lng.begin.ex:lng.end.ex],
               sub.group = "precip",
               type = "precip",
               group = "precip",
               src = "obs"),
    data.frame(date = tmp.dates,
               value = obs.flow[lng.begin.ex:lng.end.ex],
               sub.group = "obs-flow",
               type = "flow",
               group = "flow",
               src = "obs"),
    data.frame(date = tmp.dates,
               value = mod.flow[lng.begin.ex:lng.end.ex],
               sub.group = "mod-flow",
               type = "flow",
               group = "flow",
               src = "mod"),
    data.frame(date = tmp.dates,
               value = obs.bflow[lng.begin.ex:lng.end.ex],
               sub.group = "obs-baseflow",
               type = "baseflow",
               group = "flow",
               src = "obs"),
    data.frame(date = tmp.dates,
               value = mod.bflow[lng.begin.ex:lng.end.ex],
               sub.group = "mod-baseflow",
               type = "baseflow",
               group = "flow",
               src = "mod"),
    data.frame(date = min(tmp.dte.stm.ts[df.mflow$Measured[lng.begin:lng.end] == storms.peak$Measured[lng.stm]]),
               value = storms.peak$Measured[lng.stm],
               sub.group = "obs-peak-flow",
               type = "peak-flow",
               group = "flow",
               src = "obs"),
    data.frame(date = min(tmp.dte.stm.ts[df.mflow$Modelled[lng.begin:lng.end] == storms.peak$Modelled[lng.stm]]),
               value = storms.peak$Modelled[lng.stm],
               sub.group = "mod-peak-flow",
               type = "peak-flow",
               group = "flow",
               src = "mod"))

tmp.data$sub.group <- factor(tmp.data$sub.group, 
                         levels = c("precip", "obs-flow","mod-flow",
                         "obs-baseflow", "mod-baseflow",
                         "obs-peak-flow", "mod-peak-flow"))
tmp.data$group <- factor(tmp.data$group, 
                             levels = c("precip", "flow"))
tmp.data$type <- factor(tmp.data$type,
                        levels = c("precip", "flow", "baseflow", "peak-flow")) 
tmp.data$src <- factor(tmp.data$src, 
                         levels = c("obs", "mod"))

  # calculate pars for plot
  # set y-limits for storm
  tmp.ylims <- 
    c(10 ^ 
        (
          floor(
            log10(
              min(
                tmp.data$value[with(tmp.data, group == "flow")])))),
      10 ^ (
        ceiling(
          log10(
            max(tmp.data$value[with(tmp.data, group == "flow")])))))
  # set x-limits for current storm
  tmp.xlims <- c(dte.flows[lng.begin.ex], dte.flows[lng.end.ex])

  # precip plot
  p.precip <- ggplot(data = tmp.data[tmp.data$sub.group == "precip", ],
                         aes(x = date, y = value)) +
    labs(title = tmp.title) +
    xlim(tmp.xlims) + 
    ylim(c(0, max(c(tmp.data$value[tmp.data$group == "precip"], 0.1)))) +
    ylab("precip in") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.margin=unit(c(1,1,0,1), "lines")) + 
    geom_segment(aes(xend = date, 
                     y = rep(0, length(date)), yend=value))

 ## flow plot
  tmp.flow <- tmp.data[tmp.data$group == "flow", ]
  p.flow <- ggplot() + 
    geom_line(data = tmp.flow[tmp.flow$type != "peak-flow", ], 
              aes(x = date, y = value, 
                  breaks = sub.group, color = src, linetype = type)) +
    geom_point(data = tmp.flow[tmp.flow$type == "peak-flow", ], 
              aes(x = date, y = value, 
                  breaks = sub.group, color = src, shape = src), size = 5) +
    xlab("") + ylab("flow (cfs)") +
    scale_shape_manual(name = "peak-flow", values = c(1,2)) +
    scale_linetype_manual(name = "type of flow", values = c(1, 2)) +
    scale_color_manual(name = "flow source", values = c("blue", "red")) +
    guides(color = FALSE, linetype = FALSE, shape = FALSE) +
    theme(plot.margin=unit(c(0,1,1,1), "lines"))

  p.storm <- list()
  p.storm[[1]] <- p.precip
  p.storm[[2]] <- p.flow
  
  return(p.storm)
}

storm_grid <- function(p.storm) {
  ## wrapper function to grid.arrange to plot storm ggplot
  ## objects in p.storm (precip and flow) output from storm_plot function
  
  grid.thm <- theme(plot.margin=unit(c(1,1,-0.5,1),"lines"))
  
  g.precip <- ggplotGrob(p.storm[[1]] + grid.thm)
  g.flow <- ggplotGrob(p.storm[[2]])
  
  maxWidth = grid::unit.pmax(g.precip$widths[2:5], g.flow$widths[2:5])
  g.precip$widths[2:5] <- as.list(maxWidth)
  g.flow$widths[2:5] <- as.list(maxWidth)
  grid.arrange(g.precip, g.flow, ncol=1, heights = c(0.3, 1))

}


