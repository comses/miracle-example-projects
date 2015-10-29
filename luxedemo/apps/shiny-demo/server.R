# server end of LUXE demo app

library(shiny)
library(sqldf)
library(data.table)
library(raster)
library(rasterVis)
library(grid)

#Function to convert from agents to raster
agentsToRaster = function(x, y, tp) {
    world_size = 61
    #Create a new empty raster layer
    raster_data = raster(nrows=world_size, ncols=world_size, xmn=0, xmx=61, ymn=0, ymx=61, crs=NA, resolution=c(1, 1), vals=0)
    values = rep(0, world_size * world_size)
    #Fill raster space with values from transaction file
    for (i in 1:length(x)) {
        #pos = (60 - floor(y[i])) * 61 + floor(x[i]) + 1
        pos = (floor(y[i])) * 61 + floor(x[i]) + 1
        values[pos] = tp[i]
    }
    #Note that setValues does not change the original raster, need assignment!
    raster_data = setValues(raster_data, values)
    return(raster_data)
}

#Function to generate column names
xnamefun = function(df) {
    xname=array()
    for (i in 1:3) {
        xname[i] = paste0("sdp=",df[i,1],",sdbudget=",df[i,2])
    }
    return(xname)
}

# first, unzip and read files
fpath = "../../data"
runLog = read.csv(paste(fpath, "log", "runLog.csv", sep="/"), header=TRUE)

# prepare data structures
lu = matrix(0, 61, 61)
tp = matrix(0, 61, 61)
agents= matrix(0, 11, 1)

# flag
count = 1

# prepare figure data
s = c("L0", "L0.5", "L1", "L2")
sddf = data.frame(sdp=c(0,0.3,0),sdbudget=c(0,0,30))

shinyServer(function(input, output, session) {
    
  values = reactiveValues(starting=TRUE)
  session$onFlushed(function() {
      values$starting=FALSE
  })
    
  usdp = unique(runLog$sdp)
  usdbudget = unique(runLog$sdbudget)
    
  output$column1 = renderUI({
     inputPanel(
        selectInput("column11", "sdp", usdp, selected=sddf$sdp[1]),
        selectInput("column12", "sdbudget", usdbudget, selected=sddf$sdbudget[1]) 
     )
  })
  
  output$column2 = renderUI({
      inputPanel(
        selectInput("column21", "sdp", usdp, selected=sddf$sdp[2]),
        selectInput("column22", "sdbudget", usdbudget, selected=sddf$sdbudget[2]) 
      )
  })
  
  output$column3 = renderUI({
      inputPanel(
        selectInput("column31", "sdp", usdp, selected=sddf$sdp[3]),
        selectInput("column32", "sdbudget", usdbudget, selected=sddf$sdbudget[3]) 
      )
  })
  
  sddfgen = reactive({
      data.frame(sdp=c(input$column11, input$column21, input$column31), sdbudget=c(input$column12, input$column22, input$column32))
  })
  
  output$distPlot <- renderPlot({
      #Delays the start of renderPlot until page completely loaded. This solves a timing problem
      if (values$starting) return(NULL)
      
      sddf = sddfgen()
      
      x = NULL
      
      for (level in 1:4) { #note this refers to L0, L0.5, L1, L2
          agr = (ceiling(level/2) - 1) * 100
          bid = ifelse(level %% 2, 1, 400)
          for (i in 1:nrow(sddf)) {
              sdp = sddf[i,]$sdp
              sdbudget = sddf[i,]$sdbudget
              data = sqldf(paste("SELECT * FROM runLog WHERE random=51 AND agr=", agr, " AND bid=", bid, " AND sdp=", sdp, " AND sdbudget=", sdbudget, sep=""))
              agents = fread(paste(fpath, "sample", data$agents, sep="/"), header=TRUE)
              map = agentsToRaster(agents$parcel_location_x, agents$parcel_location_y, agents$transaction_price)
              e = extent(10, 52, 10, 52)
              map = crop(map, e)
              if (level == 1 & i == 1) {
                  x = map
              } else {
                  x = stack(x, map)
              }
          }		
      }
      
      old.par = par(no.readonly=TRUE)
      
      colors = colorRampPalette(c("white", "yellow", "red"))

      xnames = c(xnamefun(sddf), rep("", 9))
      
      print(levelplot(x, col.regions=colors(100), layout=c(3,4), names.attr=xnames))
      trellis.focus("toplevel")
      panel.text(0.02, 0.18, "L2", cex=1, font=1)
      panel.text(0.02, 0.4, "L1", cex=1, font=1)
      panel.text(0.02, 0.62, "L0.5", cex=1, font=1)
      panel.text(0.02, 0.84, "L0", cex=1, font=1)
      trellis.unfocus()
      
      par(old.par)

  }, height=750)

})
