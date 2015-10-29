# This produces a parameterized version of Figure 1 in Huang et al. (2013)
# 
###############################################################################

deployrUtils::deployrPackage("data.table")
deployrUtils::deployrPackage("raster")
deployrUtils::deployrPackage("rasterVis")
deployrUtils::deployrPackage("sqldf")
deployrUtils::deployrPackage("lattice")
deployrUtils::deployrPackage("latticeExtra")

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

#Read the CSV file
#fpath = "."
runLog = read.csv(file=deployrUtils::deployrExternal("luxedemo/data/log/runLog.csv", isPublic=TRUE), header=TRUE)
#Land use matrix
lu = matrix(0, 61, 61)
#Transaction price matrix
tp = matrix(0, 61, 61)
#An 11 col 1 row matrix of agents (corresponding to one line in transRecords)
agents = matrix(0, 11, 1)

#Parameters
deployrUtils::deployrInput('{"name": "sdp2", "label": "Standard deviation of preference to proximity", "render": "numeric", "default": 0.4, "valueList": [0, 0.1, 0.2, 0.3, 0.4, 0.5]}')
deployrUtils::deployrInput('{"name": "sdb3", "label": "Standard deviation of budget", "render": "integer", "default": 30, "valueRange": [0, 50, 10]}')

sddf = data.frame(sdp = c(0, sdp2, 0), sdbudget = c(0, 0, sdb3))

par(mfrow=c(4,3))

#Create trellis graphs
x = NULL

for (level in 1:4) { #note this refers to L0, L0.5, L1, L2
    agr = (ceiling(level/2) - 1) * 100
    bid = ifelse(level %% 2, 1, 400)
    for (i in 1:nrow(sddf)) {
        sdp = sddf[i,]$sdp
        sdbudget = sddf[i,]$sdbudget
        data = sqldf(paste("SELECT * FROM runLog WHERE random=51 AND agr=", agr, " AND bid=", bid, " AND sdp=", sdp, " AND sdbudget=", sdbudget, sep=""))
        agents = fread(deployrUtils::deployrExternal(paste("luxedemo/data/sample", data$agents, sep="/"), isPublic=TRUE), header=TRUE)
        setnames(agents, 1:11, c("ID", "step", "x", "y", "tp", "pre", "dist", "osa", "distc", "budge", "util"))
        map = agentsToRaster(agents$x, agents$y, agents$tp)
        e = extent(10, 52, 10, 52)
        map = crop(map, e)
        #plot(map, legend=FALSE, xaxt='n', yaxt='n', bty='n', box=FALSE)
        if (level == 1 & i == 1) {
            x = map
        } else {
            x = stack(x, map)
        }
    }		
}

# Print graphs and add legend and text
colors = colorRampPalette(c("white", "yellow", "red"))
xnames = c("homogeneous", paste0("SDP=",sdp2), paste0("SDB=",sdb3), rep("", 9))

trellis.device(device="png", filename="figure1.png")
print(levelplot(x, col.regions=colors(100), layout=c(3,4), names.attr=xnames))
trellis.focus("toplevel")
panel.text(0.08, 0.18, "L2", cex=1, font=1)
panel.text(0.08, 0.4, "L1", cex=1, font=1)
panel.text(0.08, 0.62, "L0.5", cex=1, font=1)
panel.text(0.08, 0.84, "L0", cex=1, font=1)

trellis.unfocus()
dev.off()