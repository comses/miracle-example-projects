# This produces one sub-plot of Figure 8 in Huang et al. (2013)
# 
###############################################################################

deployrUtils::deployrPackage("data.table")
deployrUtils::deployrPackage("sqldf")
deployrUtils::deployrPackage("lattice")
deployrUtils::deployrPackage("latticeExtra")

#Parameters
deployrUtils::deployrInput('{"name": "sdp2", "label": "Standard deviation of preference to proximity", "render": "numeric", "default": 0.4, "valueList": [0, 0.1, 0.2, 0.3, 0.4, 0.5]}')
deployrUtils::deployrInput('{"name": "sdb3", "label": "Standard deviation of budget", "render": "integer", "default": 30, "valueRange": [0, 50, 10]}')

runLog = read.csv(file=deployrUtils::deployrExternal("luxedemo/data/log/runLog.csv", isPublic=TRUE), header=TRUE)

v = c(1, 400, 101, 500)
s = c("L0", "L0.5", "L1", "L2")

#Add a new column "Level"
l = rep(0, 5760)
bid = runLog$bid
agr = runLog$agr
for (i in 1:nrow(runLog)) {
    v1 = agr[i] + bid[i]
    l[i] = s[match(v1, v)]
}
runLog = cbind(runLog, level=l)
runLog = runLog[c(28,1:27)]

#Add a new column hhh
h = rep(NA, 5760)
hv = c(0, sdp2, sdb3)
hs = c("SDB=0, SDP=0", paste0("SDB=0, SDP=",hv[2]), paste0("SDB=", hv[3], " SDP=0"))
sdp = runLog$sdp
sdbudget = runLog$sdbudget
for(i in 1:nrow(runLog)){
    hv1 = sdp[i] + sdbudget[i]
    h[i] = hs[match(hv1, hv)]
}
runLog = cbind(runLog, hhh=h)

#Create data for the graph
tpdata = NULL
for (level in s) {
    for (hhh in hs) {
        data=sqldf(paste("SELECT * FROM runLog WHERE random=51 AND level='", level, "' AND hhh='", hhh, "'", sep=""))
        agents = fread(deployrUtils::deployrExternal(paste("luxedemo/data/sample", data$agents, sep="/"), isPublic=TRUE), header=TRUE)
        agents = cbind(agents, c(1:nrow(agents)))
        setnames(agents, 1:12, c("ID", "step", "x", "y", "tp", "pre", "dist", "osa", "distc", "budge", "util", "seq"))
        agentsdata = data.frame(level=level, hhh=hhh, x=agents$x, y=agents$y, tp=agents$tp, seq=agents$seq)
        tpdata = rbind(tpdata, agentsdata)
    }
}

tpdata = data.frame(tpdata, index=1:nrow(tpdata))
colors = colorRampPalette(c("red", "green", "blue"))

test.data=tpdata$seq
seqcolors = level.colors(test.data, at=do.breaks(range(test.data), 20), col.regions=colors, colors=TRUE)
#tpdata = cbind(tpdata, stepcolors=stepcolors)
figure8 = cloud(tp ~ x + y | hhh + level, data=tpdata, panel.3d.cloud=panel.3dbars, as.table=TRUE, xbase=1, ybase=1, col.facet = seqcolors[tpdata$index], colorkey = list(col=colors, at=do.breaks(range(test.data),20)))
figure8 = useOuterStrips(figure8)
trellis.device(device="png", filename="figure8.png", width=600, height=750)
print(figure8)
dev.off()