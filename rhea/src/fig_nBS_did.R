## fig_nBS_did.R --

## This script produces a grid of plots containing two variables over time
## with each plot representing results from a distinct seed.  Sequences of
## increasing values are colored in green.

## The figures produced by this script appear in [Filatova et
## al. 2015] (Figure 4).

## The following parameters may be altered:

## 'yvar' is a pair of any of outcome variables (as quoted strings):
##   default: c("nB","nS")
##   valid: any pair of:
##     avPriceFP0, avPriceFP100, avPriceFP500, avPriceCoastFront, nB,
##     nS, nNC, nOS, nOB, nBL, nT, nD

## Refer to UT DID documentation for description of outcome fields

## 'S' is the seeds of the runs to plot:
##   default: 0:24
##   valid:   any subset of 0:24, e.g., 0:5 or c(1,3,5,9) etc.
## **An 'S' value of a single seed will produce a non-matrix plot

## 'abs', if TRUE sets plots to the same absolute scale in matrix
##   default: TRUE
##   valid:   any one of c(TRUE,FALSE)

## 'lwd' is line width
##   default: .5
##   valid:   reasonable values between c(.5,2.0)

## 'cex' is the green points size
##   default: .75
##   valid:   reasonable values between c(.5,2.0)


source("rhea/src/lib/ut_did.R"); ## Read in function library

## Read in data if not already loaded
if (!exists("x.bs.200")) x.bs.200 = read.csv(deployrUtils::deployrExternal("rhea/data/x_bs_200.csv", isPublic = TRUE));

## Assign defaults to parameters, if necessary
if (!exists2("data")) data = x.bs.200[x.bs.200$RH==1&x.bs.200$Ins==0,]; ## Dataset containing variables in yvar
#if (!exists2("S")) S = 0:24; ## Value(s)/range of seeds to display
deployrUtils::deployrInput('{"name": "S", "label": "Random seed for the selected model run", "render": "numeric", "default": 3, "valueList": [0, 24, 1]}')
if (!exists2("yvar")) yvar = c("nB","nS"); ## Two variables in data to plot
if (!exists2("svar")) svar = "s"; ## Variable with the seeds
if (!exists2("tvar")) tvar = "t"; ## Variable with time step values
if (!exists2("scale")) scale = 1000; ## Dividing factor for yvar (for purposes of displaying shorter numerical labels)
if (!exists2("abs")) abs = TRUE; ## Is the y -axis scale fixed and absolute across all plots?
if (!exists2("lwd")) lwd = .5; ## Line width
if (!exists2("cex")) cex = .75; ## Green point size
if (!exists2("outfile")) outfile = ""; ## PDF output file name

## Generate the plot/results
fig.nBS.did(data=data,S=S,yvar=yvar,svar=svar,tvar=tvar,scale=scale,abs=abs,lwd=lwd,cex=cex,outfile=outfile);
