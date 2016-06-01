## fig_nB_did.R --

## This function can plot trajectories of two outcomes for individual runs
## and also graphically highlight those values which exceed a threshold.

## The figures produced by this script are similar to the ones in
## [Filatova et al. 2015] (Figure 5).

## The following parameters may be altered:
## 'yvar' is a pair of any of outcome variables
##   default: c("nB","nS")
##   valid: any pair of (as quoted strings):
##     avPriceFP0, avPriceFP100, avPriceFP500, avPriceCoastFront, nB,
##     nS, nNC, nOS, nOB, nBL, nT, nD

## Refer to UT DID documentation for description of outcome fields

## 'crit' can be any value (dependent on outcome variable) above which
## trails are colored:
##   default: 1300
##   valid:   <any numeric value>

## 'yM' (simple vertical line):
##   default: 30
##   valid:   any value between 1 and 200 

## 'xlab' and 'ylab' may be set to any string

source(deployrUtils::deployrExternal("rhea/src/lib/ut_did.R", isPublic = TRUE)); ## Read in function library

## Read in data if not already loaded
if (!exists("x.bs.200")) x.bs.200 = read.csv(deployrUtils::deployrExternal("rhea/data/x_bs_200.csv", isPublic = TRUE));

## Assign defaults to parameters, if necessary
if (!exists2("data")) data = x.bs.200;
## Dataset/data.frame containing variables in yvar

if (!exists2("yvar")) yvar = c("nB","nS");
## Two outcome variables of data to plot.  Other 

if (!exists2("S")) S = 0:24; ## Values/range of the seeds to consider
if (!exists2("svar")) svar = "s"; ## Variable with the seeds
if (!exists2("tvar")) tvar = "t"; ## Variable with time step values
if (!exists2("yM")) yM = 30; ## Vertical dashed line denoting mortgage years (assumes each tick is semi-annual)
if (!exists2("crit")) crit = 1300; ## Critical value of yvar[1], above which trajectories are highlighted by different colors
if (!exists2("sd.mult")) sd.mult = 0; ## If sd.mult = 0, then use crit directly. If sd.mult > 0, then crit = µ + ·sd.mult (of yvar[1]). If sd.mult = -1, then crit = max(yvar[2])
if (!exists2("xlab")) xlab = "Ticks"; ## Label of x -axis
if (!exists2("ylab")) ylab = "Buyers/Sellers"; ## Label of y -axis
if (!exists2("y0")) y0 = 500; ## Minimum of y-axis values displayed
if (!exists2("outfile")) outfile = ""; ## PDF output file name

## Generate the plot/results
fig.nB.did(data=data,yvar=yvar,S=S,svar=svar,tvar=tvar,yM=yM,crit=crit,sd.mult=sd.mult,xlab=xlab,ylab=ylab,y0=y0,outfile=outfile);
