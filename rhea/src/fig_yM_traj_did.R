## fig_yM_traj_did.R -

## This function plots the trajectory of model behavior in a 2D space
## defined two of the outcomes.

## The figures produced by this script are similar to the ones in
## [Filatova et al. 2015] (Figures 6-8).

## The following parameters may be altered:

## 'x' and 'y' can be any of the outcome variables:
##   default 'x': "nT"
##   default 'y': "nB"
##   valid: any of (as quoted string):
##     avPriceFP0, avPriceFP100, avPriceFP500, avPriceCoastFront, nB,
##     nS, nNC, nOS, nOB, nBL, nT, dBS"

## Refer to UT DID documentation for description of outcome fields

## 'fx' and 'fy' are transformation functions:
##   default for both: I
##   valid: any one of c(log,logp1,sqrt) (not as quoted string)

## 's' is the seed of the plot
##   default: 1
##   valid: any one of 0:24

source("rhea/src/lib/ut_did.R"); ## Read in function library

## Read in data if not already loaded
if (!exists("x.yM.10.200")) x.yM.10.200 = read.csv(deployrUtils::deployrExternal("rhea/data/x_yM_10_200.csv", isPublic = TRUE));

## Assign defaults to parameters, if necessary
if (!exists2("data")) data = x.yM.10.200; ## Dataset containing variables x and y
if (!exists2("x")) x = "nT"; ## Variable for x -axis (default = number of trades)
if (!exists2("y")) y = "nB"; ## Variable for y -axis (default = number of buyers)
#if (!exists2("s")) s = 1; ## Single seed value from data to plot (assuming data contains outcomes from multiple seeds)
deployrUtils::deployrInput('{"name": "s", "label": "Random seed for the selected run", "render": "numeric", "default": 0, "valueRange": [0, 24, 1]}')
if (!exists2("fx")) fx = I; ## Transformation function for x variable (default = indicator function)
if (!exists2("fy")) fy = I; ## Transformation function for y variable (default = indicator function)
if (!exists2("xlab")) xlab = "Number of Trades"; ## Label for x -axis
if (!exists2("ylab")) ylab = "Number of Buyers"; ## Label for y -axis
if (!exists2("outfile")) outfile = ""; ## PDF output file name

## Generate the plot/results
fig.yM.traj.did(data=data,x=x,y=y,s=s,fx=fx,fy=fy,xlab=xlab,ylab=ylab,outfile=outfile);
