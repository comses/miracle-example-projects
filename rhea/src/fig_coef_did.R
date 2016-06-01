## fig_coef.R --

## This function shows separate lines for three conditions as defined by
## yvar. Note that the x-axis values are divided by 2 to reflect the
## semi-annual frequency of each tick.

## The following parameters may be altered:

## 'Y' can be any of outcome variables:
##   default: "dBS"
##   valid: any one of (as quoted string): 
##     avPriceFP0, avPriceFP100, avPriceFP500, avPriceCoastFront, nB, nS,
##     nNC, nOS, nOB, nBL, nT, dBS

## Refer to UT DID documentation for description of outcome fields

## 'nbcs' can be any subset of the buyer coefficients values:
##   default: c(0.60,0.65,0.70,0.75,0.80)
##   valid:   any subset of c(0.60,0.65,0.70,0.75,0.80)

source(deployrUtils::deployrExternal("rhea/src/lib/ut_did.R", isPublic = TRUE)); ## Read in function library

## Read in data if not already loaded
if (!exists("x.bs.coef.200")) x.bs.coef.200 = read.csv(deployrUtils::deployrExternal("rhea/data/x_bs_coef_200.csv", isPublic = TRUE));

## Assign defaults to parameters, if necessary
if (!exists2("data")) data = x.bs.coef.200; ## Dataset containing variables of yvar
#if (!exists2("Y")) Y = "dBS"; ## Variable for y-axis
deployrUtils::deployrInput('{"name": "Y", "label": "Variable for the y-axis", "render": "character", "default": "dBS", "valueList": ["avPriceFP0", "avPriceFP100", "avPriceFP500", "avPriceCoastFront", "nB", "nS", "nNC", "nOS", "nOB", "nBL", "nT", "dBS"]}')
if (!exists2("yvar")) yvar = c("RH","Ins","nBC"); ## Variables of the experimental conditions
if (!exists2("RH")) RH = 0:1; ## 1st variables' values/range (default for realtor hedonics)
if (!exists2("Ins")) Ins = 0:1; ## 2nd variables' values/range (default for insurance)
#if (!exists2("nbcs")) nbcs = c(0.60,0.65,0.70,0.75,0.80); ## 3rd variables' values/range (default for new buyer coefficient)
deployrUtils::deployrInput('{"name": "nbcs", "label": "3rd variable values and range (default for new buyer coefficent)", "render": "vector", "default": [0.60,0.65,0.70,0.75,0.80], "valueSet": [0.60,0.65,0.70,0.75,0.80]}')
if (!exists2("len")) len = 15; ## Number of points to plot on top of lines
if (!exists2("ylab")) ylab = "Buyers - Sellers"; ## Label for y -axis
if (!exists2("xlab")) xlab = "Year"; ## Label for x -axis
if (!exists2("legend")) legend = c("-A -I","-A +I","+A -I","+A +I"); ## Legend labels for all combinations of first two conditions
if (!exists2("off.div")) off.div = 0.75; ## Offset divisor for x -axis position to prevent overlapping points
if (!exists2("rel")) rel = FALSE; ## Are outcomes are relative to first value of Y (as a percentage)?
if (!exists2("outfile")) outfile = ""; ## PDF output file name

## Generate the plot/results
fig.coef.did(data=data,Y=Y,yvar=yvar,RH=RH,Ins=Ins,nbcs=nbcs,len=len,ylab=ylab,xlab=xlab,legend=legend,off.div=off.div,rel=rel,outfile=outfile);
