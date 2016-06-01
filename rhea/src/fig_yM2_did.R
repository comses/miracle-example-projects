## fig_yM2_did.R --

## This function implicitly seeks all combinations of the 1st three
## variables (enumerated in the plotted rows) against a fourth variables
## (enumerated in the plotted columns).  These four variables are specified
## in zvar. In each plot, the outcomes of yvar are shown for all seeds in
## S.

## The following parameters may be altered:

## 'yvar' is a pair of any of outcome variables:
##   default: c("nB","nS")
##   valid: any pair of (as quoted strings):
##     avPriceFP0, avPriceFP100, avPriceFP500, avPriceCoastFront,
##     avPriceNowFP0, avPriceNowFP100, avPriceNowFP500, avPriceNowCF,
##     nTradedNowFP0, nTradedNowFP100, nTradedNowFP500, nTradedNowCF,
##     nB, nS, nNC, nOS, nOB, nBL, nT, dBS"

## Refer to UT DID documentation for description of outcome fields

## 'RH','Ins','G' are the ranges of hedonic, insurance, and income
## growth flags.
## **If each of {RH,Is,G,YM} is set to single value, then a single plot
## is produced.
##   default: 0:1
##   valid:   any subset of 0:1

## 'YM' is the set of mortgage years values
##   default: c(10,20,30,40)
##   valid: default or any one of c(10,20,30,40)

## 'S' can be any subset of seeds:
##   default: 0:24
##   valid: any subset of 0:24, e.g., 0:5 or c(1,3,5,9) etc.

source(deployrUtils::deployrExternal("rhea/src/lib/ut_did.R", isPublic = TRUE)); ## Read in function library

## Read in data if not already loaded
if (!exists("x.yM2.200")) x.yM2.200 = read.csv(deployrUtils::deployrExternal("rhea/data/x_yM2_200.csv", isPublic = TRUE));

## Assign defaults to parameters, if necessary
if (!exists2("data")) data = x.yM2.200; ## Dataset containing variables of yvar and zvar
#if (!exists2("S")) S = 0; ## Range or vector of seeds to plot
deployrUtils::deployrInput('{"name": "S", "label": "Seed(s) for the selected model run(s)", "render": "numeric", "default": 0, "valueRange": [0, 24, 1]}')
if (!exists2("mfrow")) mfrow = c(8,4); ## Dimensional arrangement of subplots
if (!exists2("yvar")) yvar = c("nB","nS"); ## Two variables of data to plot
if (!exists2("svar")) svar = "s"; ## Variable with seeds
if (!exists2("zvar")) zvar = c("RH","Ins","G","yM"); ## Vector specifying the experimental conditions/variables
if (!exists2("RH")) RH = 0:1; ## 1st condition values/range (default: realtor hedonics)
if (!exists2("Ins")) Ins = 0:1; ## 2nd condition values/range (default: insurance)
if (!exists2("G")) G = 0:1; ## 3rd condition values/range (default: income growth)
if (!exists2("YM")) YM = c(10,20,30,40); ## 4th condition values/range (default: mortgage years)
if (!exists2("outfile")) outfile = ""; ## PDF output file name

## Generate the plot/results
fig.yM2.did(data=data,S=S,mfrow=mfrow,yvar=yvar,svar=svar,zvar=zvar,RH=RH,Ins=Ins,G=G,YM=YM,outfile=outfile);
