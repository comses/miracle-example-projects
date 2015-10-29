## fig_deprec_rdist_did.R

## This function displays an MDS of distances among an outcome set (over
## time) from values of a single experimental variable over selected seeds.

## The following parameters may be altered:

## 'D' is the depreciation mode
##   default: 0
##   valid:   any one of c(0,1,2)

## 'S' is the single run seed associated with the output
##   default: 0
##   valid:   any one of 0:24

## 'ignore' is the data columns to ignore in the MDS
##   default: c(1,2,3,26-0:4)
##   valid:   1:26
##   warning: in order for the MDS to be valid, ignore must contain
##            1:3 and 22:26 (or 26-0:4)

## 'scalingfn' is the MDS scaling function
##   default: cmdscale
##   valid:   cmdscale, monoMDS, isoMDS
##   warning: isoMDS requires library(MASS)
##   warning: monoMDS requires library(vegan)

deployrUtils::deployrPackage("MASS")
deployrUtils::deployrPackage("vegan")

source("rhea/src/lib/ut_did.R"); ## Read in function library

## Read in data if not already loaded
if (!exists("x.deprec.60")) x.deprec.60 = read.csv(deployrUtils::deployrExternal("rhea/data/x_deprec_60.csv", isPublic = TRUE));

## Assign defaults to parameters, if necessary
if (!exists2("data")) data = x.deprec.60; ## Data containing outcomes
#if (!exists2("D")) D = 0; ## Experimental conditions (default: depreciation modes)
deployrUtils::deployrInput('{"name": "DM", "label": "Depreciation mode", "render": "numeric", "default": 2, "valueList": [0, 1, 2]}')
if (!exists2("S")) S = 2; ## Selected seeds to visualize
if (!exists2("ignore")) ignore = c(1,2,3,26-0:4); ## Column ids of variables in the data to ignore in distance matrix calculation
if (!exists2("mode")) mode = "rel"; ## Transformation of the outcomes used in the MDS (default: relative to values at time = 0; alternative: "scale" for standardization)
#if (!exists2("scalingfn")) scalingfn = "cmdscale"; ## Name of MDS scaling function (default: classical MDS, alternatives are "isoMDS" and "monoMDS" which require the MASS and vegan libraries respectively)
deployrUtils::deployrInput('{"name": "scalingfn", "label": "Name of MDS scaling function", "render": "character", "default": "isoMDS", "valueList": ["cmdscale", "isoMDS", "monoMDS", "sammon"]}')
if (!exists2("outfile")) outfile = ""; ## PDF output file name

## Generate the plot/results
fig.deprec.rdist.did(data=data,D=DM,S=S,ignore=ignore,mode=mode,scalingfn=scalingfn,outfile=outfile);
