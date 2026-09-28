library(sp)
library(raster)
library(terra)
library(sf)

library(rjson)
library(geojsonR)

library(doMC)
library(doParallel)


###############################
args <- commandArgs(trailingOnly = TRUE)
print(args)

numSite <- 1
cc      <- as.integer(args[1])
print(numSite)
print(cc)
# numSite <- 1; cc <- 1



###############################
params <- fromJSON(file='/home/hjlee/project/mangrove/code/PLCM_Parameters.json')
source(params$setup$rFunctions)



########################################
## Get site name, image directory and coordinate
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

ckDir <- paste0(params$setup$outDir,strSite,'/chunk_1721')     #_1721, _2024, _2326
print(ckDir)

ckNum <- sprintf('%03d',cc)
file <- list.files(path=ckDir,pattern=glob2rx(paste0('*',ckNum,'.rda')),full.names=T)

load(file)


##########################################
numPix <- dim(band1)[1]
# numPix <- 1
phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr

f_mat <- matrix(NA,numPix,58*length(phenYrs))

for (i in 1:numPix){

  f_mat[i,] <- GetFeatures(band1[i,],band2[i,],band3[i,],band4[i,],dates,phenYrs,params)
  
  if(i%%10000==0) print(i)
}


# Save
ckPheDir <- paste0(params$setup$outDir,strSite,'/chunk_feat_1721')     #_1721, _2024, _2326
if (!dir.exists(ckPheDir)) {dir.create(ckPheDir)}

save(f_mat,file=paste0(ckPheDir,'/chunk_feat_',ckNum,'.rda'))
