library(sp)
library(raster)
library(terra)
library(sf)

library(rjson)
library(geojsonR)

library(doMC)
library(doParallel)

###############################
args <- commandArgs()
print(args)

numSite <- 1

###############################
params <- rjson::fromJSON(file='/home/hjlee/project/mangrove/code/PLCM_Parameters.json')
source(params$setup$rFunctions)


########################################
## Get site name, image directory and coordinate
geojsonDir <- params$setup$geojsonDir

siteInfo <- GetSiteInfo(numSite,geojsonDir,params)

imgDir <- siteInfo[[1]]; strSite <- siteInfo[[2]]
print(paste(strSite,';',imgDir))

cLong <- siteInfo[[3]];cLat <- siteInfo[[4]]
print(paste(cLong,';',cLat))

########################################
## List of files
dfileSR   <- list.files(path=imgDir,pattern=glob2rx('*MS_SR*.tif'),recursive=T)
dfileUDM2 <- list.files(path=imgDir,pattern=glob2rx('*_udm2*.tif'),recursive=T)

fileSR   <- list.files(path=imgDir,pattern=glob2rx('*MS_SR*.tif'),recursive=T,full.names=T)
fileUDM2 <- list.files(path=imgDir,pattern=glob2rx('*_udm2*.tif'),recursive=T,full.names=T)


# Dates
yy <- c(); mm <- c(); dd <- c()
for(i in 1:length(dfileSR)){
  yy[i] <- substr(unlist(strsplit(dfileSR[i],'/'))[4],3,4)
  mm[i] <- substr(unlist(strsplit(dfileSR[i],'/'))[4],5,6)
  dd[i] <- substr(unlist(strsplit(dfileSR[i],'/'))[4],7,8)
  
  dfileSR[i] <- unlist(strsplit(dfileSR[i],'/'))[4]
  dfileUDM2[i] <- unlist(strsplit(dfileUDM2[i],'/'))[4]
}
dates_all <- as.Date(paste(mm,'/',dd,'/',yy,sep=''),'%m/%d/%y')
dates <- unique(dates_all)

#
datesod <- order(dates)
setwd('/home/hjlee/project/mangrove/figures/')
png(filename='numofimage.png',width=7.5,height=6.5,unit='in',res=300)
par(oma=c(1,1,1,1),mar=c(4,4,1,1),mgp=c(2.5,1,0))
plot(dates[datesod],
     xlab='Number of Image',ylab='Dates',cex.axis=1.2,cex.lab=1.5)
dev.off()
#

print(length(dates))



#
# Image process
#
outDir <- paste0(params$setup$outDir,strSite)
if (!dir.exists(outDir)) {dir.create(outDir)}  

## Get Site Shapefile and base image
siteWin <- GetSiteShp(fileSR,cLong,cLat)
# if(numSite==88 | numSite==89 | numSite==97 | numSite==101 | numSite==102){
imgBase <- GetBaseImg(fileSR,siteWin,outDir,save=T)
# }else{
imgBase <- raster(paste0(outDir,'/base_image.tif'))
# }


##
registerDoMC(params$setup$numCores)
registerDoMC()

outDir <- paste0(params$setup$outDir,strSite,'/mosaic')
if (!dir.exists(outDir)) {dir.create(outDir)}  

foreach(dd=1:length(dates)) %dopar% {
  
# for(dd in 1:length(dates)){
  
  imgMulti <- which(substr(dfileSR,1,8)==paste0(substr(dates[dd],1,4),substr(dates[dd],6,7),substr(dates[dd],9,10)))
  
  imgVaild <- c()
  for(mm in 1:length(imgMulti)){
    log <- try({
      img <- raster(fileSR[imgMulti[mm]])
      img <- crop(img,siteWin)
      },silent=T)
    if(inherits(log,'try-error')){ 
      next  
    }else{
      numBand <- nbands(raster(fileSR[imgMulti[mm]]))
      if(numBand==4){
        imgVaild <- c(imgVaild,imgMulti[mm])   
      }
    }
  }
  
  if(length(imgVaild)>0){
    
    imgB <- vector('list',length(imgVaild))
    
    for(mm in 1:length(imgVaild)){
      ii <- imgVaild[mm]
      
      img <- raster(fileSR[ii])
      numBand <- nbands(img)
      
      str <- substr(dfileSR[ii],1,23)
      
      imgP <- vector('list',numBand)
      for(i in 1:numBand){
        imgT <- raster(fileSR[ii],band=i)
        imgT <- crop(imgT,siteWin)
        
        if(length(which(substr(dfileUDM2,1,23)==str))==1){
          log <- try({
            udmT <- raster(fileUDM2[which(substr(dfileUDM2,1,23)==str)],band=8)
            udmT <- crop(udmT,siteWin)
          },silent=T)
          if(inherits(log,'try-error')){ 
            next
          }else{
            imgT[udmT>0] <- NA  
          }
          log <- try({
            udm2T <- raster(fileUDM2[which(substr(dfileUDM2,1,23)==str)])
            udm2T <- crop(udm2T,siteWin)  
          },silent=T)
          if(inherits(log,'try-error')){ 
            imgT <- imgT
          }else{
            imgT[udm2T!=1] <- NA
          }
        }
        imgP[[i]] <- imgT
      }
      
      imgB[[mm]] <- brick(imgP)
    }
    
    temp1 <- vector('list',(length(imgVaild)+1))
    temp2 <- vector('list',(length(imgVaild)+1))
    temp3 <- vector('list',(length(imgVaild)+1))
    temp4 <- vector('list',(length(imgVaild)+1))
    for(i in 1:length(imgVaild)){
      temp1[[i]] <- raster(imgB[[i]],1)
      temp2[[i]] <- raster(imgB[[i]],2)
      temp3[[i]] <- raster(imgB[[i]],3)
      temp4[[i]] <- raster(imgB[[i]],4)
    }
    temp1[[(length(imgVaild)+1)]] <- imgBase
    temp2[[(length(imgVaild)+1)]] <- imgBase
    temp3[[(length(imgVaild)+1)]] <- imgBase
    temp4[[(length(imgVaild)+1)]] <- imgBase
    
    for(i in 1:length(imgVaild)){
      log <- try(compareRaster(temp1[[i]],imgBase,extent=F,rowcol=F),silent=T)
      if(inherits(log,'try-error')){
        temp1[[i]] <- projectRaster(temp1[[i]],imgBase)    
      }
      log <- try(compareRaster(temp2[[i]],imgBase,extent=F,rowcol=F),silent=T)
      if(inherits(log,'try-error')){
        temp2[[i]] <- projectRaster(temp2[[i]],imgBase)    
      }
      log <- try(compareRaster(temp3[[i]],imgBase,extent=F,rowcol=F),silent=T)
      if(inherits(log,'try-error')){
        temp3[[i]] <- projectRaster(temp3[[i]],imgBase)    
      }
      log <- try(compareRaster(temp4[[i]],imgBase,extent=F,rowcol=F),silent=T)
      if(inherits(log,'try-error')){
        temp4[[i]] <- projectRaster(temp4[[i]],imgBase)    
      }
    }
    temp1$fun <- mean; temp2$fun <- mean; temp3$fun <- mean; temp4$fun <- mean
    temp1$na.rm <- T; temp2$na.rm <- T; temp3$na.rm <- T; temp4$na.rm <- T
    rast1 <- do.call(mosaic,temp1)  
    rast2 <- do.call(mosaic,temp2)  
    rast3 <- do.call(mosaic,temp3)  
    rast4 <- do.call(mosaic,temp4)  
    
    # Brick bands
    Rast <- brick(rast1,rast2,rast3,rast4)
    
    # Save
    
    outFile    <- paste0(outDir,'/',substr(dates[dd],1,4),substr(dates[dd],6,7),substr(dates[dd],9,10),'_cliped_mosaic.tif')
    writeRaster(Rast, filename=outFile, format="GTiff", overwrite=TRUE)

    print(outFile)
  } 
} 


print(length(list.files(path=outDir)))

print(length(dates))

