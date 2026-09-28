###############################
library(rjson)
params <- rjson::fromJSON(file='/home/hjlee/project/mangrove/code/PLCM_Parameters.json')

geojsonList <- list.files(params$setup$geojsonDir)
sites <- 1:length(geojsonList)

gc()
###############################
# 01_Make image mosaic
setwd(paste0(params$setup$logDir,'01'))
system(paste0('/opt/pbs/bin/qsub -V -l select=1:ncpus=36:mem=200gb:host=gel00 -l walltime=12:00:00 -- /home/hjlee/project/mangrove/code/run_script_01.sh 2'))


# 02_Make image chunks
setwd(paste0(params$setup$logDir,'02'))
numSite = 1
nodes <- c("gel00","gel01","gel02")

for(numSite in numSite){
  nn <- sprintf('%03d',numSite)
  for(i in 1:params$setup$numChunks){
    node <- nodes[(i-1) %% 3 + 1]
    system(paste(
      paste0('/opt/pbs/bin/qsub -l select=1:ncpus=1:mem=8gb:host=', node),
      '-l walltime=06:00:00 --',
      paste0(params$setup$rScripts,'run_script_02.sh'),
      nn, i
    ))
  }
}


# 03_Make features chunks
setwd(paste0(params$setup$logDir,'03'))
numSite <- 1
nodes <- c("gel00", "gel01", "gel02")

for(numSite in numSite){
  
  for(cc in 1:params$setup$numChunks){
    
    node <- nodes[(cc - 1) %% length(nodes) + 1]
    
    system(paste(
      paste0('/opt/pbs/bin/qsub -l select=1:ncpus=1:mem=8gb:host=', node),
      '-l walltime=12:00:00 --',
      paste0(params$setup$rScripts, 'run_script_03.sh'),
      cc
    ))
  }
}


# 04_Fit model
setwd(paste0(params$setup$logDir,'04'))
for(numSite in 1){
  nn <- sprintf('%03d',numSite)
  for(yy in 1:6){
    system(paste('qsub -V -pe omp 4 -l h_rt=12:00:00 ',params$setup$rScripts,'run_script_04.sh ',nn,yy,sep=''))
  }
}


# 04_1_apply model
setwd(paste0(params$setup$logDir,'04_1'))
for(numSite in sites){
  nn <- sprintf('%03d',numSite)
  numYear <- params$setup$phenEndYr - params$setup$phenStartYr + 1
  for(yy in 1:numYear){
    system(paste('qsub -V -pe omp 4 -l h_rt=12:00:00 -m n ',params$setup$rScripts,'run_script_04_1.sh ',nn,yy,sep=''))
  }
}


setwd(paste0(params$setup$logDir, '04'))
for (y in 1:7) {  # yToDo: 1~7 → 2018~2024
  siteCode <- sprintf('%03d', 1)
  yy       <- sprintf('%01d', y)
  jobArg   <- paste0(siteCode, yy)

  system(paste0("qsub -V -pe omp 4 -l h_rt=12:00:00 ",
                params$setup$rScripts, "run_script_model.sh ", jobArg))
}
