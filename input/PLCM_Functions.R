#Functions for Planet LSP - PlanetScope Land Surface Phenology

#Functions written by Josh Gray, Douglas Bolton, and Minkyu
###############


#---------------------------------------------------------------------
#Despike time-series
#Uses the Maja approach to get positive time-series spikes associated with clouds
#See Section 2.1.4 in this document: https://www.theia-land.fr/sites/default/files/imce/produits/atbd_maja_071217.pdf
#The maja approach relies on both the blue and red bands (If change is red is 1.5x the change in blue, it is likely a land condition change)
#In addition, a 3 point method is used to look for negative spikes in EVI2. This step is added because the Maja approach won't capture
#spikes associated with cloud shadow
#Vectorized to avoid looping through dates
#Douglas Bolton
#
# Added the outlier detection for VI 
# that is based on the double-differenced time series, 
# using the median of absolute deviation about the median (MAD) (Sachs, 1996).
# Minkyu Moon
#---------------------------------------------------------------------
CheckSpike_MultiBand <- function(blue,red,vi, dates, pheno_pars){
  
  blue_og <- blue # preserve original vector
  red_og <- red # preserve original vector
  vi_og <- vi # preserve original vector
  dates_og <- as.numeric(dates)
  
  good <- !is.na(blue_og) & !is.na(red_og) & !is.na(vi_og)
  
  x_outs <- matrix(F, length(blue_og)) # create the outlier output vector
  count <- 0 
  while (count < pheno_pars$maxDespikeIterations) {
    count <- count + 1
    
    bS <- blue_og[good] # subset to non missing values
    rS <- red_og[good]
    eS <- vi_og[good] 
    dS <- dates_og[good] #subset date vector 
    
    ind1 <- 1:(length(dS)-2)  #Get indices for first, second, and third images
    ind2 <- 2:(length(dS)-1)
    ind3 <- 3:length(dS)
    
    dDiff1 <- dS[ind2] - dS[ind1]
    bDiff1 <- bS[ind2] - bS[ind1]
    rDiff1 <- rS[ind2] - rS[ind1]
    bTest1 <- bDiff1 > (0.03 * (1 + dDiff1/30))
    rTest1 <- rDiff1 < (1.5 * bDiff1)
    
    dDiff2 <- dS[ind3] - dS[ind2]
    bDiff2 <- bS[ind2] - bS[ind3]    #2 minus 3 because we are investigating 2 as a peak
    rDiff2 <- rS[ind2] - rS[ind3]
    bTest2 <- bDiff2 > (0.03 * (1 + dDiff2/30))
    rTest2 <- rDiff2 < (1.5 * bDiff2)
    
    majaTest <- bTest1 & rTest1 & bTest2 & rTest2
    
    dayFrac <- (dS[ind2]-dS[ind1]) / (dS[ind3]-dS[ind1])   #Calculate time fraction of date 1 to 2 compared to date 1 to 3
    fitVal <- eS[ind1] + (eS[ind3] - eS[ind1]) * dayFrac   #Calculate value at point 2 if a straight line is drawn from point 1 to 3.  
    dev1 <- eS[ind2] - eS[ind1]
    dev2 <- eS[ind2] - eS[ind3]
    dev <- fitVal - eS[ind2]
    devRatio <- dev / (eS[ind3] - eS[ind1])
    dDiff <- dS[ind3] - dS[ind1]
    
    #look for negative spikes in vi
    eTest <- (dev > pheno_pars$minResid) & (abs(devRatio) > pheno_pars$spikeThresh) & (dDiff < pheno_pars$maxDistance)   
    
    
    # Spikes in vi based on the double-differenced time series
    eDiff <- (eS[ind2] - eS[ind1]) - (eS[ind3] - eS[ind2])

    z   <- pheno_pars$MADspikeThresh

    Md  <- median(eDiff)
    MAD <- median(abs(eS[ind2]-Md))

    madTest <- eDiff < (Md-(z*MAD/0.6745)) | eDiff > (Md+(z*MAD/0.6745))
    
    
    ##
    check <- majaTest | eTest | madTest
    # check <- majaTest | eTest
    check <- c(FALSE,check,FALSE)
    check[is.na(check) | is.infinite(check)] <- FALSE
    if (sum(check) == 0) {break}    #Break if no observations have been despiked
    x_outs[good] <- check   #expand to size of original x, accounting for missing values
    good[x_outs] <- FALSE              #remove the despiked values from the pixels of interest and try again 
  }
  return(x_outs)
}




#----------------------------------------------------------
# Fit a cubic spline to the Planet time-series
# Written by Josh Gray, Douglas Bolton, and Minkyu Moon
#----------------------------------------------------------
Smooth_VI <- function(x, dates, pred_dates, pheno_pars, dormant_value) {
    #Get index of pixels with good values
    ind <- !is.na(x)  
    # smooth with a spline to get continuous daily series
    spl <- smooth.spline(dates[ind], x[ind], spar=pheno_pars$splineSpar)
    # weighted version
    xSmooth <- predict(spl, as.numeric(pred_dates))$y
    # screen and fill values less than the the dormant value
    xSmooth[xSmooth < dormant_value] <- dormant_value
    return(xSmooth)
}


Smooth_Bands <- function(x, dates, pred_dates, pheno_pars){
  #Get index of pixels with good values
  ind <- !is.na(x)  
  # smooth with a spline to get continuous daily series
  spl <- smooth.spline(dates[ind], x[ind], spar=pheno_pars$splineSpar)
  # weighted version
  xSmooth <- predict(spl, as.numeric(pred_dates))$y
  
  # determine upper and lower bound of good data (set in json file)
  lowBound <- quantile(x, probs=pheno_pars$bandLimits[1],na.rm=T) 
  upperBound <- quantile(x, probs=pheno_pars$bandLimits[2],na.rm=T) 
  
  xSmooth[xSmooth < lowBound] <- lowBound
  xSmooth[xSmooth > upperBound] <- upperBound
  return(xSmooth)
}



#----------------------------------------------------------
#New from Josh - 2018-10-31
#Finds time-series peaks
#Josh Gray
#----------------------------------------------------------
FindPeaks <- function(x, mag_order=T){
  # Function to identify peaks in time series x (or troughs if x=-x), supports "flat top" peaks
  # if mag_order is TRUE, peaks are returned in order of increasing magnitude (of x)
  d <- diff(x)
  d_code <- (d > 0) + (2 * (d < 0)) # 0=no change, 1=inc, 2=dec
  peaks <- unlist(gregexpr("12", paste(d_code, collapse=""))) # no match is -1
  if(peaks[1] == -1) peaks <- NULL
  flat_peaks <- unlist(gregexpr("10+2", paste(d_code, collapse=""))) # no match is -1
  if(flat_peaks[1] == -1) flat_peaks <- NULL
  d_code_rle <- rle(d_code)
  flat_peaks <- flat_peaks + round(d_code_rle$l[match(flat_peaks, cumsum(d_code_rle$l)) + 1] / 2)
  # all_peaks <- c(ifelse(peaks[1] == -1, NULL, peaks + 1), ifelse(flat_peaks[1] == -1, NULL, flat_peaks + 1))
  peaks <- sort(c(peaks + 1, flat_peaks + 1))
  if(mag_order) return(peaks[order(x[peaks])])
  return(peaks)
}


#----------------------------------------------------------
#New from Josh - 2018-10-31
#Determines valid segments 
#Josh Gray
#----------------------------------------------------------
GetSegs <- function(peaks, x, pars, peak=NA){
  # identifies valid increasing-decreasing segments in x subject to the parameters in pars
  # returns a list of segments: c(start, peak, end). DON'T call directly w/ peak!=NA
  # NOTE: returned segments will not necessarily be in order, and may not completely partition x
  
  # ensure that peaks are in increasing order of x's magnitude
  tmp_peaks <- peaks[order(x[peaks])] # so we only have to sort once if they're in the wrong order
  if(!identical(tmp_peaks, peaks)) peaks <- tmp_peaks
  
  # if no peak is specified, we start at the beginning
  if(is.na(peak)) peak <- peaks[1]
  
  # get the next largest peak; will be NA if this peak is the highest (last one to do)
  next_highest_peak <- peaks[which(peaks == peak) + 1]
  
  # check if we're doing C5-style relative amplitude and peak identification
  # if(!is.na(pars$rel_amp_frac) & !is.na(pars$rel_peak_frac)){
  #   global_max <- max(x, na.rm=T)
  #   seg_thresh <- (global_max - min(x, na.rm=T)) * pars$rel_amp_frac
  #   peak_thresh <- global_max * pars$rel_peak_frac
  # }else{
  #   seg_thresh <- pars$min_seg_amplitude
  #   peak_thresh <- 0
  # }
  
  # we could have any combination of rel_amp_frac, rel_peak_frac, and min_seg_amplitude specified
  # initialize seg_thresh and peak_thresh to zero
  # determine the "global max/min", if peak_frac is specified, set it, if amp_frac is specified, set it
  # if min_seg_amplitude is set, choose the max of that and amp_frac
  seg_thresh <- peak_thresh <- 0
  global_max <- max(x[(pars$splineBuffer+1):(pars$splineBuffer+365)], na.rm=T) #find gobal min/max within a target year
  global_min <- min(x[(pars$splineBuffer+1):(pars$splineBuffer+365)], na.rm=T)
  if(!is.na(pars$rel_amp_frac)) seg_thresh <- (global_max - global_min) * pars$rel_amp_frac
  #if(!is.na(pars$rel_peak_frac)) peak_thresh <- global_max * pars$rel_peak_frac
  if(!is.na(pars$min_seg_amplitude)) seg_thresh <- max(pars$min_seg_amplitude, seg_thresh)
  
  # checks if the period preceding the peak covers enough amplitude
  # search before the peak up to the maximum of: previous peak, the head of x, or the peak - max_increase_length
  previous_peaks <- peaks[peaks - peak < 0]
  previous_peak <- NA
  if(length(previous_peaks) > 0) previous_peak <- max(previous_peaks)
  search_start <- max(1, peak - pars$max_increase_length, previous_peak, na.rm=T)
  search_end <- peak
  # get the index of the closest minimum value within the search window
  # NOTE: should maybe retrieve the troughs here with FindPeaks(-x) instead
  # in the event of repeated minimum values, we take the closest one here
  inc_min_ind <- max(which(x[search_start:search_end] == min(x[search_start:search_end], na.rm=T)) + search_start - 1, na.rm=T)
  seg_amp <- x[peak] - x[inc_min_ind] # get the increasing segment amplitude
  # if(seg_amp > pars$min_seg_amplitude){
  if((seg_amp >= seg_thresh) & (x[peak] >= peak_thresh)){
    # check for a valid decreasing segment
    next_peaks <- peaks[peaks - peak > 0]
    next_peak <- NA
    if(length(next_peaks) > 0) next_peak <- min(next_peaks)
    # search after the peak up to the minimum of: next peak, the tail of x, or the max_decrease_length
    search_start <- peak
    search_end <- min(length(x), peak + pars$max_decrease_length, next_peak, na.rm=T)
    # get the index of the closest minimum value within the search window
    # NOTE: see above note about finding troughs instead
    dec_min_ind <- min(which(x[search_start:search_end] == min(x[search_start:search_end], na.rm=T)) + search_start - 1, na.rm=T)
    seg_amp <- x[peak] - x[dec_min_ind] # get the decreasing segment amplitude
    # if(seg_amp > pars$min_seg_amplitude){
    if(seg_amp >= seg_thresh){
      # we found a valid segment, store it as a list with a single vector: c(start, peak, end)
      tmp_seg <- list(c(inc_min_ind, peak, dec_min_ind))
      # if this isn't the last peak, then call CheckSegRec again w/ next highest peak
      if(!is.na(next_highest_peak)){
        return(c(tmp_seg, GetSegs(peaks, x, pars, peak=next_highest_peak)))
      }else{
        # that was the last peak, and it was valid
        return(tmp_seg) # covers the case where there's only one valid peak
      }
    }else{
      # increase was valid, but decrease was not
      peaks <- peaks[-which(peaks == peak)] # remove peak from peaks list
      # if this isn't the last peak, then call CheckSegRec again w/ next highest peak
      if(!is.na(next_highest_peak)){
        return(GetSegs(peaks, x, pars, peak=next_highest_peak))
      }else{
        # that was the last peak, and it was invalid
        return(NULL)
      }
    }
  }else{
    # increase segment not valid
    peaks <- peaks[-which(peaks == peak)] # remove peak from peaks list
    # if this isn't the last peak, then call CheckSegRec again w/ next highest peak
    if(!is.na(next_highest_peak)){
      return(GetSegs(peaks, x, pars, peak=next_highest_peak))
    }else{
      # that was the last peak, and it was invalid
      return(NULL)
    }
  }
}


#----------------------------------------------------------
#Get phenology dates from segments. Also pull the peak date
#Josh Gray. Updated by Douglas Bolton to include peak date and cleaned
#----------------------------------------------------------
GetPhenoDates <- function(segs, x, dates, pheno_pars){
  pheno_dates <- list()
  
  #Pull greenup dates
  for(gup_thresh in pheno_pars$gup_threshes){
    pheno_dates <- c(pheno_dates, list(dates[unlist(lapply(segs, GetSegThresh, x, gup_thresh, gup=T), use.names=F)]))
  }
  
  #Pull peak dates
  pheno_dates <- c(pheno_dates, list(dates[sapply(segs, "[[", 2)]))
  
  #Pull greendown dates
  for(gdown_thresh in pheno_pars$gdown_threshes){
    pheno_dates <- c(pheno_dates, list(dates[unlist(lapply(segs, GetSegThresh, x, gdown_thresh, gup=F), use.names=F)]))
  }
  return(pheno_dates)
}


#----------------------------------------------------------
#Josh Gray
#----------------------------------------------------------
GetThresh <- function(thresh_value, x, first_greater=T, gup=T){
  # returns the index of the first/last value  of x that is greater/less than the value of thresh.
  # If gup is False (greendown) then it returns the first/last value of x that is less/greater than
  # the value of thresh. first/last and greater/less determined by first_greater
  # NOTE: if thresh is 1 or 0, rounding error can be a problem. Now we round the threshold and each
  # of the evi values to 6 decimal places to compensate
  
  if(gup){
    if(first_greater){
      return(min(which(round(x, 6) >= round(thresh_value, 6))))
    }else{
      return(max(which(round(x, 6) <= round(thresh_value, 6))))
    }
  }else{
    if(first_greater){
      return(min(which(round(x, 6) <= round(thresh_value, 6))))
    }else{
      return(max(which(round(x, 6) >= round(thresh_value, 6))))
    }
  }
}

#----------------------------------------------------------
#Josh Gray
#----------------------------------------------------------
GetSegThresh <- function(seg, x, thresh, gup=T){
  if(gup){
    # check for valid greenup segment
    if(!is.na(seg[1]) & !is.na(seg[2])){
      gup_thresh <- x[seg[1]] + ((x[seg[2]] - x[seg[1]]) * thresh)
      gup_thresh_index <- GetThresh(gup_thresh, x[seg[1]:seg[2]], first_greater=T, gup=T)
      return(gup_thresh_index + seg[1] - 1)
    }else{
      return(NA)
    }
  }else{
    # check for valid greendown segment
    if(!is.na(seg[2]) & !is.na(seg[3])){
      gdown_thresh <- x[seg[3]] + ((x[seg[2]] - x[seg[3]]) * thresh)
      gdown_thresh_index <- GetThresh(gdown_thresh, x[seg[2]:seg[3]], first_greater=F, gup=F)
      return(gdown_thresh_index + seg[2] - 1)
    }else{
      return(NA)
    }
  }
}



#----------------------------------------------------------
#Developed by Josh Gray, updated by Douglas Bolton for HLS
#----------------------------------------------------------
GetSegMetrics <- function(seg, x_smooth, x_raw, smooth_dates, raw_dates){
  if(any(is.na(seg))){return(NA)}
  # get the subset of the smoothed and original time series
  tmp_seg_smooth <- x_smooth[seg[1]:seg[3]]
  tmp_gup_smooth <- x_smooth[seg[1]:seg[2]]
  tmp_gdown_smooth <- x_smooth[seg[2]:seg[3]]
  
  # get the full segment minimum/maximum SVI
  seg_min <- min(tmp_seg_smooth, na.rm=T)
  seg_max <- max(tmp_seg_smooth, na.rm=T)
  seg_amp <- seg_max - seg_min
  
  # get the segment integrated SVI: the sum of values.
  #For MODIS C6, this is the sum of values above the minimum evi. 
  seg_int <- sum(tmp_seg_smooth)
  
  
  # organize greenup segment
  ######################################
  gup_raw_date_inds <- which(raw_dates >= smooth_dates[seg[1]] & raw_dates <= smooth_dates[seg[2]]) # indices in raw data of gup segment
  gup_smooth_date_inds <- match(raw_dates[gup_raw_date_inds], smooth_dates) # indices of raw dates in smooth dates
  
  raw_dates_gup <- raw_dates[gup_raw_date_inds]
  gup_raw_data <- x_raw[gup_raw_date_inds] # get the raw data associated with the gup segment (this is the pre-filled, despiked version)
  gup_smooth_data <- x_smooth[gup_smooth_date_inds] # get the smoothed values associated with each raw data value
  
  gup_numObs <- sum(!is.na(gup_raw_data)) 
  
  
  # organize greendown segment
  ######################################
  gdown_raw_date_inds <- which(raw_dates >= smooth_dates[seg[2]] & raw_dates <= smooth_dates[seg[3]]) # indices in raw data of gdown segment
  gdown_smooth_date_inds <- match(raw_dates[gdown_raw_date_inds], smooth_dates) # indices of raw dates in smooth dates
  
  raw_dates_gdown <- raw_dates[gdown_raw_date_inds]
  gdown_raw_data <- x_raw[gdown_raw_date_inds] # get the raw data associated with the gdown segment (this is the pre-filled, despiked version)
  gdown_smooth_data <- x_smooth[gdown_smooth_date_inds] # get the smoothed values associated with each raw data value
  
  gdown_numObs <- sum(!is.na(gdown_raw_data)) 
  
  
  if (gup_numObs == 0 | gdown_numObs == 0) {return(rep(NA,9))}
  
  
  ###Get the observation density for each period
  #This approach counts snow filled values as good values, since snow images are valuable for pinning down dormant period
  ###ind the biggest gap between images
  gup_seg_rsquared <- 1 - (sum((gup_raw_data - gup_smooth_data)^2, na.rm=T) / sum((gup_raw_data - mean(gup_raw_data, na.rm=T))^2, na.rm=T))
  gup_seg_rsquared[is.infinite(gup_seg_rsquared)] <- NA
  gup_maxgap <- max(diff(c(smooth_dates[seg[1]],raw_dates_gup[!is.na(gup_raw_data)],smooth_dates[seg[2]])))
  # gup_maxgap_frac <- gup_maxgap / (seg[2] - seg[1]) 
    
  gdown_seg_rsquared <- 1 - (sum((gdown_raw_data - gdown_smooth_data)^2, na.rm=T) / sum((gdown_raw_data - mean(gdown_raw_data, na.rm=T))^2, na.rm=T))
  gdown_seg_rsquared[is.infinite(gdown_seg_rsquared)] <- NA
  gdown_maxgap <- max(diff(c(smooth_dates[seg[2]],raw_dates_gdown[!is.na(gdown_raw_data)],smooth_dates[seg[3]])))
  # gdown_maxgap_frac <- gdown_maxgap / (seg[3] - seg[2]) 
  
  return(c(seg_amp, seg_max, seg_int, 
           gup_seg_rsquared, gup_numObs, gup_maxgap, 
           gdown_seg_rsquared, gdown_numObs, gdown_maxgap))
}



#----------------------------------------------------------
#When a cycle is not detected, return a subset of metrics for the calendar year
#For now, returning evi amplitude, evi maximum, and evi area
#Douglas Bolton
#adapted by Minkyu Moon  
#----------------------------------------------------------
annualMetrics <- function(viSub, dateSub, smoothed_bns, pred_dates, yr, pheno_pars, vi_dorm) {
  
  out <- c(rep(NA,40))
  try({
    inyear <- as.numeric(format(dateSub,'%Y')) == yr
    viObs_inyear  <- viSub[inyear]
    numObs <- sum(!is.na(viObs_inyear))
    
    inyear <- as.numeric(format(pred_dates,'%Y')) == yr
    
    for(i in 1:8){
      vi_inyear  <- smoothed_bns[inyear,i]
      seg_min <- min(vi_inyear,na.rm=T) * 10000  
      seg_max <- max(vi_inyear,na.rm=T) * 10000  
      seg_med <- median(vi_inyear,na.rm=T) * 10000  
      seg_int <- sum(vi_inyear,na.rm=T) * 100
      seg_amp <- seg_max-seg_min  
      
      out[(5*(i-1)+1):(5*i)] <- c(seg_min, seg_max, seg_med, seg_int, seg_amp)
    }
    
    if((seg_max > 10000 | seg_max < 0 | seg_amp > 10000 | seg_amp < 0) ){
      out <- c(rep(NA,40))
      return(out)
    }
    
  },silent=T)
  
  return(out)
}


#----------------------------------------------------------
#Make image composites for specific dates from a smoothed time-series
#Does all bands at the same time
#Douglas Bolton
#----------------------------------------------------------
MakeComposites <- function(compositeDates,pred_dates,smooth_b2, smooth_b3, smooth_b4, smooth_b5, smooth_b6, smooth_b7) {
  
  days <- match(compositeDates,pred_dates)
  
  #Ensure that all days were found. If not, don't return any, as this means something is incorrect
  if (length(days) != length(compositeDates)) {return(NA)}
  
  out <- rbind(smooth_b2[days],smooth_b3[days],smooth_b4[days],
               smooth_b5[days],smooth_b6[days],smooth_b7[days])
  
  return(out)
  
}


#---------------------------------------------------------------------
#Calculate QA 
#Code adapted by Minkyu Moon  
#---------------------------------------------------------------------
GetQAs <- function(gup_rsq, gdown_rsq, gup_maxgap, gdown_maxgap, theOrd, qa_pars){

  if(length(theOrd)==1){
    Rsq    <- min(c(gup_rsq/10000,gdown_rsq/10000))
    maxGap <- max(c(gup_maxgap,gdown_maxgap))
    
    if(maxGap <= qa_pars$maxGap_high_quality){
      if(Rsq >= qa_pars$min_r2_high_quality){
        qual_1 <- 1
      }else{
        qual_1 <- 2
      }
    }else{
      if(Rsq >= qa_pars$min_r2_high_quality){
        qual_1 <- 2
      }else{
        qual_1 <- 3
      }
    }
    
    return(qual_1)
    
  }else{
    # First cycle
    Rsq    <- min(c(gup_rsq[theOrd[1]]/10000,gdown_rsq[theOrd[1]]/10000))
    maxGap <- max(c(gup_maxgap[theOrd[1]],gdown_maxgap[theOrd[1]]))
    
    if(maxGap <= qa_pars$maxGap_high_quality){
      if(Rsq >= qa_pars$min_r2_high_quality){
        qual_1 <- 1
      }else{
        qual_1 <- 2
      }
    }else{
      if(Rsq >= qa_pars$min_r2_high_quality){
        qual_1 <- 2
      }else{
        qual_1 <- 3
      }
    }
    # Second cycle
    Rsq    <- min(c(gup_rsq[theOrd[2]]/10000,gdown_rsq[theOrd[2]]/10000))
    maxGap <- max(c(gup_maxgap[theOrd[2]],gdown_maxgap[theOrd[2]]))
    
    if(maxGap <= qa_pars$maxGap_high_quality){
      if(Rsq >= qa_pars$min_r2_high_quality){
        qual_2 <- 1
      }else{
        qual_2 <- 2
      }
    }else{
      if(Rsq >= qa_pars$min_r2_high_quality){
        qual_2 <- 2
      }else{
        qual_2 <- 3
      }
    }
  }

  return(list(qual_1,qual_2))
  
}



#---------------------------------------------------------------------
#Calculate weights
#Douglas Bolton
#---------------------------------------------------------------------
calculateWeights <- function(smoothMat_Masked, numDaysFit, numYrs, pheno_pars) {
  outWeights <- array(0, dim = c(numDaysFit, numYrs, numYrs))
  for (y in 1:numYrs) {
    #Only compare on dates that have splined data in target year
    ind <- !is.na(smoothMat_Masked[,y])
    sub_vi <- smoothMat_Masked[ind,]
    numGoodDays <- colSums(!is.na(sub_vi))  #how many days actually have splined data?
    
    #What approach to use for weighting
    ######
    #Calculate euclidean distance
    eucl <- colSums((sub_vi - sub_vi[,y])^2,na.rm=T)^0.5 
    
    #Now calculate euclidean distance assuming the average through the year
    #Scale euculidean distances between this value and a perfect fit (0 to 1)
    theAvg <- matrix(mean(sub_vi[,y],na.rm=T),length(sub_vi[,y]),numYrs,byrow=T)
    theAvg[is.na(sub_vi)] <- NA   #only calculate for days that have data
    
    max_eucl <- colSums((theAvg - sub_vi[,y])^2,na.rm=T)^0.5    #calculate eucidean distance for this case
    scaled_eucl <- 1 - (eucl / max_eucl)
    scaled_eucl[scaled_eucl < 0] <- 0
    
    #Weigh as the scaled euclidean distance (0 = same/worse than assuming average, 1 = perfect fit)
    weight <- pheno_pars$maxWeight * scaled_eucl
    weight[numGoodDays < pheno_pars$minDaysForSplineComparison]  <- 0
    weight[is.na(weight)] <- 0
    weight[is.infinite(weight)] <- 0
    outWeights[,,y] <- matrix(weight,numDaysFit,numYrs,byrow=T)
  }
  return(outWeights)
}




#---------------------------------------------------------------------
# Get feature metrics for each pixel
# Code by Minkyu Moon  
#---------------------------------------------------------------------
GetFeatures <- function(blue, green, red, nir, dates, phenYrs, params){
  
  # Despike, calculate dormant value, fill negative VI values with dormant value
  log <- try({    
    
    pheno_pars <- params$phenology_parameters
    qa_pars    <- params$qa_parameters
    
    b2 <- blue/10000; b3 <- green/10000; b4 <- red/10000; b5 <- nir/10000
    
    i1   <- (b5 - b4) / (b5 + b4) # NDVI
    i2   <- (b5 - b3) / (b5 + b3) # GNDVI
    i3   <- b5 * i1               # NIRv
    vi   <- 2.5*(b5 - b4) / (b5 + 2.4*b4 + 1) # EVI2
    
    
    # Spikes check, and remove
    spikes     <- CheckSpike_MultiBand(b2, b4, vi, dates, pheno_pars)
    b2[spikes] <- NA; b3[spikes] <- NA; b4[spikes] <- NA; b5[spikes] <- NA
    i1[spikes] <- NA; i2[spikes] <- NA; i3[spikes] <- NA; vi[spikes] <- NA
    
    # Replace negative VIs with dormant value
    dormIms <- dates >= pheno_pars$dormStart & dates <= pheno_pars$dormEnd
    vi_dorm <- quantile(vi[dormIms & vi>0],probs=pheno_pars$dormantQuantile,na.rm=T)   # Calc vi dormant value using non-negative VIs
    
    
    #now calculate dormancy values and fill individual bands
    dormObs <- dormIms & vi < vi_dorm    #Defining dormant observations for bands as median on dates when vi < vi_dorm
    b2_dorm <- median(b2[dormObs], na.rm=T); b2[dormObs] <- b2_dorm
    b3_dorm <- median(b3[dormObs], na.rm=T); b3[dormObs] <- b3_dorm
    b4_dorm <- median(b4[dormObs], na.rm=T); b4[dormObs] <- b4_dorm
    b5_dorm <- median(b5[dormObs], na.rm=T); b5[dormObs] <- b5_dorm
    i1_dorm <- median(i1[dormObs], na.rm=T); i1[dormObs] <- i1_dorm
    i2_dorm <- median(i1[dormObs], na.rm=T); i2[dormObs] <- i2_dorm
    i3_dorm <- median(i1[dormObs], na.rm=T); i3[dormObs] <- i3_dorm
    vi[vi < vi_dorm] <- vi_dorm
    
    
    #
    splineStart <- as.Date(as.Date(paste0(phenYrs,'-01-01')) - pheno_pars$splineBuffer) 
    numDaysFit  <- 365 + (pheno_pars$splineBuffer * 2)    
    splineEnd   <- splineStart+(numDaysFit-1)
    all_dates   <- seq(min(splineStart), max(splineEnd), by="day")
    
    daysVec <- 1:numDaysFit
    inYear  <- daysVec > pheno_pars$splineBuffer & daysVec <= (pheno_pars$splineBuffer+365)
    prevYear <- daysVec <= pheno_pars$splineBuffer
    nextYear <- daysVec > (pheno_pars$splineBuffer+365)
    
    numYrs <- length(phenYrs)
    vecLength <- numDaysFit*numYrs
    
    #
    smoothMat <- matrix(NA, numDaysFit, numYrs)
    b2Mat <- smoothMat; b3Mat <- smoothMat; b4Mat <- smoothMat; b5Mat <- smoothMat
    i1Mat <- smoothMat; i2Mat <- smoothMat; i3Mat <- smoothMat
    
  },silent=TRUE)
  #If there is an error despiking or other initial steps, return NAs
  if(inherits(log, "try-error")){return(matrix(NA,58*length(phenYrs)))}   
  
  
  outAll <- c()
  for(y in 1:numYrs){
    log <- try({    
      
      pred_dates <- seq(splineStart[y], splineEnd[y], by="day")
      
      dateRange <- dates >= splineStart[y] & dates <= splineEnd[y] & !is.na(vi)   
      dateSub <- dates[dateRange]
      viSub <- vi[dateRange]
      
      smoothed_vi <- Smooth_VI(viSub, dateSub, pred_dates, pheno_pars, vi_dorm)  #Fit spline
      
      # Number of clear observation
      numObs <- sum(dateSub > splineStart[y]+184 & dateSub < splineEnd[y]-184)   #Number of observations in year
      
      
      # Fit spline to individual bands
      smooth_b2 <- Smooth_Bands(b2[dateRange], dateSub, pred_dates, pheno_pars)
      smooth_b3 <- Smooth_Bands(b3[dateRange], dateSub, pred_dates, pheno_pars)
      smooth_b4 <- Smooth_Bands(b4[dateRange], dateSub, pred_dates, pheno_pars)
      smooth_b5 <- Smooth_Bands(b5[dateRange], dateSub, pred_dates, pheno_pars)
      smooth_i1 <- Smooth_Bands(i1[dateRange], dateSub, pred_dates, pheno_pars)
      smooth_i2 <- Smooth_Bands(i2[dateRange], dateSub, pred_dates, pheno_pars)
      smooth_i3 <- Smooth_Bands(i3[dateRange], dateSub, pred_dates, pheno_pars)
      
      smoothed_bns <- cbind(smooth_b2, smooth_b3, smooth_b4, smooth_b5, smooth_i1, smooth_i2, smooth_i3, smoothed_vi)
      
      
      
      ################################################
      #Fit phenology
      peaks <- FindPeaks(smoothed_vi)
      if (all(is.na(peaks))) {outAll <- c(outAll,annualMetrics(viSub,dateSub,smoothed_bns,pred_dates,phenYrs[y],pheno_pars,vi_dorm),rep(NA,18));next}
      
      #Find full segments
      full_segs <- GetSegs(peaks, smoothed_vi, pheno_pars)
      if (is.null(full_segs)) {outAll <- c(outAll,annualMetrics(viSub,dateSub,smoothed_bns,pred_dates,phenYrs[y],pheno_pars,vi_dorm),rep(NA,18));next}
      
      #Only keep segments with peaks within year *****
      full_segs <- full_segs[inYear[sapply(full_segs, "[[", 2)] ]  #check if peaks are in the year
      if (length(full_segs)==0) {outAll <- c(outAll,annualMetrics(viSub,dateSub,smoothed_bns,pred_dates,phenYrs[y],pheno_pars,vi_dorm),rep(NA,18));next}
      
      #Get PhenoDates
      pheno_dates <- GetPhenoDates(full_segs, smoothed_vi, pred_dates, pheno_pars)
      phen <- unlist(pheno_dates, use.names=F)
      phen <- phen - as.numeric(as.Date(paste0((as.numeric(phenYrs[y])-1),'-12-31')))
      if (all(is.na(phen))) {outAll <- c(outAll,annualMetrics(viSub,dateSub,smoothed_bns,pred_dates,phenYrs[y],pheno_pars,vi_dorm),rep(NA,18));next}
      
      #EVI layers
      seg_metrics <- lapply(full_segs, GetSegMetrics,smoothed_vi,viSub,pred_dates,dateSub) #full segment metrics
      un <- unlist(seg_metrics, use.names=F)
      ln <- length(un)
      seg_amp <- un[seq(1, ln, by=9)] * 10000
      seg_max <- un[seq(2, ln, by=9)] * 10000
      seg_int <- un[seq(3, ln, by=9)] * 100  
      gup_rsq <- un[seq(4, ln, by=9)] * 10000
      gup_maxgap <- un[seq(6, ln, by=9)] 
      gdown_rsq <- un[seq(7, ln, by=9)] * 10000
      gdown_maxgap <- un[seq(9, ln, by=9)]
      
      
      ##
      theOrd <- order(seg_amp,decreasing=T)   
      
      # Filter for bad EVI layers
      if(seg_max[theOrd[1]] > 10000 | seg_max[theOrd[1]] < 0 | seg_amp[theOrd[1]] > 10000 | seg_amp[theOrd[1]] < 0){
        outAll <- c(outAll,rep(NA,58));next}
      
      numRecords <- length(seg_amp)  #how many cycles were recorded
      naCheck <- is.na(seg_amp)
      numCyc <- sum(naCheck == 0)  #how many cycles have good data (seg metrics has valid observations)
      
      
      # QA
      if(length(full_segs)==1){
        qual_1 <- GetQAs(gup_rsq, gdown_rsq, gup_maxgap, gdown_maxgap, theOrd, qa_pars)
      }else{
        qual_1 <- GetQAs(gup_rsq, gdown_rsq, gup_maxgap, gdown_maxgap, theOrd, qa_pars)[[1]][1]
        qual_2 <- GetQAs(gup_rsq, gdown_rsq, gup_maxgap, gdown_maxgap, theOrd, qa_pars)[[2]][1]
      } 
      
      
      ################################################
      if(numCyc == 0){outAll <- c(outAll,annualMetrics(viSub,dateSub,smoothed_bns,pred_dates,phenYrs[y],pheno_pars,vi_dorm),rep(NA,18));next}
      
      if(numRecords == 1) {
        out <- c(annualMetrics(viSub,dateSub,smoothed_bns,pred_dates,phenYrs[y],pheno_pars,vi_dorm),
                 1,phen,qual_1,c(rep(NA,7),4),numObs)
      }else{
        phen1 <- phen[seq(theOrd[1], length(phen), by = numRecords)]
        phen2 <- phen[seq(theOrd[2], length(phen), by = numRecords)]
        if(naCheck[theOrd[2]]){
          out <- c(annualMetrics(viSub,dateSub,smoothed_bns,pred_dates,phenYrs[y],pheno_pars,vi_dorm),
                   numCyc,phen1,qual_1,c(rep(NA,7),4),numObs)  
        }else{
          out <- c(annualMetrics(viSub,dateSub,smoothed_bns,pred_dates,phenYrs[y],pheno_pars,vi_dorm),
                   numCyc,phen1,qual_1,phen2,qual_2,numObs)  
        }
      }
      
    },silent=TRUE) #End of the try block
    if(inherits(log, "try-error")){outAll <- c(outAll,rep(NA,58))
    }else{outAll <- c(outAll,out);remove(out)}
    
  }
  
  return(outAll)
}




#---------------------------------------------------------------------
# Get site name, image directory and coordinate
# Minkyu Moon 
#---------------------------------------------------------------------
GetSiteInfo <- function(numSite, geojsonDir, params){
  
    gjList <- list.files(path=paste0(geojsonDir),pattern=glob2rx('*.geojson'))
    gjListFull <- list.files(path=paste0(geojsonDir),pattern=glob2rx('*.geojson'),full.names=T)
    
    gparams <- FROM_GeoJson(gjListFull[numSite])
    
    strSite <- paste(strsplit(strsplit(gjList[numSite],'[.]')[[1]][1],' ')[[1]],collapse='_')
    
    imgDir <- dir(path=params$setup$dataDir,pattern=glob2rx(paste0('*',strSite,'*')),full.names=T)
    
  
  temp <- unlist(strsplit(imgDir,'/'))
  strSite <- temp[length(temp)]
  
  # Longitude and Latitude
  cLong <- (min(gparams$features[[1]]$geometry$coordinates[,1])+max(gparams$features[[1]]$geometry$coordinates[,1]))/2
  cLat  <- (min(gparams$features[[1]]$geometry$coordinates[,2])+max(gparams$features[[1]]$geometry$coordinates[,2]))/2
  
  return(list(imgDir=imgDir,strSite=strSite,cLong=cLong,cLat=cLat))
}



#---------------------------------------------------------------------
# Get site shape file
# Minkyu Moon 
#---------------------------------------------------------------------
GetSiteShp <- function(fileSR, cLong, cLat){
  
  # Shape file for 10 by 10 km window
  geog_crs = CRS("+proj=longlat +datum=WGS84")
  r <- raster(fileSR[1])
  utm_crs <- crs(r) 
  # utm_crs = raster(fileSR[1])@crs
  site <- data.frame(1,cLong,cLat)
  colnames(site) <- c('id','lon','lat')
  xy   <- site[,c(2,3)]
  bb   <- SpatialPointsDataFrame(coords=xy,data=site,proj4string=geog_crs)
  bb   <- spTransform(bb,utm_crs)
  
  x1 <- bb@coords[1] - 3000; x2 <- bb@coords[1] + 3000
  y1 <- bb@coords[2] - 3000; y2 <- bb@coords[2] + 3000
  xCoor <- c(x1,x2,x2,x1); yCoor <- c(y1,y1,y2,y2); xym <- cbind(xCoor,yCoor)
  p   <- Polygon(xym); ps  <- Polygons(list(p),1); sps <- SpatialPolygons(list(ps))
  proj4string(sps) <- utm_crs; data <- data.frame(f=99.9)
  
  siteWin <- SpatialPolygonsDataFrame(sps,data)
  
  return(siteWin)
}


#---------------------------------------------------------------------
# Get base image
# Minkyu Moon 
#---------------------------------------------------------------------
GetBaseImg <- function(fileSR, siteWin, outDir, save=TRUE){
  
  # Base Image
  for(i in 1:1000){
    log <- try({img1 <- raster(fileSR[i])},silen=T)
    if(inherits(log,'try-error')) next
    
    try(temp <- intersect(img1,siteWin),silent=T)
    if(temp@extent[1]>0) break
  }
  img1 <- crop(img1,siteWin)
  
  numImg <- 50
  imgBase <- vector('list',numImg)
  set.seed(456123)
  sam <- sample(1:length(fileSR),numImg)
  for(i in 1:numImg){
    log <- try({temp <- raster(fileSR[sam[i]])},silent=T)
    if(inherits(log,'try-error')){
      imgBase[[i]] <- img1
    }else{
      imgBase[[i]] <- raster(fileSR[sam[i]])
    }
  }
  
  for(i in 1:numImg){
    log <- try(compareRaster(imgBase[[i]],img1,extent=F,rowcol=F),silent=T)
    if(inherits(log,'try-error')){
      imgBase[[i]] <- projectRaster(imgBase[[i]],img1)    
    }
    log <- try(imgBase[[i]] <- crop(imgBase[[i]],siteWin),silent=T)
    if(inherits(log,'try-error')){
      imgBase[[i]] <- img1
    }
  }
  imgBase$fun <- mean; 
  imgBase$na.rm <- T
  imgBase <- do.call(mosaic,imgBase)  
  values(imgBase) <- NA
  
  if(save==TRUE){
    # Save Base Image
    writeRaster(imgBase, filename=paste0(outDir,'/base_image.tif'), format="GTiff", overwrite=TRUE)  
  }
  
  
  
  
  return(imgBase)
}







