PELT.trendARpJOIN=function(data,p=p,pen=0,minseglen=1,verbose=FALSE){
  # data is a n length vector
  # pen is a positive value for the penalty
  trendARpsegfit=function(data,start,previousbeta,p){
    # assumes that data is the data for the segment
    n=length(data)
    t=start:(start+n-1)
    filtered=data-previousbeta+previousbeta*(t-start)/n
    X=(t-start)/n
    trendfit=lm(filtered~-1+X)
    arfit=arima(resid(trendfit),order=c(p,0,0),include.mean=FALSE, method="ML")
    #loglik=n*log(arfit$sigma2)-log(1-arfit$coef^2)+
    #  (1-arfit$coef^2)*(resid(trendfit)[1])^2/arfit$sigma2+ # first obs
    #  (1/arfit$sigma2)*sum((resid(trendfit)[-1]-arfit$coef*resid(trendfit)[-n])^2) # remaining obs
    return(c(-2*logLik(arfit),coef(trendfit)[1])) # 1 is the trend estimate as no intercept
  }
  
  FIRSTtrendARpsegfit=function(data,start,p){
    # assumes that data is the data for the segment
    n=length(data)
    t=start:(start+n-1)
    X=(t-start)/n
    trendfit=lm(data~X)
    arfit=arima(resid(trendfit),order=c(p,0,0),include.mean=FALSE, method="ML")
    return(c(-2*logLik(arfit),coef(trendfit)[2]+coef(trendfit)[1])) # 2 is the trend estimate
  }
  
  n=length(data)
  
  lastchangecpts=rep(0,n+1)
  lastchangelike=matrix(NA,ncol=2,nrow=n+1)
  lastchangelike[1,]=c(-pen,0)
  checklist=NULL
  for(i in minseglen:(2*minseglen-1)){
    lastchangelike[i+1,]=FIRSTtrendARpsegfit(data[1:i], start=0,p=p)
    lastchangecpts[i+1]=0
  }
  checklist=c(0,minseglen)
  for(tstar in (2*minseglen):n){
    tmplike=unlist(lapply(checklist,FUN=function(tmpt){
      if(tmpt==0){
        return(FIRSTtrendARpsegfit(data[(tmpt+1):tstar], start=tmpt,p=p)[1])
      }
      return(lastchangelike[tmpt+1,1]+trendARpsegfit(data[(tmpt+1):tstar], start=tmpt,previousbeta=lastchangelike[tmpt+1,2],p=p)[1]+pen)
    }))
    lastchangelike[tstar+1,1]=min(tmplike,na.rm=TRUE)
    lastchangecpts[tstar+1]=checklist[which.min(tmplike)[1]]
    lastchangelike[tstar+1,2]=ifelse(checklist[which.min(tmplike)[1]]==0,FIRSTtrendARpsegfit(data[1:tstar], start=0,p=p)[2],
                                     trendARpsegfit(data[(lastchangecpts[tstar+1]+1):tstar], start=lastchangecpts[tstar+1],previousbeta=lastchangelike[lastchangecpts[tstar+1]+1,2],p=p)[2])
    checklist=checklist[tmplike<=(lastchangelike[tstar+1,1]+pen)]
    checklist=c(checklist,tstar-minseglen+1)
    if(verbose){if(tstar%%10==0){print(paste("Finished",tstar))}}
  }
  fcpt=NULL
  last=n
  while(last!=0){
    fcpt=c(fcpt,lastchangecpts[last+1])
    last=lastchangecpts[last+1]
  }
  return(sort(fcpt)[-1])
}

# ---------------- FIT with continuity (match FIXARp logic) ----------------
fit.trendARpJOIN <- function(data, cpts, p, dates=NULL, plot=TRUE, fit=TRUE,
                             add.ar=FALSE, title="Data", pred=FALSE){
  library(forecast)

  ## joined segment: no intercept; continuity enforced by previous level
  trendARpsegfit <- function(data, start, prev_level, p, pred=FALSE){
    n <- length(data)
    t <- start:(start+n-1)
    X <- (t - start+1) / n

    filtered  <- data - prev_level + prev_level * X
    trendfit  <- lm(filtered ~ -1 + X)
    arfit     <- arima(resid(trendfit), order=c(p,0,0),
                       include.mean=FALSE, method="ML")

    if(!pred){
      beta_per_year <- as.numeric(coef(trendfit)["X"])/n
      fitted_trend  <- as.numeric(fitted(trendfit)) + prev_level - prev_level*X
      return(list(
        coef      = c(beta_per_year, coef(arfit)),   # length = 1+p
        arfit     = as.numeric(fitted(arfit)),
        trendfit  = fitted_trend
      ))
    }else{
      fc <- predict(arfit)
      fc$pred <- fc$pred + coef(trendfit)
      return(fc)
    }
  }

  ## first segment: free intercept; also return a "carry level" for the join
  FIRSTtrendARpsegfit <- function(data, start, p){
    n <- length(data)
    t <- start:(start+n-1)
    X <- (t - start+1) / n

    trendfit <- lm(data ~ X)
    arfit    <- arima(resid(trendfit), order=c(p,0,0),
                      include.mean=FALSE, method="ML")

    intercept      <- as.numeric(coef(trendfit)[1])
    slope_segment  <- as.numeric(coef(trendfit)["X"])      # per-X-unit slope
    carry_level    <- intercept + slope_segment            # level to pass at knot
    beta_per_year  <- slope_segment / n

    return(list(
      # coef columns: [1]=carry_level (for continuity), [2]=beta_per_year, [3..]=AR
      coef      = c(intercept, beta_per_year, coef(arfit)),  # length = 2+p
      arfit     = as.numeric(fitted(arfit)),
      trendfit  = as.numeric(fitted(trendfit))
    ))
  }

  # ---- set up and fit all segments ----
  cpts <- c(0, cpts, length(data))
  nseg <- length(cpts) - 1

  # coeffs: (CarryLevel, Beta_per_year, AR1..ARp)
  coeffs <- matrix(NA_real_, nrow=nseg, ncol=2+p)
  colnames(coeffs) <- c("CarryLevel", "Beta", paste0("AR", seq_len(p)))

  segments <- vector("list", nseg)

  # first segment
  segments[[1]] <- FIRSTtrendARpsegfit(data[(cpts[1]+1):cpts[2]], start=cpts[1], p=p)
  coeffs[1, ]   <- segments[[1]]$coef

  # pass the actual last fitted level to the next segment
  prev_level <- tail(segments[[1]]$trendfit, 1)

  if(nseg >= 2){
    for(i in 2:nseg){
      seg <- trendARpsegfit(data[(cpts[i]+1):cpts[i+1]],
                            start=cpts[i],
                            prev_level=prev_level,
                            p=p,
                            pred=FALSE)
      segments[[i]] <- seg
      # pad to (2+p): NA carry level (not used for joined), then beta + ARs
      coeffs[i, ] <- c(NA_real_, seg$coef)
      prev_level  <- tail(seg$trendfit, 1)
    }
  }

  # ---- assemble fitted trend ----
  fit_vec <- unlist(lapply(segments, function(x) x$trendfit))
  if(add.ar){
    fit_vec <- fit_vec + unlist(lapply(segments, function(x) x$arfit))
  }

  # ---- optional save/plot ----
  if(fit){
    if(length(fit_vec) > 55){
      save(dates, fit_vec, data, file=paste0(title,"_TrendAR",p,"JOIN_fit",".RData"))
    } else {
      save(dates, fit_vec, data, file=paste0(title,"_TrendAR",p,"JOINshort_fit",".RData"))
    }
  }

  if(plot){
    if(!is.null(dates)){
      if(length(dates)!=length(data)) stop("Length of dates and data must be the same")
      plot(dates, data, type='l', main=title, xlab="Year", ylab="Anomaly (°C)")
      lines(dates, fit_vec, col='blue')
    } else {
      ts.plot(data, main=title)
      lines(seq_along(data), fit_vec, col='blue')
    }
  }

  return(list(dates=dates, fit=fit_vec, coeffs=coeffs))
}


estimate_segment_sd <- function(data, fitobj, cpts){
  # cpts are changepoints in indices (no 0, no n yet)
  cpts <- c(0, cpts, length(data))
  seg_sd <- numeric(length(cpts)-1)

  for(i in 1:(length(cpts)-1)){
    idx <- (cpts[i]+1):cpts[i+1]
    resid <- data[idx] - fitobj$fit[idx]

    # AR(1) filtered residuals for this segment
    phi <- fitobj$coeffs[i, "AR1"]
    
    # one-step filtered noise
    eps <- resid[-1] - phi*head(resid,-1)
    seg_sd[i] <- sd(eps, na.rm=TRUE)
  }
  return(seg_sd)
}


simulate_trendARpJOIN <- function(y, fitobj, cpts){
  coeffs <- fitobj$coeffs
  yfit   <- fitobj$fit
  n      <- length(yfit)

  # ---- estimate AR-noise SD per segment ----  
  sds <- estimate_segment_sd(y, fitobj, cpts)

  sim <- numeric(n)

  ### ======================
  ### Segment 1 (with intercept)
  ### ======================
  intercept <- coeffs[1,"CarryLevel"]
  beta1 <- coeffs[1,"Beta"]
  phi1  <- coeffs[1,"AR1"]
  sd1   <- sds[1]

  cpts_full <- c(0, cpts, length(data))
  end1 <- cpts_full[2]   # first segment ends at this index
  eps1 <- as.numeric(arima.sim(list(ar=phi1), n=end1, sd=sd1))
  sim[1:end1] <- intercept + beta1*(1:end1) + eps1
  end_seg_1 = intercept + beta1*end1
  ### ======================
  ### Joined continuation segments (continuous)
  ### ======================
  if(nrow(coeffs) > 1){
    end_seg_i = end_seg_1
    for(seg in 2:nrow(coeffs)){
      beta <- coeffs[seg,"Beta"]
      phi  <- coeffs[seg,"AR1"]
      sd   <- sds[seg]

      start <- ifelse(seg==2, end1+1, which(is.na(coeffs[,1]))[seg-2] + 1)
      end   <- ifelse(seg==nrow(coeffs), n, which(is.na(coeffs[,1]))[seg-1])
      len   <- end - start + 1

      eps <- as.numeric(arima.sim(list(ar=phi), n=len, sd=sd))

      for(t in start:end){
        sim[start:end] <- end_seg_i + beta*(1:len) + eps
      }
      end_seg_i = end_seg_i + beta*len
    }
  }
  return(sim)
}

PELT.trendFIXARpJOIN=function(data,pen=0,arp=c(0),minseglen=1,verbose=FALSE){
  # data is a n length vector
  # pen is a positive value for the penalty
  # arp is the vector of the fixed AR(p) coefficient
  trendARpsegfit=function(data,start,arp,previousbeta){
    # assumes that data is the data for the segment
    n=length(data)
    t=start:(start+n-1)
    filtered=data-previousbeta+previousbeta*(t-start)/n
    X=(t-start)/n
    trendfit=lm(filtered~-1+X)
    arfit=arima(resid(trendfit),order=c(length(arp),0,0),fixed=arp,include.mean=FALSE, method="ML")
    return(c(-2*logLik(arfit),coef(trendfit)[1])) # 1 is the trend estimate as no intercept
  }
  
  FIRSTtrendARpsegfit=function(data,start,arp){
    # assumes that data is the data for the segment
    n=length(data)
    t=start:(start+n-1)
    X=(t-start)/n
    trendfit=lm(data~X)
    arfit=arima(resid(trendfit),order=c(length(arp),0,0),fixed=arp,include.mean=FALSE, method="ML")
    return(c(-2*logLik(arfit),coef(trendfit)[2]+coef(trendfit)[1])) # 2 is the trend estimate
  }
  
  n=length(data)
  
  lastchangecpts=rep(0,n+1)
  lastchangelike=matrix(NA,ncol=2,nrow=n+1)
  lastchangelike[1,]=c(-pen,0)
  checklist=NULL
  for(i in minseglen:(2*minseglen-1)){
    lastchangelike[i+1,]=FIRSTtrendARpsegfit(data[1:i], start=0,arp=arp)
    lastchangecpts[i+1]=0
  }
  checklist=c(0,minseglen)
  for(tstar in (2*minseglen):n){
    tmplike=unlist(lapply(checklist,FUN=function(tmpt){
      if(tmpt==0){
        return(FIRSTtrendARpsegfit(data[(tmpt+1):tstar], start=tmpt,arp=arp)[1])
      }
      return(lastchangelike[tmpt+1,1]+trendARpsegfit(data[(tmpt+1):tstar], start=tmpt,arp=arp,previousbeta=lastchangelike[tmpt+1,2])[1]+pen)
    }))
    lastchangelike[tstar+1,1]=min(tmplike,na.rm=TRUE)
    lastchangecpts[tstar+1]=checklist[which.min(tmplike)[1]]
    lastchangelike[tstar+1,2]=ifelse(checklist[which.min(tmplike)[1]]==0,FIRSTtrendARpsegfit(data[1:tstar], start=0,arp=arp)[2],
                                     trendARpsegfit(data[(lastchangecpts[tstar+1]+1):tstar], start=lastchangecpts[tstar+1],arp=arp,previousbeta=lastchangelike[lastchangecpts[tstar+1]+1,2])[2])
    checklist=checklist[tmplike<=(lastchangelike[tstar+1,1]+pen)]
    checklist=c(checklist,tstar-minseglen+1)
    if(verbose){if(tstar%%10==0){print(paste("Finished",tstar))}}
  }
  fcpt=NULL
  last=n
  while(last!=0){
    fcpt=c(fcpt,lastchangecpts[last+1])
    last=lastchangecpts[last+1]
  }
  return(sort(fcpt)[-1])
}



# Now for a given output from the above (set of changepoints) we want to get
# the final fit for each segment and a plot
fit.trendFIXARpJOIN=function(data,cpts,arp,dates=NULL,plot=T,fit=T,add.ar=F,title="Data"){
  library(forecast)
  trendFIXARpsegfit=function(data,start,arp,previousbeta){
    # assumes that data is the data for the segment
    n=length(data)
    t=start:(start+n-1)
    filtered=data-previousbeta+previousbeta*(t-start)/n
    X=(t-start)/n
    trendfit=lm(filtered~-1+X)
    arfit=arima(resid(trendfit),order=c(length(arp),0,0),fixed=arp,include.mean=FALSE, method="ML")
    return(list(coef=c(coef(trendfit),coef(arfit)),arfit=fitted(arfit),
                trendfit=fitted(trendfit)+previousbeta-previousbeta*(t-start)/n))
  }
  
  FIRSTtrendFIXARpsegfit=function(data,start,arp){
    # assumes that data is the data for the segment
    n=length(data)
    t=start:(start+n-1)
    X=(t-start)/n
    trendfit=lm(data~X)
    arfit=arima(resid(trendfit),order=c(length(arp),0,0),fixed=arp,include.mean=FALSE, method="ML")
    return(list(coef=c(coef(trendfit)[2]+coef(trendfit)[1],coef(arfit)),arfit=fitted(arfit),trendfit=fitted(trendfit)))
  }
  
  
  cpts=c(0,cpts,length(data))
  segments=list()
  coeffs=matrix(NA,nrow=length(cpts)-1,ncol=length(arp)+1)
  
  segments[[1]]=FIRSTtrendFIXARpsegfit(data[(cpts[1]+1):cpts[2]],start=cpts[1],arp=arp)
  coeffs[1,]=segments[[1]]$coef
  if(length(cpts)>2){
    for(i in 2:(length(cpts)-1)){
      segments[[i]]=trendFIXARpsegfit(data[(cpts[i]+1):cpts[i+1]],start=cpts[i],arp=arp,previousbeta=coeffs[i-1,1])
      coeffs[i,]=segments[[i]]$coef
    }
  }
  
  if(plot|fit){
    fit=unlist(lapply(segments,FUN=function(x){x$trendfit}))
    print(paste(title,"Resid ar:",coef(arima(data-fit,order=c(length(arp),0,0),include.mean=FALSE,method="ML"))))
    
    if(add.ar){
      fit=fit+unlist(lapply(segments,FUN=function(x){x$arfit}))
    }
    return(list(dates=dates,fit=fit,coeffs=coeffs))
    
    # if (length(fit) > 55){#full time series
    #   
    #   myfile = paste0(title, "_TrendFIXAR",length(arp),"JOIN_fit", ".RData")
    #   save(dates,fit,data,file=myfile)
    # }else{#short
    #   myfile = paste0(title, "_TrendFIXAR",length(arp),"JOINshort_fit", ".RData")
    #   save(dates,fit,data,file=myfile)
    # }
  }
  if(plot){
    if(!is.null(dates)){
      if(length(dates)!=length(data)){stop("Length of dates and data must be the same")}
      plot(dates,data,type='l',main=title,xlab="Year",ylab="Anomaly (°C)")
      lines(dates,fit,col='blue')
      #abline(v=dates[cpts[-c(1,length(cpts))]],col='blue')
    }
    else{
      ts.plot(data,main=title)
      lines(1:length(data),fit,col='blue')
      #abline(v=cpts[-c(1,length(cpts))],col='blue')
    }
  }
  
  return(coeffs)
}



acfpacfcpts=function(path,acf=T,pacf=T, title='Data'){
  load(path)
  resid=data-fit
  if(acf){print(acf(resid,main=paste(title, "ACF")))}
  if(pacf){print(pacf(resid,main=paste(title, "PACF")))}
}
