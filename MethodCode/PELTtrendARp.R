
PELT.trendARp=function(data,p=1, pen=0,minseglen=1,verbose=FALSE){
  # data is a n length vector
  # pen is a positive value for the penalty
  trendARpsegfit=function(data,start,p){
    # assumes that data is the data for the segment
    n=length(data)
    X=start:(start+n-1)
    trendfit=lm(data~X)
    arfit=arima(resid(trendfit),order=c(p,0,0),include.mean=FALSE, method="ML")
    return(-2*logLik(arfit))
  }
  
  n=length(data)
  
  lastchangecpts=rep(0,n+1)
  lastchangelike=-pen
  checklist=NULL
  for(i in minseglen:(2*minseglen-1)){
    lastchangelike[i+1]=trendARpsegfit(data[1:i], start=0,p=p)
    lastchangecpts[i+1]=0
  }
  checklist=c(0,minseglen)
  for(tstar in (2*minseglen):n){
    tmplike=unlist(lapply(checklist,FUN=function(tmpt){return(lastchangelike[tmpt+1]+trendARpsegfit(data[(tmpt+1):tstar], start=tmpt,p=p)+pen)}))
    lastchangelike[tstar+1]=min(tmplike,na.rm=TRUE)
    lastchangecpts[tstar+1]=checklist[which.min(tmplike)[1]]
    checklist=checklist[tmplike<=(lastchangelike[tstar+1]+pen)]
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




fit.trendARp = function(data, cpts, p=1, dates=NULL, plot=T, fit=T, add.ar=F, title="Data") {
  library(forecast)
  
  trendARpsegfit = function(data, start, p) {
    n = length(data)
    X = start:(start + n - 1)
    trendfit = lm(data ~ X)
    arfit = arima(resid(trendfit), order=c(p,0,0), include.mean=FALSE, method="ML")
    return(list(
      coef     = c(coef(trendfit), coef(arfit)),
      arfit    = fitted(arfit),       # AR fitted values (on residuals)
      trendfit = fitted(trendfit)     # trend fitted values
    ))
  }
  
  cpts = c(0, cpts, length(data))
  segments = list()
  coeffs = matrix(NA, nrow=length(cpts)-1, ncol=p+2)
  
  for (i in 1:(length(cpts)-1)) {
    segments[[i]] = trendARpsegfit(data[(cpts[i]+1):cpts[i+1]], start=cpts[i], p=p)
    coeffs[i,] = segments[[i]]$coef
  }
  
  # ---- Combine fitted values across segments ----
  fit_trend = unlist(lapply(segments, function(x) x$trendfit))
  fit_ar    = unlist(lapply(segments, function(x) x$arfit))    # AR residual fit
  fit_total = fit_trend + fit_ar                                # trend + AR combined

  # ---- Plot ----
  if (plot) {
    if (!is.null(dates)) {
      if (length(dates) != length(data)) stop("Length of dates and data must be the same")
      plot(dates, data, type="l", main=title, xlab="Year", ylab="Anomaly (°C)")
      lines(dates, fit_trend, col="blue",  lwd=2)
      lines(dates, fit_total, col="green", lwd=2, lty=2)   # optional: total fit
      abline(v=dates[cpts[-c(1, length(cpts))]], col="red", lty=2)
    } else {
      plot(data, type="l", main=title)
      lines(fit_trend, col="blue",  lwd=2)
      lines(fit_total, col="green", lwd=2, lty=2)
      abline(v=cpts[-c(1, length(cpts))], col="red", lty=2)
    }
  }
  
  # ---- Single clean return ----
  return(list(
    dates     = dates,
    fit_trend = fit_trend,   # trend component only
    fit_ar    = fit_ar,      # AR component (on residuals)
    fit_total = fit_total,   # trend + AR
    coeffs    = coeffs
  ))
}

# Now for a given output from the above (set of changepoints) we want to get
# the final fit for each segment and a plot
fit.trendARpold=function(data,cpts,p=1, dates=NULL,plot=T,fit=T,add.ar=F,title="Data"){
  library(forecast)
  trendARpsegfit=function(data,start,p){
    # assumes that data is the data for the segment
    n=length(data)
    X=start:(start+n-1)
    trendfit=lm(data~X)
    arfit=arima(resid(trendfit),order=c(p,0,0),include.mean=FALSE, method="ML")
    return(list(coef=c(coef(trendfit),coef(arfit)),arfit=fitted(arfit),trendfit=fitted(trendfit)))
  }
  
  cpts=c(0,cpts,length(data))
  segments=list()
  coeffs=matrix(NA,nrow=length(cpts)-1,ncol=p+2)
  arfit=list()

  for(i in 1:(length(cpts)-1)){
    segments[[i]]=trendARpsegfit(data[(cpts[i]+1):cpts[i+1]],start=cpts[i],p=p)
    coeffs[i,]=segments[[i]]$coef
    #arfit[i,]=segments[[i]]$arfit
  }
# ---- Combine fitted trend only ----
fit_trend <- unlist(lapply(segments, function(x) x$trendfit))

  # ---- Plot the trend ----
  if (plot) {
    if (!is.null(dates)) {
      if (length(dates) != length(data)) stop("Length of dates and data must be the same")
      plot(dates, data, type="l", main=title, xlab="Year", ylab="Anomaly (°C)")
      lines(dates, fit_trend, col="blue", lwd=2)
      abline(v = dates[cpts[-c(1, length(cpts))]], col="red", lty=2)  # mark breakpoints
    } else {
      plot(data, type="l", main=title)
      lines(fit_trend, col="blue", lwd=2)
      abline(v = cpts[-c(1, length(cpts))], col="red", lty=2)
    }
  }

# ---- Return ----
return(list(dates = dates,
            fit_trend = fit_trend,
            coeffs = coeffs))

  # ✅ Now return AFTER plotting
  if(fit){
    return(list(dates=dates, fit=fit_vals, coeffs=coeffs))
  } else {
    return(coeffs)
  }
  
  return(coeffs)
}





