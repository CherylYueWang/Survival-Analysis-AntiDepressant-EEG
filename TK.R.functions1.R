
g.emphaz<-function(data = NULL,type="hhat",legend = NULL, main = NULL)
{
  ## Purpose:  1. draw empirical hazards
  ##           2. prints out actual hazard values in the order hitilde, hihat 
  ## ------------------------------------------------------------------------
  ## Arguments:
  ##   data     a Surv object or a list of Surv objects 
  ##   type     what should be drawn? "ht" for hitilde or "hhat" for hihat
  ##   main     main title
  ##   legend   
  ## ------------------------------------------------------------------------
  ## Author: Werner Stahel, Date:  2 Aug 2004, 17:21 / WSt Aug 04
  ## Second Author:  Mara Tableman, Date: 30 October 2004
  f.emphaz<-function(w)
  { ## creates a list of data.frames with "time","ht" for hitilde,"hhat" for hihat
    ##    for two different treatment groups.
    f.f <- function(s)
    {
      sf <- summary(survfit(Surv(s[, 1], s[, 2])~1, type = "kaplan-meier"))
      time <- sf$time
      hitilde <- round(sf$n.event/sf$n.risk,digits=3)
      hihat <- round(hitilde/c(diff(time), NA),digits=3)
      hihat[length(time)] <- hihat[(length(time) - 1)]
      data.frame(time = sf$time, ht = hitilde, hhat = hihat)
    }
    if(length(dim(data)) == 2)
      f.f(w)
    else lapply(w, f.f)
  }
  if(!is.null(data)) emphaz <- f.emphaz(data)
  print(emphaz)
  if(is.data.frame(emphaz))
    emphaz <- list(emphaz)
  xrg <- range(sapply(emphaz, function(x)range(x[, "time"])))
  yrg <- range(sapply(emphaz, function(x)range(x[, "ht"])))
  if (type =="hhat"){ 
    yrg <- range(sapply(emphaz, function(x)range(x[,  "hhat" ])))}
  
  plot(xrg, yrg, type = "n", xlab = "observed failure times", ylab = "hazard at time i")
  nds <- length(emphaz)
  for(li in 1:nds) {
    ldat <- emphaz[[li]]
    lines(ldat[, "time"], ldat[, type],type="s", lty = li)
    points(ldat[, "time"], ldat[, type], pch = 16, cex = 0.5)
  }
  if(!is.null(legend)) {
    usr <- par("usr")
    legend(0.9 * usr[1] + 0.1 * usr[2], 0.05 * usr[3] + 0.95 *
             usr[4], rep(legend, length = nds), lty = 1:2)
  }
  title(main)
  box()
}
##=============================================================
hazard.km <- function(data){
  ## Author: Mara Tableman  Date: 20 November 2002 
  ## Purpose:  To compute the two types of empirical hazards,
  ##           and the cumulative hazards along their s.e.'s
  ## Arguments:  data is survfit object
  time <- summary(data)$time
  ni <- summary(data)$n.risk
  di <- summary(data)$n.event
  surv <- summary(data)$surv
  stderr <- summary(data)$std.err
  hitilde <- di/ni
  tau <- diff(time, lag = 1)	#length of interval to right of ti
  tau[length(tau) + 1] <- NA
  hihat <- hitilde/tau
  Hhat <-  -log(surv)
  Htilde <- cumsum(hitilde)
  sqri <- di/(ni^2)
  se.Hhat <- stderr/surv
  se.Htilde <- sqrt(cumsum(sqri))
  hazardtable <- round(data.frame(time, ni, di, hihat, hitilde, Hhat,se.Hhat, Htilde, se.Htilde),4)
  print(hazardtable)
  on.exit()
  "hazard.km:done"
}
##=============================================================
plot.logt.x<-function(time,status,x,xlab="x",ylab="ordered log data")
{  ##Purpose:  To produce a plot of log(t) against each predictor variable x.
  ##     A least squares line is drawn through the point cloud.
  ##     Note that only the uncensored times and their corresponding
  ##     values are plotted.
  ##-----------------------------------------------------------------
  ##Arguments:  time is the survival time varibale
  ##            status is the variable use for censored or not
  ##            x is one of the j=1,...,m predicotr variables
  ##-----------------------------------------------------------------
  ##Author: Mara Tableman, Date: 3.July 2004
  time.u<-time[status==1]
  x.u<-x[status==1]
  plot(x.u,log(time.u),type="p",xlab=xlab,ylab=ylab)
  regr<-lsfit(x.u,log(time.u))
  abline(regr$coef)
}
##===================================================================
qq.gamma <- function(data, time, status, xlab = "gamma quantiles \n based on MLE's", ylab = "ordered data", pch = NULL,lty = NULL, col = NULL)
{
  ## Purpose: qqplot for gamma distribution for one sample
  ##-------------------------------------------------------------------
  ## Arguments: ##  data  A Surv object, or a list of Surv object 
  ##         time, and status, are each a list or vector
  ##                 
  ##-------------------------------------------------------------------
  ## Author: Mara Tableman, Date: 2 Janury 2003, 01:05
  on.exit(browser())
  t.c <- class(data)
  if((!is.null(t.c)) && t.c == "Surv")
    data <- list(data)
  t.k <- length(data)
  #t.mt <- class(time)
  #if((!is.null(t.mt)) & (class(t.mt)!="list"))
  #	time <- list(time)
  #t.ms <- class(status)
  #if((!is.null(t.ms)) & (class(t.ms)!="list"))
  #status <- list(status)
  if(is.null(pch))
    pch <- 1:t.k
  if(is.null(lty))
    lty <- 1:t.k
  if(is.null(col))
    col <- 1 + 1:t.k
  t.sf <- lapply(data, function(s)
  {
    t.s <- summary(survfit(s, type = "kaplan-meier"))
    t.ss <- t.s$surv
    t.s$surv <- (t.ss + c(1, t.ss[ - length(t.ss)]))/2
    t.s
  }
  )
  t.rgs <- range(sapply(t.sf, function(s)
    range(s$surv)))
  t.rgt <- range(sapply(t.sf, function(s)
    range(s$time)))
  for(t.g in 1:t.k) {
    sr <- vector(mode = "numeric", length = 2)
    status1 <- status[[t.g]]
    time1 <- time[[t.g]]
    minusloglik.gamma <- function(sr, time1, status1)
    {
      - log(prod(dgamma(time1, sr[1], sr[2])^(status1) * (1 - pgamma(time1, sr[1], sr[2]))^(1 - status1)))
    }
    mlegamma <- nlminb(start = c(0.05, 0.05), obj = minusloglik.gamma, lower = c(0.001, 0.0001), upper = c(Inf, Inf), time1 = time1, status1 = status1)
    mles <- mlegamma$par
    names(mles) <- c("shape", "rate")
    print(mles)
    par(new = T)
    xmax <- max(qgamma(1 - t.rgs, mles[1], mles[2]))
    ymax <- max(t.rgt)
    plot(qgamma(1 - t.rgs, mles[1], mles[2]), t.rgt, type = "n", xlab = " ", ylab = "  ", xlim = c(0, xmax), ylim = c(0,ymax), axes = F)
    points(qgamma(1 - t.sf[[t.g]]$surv, mles[1], mles[2]), t.sf[[t.g]]$time, pch = pch[t.g], col = 1, lwd = 2)
    reg <- lsfit(qgamma(1 - t.sf[[t.g]]$surv, mles[1], mles[2]), t.sf[[t.g]]$time)
    abline(reg, lty = lty[t.g], col = 1, lwd = 2)
  }
  abline(h = 0, v = 0)
  abline(a=0,b=1,lwd=3)
  mtext(side = 1, "gamma quantiles \n based on MLE's", line=1.3,cex = 1.7)
  mtext(side = 2, "ordered data", cex = 1.8, line = -0.3)
  on.exit()
  "qq.gamma:done"
}

##==========================================================
qq.loglogisreg<-function(data, loglogis.fit, xlab = "standard logistic quantiles", ylab = "ordered log data", pch =NULL,lty= NULL, col = NULL)
{
  ## Purpose: qqplot for loglogistic distribution  
  ##          Fit with one covariate, several distinct levels.  Same slope,
  ##          but different intercepts.
  ##-------------------------------------------------------------------
  ## Arguments: ##  data  A Surv object, or a list of such objects
  ##            ##  loglogis.fit is a survreg object with dist="loglogistic"  
  ##-------------------------------------------------------------------
  ## Author: Mara Tableman, Date: 2 December, 2002
  on.exit(browser())
  t.c <- class(data)
  if((!is.null(t.c)) && t.c == "Surv")
    data <- list(data)
  t.k <- length(data)
  
  # if(is.null(pch))
  pch <- 1:t.k
  # if(is.null(lty))
  lty <- 1:t.k
  # if(is.null(col))
  col <- 1 + 1:t.k
  t.sf <- lapply(data, function(s)
  {
    t.s <- summary(survfit(s~1, type = "kaplan-meier"))
    t.ss <- t.s$surv
    t.s$surv <- (t.ss + c(1, t.ss[ - length(t.ss)]))/2
    t.s
  }
  )
  t.rgs <- range(sapply(t.sf, function(s)
    range(s$surv)))
  t.rgt <- range(sapply(t.sf, function(s)
    range(s$time)))
  plot(qlogis(1 - t.rgs), log(t.rgt), type = "n", ylim = c(min(log(t.rgt)), max(log(t.rgt))), xlab = xlab, 
       ylab = ylab)
  for(t.g in 1:t.k) {
    points(qlogis(1 - t.sf[[t.g]]$surv), log(t.sf[[t.g]]$time), pch = pch[t.g], col = 1)
    mutilde <- as.numeric(levels(factor(loglogis.fit$linear.predictors)))
    scale <- rep(loglogis.fit$scale, length(mutilde))
    abline(mutilde[t.g], scale[t.g], lty = lty[t.g], col = 1)
  }
  on.exit()
  "qq.loglogisreg:done"
}
##===========================================================
qq.loglogistic <- function(data, xlab = "standard logistic quantiles", ylab = 
                             "ordered log data", pch = NULL, lty = NULL, col = NULL)
{
  ## Purpose: qqplot for loglogistic distribution for one or several samples 
  ##-------------------------------------------------------------------
  ## Arguments: ##  data  A Surv object, or a list of such objects
  ##    
  ##--------------------------------------------------------------------
  ## Author: Mara Tableman, Date: 24 December 2002, 18:45
  on.exit(browser())
  t.c <- class(data)
  if((!is.null(t.c)) && t.c == "Surv")
    data <- list(data)
  t.k <- length(data)
  if(is.null(pch))
    pch <- 1:t.k
  if(is.null(lty))
    lty <- 1:t.k
  if(is.null(col))
    col <- 1 + 1:t.k
  t.sf <- lapply(data, function(s)
  {
    t.s <- summary(survfit(s~1, type = "kaplan-meier"))
    t.ss <- t.s$surv
    t.s$surv <- (t.ss + c(1, t.ss[ - length(t.ss)]))/2
    t.s
  }
  )
  t.rgs <- range(sapply(t.sf, function(s)
    range(s$surv)))
  t.rgt <- range(sapply(t.sf, function(s)
    range(s$time)))
  plot(qlogis(1 - t.rgs, 0, 1), log(t.rgt), type = "n", xlab = xlab, ylab = ylab)
  for(t.g in 1:t.k) {
    points(qlogis(1 - t.sf[[t.g]]$surv, 0, 1), log(t.sf[[t.g]]$time), pch = pch[t.g], col = 1)
    t.r <- survreg(data[[t.g]] ~ 1, dist = "loglogistic")
    abline(t.r$coef, t.r$scale, lty = lty[t.g], col = 1)
  }
  on.exit()
  "qq.loglogistic:done"
}

##==========================================================
qq.lognormreg<-function(data, lognorm.fit, xlab = "standard normal quantiles", ylab = "ordered log data", pch = NULL, lty = NULL, col = NULL)
{
  ## Purpose: qqplot for lognormal regression with one covariate which has a finite number 
  ##          of distinct values which can be coverted into levels. It draws lines with same 
  ##          slope sigma, but different intercepts, whose values are the MLE's
  ##-------------------------------------------------------------------
  ## Arguments: ##  data  A Surv object, or a list of such objects
  ##            ##  lognorm.fit is a survreg object with dist="lognormal"    
  ##-------------------------------------------------------------------
  ## Author: Mara Tableman, Date: 2 December, 2002
  on.exit(browser())
  t.c <- class(data)
  if((!is.null(t.c)) && t.c == "Surv")
    data <- list(data)
  t.k <- length(data)
  if(is.null(pch))
    pch <- 1:t.k
  if(is.null(lty))
    lty <- 1:t.k
  if(is.null(col))
    col <- 1 + 1:t.k
  t.sf <- lapply(data, function(s)
  {
    t.s <- summary(survfit(s~1, type = "kaplan-meier"))
    t.ss <- t.s$surv
    t.s$surv <- (t.ss + c(1, t.ss[ - length(t.ss)]))/2
    t.s
  }
  )
  t.rgs <- range(sapply(t.sf, function(s)
    range(s$surv)))
  t.rgt <- range(sapply(t.sf, function(s)
    range(s$time)))
  plot(qnorm(1 - t.rgs), log(t.rgt), type = "n", ylim = c(min(log(t.rgt)), max(log(t.rgt))), xlab = xlab, ylab
       = ylab)
  for(t.g in 1:t.k) {
    points(qnorm(1 - t.sf[[t.g]]$surv), log(t.sf[[t.g]]$time), pch = pch[t.g], col = 1)
    mutilde <- as.numeric(levels(factor(lognorm.fit$linear.predictors)))
    scale <- rep(lognorm.fit$scale, length(mutilde))
    abline(mutilde[t.g], scale[t.g], lty = lty[t.g], col = 1)
  }
  on.exit()
  "qq.lognormreg:done"
}

qq.reg.resid.r<-function(data,time,status,fit,quantile,xlab){
  ##=================================
  ## Purpose: For parametric regression models, this constructs a
  ##          qq-plot of ordered residuals e_i=(y_i-yhat_i)/sigmahat against
  ##          the log-parametric standard quantiles z_i of either the
  ##          "weibull", "lognormal", or "loglogistic" distribution.
  ##--------------------------------------------------------------------------
  ## NOTE:  This can also be used for fitting a single sample of survival
  ##        times to a parametric model. Since there are no covariates
  ##        remember to type survReg(Surv(...,...)~1,dist="...",data=...) 
  ##        in order to estimate the intercept mu.
  ##---------------------------------------------------------------------------
  ## Arguments:   data = data.frame
  ##              time = survival time variable name in data.frame
  ##              status = name of status variable in data.frame
  ##              fit = a survreg object
  ##              quantile = "qweibull" or "qnorm" or "qlogis"  
  ##              xlab = "type your label" E.g., "standard extreme value quantiles"
  ##------------------------------------------------------------------------------
  ## Author: Mara Tableman, Date: 2 February, 2005:  repaired on 24 August 2012
  ##==============================================================================
  temp<-data
  temp$time<-time
  temp$status<-status
  temp$ei<-(log(temp$time)-predict(fit,type="lp"))/fit$scale
  temp<-temp[order(temp$ei), ]
  con<-abs(min(temp$ei))+.00001
  temp$ei<-temp$ei+con 
  km.fit<-survfit(Surv(ei,status)~1,data=temp)
  temp$km.surv<-summary(km.fit,censor=T)$surv
  if (quantile == "qweibull") {
    zi<-as.numeric(qweibull(1-temp$km.surv,1,1))
    k<-nrow(temp)
    for (i in 1:k){
      if (zi[i]!=-Inf && zi[i]!=Inf ) zi[i]<-log(zi[i])
    }
    temp$zi<-zi
    for(i in 1:k){ if (temp$zi[i]==-Inf)
    {	surv.max.1<-max(temp$km.surv[temp$status==1])
    d<-1-surv.max.1
    surv.pu<-1-d/2	
    temp$zi[i]<-log(qweibull(1-surv.pu,1,1))}
    }
    for (i in 1:k){ if (temp$zi[i]==Inf)
    {  d<-min(temp$km.surv[temp$km.surv > 0])
    surv.pl<-d/2
    temp$zi[i]<-log(qweibull(1-surv.pl,1,1))}
    }
  }
  if (quantile == "qlogis") {
    zi<-as.numeric(qlogis(1-temp$km.surv,0,1))
    k<-nrow(temp)
    for (i in 1:k){
      if (zi[i]!=-Inf && zi[i]!=Inf ) zi[i]<-zi[i]
    }
    temp$zi<-zi
    for(i in 1:k){ if (temp$zi[i]==-Inf)
    {	surv.max.1<-max(temp$km.surv[temp$status==1])
    d<-1-surv.max.1
    surv.pu<-1-d/2	
    temp$zi[i]<-qlogis(1-surv.pu,0,1)}
    }
    for (i in 1:k){ if (temp$zi[i]==Inf)
    {  d<-min(temp$km.surv[temp$km.surv > 0])
    surv.pl<-d/2
    temp$zi[i]<-qlogis(1-surv.pl,0,1)}
    }
  }
  if (quantile == "qnorm"){ 
    zi<-as.numeric(qnorm(1-temp$km.surv,0,1))
    k<-nrow(temp)
    for (i in 1:k){
      if (zi[i]!=-Inf && zi[i]!=Inf ) zi[i]<-zi[i]
    }
    temp$zi<-zi
    for(i in 1:k){ if (temp$zi[i]==-Inf)
    {	surv.max.1<-max(temp$km.surv[temp$status==1])
    d<-1-surv.max.1
    surv.pu<-1-d/2	
    temp$zi[i]<-qnorm(1-surv.pu,0,1)}
    }
    for (i in 1:k){ if (temp$zi[i]==Inf)
    {  d<-min(temp$km.surv[temp$km.surv > 0])
    surv.pl<-d/2
    temp$zi[i]<-qnorm(1-surv.pl,0,1)}
    }
  }
  temp$ei<-temp$ei-con
  ##print(temp)
  plot(temp$zi,temp$ei,xlab=xlab,ylab="ordered ei residuals",type="n",xlim=c(min(temp$zi),max(temp$zi)),ylim=c(min(temp$ei),max(temp$ei)+.15))
  points(temp$zi[temp$status==0],temp$ei[temp$status==0],pch=".",cex=3)
  points(temp$zi[temp$status==1],temp$ei[temp$status==1],pch="o",cex=1)
  lines(temp$zi[temp$status==1],temp$ei[temp$status==1],lty=1,lwd=1)
  k<-nrow(temp)
  for(i in 1:k)
  { if (temp$status[i]==0)
    
    arrows(temp$zi[i],temp$ei[i],temp$zi[i],temp$ei[i]+.15,code=2,length=.12,lwd=2.01)}
  abline(a=0,b=1,lty=4,lwd=2)
  usr<-par("usr")
  arrows(.9*usr[1]+.1*usr[2],.07*usr[3]+.93*usr[4],.9*usr[1]+.1*usr[2],.07*usr[3]+.93*usr[4]+.15,code=2,length=.12,lwd=2.01)
  text(.8*usr[1]+.2*usr[2],.05*usr[3]+.95*usr[4],"= censored")
  points(.9*usr[1]+.1*usr[2],.11*usr[3]+.89*usr[4],pch="o")
  text(.79*usr[1]+.21*usr[2],.1*usr[3]+.90*usr[4], "= uncensored") 
  ##===============================================
  ##Eample 1:  
  ##fit.weib<-survreg(Surv(time,status)~x,dist="weibull",data=motorette)
  ##qq.reg.resid.r(motorette,motorette$time,motorette$status,fit.weib,"qweibull","standard extreme value quantiles")
  ##Example 2: 
  ##fit.weib<-survreg(Surv(weeks,status)~1,dist="weibull",data=aml1)
  ##qq.reg.resid.r(aml1,aml1$weeks,aml1$status,fit.weib,"qweibull","standard extreme value quantiles")
  ##==================================================
  on.exit()
  "qq.reg.resid:done"
}

# Weibull qq plot for right censored data

qq.weibreg <- function(data,weib.fit,xlab = "standard extreme value quantiles", 
                       ylab = "ordered log data", pch = NULL, lty = NULL, col= NULL) 
{
  ## Purpose: qqplot for Weibull regression with one covariate which has a finite number 
  ##          of distinct values which can be coverted into levels. It draws lines with same 
  ##          slope sigma, but different intercepts, whose values are the MLE's
  ##-------------------------------------------------------------------
  ## Arguments: ##  data  A Surv object, or a list of such objects
  ##            ##  weib.fit is a survreg object with dist="weibull"    
  ##-------------------------------------------------------------------
  ## Author: Mara Tableman  Date: 18 December 2002, 20:30
  on.exit(browser())
  t.c <- class(data)
  if((!is.null(t.c)) && t.c == "Surv")
    data <- list(data)
  t.k <- length(data)
  if(is.null(pch))
    pch <- 1:t.k
  if(is.null(lty))
    lty <- 1:t.k
  if(is.null(col))
    col <- 1 + 1:t.k
  t.sf <- lapply(data, function(s)
  {
    t.s <- summary(survfit(s~1, type = "kaplan-meier"))
    t.ss <- t.s$surv
    t.s$surv <- (t.ss + c(1, t.ss[ - length(t.ss)]))/2
    t.s
  }
  )
  t.rgs <- range(sapply(t.sf, function(s)
    range(s$surv)))
  t.rgt <- range(sapply(t.sf, function(s)
    range(s$time)))
  plot(log(qweibull(1 - t.rgs, 1)), log(t.rgt), type = "n",ylim=c(min(log(t.rgt)), max(log(t.rgt))),
       xlab = xlab, ylab = ylab)
  for(t.g in 1:t.k) {
    points(log(qweibull(1 - t.sf[[t.g]]$surv, 1)), log(t.sf[[t.g]]$time), pch=pch[t.g], col=1)
    mutilde<-as.numeric(levels(factor(weib.fit$linear.predictors))) 
    scale<-rep(weib.fit$scale,length(mutilde))
    abline(mutilde[t.g],scale[t.g],lty=lty[t.g],col=1)   
  }
  on.exit()
  "qq.weibreg:done"
}

#xln_levels(factor(motorette$x))
#ts.1_Surv(motorette$time[as.factor(motorette$x)==xln[1]],motorette$status[as.factor(motorette$x)==xln[1]]) 
#ts.2_Surv(motorette$time[as.factor(motorette$x)==xln[2]],motorette$status[as.factor(motorette$x)==xln[2]])
#ts.3_Surv(motorette$time[as.factor(motorette$x)==xln[3]],motorette$status[as.factor(motorette$x)==xln[3]])
#weib.fit_survreg(Surv(time,status)~x, data=motorette,dist="weib")
#qq.weibreg(list(ts.1,ts.2,ts.3),weib.fit)

##==========================================================
qq.weibull<-
  function(data,scale=0, xlab = "standard extreme value quantiles", ylab = 
             "ordered log data", pch = NULL, lty = NULL, col = NULL)
  {
    ## Purpose: qqplot for Weibull distribution for one or several samples 
    ## It fits each sample with own intercept and slope (location and scale).
    ##-------------------------------------------------------------------
    ## Arguments: ##  data  A Surv object, or a list of such objects
    ##       options:  scale=0 is the default. This estimates the scale.
    ##                 scale=1 fits the exponential model. 
    ##-------------------------------------------------------------------
    ## Author: Werner Stahel, Date: 15 Aug 2002, 18:41
    on.exit(browser())
    t.c <- class(data)
    if((!is.null(t.c)) && t.c == "Surv")
      data <- list(data)
    t.k <- length(data)
    if(is.null(pch))
      pch <- 1:t.k
    if(is.null(lty))
      lty <- 1:t.k
    if(is.null(col))
      col <- 1 + 1:t.k
    t.sf <- lapply(data, function(s)
    {
      t.s <- summary(survfit(s~1, type = "kaplan-meier"))
      t.ss <- t.s$surv
      t.s$surv <- (t.ss + c(1, t.ss[ - length(t.ss)]))/2
      t.s
    }
    )
    t.rgs <- range(sapply(t.sf, function(s)
      range(s$surv)))
    t.rgt <- range(sapply(t.sf, function(s)
      range(s$time)))
    plot(log(qweibull(1 - t.rgs, 1)), log(t.rgt), type = "n", xlab
         = xlab, ylab = ylab)
    for(t.g in 1:t.k) {
      points(log(qweibull(1 - t.sf[[t.g]]$surv, 1)), 
             log(t.sf[[t.g]]$time), pch = pch[t.g], col = 1)
      t.r <- survreg(data[[t.g]] ~ 1, dist = "weibull",scale=scale)
      abline(t.r$coef, t.r$scale, lty = lty[t.g], col = 1)
    }
    on.exit()
    "qq.weibull:done"
  }

#xln_levels(factor(motorette$x))
#ts.1_Surv(motorette$time[as.factor(motorette$x)==xln[1]],motorette$status[as.factor(motorette$x)==xln[1]]) 
#ts.2_Surv(motorette$time[as.factor(motorette$x)==xln[2]],motorette$status[as.factor(motorette$x)==xln[2]])
#ts.3_Surv(motorette$time[as.factor(motorette$x)==xln[3]],motorette$status[as.factor(motorette$x)==xln[3]])
#qq.weibull(list(ts.1,ts.2,ts.3))

##======================================================
quantile.km <- function(data, p, eps, z)
{
  ## data is survfit object, p is between 0 and 1
  ## eps is epsilon 0f 0.05 or bigger, 
  ## z iz z-score for confidence coefficient
  time <- summary(data)$time
  ni <- summary(data)$n.risk
  di <- summary(data)$n.event
  surv <- summary(data)$surv
  stderr <- summary(data)$std.err
  qp <- min(time[surv <= 1 - p])	
  ## The point estimate of pth-quantile 
  se.S.qp <- stderr[surv == max(surv[surv <= 1 - p])]	
  ## S.qp is the standard error of the estimated survival
  ## probability at qp.
  u.p <- max(time[surv >= 1 - p + eps])	
  # the largest time at which surv >= 1-p+eps
  l.p <- min(time[surv <= 1 - p - eps])	
  # the smallest time at which surv <=1-p-eps
  S.u.p <- surv[time == u.p]	# survival probability at u.p
  S.l.p <- surv[time == l.p]	# survival probability at l.p
  f.qp <- (S.u.p - S.l.p)/(l.p - u.p)	
  ## estimated probability density at pth-quantile
  se.qp <- se.S.qp/f.qp	
  ## estimated standard error of the sample pth-quantile
  LCL <- qp - z * se.qp
  UCL <- qp + z * se.qp
  out <- round(data.frame(qp, se.S.qp, f.qp, se.qp, LCL, UCL), 4)
  print("summary")
  print(out, invisible(1))	
  ##print("An approximate 1-alpha confidence interval for the true pth-quantile")
  on.exit()
  "quantile.km:done"
}

##==============================================
## Produces the log of the Weibull likelihood function evaluated 
## at a parameter value specified under the null hypothesis
## Arguements:  time, status, shape=alpha (=1 for exponential), 
##              theta=1/lambda, where lambda is the scale for Weibull.
## Author:  Mara Tableman  23 December 2002  23:50.
weib.loglik.theta <- function(time,status,shape,theta)
{out<-log(prod(dweibull(time, shape, theta)^(status)*(1-pweibull(time, shape, theta))^(1 - status)))
return(out)
on.exit()
"weib.loglik.theta:done"}

##======================================================
emphazplot<-function(data,text=" "){ 
  ##=============================
  ## Purpose:   To plot the empirical hazards, hitilde and hihat,
  ##            over time of one or two samples. This also
  ##            prints out table of times, hitilde's, and hihat's for
  ##            each group.
  ## Argument: The data argument is a list of Surv objects
  ##           The text argument; e.g., text="solid line is maintained group",
  ##           is meant to identify one of the two lines if there are two groups.
  ##           The solid line identifies the second group in the list.
  ##           The default is no text.
  ## Author:  Mara Tableman  January 29, 2003  20:00
  ##=============================
  my<-lapply(data,function(s){
    ni<-summary(survfit(s,type="kaplan-meier"))$n.risk
    di<-summary(survfit(s,type="kaplan-meier"))$n.event
    time<-summary(survfit(s,type="kaplan-meier"))$time
    surv<-summary(survfit(s,type="kaplan-meier"))$surv
    hitilde<-di/ni
    tau<-diff(time,lag=1) #Length of interval to right of ti
    tau[length(tau)+1]<-NA
    hihat<-hitilde/tau
    hihat[length(time)]<-hihat[(length(time)-1)]
    survtable<-data.frame(time,hitilde,hihat)
    print(round(survtable,3))
    hitilde<-round(hitilde,3)
    hihat<-round(hihat,3)
    my.try<-list(data.frame(time,hitilde,hihat))
    
  }
  )
  t.k<-length(data)
  par(mfrow=c(2,2))
  lty<-c(5,1)
  t.rg<-range(sapply(my, function(s) range(s[[1]]$time)))
  h.rg<-range(sapply(my, function(s) range(s[[1]]$hitilde)))
  plot(t.rg,h.rg,type="n",xlab="observed failure times ",ylab=" ")
  for(t.g in 1:t.k){
    lines(my[[t.g]][[1]]$time,my[[t.g]][[1]]$hitilde,type="s",lty=lty[t.g],col=1)
    points(my[[t.g]][[1]]$time,my[[t.g]][[1]]$hitilde,pch=16,cex=.5,col=1)
    #axis(1,at=my[[t.g]][[1]]$time,outer=F)
    #axis(2,at=my[[t.g]]$hitilde,outer=F,lab=paste(my[[t.g]]$hitilde),tck=.02,cex=.5)
    mtext("hazard at time i",side=2,line=2)
    mtext("hitilde at ti",side=3,line=1)
    #mtext(text,side=3,line=-1.3)
  }
  h.rg<-range(sapply(my, function(s) range(s[[1]]$hihat)))
  plot(t.rg,h.rg,type="n",xlab="observed failure times ",ylab=" ")
  for(t.g in 1:t.k){
    lines(my[[t.g]][[1]]$time,my[[t.g]][[1]]$hihat,type="s",lty=lty[t.g],col=1)
    points(my[[t.g]][[1]]$time,my[[t.g]][[1]]$hihat,pch=16,cex=.5,col=1)
    #axis(1,at=my[[t.g]][[1]]$time,outer=F,lab=paste(my[[t.g]][[1]]$time))
    #axis(2,at=my[[t.g]]$hihat,outer=F,lab=paste(my[[t.g]]$hihat),tck=.02,cex=.5)
    mtext("hazard over each interval",side=2,line=2)
    mtext("hihat", side=3,line=1)
    #mtext(text,side=3,line=-1.3)
  }
  on.exit()
  "emphazplot:done"   
} # The function ends here.

##====== 
## Example
##t.sa<-Surv(aml$weeks[aml$group==0],aml$status[aml$group==0])
##t.sb<-Surv(aml$weeks[aml$group==1],aml$status[aml$group==1]) # maintained group
##l.sab<-list(t.sa,t.sb)
##data<-l.sab
##emphazplot(data,text="solid line is maintained group")


#=============optimal.change.point function
optimal.change.point<-function(data,end,status,trt,...){
  #==========================
  #Purpose is to put the dataset in the Andersen-Gill counting process format,
  #then to obtain the optimal change point (time); that is, the distinct time point
  #which maximizes the profile likelihood. 
  #This is the change point value to used to fit a piecewise Cox PH model over two intervals .
  #Arguments
  # data = data frame
  # end  = time of event variable
  # status = censor variable: 1 if dead, 0 if alive
  # trt = the treatment or exposure variable consisiting of two groups.  Best to code 1 and 0. 
  #=========================
  # Author:  Mara Tableman, August 5, 2010
  #==========================
  cpt<-survfit(Surv(end,status)~1,conf.type="n",data=data)$time
  data<-cbind(data,end,status,trt)
  profile<-rep(1,length(cpt))
  for (i in 1:length(cpt)){
    AG<-survSplit(data,cut=cpt[i],end="end",start="start",event="status")
    AG$ET1<-AG$trt
    AG$ET1[AG$start!=0] <-0
    AG$ET2<-AG$trt
    AG$ET2[AG$start==0] <-0
    profile[i]<-coxph(Surv(start,end,status)~ET1+ET2,data=AG)$loglik[2]
  }
  out<-data.frame(cpt,profile)
  names(out)<-c("changepoint","loglik")
  optimal.changepoint<-out[out$loglik==max(out$loglik), ]
  print(optimal.changepoint)
  plot(out,type="l",lty=1,xlab="change point (distinct survival times)",ylab="log-likelihood",...)
  title(main="Profile of the log-partial likelihood \n for a piecewise Cox PH model")
}


#=============extcox.Et function
extcox.1Et<-function(data,end,status,trt,cut){
  #==========================
  #purpose is to put the dataset in the Andersen-Gill counting process format
  #so that a piecewise Cox PH model over two intervals can be fit.
  #Arguments
  # data = data frame
  # end  = time of event variable
  # status = censor variable: 1 if dead, 0 if alive
  # trt = the treatment or exposure variable consisitng of two groups.  Best to code 1 and 0. 
  # cut = the change point t0
  #=========================
  # Author:  Mara Tableman, August 5, 2010
  #==========================
  
  data<-cbind(data,end,status,trt)
  AG<-survSplit(data,cut=cut,end="end",start="start",event="status")
  AG$ET1<-AG$trt
  AG$ET1[AG$start!=0] <-0
  AG$ET2<-AG$trt
  AG$ET2[AG$start==0] <-0
  return(AG)
}
#=============