#==================================================
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
      ni<-summary(survfit(s~1,type="kaplan-meier"))$n.risk
      di<-summary(survfit(s~1,type="kaplan-meier"))$n.event
      time<-summary(survfit(s~1,type="kaplan-meier"))$time
      surv<-summary(survfit(s~1,type="kaplan-meier"))$surv
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
    })
    
    
    
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
      mtext(text,side=1,line=4)
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
      mtext(text,side=1,line=4)
    }
    on.exit()
    "emphazplot:done"   
  } # The function ends here.
  
  ##====== 
  ## Example
  #Surv0 <- Surv(LowFat$X1,LowFat$X2)
  #Surv1 <- Surv(SatFat$X1,SatFat$X2)  # Saturated Fat group
  #data<-list(Surv0,Surv1)
  
  #emphazplot(data,text="Solid line is saturated group")