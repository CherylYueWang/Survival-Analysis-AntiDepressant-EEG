rm(list = ls()) # clear memory
source("TK.R.functions1.R")
# read data
eeg = read.csv("eeg_data_cleaned.csv",sep = ",")
new_eeg = na.omit(eeg)
data = new_eeg[,-1]
library(survival)
library(rTensor)

X_covariate = data[,4:633]

n = dim(X_covariate)[1]
eeg_tnsr=array(NA,dim=c(14,45,n))
for (obs in 1:n){
  for (q in 1:45){
    for (p in 1:14){
      eeg_tnsr[p,q,obs] = X_covariate[obs,p+14*(q-1)]
    }
  }
}
dim(eeg_tnsr)  
mpca_x=mpca(as.tensor(eeg_tnsr[,,1:n]),ranks=c(2,3),max_iter = 1000, tol=1e-1)
mpca_x$conv
attributes(mpca_x)
temp_A=mpca_x$U[[1]]
temp_B=mpca_x$U[[2]]

# percentage of variation explained by MPCA
eeg_mean=apply(eeg_tnsr,c(1,2),mean)
fnorm_den=rep(NA,n)
for(j in 1:n){
  tmp=eeg_tnsr[,,j]-eeg_mean
  fnorm_den[j]=sum(tmp^2)
  
}
temp_U = NULL
temp_u_vec=NULL
fnorm_num_MPCA=rep(NA,n)

for (i in 1:n){
  
  # tempU = t(temp_A)%*%attributes(mpca_x$est)$data[,,i]%*%temp_B # the same as 
  
  tempU = t(temp_A)%*%eeg_tnsr[,,i]%*%temp_B		
  temp_U = rbind(temp_U,as.matrix(tempU))
  
  temp_u_vec = rbind(temp_u_vec,as.vector(tempU))
  if(i<=n){
    fnorm_num_MPCA[i]=sum((attributes(mpca_x$est)$data[,,i]-eeg_mean)^2)
  }
  
}
sum(fnorm_num_MPCA)/sum(fnorm_den)
# 82% of variantion is explianed by MPCA.

# cox model
# reduced data

# define the new covariates
U_covariate = matrix(NA,n,6)
for (i in 1:n){
  temp = t(temp_A)%*%eeg_tnsr[,,i]%*%temp_B
  U_covariate[i,] = as.vector(temp)
  
}

reduced_data = cbind(data[,2:3],U_covariate) # 68*(3+3*3)
colnames(reduced_data) = c("time","status","V1","V2","V3","V4","V5","V6")

attach(reduced_data)

# parametric
weibull.fit = survreg(Surv(time,status)~., data=reduced_data, dist="weib")
log_logistic.fit =survreg(Surv(time,status)~., data=reduced_data, dist="loglogistic")
log_normal.fit = survreg(Surv(time,status)~., data=reduced_data, dist="lognormal")
exp.fit = survreg(Surv(time,status)~., data=reduced_data, dist="weib",scale = 1)
par(mfrow = c(2,2))
qq.weibull(Surv(time,status)) # weibull
qq.loglogistic(Surv(time,status)) # loglogistic
qq.lognormreg(Surv(time,status),log_normal.fit) # lognormal
qq.weibull(Surv(time,status),scale = 1) # exponential

summary(exp.fit)
summary(log_normal.fit)
summary(log_logistic.fit)
summary(weibull.fit)

extractAIC(exp.fit)
extractAIC(log_normal.fit)
extractAIC(log_logistic.fit)
extractAIC(weibull.fit)
aicstat=c(extractAIC(exp.fit),extractAIC(log_normal.fit),extractAIC(log_logistic.fit),extractAIC(weibull.fit))

aicstat=aicstat[c(2,4,6,8)]
names(aicstat)=c("Exp","Lognormal","loglogistic","Weib")
aicstat

# exponential
summary(exp.fit)
exp_reduce_coeff = exp.fit$coefficients

# point estimator of hazard ratio
exp(-exp_reduce_coeff[-1])

# 95% CI for reduced hazard ratio e^betah


lowerbd = c()
uperbd = c()
for (i in 1:6){
  deri = -exp(-exp_reduce_coeff[i+1])
  cov_be = deri*exp.fit$var[i+1,i+1]*deri
  point = exp(-exp_reduce_coeff[i+1])
  lbd_h = point-1.96*sqrt(cov_be)
  ubd_h = point+1.96*sqrt(cov_be)
  lowerbd[i] = lbd_h
  uperbd[i] = ubd_h
}

uperbd = as.matrix(uperbd)
lowerbd = as.matrix(lowerbd)
point = as.matrix(exp(-exp_reduce_coeff[-1]))
CI = cbind(lowerbd,point,uperbd)
rownames(CI) = c("V1","V2","V3","V4","V5","V6")
colnames(CI) = c("lower bound","point est","upper bound")
CI

# transform back to original 14*45
reduce_coeff_matrix = matrix(exp.fit$coefficients[-1],2,3)
coeff_matrix = temp_A%*%reduce_coeff_matrix%*%t(temp_B)

original_cov_matrix = kronecker(temp_B,temp_A)%*%exp.fit$var[2:7,2:7]%*%kronecker(t(temp_B),t(temp_A)) # 630*630
original_coeff = as.vector(coeff_matrix)
# Point estimator of hazard ratio in original scale
write.csv(round(exp(-coeff_matrix),5),"exp_original_HR.csv") # hazard ratio matrix

# confidence interval for original hazard ratio
lbd_ori = rep(0,630)
ubd_ori = rep(0,630)


for (i in 1:630){
  deri = -exp(-original_coeff[i])
  cov_be = deri*original_cov_matrix[i,i]*deri
  point = exp(-original_coeff[i])
  lbd_h = point-1.96*sqrt(cov_be)
  ubd_h = point+1.96*sqrt(cov_be)
  lbd_ori[i] = lbd_h
  ubd_ori[i] = ubd_h
  
}

write.csv(round(matrix(lbd_ori,14,45),5),"exp_original_lower.csv")
write.csv(round(matrix(ubd_ori,14,45),5),"exp_original_upper.csv")

# cox snell residual plots
summary(exp.fit)
scale=exp(exp_reduce_coeff[1]+exp_reduce_coeff[2]*reduced_data$V1+exp_reduce_coeff[3]*reduced_data$V2+
            exp_reduce_coeff[4]*reduced_data$V3+exp_reduce_coeff[5]*reduced_data$V4+
            exp_reduce_coeff[6]*reduced_data$V5+exp_reduce_coeff[7]*reduced_data$V6)
s=exp(-(time/scale))
H=-log(s)
r <- (H + 0.693*(reduced_data$status==0)) # adjusting for censored data and this is the residuals.
KM.r = survfit(Surv(r,reduced_data$status)~1) #survival analysis on residuals, r, as if it is the time variable.
LSr = -log(KM.r$surv) 
r = sort(unique(r)) 
par(mfrow=c(1,1))
plot(r,LSr,main="Cox-Snell residual plot from Exponential model")#this is the LS plot.
abline(0,1,col="red")

# cox model

cox_full<-coxph(Surv(time,status)~.,data = reduced_data)
cox_full

# confidence interval for hazard ratio
reduce_coeff = cox_full$coefficients
reduce_cov = cox_full$var
point = exp(reduce_coeff)
lbd = rep(0,6)
ubd = rep(0,6)
for (i in 1:6){
  derivat = exp(reduce_coeff[i])
  var = derivat*reduce_cov[i,i]*derivat
  lower = derivat - 1.96*sqrt(var)
  uper = derivat + 1.96*sqrt(var)
  lbd[i] = lower
  ubd[i] = uper
}
point = as.matrix(point)
lbd = as.matrix(lbd)
ubd = as.matrix(ubd)
CI = cbind(lbd,point,ubd)
colnames(CI) = c("lower bound", "point est", "upper bound")
CI
# transform back to the original 14*45 scale
reduce_coeff_matrix = matrix(reduce_coeff,2,3)

coeff_matrix = temp_A%*%reduce_coeff_matrix%*%t(temp_B)

write.csv(round(exp(coeff_matrix),5),"cox_original_HR.csv") # hazard ratio matrix


# 95% CI
# use transformed covariance matrix and delta method
original_cov_matrix = kronecker(temp_B,temp_A)%*%reduce_cov%*%kronecker(t(temp_B),t(temp_A)) # 630*630
original_coeff = as.vector(coeff_matrix) # 630 length
# confidence interval for original hazard ratio
lbd_ori = rep(0,630)
ubd_ori = rep(0,630)





for (i in 1:630){
  derivat = exp(original_coeff[i])
  var = derivat*original_cov_matrix[i,i]*derivat
  lower = derivat - 1.96*sqrt(var)
  uper = derivat + 1.96*sqrt(var)
  lbd_ori[i] = lower
  ubd_ori[i] = uper
}

write.csv(round(matrix(lbd_ori,14,45),5),"cox_original_CI_lower.csv")
write.csv(round(matrix(ubd_ori,14,45),5),"cox_original_CI_upper.csv")

# no significant effect


# PH assumption
PH=cox.zph(cox_full)
PH
par(mfrow=c(2,3))
plot(PH)

# residual analysis
scatter.smooth(reduced_data$V1,resid(cox_full),type="p",pch=20, cex=1,xlab="V1",ylab="Martingale residual")
scatter.smooth(reduced_data$V2,resid(cox_full),type="p",pch=20, cex=1,xlab="V2",ylab="Martingale residual")
scatter.smooth(reduced_data$V3,resid(cox_full),type="p",pch=20, cex=1,xlab="V3",ylab="Martingale residual")
scatter.smooth(reduced_data$V4,resid(cox_full),type="p",pch=20, cex=1,xlab="V4",ylab="Martingale residual")
scatter.smooth(reduced_data$V5,resid(cox_full),type="p",pch=20, cex=1,xlab="V5",ylab="Martingale residual")
scatter.smooth(reduced_data$V6,resid(cox_full),type="p",pch=20, cex=1,xlab="V6",ylab="Martingale residual")

# The deviance residual plots to detect outliers:


dresid <- resid(cox_full,type="deviance") # deviance residual
par(mfrow=c(1,1))
plot(dresid,pch=3, cex=1, col='magenta')
abline(h=0)
par(mfrow=c(2,3))
plot(reduced_data$V1,dresid,pch=2, cex=1) 
abline(h=0)


plot(reduced_data$V2,dresid,type="p",pch=".")
abline(h=0)

plot(reduced_data$V3,dresid, pch=2, cex=1)
abline(h=0)

plot(reduced_data$V4,dresid, pch=3, cex=1)
abline(h=0)

plot(reduced_data$V5,dresid,pch=19, cex=1) 
abline(h=0)
plot(reduced_data$V6,dresid,pch=19, cex=1) 
abline(h=0)


# Schoenfeld residuals to examine fit and detect outlying covariate values

detail <- coxph.detail(cox_full) # detailed coxph object 
time <- detail$y[,2]  # ordered times including censored ones 
status <- detail$y[,3]  # censoring status
sch <- resid(cox_full,type="schoenfeld") # Schoenfeld
plot(time[status==1],sch[,1],xlab="Ordered survival time",
     ylab="Schoenfeld residual for V1") 

plot(time[status==1],sch[,2],xlab="Ordered survival time",
     ylab="Schoenfeld residual for V2") 

plot(time[status==1],sch[,3],xlab="Ordered survival time",
     ylab="Schoenfeld residual for V3") 

plot(time[status==1],sch[,4],xlab="Ordered survival time",
     ylab="Schoenfeld residual for V4") 
plot(time[status==1],sch[,5],xlab="Ordered survival time",
     ylab="Schoenfeld residual for V5") 
plot(time[status==1],sch[,6],xlab="Ordered survival time",
     ylab="Schoenfeld residual for V6") 


# df beta
bresid <- resid(cox_full,type="dfbetas")
index <- seq(1:dim(reduced_data)[1])
plot(index,bresid[,1],type="h",ylab="scaled change in coef",
     xlab="observation")
plot(index,bresid[,2],type="h",ylab="scaled change in coef",
     xlab="observation")
plot(index,bresid[,3],type="h",ylab="scaled change in coef",
     xlab="observation")
plot(index,bresid[,4],type="h",ylab="scaled change in coef",
     xlab="observation")
plot(index,bresid[,5],type="h",ylab="scaled change in coef",      
     xlab="observation")
plot(index,bresid[,6],type="h",ylab="scaled change in coef",      
     xlab="observation")

# cox snell residual plots
par(mfrow=c(1,1))
rc <- abs(reduced_data$status - cox_full$residuals) # creating Cox-Snell residuals using Martingales.

km.rc <- survfit(Surv(rc,reduced_data$status) ~ 1)
summary.km.rc <- summary(km.rc)
rcu <- summary.km.rc$time # Cox-Snell residuals of uncensored points.
surv.rc <- summary.km.rc$surv
plot(rcu,-log(surv.rc),pch=20, cex=2,col='blue',
     xlab="Cox-Snell residual rc",ylab="Cumulative hazard on rc") 
abline(a=0,b=1); abline(v=0); abline(h=0)

