###############
##### EEG #####
###############
#setwd("532project/")
rm(list = ls())
load("train83_yue_wang.Rdata")
eeg83_all=train_data_ap_yue_wang$EEGClosed # eeg[[i]], i=1,...,72 is 92 times 511 matricies
#eeg_all=train51.data$EEGOpen
dim(eeg83_all[[1]])

electrodesC=rownames(eeg83_all[[1]])

#start at 0.25 Hz with consistent 0.25 increment
#delta <4
#theta: 4-7Hz
#alpha: 8-15Hz
#beta:16-31
#gamma >32

t0=4/0.25;t1=15/0.25
freq_c=seq(4,15,length=45)
freq=c((t0+1):(t1+1))
q=length(freq)

chan_c=c("P9","P10","P7","P8","P5","P6","PO7","PO8","PO3","PO4","O1","O2","POZ","OZ")
chan=c(1:72)[electrodesC %in% chan_c]
#chan=c(1:72)
p=length(chan)

eeg <- lapply(eeg83_all, function(x) x[chan,freq])

### MPCA R page for initial
n=length(eeg)
eeg_tnsr=array(NA,dim=c(p,q,n))
for(i in 1:n){
  
  eeg_tnsr[,,i]=as.matrix(eeg[[i]])
  
}

dim(eeg_tnsr)  

# write out the data
vec_eeg = matrix(0,n,p*q)
for (i in 1:n){
  vec_eeg[i,] = as.vector(eeg_tnsr[,,i])
}
# read outcome
outcome = read.csv("out.csv")
output = cbind(outcome, vec_eeg)
names = matrix(0,p,q)
for (colname in 1:p){
  for(rowname in 1:q){
    names[colname,rowname] = paste(chan_c[colname],"_HZ",(freq[rowname]-1)/4,sep = "")
  }
}
vec_name = as.vector(names)
full_name = cbind("ID","time","status",t(vec_name))
colnames(output) = full_name
write.csv(output, "eeg_data.csv")
#eeg_tnsr=log(eeg_tnsr)