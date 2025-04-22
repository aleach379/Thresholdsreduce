
### Libraries that need to be loaded ###

library(metafor)
library(readr)
library(meta)

### We import the dataset and attach the variables ###

d0 <- read_csv("C:/Users/kboura01/OneDrive - University of Cyprus/AgStat/meta-analysis/Threshold database_March2023.csv")
View(d0)

attach(d0)
table(SOURCE) # this refers to the number of rows (inputs) corresponding to each source (study)

colnames(d0)
#######################
#### ENEMY DENSITY ####
#######################

# We select the subset of the file related to the variables under study and 
# transform the variables into the appropriate format.

d01=d0[,c("SOURCE","trt_functional","replicates","enemy_density","natdensity_SE")]

View(d01)

d01$SOURCE=as.numeric(d01$SOURCE)
d01$trt_functional=as.factor(d01$trt_functional)
d01$replicates=as.numeric(d01$replicates)
d01$enemy_density=as.numeric(d01$enemy_density)
d01$natdensity_SE=as.numeric(d01$natdensity_SE)

d01=d01[-which(is.na(d01$enemy_density)),] # we exclude the rows that do not provide an effect size
dim(d01) # the number or rows and columns of the final dataset
table(d01$SOURCE)

# In the following procedure, we create three tables for control, standard, and threshold, respectively.
# Specifically, for each study (source) that has multiple entries for the same design for trt_functional (i.e., control, standard, and threshold), 
# we calculate the weighted mean effect size and the corresponding standard error, 
# setting the final number of replicates as the sum of replicates across all entries. 
# This process is applied separately for each study and each design.

dvth=data.frame(SOURCE=NA, trt_functional=NA, replicates=NA, enemy_density=NA, natdensity_SE=NA)
dvst=data.frame(SOURCE=NA, trt_functional=NA, replicates=NA, enemy_density=NA, natdensity_SE=NA)
dvc=data.frame(SOURCE=NA, trt_functional=NA, replicates=NA, enemy_density=NA, natdensity_SE=NA)

ind=1
for (i in as.numeric(names(table(d01$SOURCE)))) {
  
  dv0th=d01[which((d01$SOURCE==i) & (d01$trt_functional=="threshold")),]
  dvth[ind,1]=i
  dvth[ind,2]="threshold"
  
  if (sum(is.na(dv0th[,3]))==0) {  
    
    dvth[ind,3]=round(sum(dv0th[,3]))
    dvth[ind,4]=round(sum(dv0th[,3]*dv0th[,4])/sum(dv0th[,3]),3)
    dvth[ind,5]=round(sum(dv0th[,3]*dv0th[,5])/sum(dv0th[,3]),3) } else {
      
      dvth[ind,3]=NA
      dvth[ind,4]=round(mean(as.numeric(unlist(dv0th[,4]))),3)
      dvth[ind,5]=round(mean(as.numeric(unlist(dv0th[,5]))),3) }

  ################################
  
  dv0st=d01[which((d01$SOURCE==i) & (d01$trt_functional=="standard")),]
  dvst[ind,1]=i
  dvst[ind,2]="standard"
  
  if (sum(is.na(dv0st[,3]))==0) {  
    
    dvst[ind,3]=round(sum(dv0st[,3]))
    dvst[ind,4]=round(sum(dv0st[,3]*dv0st[,4])/sum(dv0st[,3]),3)
    dvst[ind,5]=round(sum(dv0st[,3]*dv0st[,5])/sum(dv0st[,3]),3) } else {
      
      dvst[ind,3]=NA
      dvst[ind,4]=round(mean(as.numeric(unlist(dv0st[,4]))),3)
      dvst[ind,5]=round(mean(as.numeric(unlist(dv0st[,5]))),3) }
  
  ################################
  
  dv0c=d01[which((d01$SOURCE==i) & (d01$trt_functional=="control")),]
  dvc[ind,1]=i
  dvc[ind,2]="control"
  
  
  if (sum(is.na(dv0c[,3]))==0) {  
    
    dvc[ind,3]=round(sum(dv0c[,3]))
    dvc[ind,4]=round(sum(dv0c[,3]*dv0c[,4])/sum(dv0c[,3]),3)
    dvc[ind,5]=round(sum(dv0c[,3]*dv0c[,5])/sum(dv0c[,3]),3) } else {
      
      dvc[ind,3]=NA
      dvc[ind,4]=round(mean(as.numeric(unlist(dv0c[,4]))),3)
      dvc[ind,5]=round(mean(as.numeric(unlist(dv0c[,5]))),3)
      
    }
  
  ind=ind+1
  
}

View(dvth)
View(dvst)
View(dvc)

# To perform pairwise comparisons between different designs using the methodology of Nakagawa et al. (2023), 
# we need studies that have both designs for which the comparison is made 
# to calculate the corresponding log-ratio of effect sizes for each specific study. 
# With the following procedure, we retain studies that have common pairs of 
# standard-control, threshold-control, and threshold-standard, respectively.

### control - standard ### 

del01=union(which(is.na(dvst[,3])), which(is.na(dvc[,3])))

Dst01=dvst[-del01,]
Dc01=dvc[-del01,]

ex01=union(which(Dc01$replicates==0), which(Dst01$replicates==0))

dim(Dst01) ; dim(Dc01)

### control - threshold ### 

del02=union(which(is.na(dvth[,3])), which(is.na(dvc[,3])))

Dth02=dvth[-del02,]
Dc02=dvc[-del02,]

ex02=union(which(Dc02$replicates==0), which(Dth02$replicates==0))

dim(Dc02) ; dim(Dth02)

### standard - threshold ### 

del12=union(which(is.na(dvth[,3])), which(is.na(dvst[,3])))

Dst12=dvst[-del12,]
Dth12=dvth[-del12,]

ex12=union(which(Dth12$replicates==0), which(Dst12$replicates==0))

dim(Dst12) ; dim(Dth12)

################################################
# Geary's test
################################################
# Function to calculate Geary's "number"

geary = function(mean, sd, n){
  (1 / (sd / mean)) * ((4*n)^(3/2) / (1 + 4*n))
}

# We apply the Geary test to all log-ratios and exclude those that have values below 3.

### control - standard ### 


dst01=which(!is.na(Dst01[,5]))
dc01=which(!is.na(Dc01[,5]))

for (k in dst01) {
  
  f=geary(Dst01[k,4],Dst01[k,5],Dst01[k,3])
#  print(c(k,f))
  
  if (f<3) {print(k)} 
  
}


for (k in dc01) {
  
  f=geary(Dc01[k,4],Dc01[k,5],Dc01[k,3])
#  print(c(k,f))
  
  if (f<3) {print(k)} 
  
}


### control - threshold ### 


dth02=which(!is.na(Dth02[,5]))
dc02=which(!is.na(Dc02[,5]))

for (k in dth02) {
  
  f=geary(Dth02[k,4],Dth02[k,5],Dth02[k,3])
#  print(c(k,f))
  
  if (f<3) {print(k)} 
  
}


for (k in dc02) {
  
  f=geary(Dc02[k,4],Dc02[k,5],Dc02[k,3])
#  print(c(k,f))
  
  if (f<3) {print(k)} 
  
}


### standard - threshold ### 

dst12=which(!is.na(Dst12[,5]))
dth12=which(!is.na(Dth12[,5]))

for (k in dst12) {
  
  f=geary(Dst12[k,4],Dst12[k,5],Dst12[k,3])
#  print(c(k,f))
  
  if (f<3) {print(k)} 
  
}


for (k in dth12) {
  
  f=geary(Dth12[k,4],Dth12[k,5],Dth12[k,3])
#  print(c(k,f))
  
  if (f<3) {print(k)} 
  
}


# OK 

# Next, we apply the proposed “hybrid” method to estimate the standard error 
# of the effect size for the studies that they do not report it. 

################
##### lnRR #####
################

### control - standard ### 

### complete cases ###

cc01=intersect(dc01,dst01)
#cc01=cc01[-5]
CVc01=Dc01[cc01,5]/Dc01[cc01,4]*sqrt(Dc01[cc01,3])
CVst01=Dst01[cc01,5]/Dst01[cc01,4]*sqrt(Dst01[cc01,3])


lnRRc01 = log( Dst01[cc01,4] / Dc01[cc01,4]) + (CVst01^2/Dst01[cc01,3] - CVc01^2/Dc01[cc01,3]) / 2
vlnRRc01 = CVst01^2/Dst01[cc01,3] + CVc01^2/Dc01[cc01,3] + CVst01^4/(2*Dst01[cc01,3]^2) + CVc01^4/(2*Dc01[cc01,3])

### missing cases ### 

c101=(sum(Dc01[cc01,3]*CVc01)/sum(Dc01[cc01,3]))^2
st101=(sum(Dst01[cc01,3]*CVst01)/sum(Dst01[cc01,3]))^2

lnRRm01 = log( Dst01[-cc01,4] / Dc01[-cc01,4]) + (st101/Dst01[-cc01,3] - c101/Dc01[-cc01,3]) / 2
vlnRRm01 = (st101/Dst01[-cc01,3] + c101/Dc01[-cc01,3] +  st101^2/(2*Dst01[-cc01,3]^2) + c101^2/(2*Dc01[-cc01,3]^2) ) 


### control - threshold ### 

### complete cases ###

cc02=intersect(dc02,dth02)

CVc02=Dc02[cc02,5]/Dc02[cc02,4]*sqrt(Dc02[cc02,3])
CVth02=Dth02[cc02,5]/Dth02[cc02,4]*sqrt(Dth02[cc02,3])

lnRRc02 = log( Dth02[cc02,4] / Dc02[cc02,4]) + (CVth02^2/Dth02[cc02,3] - CVc02^2/Dc02[cc02,3]) / 2
vlnRRc02 = CVth02^2/Dth02[cc02,3] + CVc02^2/Dc02[cc02,3] + CVth02^4/(2*Dth02[cc02,3]^2) + CVc02^4/(2*Dc02[cc02,3])

### missing cases ### 

c102=(sum(Dc02[cc02,3]*CVc02)/sum(Dc02[cc02,3]))^2
th102=(sum(Dth02[cc02,3]*CVth02)/sum(Dth02[cc02,3]))^2


lnRRm02 = log( Dth02[-cc02,4] / Dc02[-cc02,4]) + (th102/Dth02[-cc02,3] - c102/Dc02[-cc02,3]) / 2
vlnRRm02 = (th102/Dth02[-cc02,3] + c102/Dc02[-cc02,3] +  th102^2/(2*Dth02[-cc02,3]^2) + c102^2/(2*Dc02[-cc02,3]^2) ) 


### standard - threshold ### 

### complete cases ###

cc12=intersect(dth12,dst12)

#cc12=cc12[-6]

CVth12=Dth12[cc12,5]/Dth12[cc12,4]*sqrt(Dth12[cc12,3])
CVst12=Dst12[cc12,5]/Dst12[cc12,4]*sqrt(Dst12[cc12,3])

lnRRc12 = log( Dth12[cc12,4] / Dst12[cc12,4]) + (CVth12^2/Dth12[cc12,3] - CVst12^2/Dst12[cc12,3]) / 2
vlnRRc12 = CVth12^2/Dth12[cc12,3] + CVst12^2/Dst12[cc12,3] + CVth12^4/(2*Dth12[cc12,3]^2) + CVst12^4/(2*Dst12[cc12,3])

### missing cases ### 

st112=(sum(Dst12[cc12,3]*CVst12)/sum(Dst12[cc12,3]))^2
th112=(sum(Dth12[cc12,3]*CVth12)/sum(Dth12[cc12,3]))^2


lnRRm12 = log( Dth12[-cc12,4] / Dst12[-cc12,4]) + (th112/Dth12[-cc12,3] - st112/Dst12[-cc12,3]) / 2
vlnRRm12 = (th112/Dth12[-cc12,3] + st112/Dst12[-cc12,3] +  th112^2/(2*Dth12[-cc12,3]^2) + st112^2/(2*Dst12[-cc12,3]^2) ) 


### We fit random-effects models using rma()

### control - standard ### 

lnRRi01=c(lnRRc01,lnRRm01)
vlnRRi01=c(vlnRRc01,vlnRRm01)

metadata01=data.frame(lnRRi01, vlnRRi01)

A01=rma(lnRRi01, vlnRRi01, data=metadata01)

# the estimated meta-analytic ratio along with the lower and upper bounds of the estimation 

exp(A01$beta)
exp(A01$ci.lb)
exp(A01$ci.ub)


### funnel plot

funnel(A01, main="Funnel plot standard/control for enemy density")


### control - threshold ### 

lnRRi02=c(lnRRc02,lnRRm02)
vlnRRi02=c(vlnRRc02,vlnRRm02)

metadata02=data.frame(lnRRi02, vlnRRi02)

A02=rma(lnRRi02, vlnRRi02, data=metadata02)

# the estimated meta-analytic ratio along with the lower and upper bounds of the estimation 

exp(A02$beta)
exp(A02$ci.lb)
exp(A02$ci.ub)

# funnel plot

funnel(A02, main="Funnel plot threshold/control for enemy density")


### standard - threshold ### 

lnRRi12=c(lnRRc12,lnRRm12)
vlnRRi12=c(vlnRRc12,vlnRRm12)

metadata12=data.frame(lnRRi12, vlnRRi12)

A12=rma(lnRRi12, vlnRRi12, data=metadata12)

# the estimated meta-analytic ratio along with the lower and upper bounds of the estimation 

exp(A12$beta)
exp(A12$ci.lb)
exp(A12$ci.ub)

# funnel plot

funnel(A12, main="Funnel plot threshold/standard for enemy density")


###########################################################
# We create a table with the results of the fitted models #
###########################################################

Re=matrix(NA,3,7)

Re[1,1]=length(lnRRi01) ; Re[1,2]=round(A01$I2,2) ; Re[1,3]=round(exp(A01$beta),2)  
Re[1,4]=round(exp(A01$beta)*A01$se,3) ; Re[1,5]=round(A01$pval,4)  
Re[1,6]=round(exp(A01$ci.lb),2) ; Re[1,7]=round(exp(A01$ci.ub),2)

Re[2,1]=length(lnRRi02) ; Re[2,2]=round(A02$I2,2) ; Re[2,3]=round(exp(A02$beta),2)  
Re[2,4]=round(exp(A02$beta)*A02$se,3) ; Re[2,5]=round(A02$pval,4)  
Re[2,6]=round(exp(A02$ci.lb),2) ; Re[2,7]=round(exp(A02$ci.ub),2)

Re[3,1]=length(lnRRi12) ; Re[3,2]=round(A12$I2,2) ; Re[3,3]=round(exp(A12$beta),2)  
Re[3,4]=round(exp(A12$beta)*A12$se,3) ; Re[3,5]=round(A12$pval,4)  
Re[3,6]=round(exp(A12$ci.lb),2) ; Re[3,7]=round(exp(A12$ci.ub),2)


colnames(Re)=c("k", "I2", "est", "se", "p.value", "ci.lb", "ci.ub")  

# we save the results in a text file 

write.table(Re, "C:/Users/kboura01/OneDrive - University of Cyprus/AgStat/meta-analysis/r-codes for metanalysis - final/re-enemy_density.txt", row.names=F)

# forest plots 

par(mfrow=c(1,3))

forest.rma(A01, cex=1.4, pch=16, cex.lab=1.4, cex.axis=1.4, shade = T, header=T)
mtext("log(st/c)",3, line=-3, cex=2.1)
forest.rma(A02, cex=1.4, pch=16, cex.lab=1.4, cex.axis=1.4, shade = T, header=T)
mtext("log(thr/c)",3, line=-3, cex=2.1)
forest.rma(A12, cex=1.4, pch=16, cex.lab=1.4, cex.axis=1.4, shade = T, header=T)
mtext("log(thr/st)",3, line=-3, cex=2.1)
mtext("enemy density",3, line=-3, cex=3, oute=TRUE)

