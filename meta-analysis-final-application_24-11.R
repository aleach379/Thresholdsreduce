
### Libraries that need to be loaded ###

library(metafor)
library(readxl)
library(meta)

### We import the dataset and attach the variables ###

d0 = read_excel("C:/Users/kboura01/OneDrive - University of Cyprus/AgStat/meta-analysis/Final_threshold_database_June2023.xlsx")
View(d0)

attach(d0)
table(SOURCE) # this refers to the number of rows (inputs) corresponding to each source (study)


#####################
#### APPLICATION ####
#####################

# We select the subset of the file related to the variables under study and 
# transform the variables into the appropriate format.

d01=d0[,c("SOURCE","trt_functional","replicates","application","app_SE")]

View(d01)

d01$SOURCE=as.numeric(d01$SOURCE)
d01$trt_functional=as.factor(d01$trt_functional)
d01$replicates=as.numeric(d01$replicates)
d01$application=as.numeric(d01$application)
d01$app_SE=as.numeric(d01$app_SE)

dim(d01) # the number or rows and columns of the final dataset
View(d01)

# In the following procedure, we create three tables for control, standard, and threshold, respectively.
# Specifically, for each study (source) that has multiple entries for the same design for trt_functional (i.e., control, standard, and threshold), 
# we calculate the weighted mean effect size and the corresponding standard error, 
# setting the final number of replicates as the sum of replicates across all entries. 
# This process is applied separately for each study and each design.

dvth=data.frame(SOURCE=NA, trt_functional=NA, replicates=NA, application=NA, app_SE=NA)
dvst=data.frame(SOURCE=NA, trt_functional=NA, replicates=NA, application=NA, app_SE=NA)
dvc=data.frame(SOURCE=NA, trt_functional=NA, replicates=NA, application=NA, app_SE=NA)

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

Dst01=dvst
Dc01=dvc

ex01=union(which(Dc01$replicates==0), which(Dst01$replicates==0))
# 
Dc01=Dc01[-ex01,]
Dst01=Dst01[-ex01,]
dim(Dst01) ; dim(Dc01)

### control - threshold ### 

del02=union(which(is.na(dvth[,3])), which(is.na(dvc[,3])))

Dth02=dvth
Dc02=dvc

ex02=union(which(Dc02$replicates==0), which(Dth02$replicates==0))

Dc02=Dc02[-ex02,]
Dth02=Dth02[-ex02,]
dim(Dc02) ; dim(Dth02)

### standard - threshold ### 

del12=union(which(is.na(dvth[,3])), which(is.na(dvst[,3])))

Dst12=dvst
Dth12=dvth

ex12=union(which(Dth12$replicates==0), which(Dst12$replicates==0))
# 
Dst12=Dst12[-ex12,]
Dth12=Dth12[-ex12,]
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

#

for (k in dc01) {
  
  f=geary(Dc01[k,4],Dc01[k,5],Dc01[k,3])
#  print(c(k,f))
  
  if (f<3) {print(k)} 
  
}


# We cannot proceed with standard-control because all effect sizes for control are 0, 
# and all standard errors are also 0.

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

# We cannot proceed with threshold-control because all effect sizes for control are 0, 
# and all standard errors are also 0.


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



### standard - threshold ### 

### complete cases ###

cc12=intersect(dth12,dst12)

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

funnel(A12, main="Funnel plot threshold/standard for application")

###########################################################
# We create a table with the results of the fitted models #
###########################################################

Re=matrix(NA,3,7)

Re[3,1]=length(lnRRi12) ; Re[3,2]=round(A12$I2,2) ; Re[3,3]=round(exp(A12$beta),2)  
Re[3,4]=round(exp(A12$beta)*A12$se,3) ; Re[3,5]=round(A12$pval,4)  
Re[3,6]=round(exp(A12$ci.lb),2) ; Re[3,7]=round(exp(A12$ci.ub),2)

colnames(Re)=c("k", "I2", "est", "se", "p.value", "ci.lb", "ci.ub")  

# we save the results in a text file 

write.table(Re, "C:/Users/kboura01/OneDrive - University of Cyprus/AgStat/meta-analysis/r-codes for metanalysis - final/re-application.txt", row.names=F)

# forest plots

forest.rma(A12, cex=1.2, pch=16, cex.lab=1.4, cex.axis=1.2, shade = T, header=T)
mtext("log(thr/st)",3, line=-4, cex=2.1)
mtext("application",3, line=-3, cex=3, oute=TRUE)

