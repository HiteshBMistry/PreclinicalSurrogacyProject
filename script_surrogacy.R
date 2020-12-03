require(survival)

all<-read.csv("novartis_volumes.csv")
all<-all[,c(1:4,6)]
colnames(all)<-c("ID","TID","TRT","VOL","DAYS")
DIAM<-2*((3*all$VOL/(4*pi))^(1/3))
all<-cbind(all,DIAM)



quantile(all$DIAM[all$DAYS==0])


# calculate RECIST response for each treatment
# treat this like you would a clinical trial - so we wouldn't have untreated data
all<-all[all$TRT!="untreated",]
# we can now fit a time-series model and ignore drop-out and use decay, 
# growth rate and see how well it does
NID<-array(NA,dim=c(dim(all)[1],1))
NID[1]<-1
for(ii in 2:length(NID)){
  if(all$DAYS[ii]>all$DAYS[ii-1]){
    NID[ii]<-NID[ii-1]
  }else{
    NID[ii]<-NID[ii-1]+1
  }
}
all<-cbind(all,NID)


idtrt<-unique(all[,c(1,3)])
final<-data.frame(array(NA, dim=c(dim(idtrt)[1],8)))
all$TRT<-as.character(all$TRT)
all$TID<-as.character(all$TID)
all$ID<-as.character(all$ID)

# calculate PFS & extras
for (ii in 1:dim(idtrt)[1]){
    dummy<-all[all$ID==idtrt[ii,1] & all$TRT==idtrt[ii,2],]
    bsl<-dummy$DIAM[dummy$DAYS==0]
    flag<-0
    if(length(bsl==1)){
      for (kk in 2:dim(dummy)[1]){
        if(dummy$DIAM[kk]>= max(1.2*min(dummy$DIAM[1:kk]),1) & flag==0){
          final[ii,1]<-dummy$ID[kk]
          final[ii,2]<-dummy$TID[kk]
          final[ii,3]<-dummy$TRT[kk]
          final[ii,4]<-100*(min(dummy$DIAM[2:kk])-bsl)/bsl
          final[ii,5]<-dummy$DAYS[kk]
          final[ii,6]<-1
          final[ii,7]<-(min(dummy$DIAM[2:kk]))/bsl
          final[ii,8]<-100*(min(dummy$DIAM[2])-bsl)/bsl
          flag<-1
        }
      }
      if (flag==0){
        final[ii,1]<-dummy$ID[kk]
        final[ii,2]<-dummy$TID[kk]
        final[ii,3]<-dummy$TRT[kk]
        final[ii,4]<-100*(min(dummy$DIAM[2:kk])-bsl)/bsl
        final[ii,5]<-dummy$DAYS[kk]
        final[ii,6]<-0
        final[ii,7]<-(min(dummy$DIAM[2:kk]))/bsl
        final[ii,8]<-100*(min(dummy$DIAM[2])-bsl)/bsl
        flag<-1
      }
      flag<-0 
    }
}
colnames(final)<-c("ID","TTYPE","TRT","BR","PDAY","PFS","FC","IPCH")  
# so we can look at PFS for each drug
unique(final$TRT)
km<-survfit(Surv(PDAY,PFS)~TRT,data=final)
test<-print(km)
length(unique(final$TRT))
plot(km,col=1:62,xlab="Time (Days)",ylab="Fraction Progression Free")
final<-final[is.na(final$TRT)==F,]
# fit time-series
pfs<-final[,c(1,3,5,6)]
all<-merge(all,pfs,by=c("ID","TRT"),all.x=T,all.y=T)
head(all)
all<-all[order(all[,7],all[,5]),]
all<-all[all$DAYS<=all$PDAY,]
all<-all[is.na(all$DIAM)==F,]

DAYS2<-all$DAYS[2:dim(all)[1]]
DAYS2<-c(DAYS2,NA)
all<-cbind(all,DAYS2)

all<-all[all$DAYS2!=0,]
# re-set PFS flag to 0 for all preceding time-points
PFS2<-array(NA,dim=c(dim(all)[1],1))
for(ii in 2:length(PFS2)){
  if(all$NID[ii]==all$NID[ii+1]){
    PFS2[ii]<-0
  }else{
    PFS2[ii]<-all$PFS[ii]
  }
}
# Error is simply an indexing issue and can be ignored
all<-cbind(all,PFS2)
all<-all[is.na(all$PFS)==F,]
# fill in 1st and last
all$PFS2[1]<-all$PFS[1]
all$PFS2[length(all$PFS2)]<-all$PFS[length(all$PFS2)]

# let's calculate relative change too, first include baseline
BSL<-all[all$DAYS==0,c("NID","DIAM")]
colnames(BSL)<-c("NID","BSL")
all<-merge(all,BSL,by=c("NID"),all.x=T,all.y=F)
all$PCH<-100*(all$DIAM-all$BSL)/(all$BSL)
# note that the fisrt % change is in fact missing
all$PCH[all$DAYS==0]<-NA

# next we can look at rate of change from visit to visit
DLDT<-array(NA,dim=c(length(all$NID),1))
for (ii in 2:length(DLDT)){
  if(all$NID[ii]==all$NID[ii-1]){
    DLDT[ii]=(all$DIAM[ii]-all$DIAM[ii-1])/(all$DAYS[ii]-all$DAYS[ii-1])
  }
}
all<-cbind(all,DLDT)
# record best change seen as we move through time



# quick general look - care about current size and how its changing:
summary(coxph(Surv(DAYS, DAYS2, PFS2) ~ cluster(NID)+DIAM,data=all))
summary(coxph(Surv(DAYS, DAYS2, PFS2) ~ cluster(NID)+c(DLDT*100),data=all))
summary(coxph(Surv(DAYS, DAYS2, PFS2) ~ cluster(NID)+c(DLDT*100)+DIAM,data=all))
# so this tells us that in general across the entire data-set using either
# one or combined correlates quite well with progression


# next we want to look at this in each drug
trt<-unique(all$TRT)
results<-data.frame(array(NA,dim=c(length(trt),5)))

for (ii in 1:length(trt)){
  results[ii,1]<-trt[ii]
  m1<-coxph(Surv(DAYS, DAYS2, PFS2) ~ cluster(NID)+DIAM,data=all[all$TRT==trt[ii],])
  results[ii,2]<-m1$coefficients[1]
  results[ii,3]<-1-pchisq(-2*m1$loglik[1]+2*m1$loglik[2],1)
  m1<-coxph(Surv(DAYS, DAYS2, PFS2) ~ cluster(NID)+c(DLDT*100),data=all[all$TRT==trt[ii],])
  results[ii,4]<-m1$coefficients[1]
  results[ii,5]<-1-pchisq(-2*m1$loglik[1]+2*m1$loglik[2],1)
}
colnames(results)<-c("TRT","DIAM","DIAM_P","DLDT","DLDT_P")
# histogram of HR values
hist((results$DIAM),xlab="log (Hazard Ratio)",main="Diameter")
abline(v=mean(results$DIAM),col=2,lwd=2)
abline(v=mean(results$DIAM)-1*sqrt(var(results$DIAM)),col=2,lty=2,lwd=2)
abline(v=mean(results$DIAM)+1*sqrt(var(results$DIAM)),col=2,lty=2,lwd=2)
sqrt(var(results$DIAM))/mean(results$DIAM)# 49% CV

hist((results$DLDT),xlab="log (Hazard Ratio)",main="DLDT")
abline(v=mean(results$DLDT),col=2,lwd=2)
abline(v=mean(results$DLDT)-1*sqrt(var(results$DLDT)),col=2,lty=2,lwd=2)
abline(v=mean(results$DLDT)+1*sqrt(var(results$DLDT)),col=2,lty=2,lwd=2)
sqrt(var(results$DLDT))/mean(results$DLDT)# 42% CV

plot(results$DIAM,results$DLDT,pch=19,ylab="DLDT log (Hazard Ratio)",xlab="Diameter log (Hazard Ratio)")
# not a useful plot

 plot(exp(results$DIAM),-log10(results$DIAM_P),pch=19,xlim=c(1,6),
      xlab="Hazard Ratio",ylab="-log10(p-value)",main="Diameter")
 abline(h = -log10(0.05),col=2)
 100*sum(results$DIAM_P>=0.05)/length(results$DIAM_P)# 10 %
 results$TRT[results$DIAM_P>=0.05]# no trend
 100*sqrt(var(results$DIAM[results$DIAM_P<0.05]))/mean(results$DIAM[results$DIAM_P<0.05])# 45% CV
 
 plot(exp(results$DLDT),-log10(results$DLDT_P),pch=19,xlab="Hazard Ratio",ylab="-log10(p-value)",
      log="x",main="d(Diameter)/dt",xlim=c(1,1.15))
 abline(h = -log10(0.05),col=2)
 100*sum(results$DLDT_P>=0.05)/length(results$DLDT_P)# 5 %
 results$TRT[results$DLDT_P>=0.05]# no trend
 100*sqrt(var(results$DLDT[results$DLDT_P<0.05]))/mean(results$DLDT[results$DLDT_P<0.05])# 39% CV


# so it is generally prognostic, now we do the does it capture treatment effect stuff...

# cycle through all possible 2 arm trials and record HR and ratio mean
# % change
idx<-unique(final$TRT)
results<-{}
count<-1
for (ii in 1:length(idx)){
  for(jj in (ii+1):length(idx)){
    dummy<-final[final$TRT==idx[ii]|final$TRT==idx[jj],]
    res1<-summary(coxph(Surv(PDAY,PFS)~c(TRT==idx[jj]),data=dummy))$coefficients[2]
    res2<-(survdiff(Surv(PDAY,PFS)~c(TRT==idx[jj]),data=dummy))$chisq
    d1<-summary(coxph(Surv(DAYS, DAYS2, PFS2) ~ cluster(NID)+DIAM,
                      data=all[all$TRT==idx[ii]|all$TRT==idx[jj],]))
    d2<-summary(coxph(Surv(DAYS, DAYS2, PFS2) ~ cluster(NID)+DIAM+TRT,
                      data=all[all$TRT==idx[ii]|all$TRT==idx[jj],]))
    res3<- d2$logtest[1] - d1$logtest[1]
    d1<-summary(coxph(Surv(DAYS, DAYS2, PFS2) ~ cluster(NID)+DLDT,
                      data=all[all$TRT==idx[ii]|all$TRT==idx[jj],]))
    d2<-summary(coxph(Surv(DAYS, DAYS2, PFS2) ~ cluster(NID)+DLDT+TRT,
                      data=all[all$TRT==idx[ii]|all$TRT==idx[jj],]))
    res4<- d2$logtest[1] - d1$logtest[1]
    d1<-summary(coxph(Surv(DAYS, DAYS2, PFS2) ~ cluster(NID)+DIAM+DLDT,
                      data=all[all$TRT==idx[ii]|all$TRT==idx[jj],]))
    d2<-summary(coxph(Surv(DAYS, DAYS2, PFS2) ~ cluster(NID)+DIAM+DLDT+TRT,
                      data=all[all$TRT==idx[ii]|all$TRT==idx[jj],]))
    res5<- d2$logtest[1] - d1$logtest[1]
    results<-rbind(results,c(res1,res2,res3,res4,res5,idx[ii],idx[jj]))
    print(results[count,])
    count<-count+1
  }
}
# An error may occur but this can be ignored

results<-data.frame(results)
results[,1]<-as.numeric(as.character(results[,1]))
results[,2]<-as.numeric(as.character(results[,2]))
results[,3]<-as.numeric(as.character(results[,3]))
results[,4]<-as.numeric(as.character(results[,4]))
results[,5]<-as.numeric(as.character(results[,5]))


colnames(results)<-c("PFS_HR","LR_CHISQ","PC_DIAM","PC_DLDT","PC_DIAM_DLDT","D1","D2")

length(results$PFS_HR[results$LR_CHISQ>3.84])# 1103
length(results$PFS_HR)# 1830
hist(log10(results$PFS_HR),breaks=25,main="",xlab="log10 (Hazard Ratio)",ylim=c(0,250))
hist(log10(results$PFS_HR[results$LR_CHISQ>3.84]),breaks=25,main="",xlab="log10 (Hazard Ratio)",ylim=c(0,250))
hist(log10(results$PFS_HR[results$LR_CHISQ>6.64]),breaks=25,main="",xlab="log10 (Hazard Ratio)",ylim=c(0,250))

hist((results$PC_DIAM[results$LR_CHISQ>3.84]),breaks=25,main="Diameter",xlab="LR test-statistic",ylim=c(0,400))
hist((results$PC_DIAM[results$LR_CHISQ>6.64]),breaks=25,main="Diameter",xlab="LR test-statistic",ylim=c(0,400))

hist((results$PC_DLDT[results$LR_CHISQ>3.84]),breaks=25,main="DLDT",xlab="LR test-statistic",ylim=c(0,400))
hist((results$PC_DLDT[results$LR_CHISQ>6.64]),breaks=25,main="DLDT",xlab="LR test-statistic",ylim=c(0,400))

hist((results$PC_DIAM_DLDT[results$LR_CHISQ>3.84]),breaks=25,main="Diameter + DLDT",xlab="LR test-statistic",ylim=c(0,400))
hist((results$PC_DIAM_DLDT[results$LR_CHISQ>6.64]),breaks=25,main="Diameter + DLDT",xlab="LR test-statistic",ylim=c(0,400))


sum(results$PC_DIAM[results$LR_CHISQ>3.84]>3.84)/
  length(results$PC_DIAM[results$LR_CHISQ>3.84])# 77%
sum(results$PC_DLDT[results$LR_CHISQ>3.84]>3.84)/
  length(results$PC_DLDT[results$LR_CHISQ>3.84])# 83%
sum(results$PC_DIAM_DLDT[results$LR_CHISQ>3.84]>3.84)/
  length(results$PC_DIAM_DLDT[results$LR_CHISQ>3.84])# 71%

sum(results$PC_DIAM[results$LR_CHISQ>3.84]>6.64)/
  length(results$PC_DIAM[results$LR_CHISQ>3.84])# 65%
sum(results$PC_DLDT[results$LR_CHISQ>3.84]>6.64)/
  length(results$PC_DLDT[results$LR_CHISQ>3.84])# 64%
sum(results$PC_DIAM_DLDT[results$LR_CHISQ>3.84]>6.64)/
  length(results$PC_DIAM_DLDT[results$LR_CHISQ>3.84])# 58%


# so we get improvements when we look at the change and direction than either one alone
# but the caveat is that the regression coefficients are more variable as we saw in the 
# single arm analysis

plot(results$PFS_HR[results$LR_CHISQ>3.84],
     -log10(1-pchisq(results$LR_CHISQ[results$LR_CHISQ>3.84],1)),log="x",
     xlab="PFS HR",ylab="-log10(p-value)",pch=19,cex=0.5,ylim=c(0,17),col="green")
abline(h = -log10(0.05),col=1)
points(results$PFS_HR[results$LR_CHISQ<=3.84],
     -log10(1-pchisq(results$LR_CHISQ[results$LR_CHISQ<=3.84],1)),pch=19,col=2,cex=0.5)


plot(results$PFS_HR[results$LR_CHISQ>3.84],
     -log10(1-pchisq(results$PC_DIAM_DLDT[results$LR_CHISQ>3.84],1)),log="x",
     xlab="PFS HR",ylab="-log10(p-value)",pch=19,cex=0.5,ylim=c(0,17),col="green")
abline(h = -log10(0.05),col=1)
points(results$PFS_HR[results$LR_CHISQ>3.84 & results$PC_DIAM_DLDT<=3.84],
     -log10(1-pchisq(results$PC_DIAM_DLDT[results$LR_CHISQ>3.84 & results$PC_DIAM_DLDT<=3.84],1)),
     col="red",pch=19,cex=0.5)

