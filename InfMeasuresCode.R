###Code used in "Using Concurrent Functional Regression to Identify Influential Functional Measurements"
#Ie ch 3
#Simulation code in separate files
library(fdaconcur)
library(fda)
library(pracma)
library(lubridate)
library(MASS) #for mvrnorm
library(corpcor) #for make positve definite
#Application: River Stage During Flood Events

#Calculates "B0t","B1t","stderr.B0t","stderr.B1t","Residuals","AdjHatMatrix","MeanHat","FullyHat"
BetaOnly<-function(xData,yData,BasisNum=11, BetaOnly=FALSE,lambda=10^(-1)){
  n<-nrow(xData)
  nPred<-ncol(xData) #Number of predictors (ie flood events)
  
  gaittime <- seq(1:n)
  gaitrange <- c(0,n)
  gaitfine = seq(0,n,.2)
  harmaccelLfd <- vec2Lfd(c(0, (2*pi/n)^2, 0), rangeval=gaitrange)
  gaitbasis <- create.fourier.basis(gaitrange, nbasis=BasisNum)
  
  mygaitExp <- array(NA, dim = c(n,nPred,2))
  mygaitExp[1:n, ,] <- seq(1:n)
  
  for(i in 1:nPred){
    mygaitExp[,i, 1] <- xData[,i]
    mygaitExp[,i, 2] <- yData[,i]
  }
  
  gaitSmooth = smooth.basisPar(gaittime, mygaitExp, gaitbasis, Lfdobj=harmaccelLfd, lambda=lambda) #Could Edit this part Original played around no different
  gaitfd = gaitSmooth$fd
  
  names(gaitfd$fdnames) = c("Normalized time", "Event", "Height")
  gaitfd$fdnames[[3]] = c("Cong", "Cedar")
  
  congfd  = gaitfd[,1] #Predictor Stuff
  cedfd = gaitfd[,2] # Response Stuff This seems like a more detailed version of what we get from the other textbook.
  
  xfdlist   = list(const=rep(1,nPred), cong=congfd)
  betafdPar = fdPar(gaitbasis, harmaccelLfd)
  betalist  = list(const=betafdPar, cong=betafdPar)
  
  
  gaitRegress= fRegress(cedfd, xfdlist, betalist)
  betaestlist = gaitRegress$betaestlist
  cedIntercept = predict(betaestlist$const$fd, gaitfine) #B0
  congCoef = predict(betaestlist$cong$fd, gaitfine) #Slope Term B1
  
  Intercept<-approx(x=1:length(cedIntercept), y=cedIntercept, n = n)$y
  Slope<-approx(x=1:length(congCoef), y=congCoef , n = n)$y
  
  if(BetaOnly==FALSE){
    cedhatfd = gaitRegress$yhatfd$fd #
    cedhatmat = eval.fd(gaittime, cedhatfd)  #1 to 2000 compare this to yHat.t.mat
    resmat. = mygaitExp[,,2] - cedhatmat #The 2 represents cedar (response)
    SigmaE = cov(t(resmat.))
    
    cedfinemat   = eval.fd(gaitfine, cedfd)
    cedmeanvec   = eval.fd(gaitfine, mean.fd(cedfd))
    cedhatfinemat= eval.fd(gaitfine, cedhatfd)
    resmat        = cedfinemat - cedhatfinemat #finer grid than the . needed for computation below
    ncurve        = dim(mygaitExp)[2] #Says the number of flood events
    resmat0 = cedfinemat - cedmeanvec %*% matrix(1,1,ncurve)
    SSE0 = apply((resmat0)^2, 1, sum)
    SSE1 = apply(resmat^2, 1, sum)
    ced.R2 = (SSE0-SSE1)/SSE0
    
    ylim2=c(0, max(congCoef, ced.R2)) #And this is what I had plotted before exactly.
    gaitbasismat = eval.basis(gaitfine, gaitbasis)
    y2cMap = gaitSmooth$y2cMap #Needed to get error
    
    cedhatfd = gaitRegress$yhatfd$fd #
    cedhatmat = eval.fd(gaittime, cedhatfd)  #1 to 2000 compare this to yHat.t.mat
    resmat. = mygaitExp[,,2] - cedhatmat #The 2 represents cedar (response)
    SigmaE = cov(t(resmat.))
    
    fRegressList1 = fRegress(cedfd, xfdlist, betalist,
                             y2cMap=y2cMap, SigmaE=SigmaE)
    
    fRegressList2 = fRegress.stderr(fRegressList1, y2cMap, SigmaE)
    betastderrlist = fRegressList2$betastderrlist
    
    error.B0<-predict(betastderrlist[[1]], gaitfine)
    error.B1<-predict(betastderrlist[[2]], gaitfine)
    #This gives the standard errors for B0t and B1t to use for DFBETAS
    stderr.B0<-approx(x=1:length(error.B0), y=error.B0, n = n)$y
    stderr.B1<-approx(x=1:length(error.B1), y=error.B1, n = n)$y
    
    
    #Adjust to be nRow Ncol
    meanHat_i<-c()
    data1<-cbind(rep(1,nPred),t(xData)) #adds the column of 1's
    Hatiit<-matrix(NA,nrow=nPred,ncol=n)
    for(t in 2:(n+1)){
      X<-data1[,c(1,t)]
      Hat<-X%*%solve(t(X)%*%X)%*%t(X)
      Hatiit[,t-1]<-diag(Hat) #Each row is an event. The corresponding column is the Hii of that event at each timepoint
    }
    meanHat_i<-c(meanHat_i,apply(Hatiit,1,mean))#will return an average Hii for each of the 10 sampled curve
    
    out<-list(Intercept,Slope,stderr.B0,stderr.B1,resmat.,Hatiit,meanHat_i,cedhatmat)
    names(out)<-c("B0t","B1t","stderr.B0t","stderr.B1t","Residuals","AdjHatMatrix","MeanHat","FullyHat")
  }
  else {out<-list(Intercept,Slope)
  names(out)=c("B0t","B1t")
  }
  return(out)
}

FullResults<-BetaOnly(FinalXtL1star,FinalYtL1star,BasisNum=11)
B0_i_<-matrix(NA,nrow=nrow(FinalXtL1star),ncol=ncol(FinalXtL1star))
B1_i_<-matrix(NA,nrow=nrow(FinalXtL1star),ncol=ncol(FinalXtL1star))
stderrB1_i_2<-matrix(NA,nrow=nrow(FinalXtL1star),ncol=ncol(FinalXtL1star))
stderrB0_i_2<-matrix(NA,nrow=nrow(FinalXtL1star),ncol=ncol(FinalXtL1star))

#Calculates "B0t","B1t","stderr.B0t","stderr.B1t","Residuals","AdjHatMatrix","MeanHat","FullyHat"
#With each of the i events left out using the same number of Fourier basis functions as full run with all 10
for(i in 1:10){
  BasisToUse<-11
  tempBetas<-BetaOnly(FinalXtL1star[,-i],FinalYtL1star[,-i],BetaOnly=FALSE,BasisNum=BasisToUse)
  B0_i_[,i]<-tempBetas$B0t
  B1_i_[,i]<-tempBetas$B1t
  stderrB0_i_2[,i]<-tempBetas$stderr.B0t
  stderrB1_i_2[,i]<-tempBetas$stderr.B1t
}

DFBETAS011<-matrix(NA,nrow=nrow(FinalXtL1star),ncol=ncol(FinalXtL1star))
DFBETAS111<-matrix(NA,nrow=nrow(FinalXtL1star),ncol=ncol(FinalXtL1star))
DFBETA011<-matrix(NA,nrow=nrow(FinalXtL1star),ncol=ncol(FinalXtL1star))
DFBETA111<-matrix(NA,nrow=nrow(FinalXtL1star),ncol=ncol(FinalXtL1star))
for(i in 1:10){
  DFBETAS011[,i]<-(FullResults$B0t-B0_i_[,i])/stderrB0_i_2[,i]
  DFBETAS111[,i]<-(FullResults$B1t-B1_i_[,i])/stderrB1_i_2[,i]
  DFBETA011[,i]<-(FullResults$B0t-B0_i_[,i])
  DFBETA111[,i]<-(FullResults$B1t-B1_i_[,i])
}

#all 10 Plots B_p(t) vs B_p(i)(t):
for(i in 1:10){
  par(mfrow=c(1,1))
  par(mar=c(5.6,4.6,4.1,2.1))
  plot(B0_i_[,i],type = "l",lwd = 2,col=2,ylab = expression(paste(beta [0](t)," Estimate")),ylim=c(-50,10),sub = paste("Event",i),main=expression(paste("Full ", hat(beta) [0](t)," vs. ", hat(beta)[0][(i)](t))))
  lines(FullResults$B0t,col=1,lwd=2)
  legend("bottomright", legend = c("Using All 10 Events", "With event i Removed"),
         lty = 1,col = c(1,2))
}
for(i in 1:10){
  par(mfrow=c(1,1))
  par(mar=c(5.6,4.6,4.1,2.1))
  plot(B1_i_[,i],type = "l",lwd = 2,col=2,ylab = expression(paste(beta [1](t)," Estimate")),ylim=c(-2,4),sub = paste("Event",i),main=expression(paste("Full ", hat(beta) [1](t)," vs. ", hat(beta)[1][(i)](t))))
  lines(FullResults$B1t,col=1,lwd = 2,)
  legend("bottomright", legend = c("Using All 10 Events", "With event i Removed"),
         lty = 1,col = c(1,2))
}
#plot for February 2020 Event
i=10
plot(B0_i_[,i],type = "l",col=2,lwd=3,lty=2,ylab = expression(paste(beta [0](t)," Estimate")),ylim=c(-50,10),sub = "February 2020",main=expression(paste("Full ", hat(beta) [0](t)," vs. ", hat(beta)[0][(i)](t))))
lines(FullResults$B0t,col=1,lwd=3)
legend("bottomright", legend = c("Using All 10 Events", "With February 2020 event Removed"),
       lty = c(1,2),col = c(1,2),lwd=c(2,2))

plot(B1_i_[,i],type = "l",col=2,lwd=3,lty = 2,ylab = expression(paste(beta [1](t)," Estimate")),ylim=c(-2,4),sub = "February 2020",main=expression(paste("Full ", hat(beta) [1](t)," vs. ", hat(beta)[1][(i)](t))))
lines(FullResults$B1t,col=1,lwd=3)
legend("bottomright", legend = c("Using All 10 Events", "With February 2020 event Removed"),
       lty = c(1,2),col = c(1,2),lwd=c(2,2))


#ALL DFBETAS(t)
plot(DFBETAS011[,1],type="l",lwd=2,yaxt="n",main = expression("All events' DFBETAS"[0]), ylab= expression("DFBETAS"[0]), ylim = c(-3,3))
for(i in 2:10){
  lines(DFBETAS011[,i],col=i,lwd=2)
}
abline(h=-2/sqrt(10), lty=2,lwd=1.5)
abline(h=2/sqrt(10), lty=2,lwd=1.5)
x<-c(-3,-2,-0.63,0,0.63,2,3)
axis(2,at=x,labels=x)

plot(DFBETAS111[,1],type="l",lwd=2,yaxt="n",main = expression("All events' DFBETAS"[1]), ylab= expression("DFBETAS"[1]), ylim = c(-3,3))
for(i in 2:10){
  lines(DFBETAS111[,i],col=i,lwd=2)
}
abline(h=-2/sqrt(10), lty=2,lwd=1.5)
abline(h=2/sqrt(10), lty=2, lwd=1.5)
x<-c(-3,-2,-0.63,0,0.63,2,3)
axis(2,at=x,labels=x)
#Taking Means of DFbetas
round(apply(abs(DFBETAS011),2,mean),3) #Changed to be the one i wanted
round(apply(abs(DFBETAS111),2,mean),3) #Changed to be the one i wanted


##Cooks Distance
MSE_t<-apply((FullResults$Residuals)^2,1,mean) #Gives MSE at each t

##Need the predicted values of Yj when i is left out.
yHat_mini<-matrix(NA,nrow=nrow(FinalXtL1star),ncol=ncol(FinalXtL1star))
for(i in 1:10){
  yHat_mini[,i]<-B0_i_[,i]+B1_i_[,i]*FinalXtL1star[,i]
}

CooksD_ti<-matrix(0:0, nrow=nrow(FinalXtL1star),ncol=ncol(FinalXtL1star))
for(i in 1:10){
  #Get predicted value of each event using B0 and B1 calculated without i
  yHatj_<-matrix(NA,nrow=nrow(FinalXtL1star),ncol=ncol(FinalXtL1star))
  for(j in 1:10){
    yHatj_[,j]<-B0_i_[,i]+B1_i_[,i]*FinalXtL1star[,j] #B0_i_[,i]+B1_i_[,i] calculation of B0 and B1 when i is not included
  }
  for(t in 1:nrow(FinalXtL1star)){
    CooksD_ti[t,i]<-sum((as.numeric(FullResults$FullyHat[t,i])-yHatj_[t,i])^2)/((1+1)*MSE_t[t])
  }
}

Cookscutoff<-qf(.5,2,8-2) #.5 percentile, and 2 is number of parameters
for(i in 1:10){
  plot(CooksD_ti[,i],lwd=2,type="l",main=paste("Cook's Distance for event", i), ylim=c(0,10),ylab = "Cook's Distance")
  abline(h=Cookscutoff,lty =2, lwd = 2 )
  # print(sum(CooksD_ti[,i]>0.7568285)/2000)
}

#Cooks D plot
plot(CooksD_ti[,10],lwd=2,type="l",main=paste("Cook's Distance for the February 2020 Event"), ylim=c(0,20),ylab = "Cook's Distance")
abline(h=Cookscutoff,lty =2,lwd=2)

plot(CooksD_ti[,1],lwd=2,type="l",main=paste("Cook's Distance for the August 1995 Event"), ylim=c(0,20),ylab = "Cook's Distance")
abline(h=Cookscutoff,lty =2,lwd=2)

#Cooks D Mean
round(apply(CooksD_ti,2,mean),3)

##DFFITS
yHat_mini<-matrix(NA,nrow=nrow(FinalXtL1star),ncol=ncol(FinalXtL1star))
for(i in 1:10){
  yHat_mini[,i]<-B0_i_[,i]+B1_i_[,i]*FinalXtL1star[,i]
}

DFFIT<-FullResults$FullyHat-yHat_mini

#For MSE_mini minus i we want the true-predicted at each t
#for all 9 events when the 10th event is left out
MSE_mini_t<-matrix(NA,nrow=nrow(FinalXtL1star),ncol=ncol(FinalXtL1star))
yHatsNoi<-matrix(NA,nrow=nrow(FinalXtL1star),ncol=ncol(FinalXtL1star)-1)
n=1:10
for(i in n){
  k=1
  for(j in n[-i]){
    yHatsNoi[,k]<-B0_i_[,i]+B1_i_[,i]*FinalXtL1star[,j] #goes through all but i using the betas with i left out
    k=k+1
  }
  MSE_mini_t[,i]<-apply((FinalYtL1star[,-i]-yHatsNoi)^2,1,mean) #Get an MSE across the other events when event i is left out
}
DFFITS<-matrix(NA,nrow=nrow(DFFIT),ncol=ncol(DFFIT))
for(t in 1:nrow(DFFIT)){
  DFFITS[t,]<-DFFIT[t,]/sqrt(MSE_mini_t[t,]*FullResults$AdjHatMatrix[,t])
}

par(mfrow=c(1,1))
for(i in 1:10){
  plot(DFFITS[,i],lwd=2,type="l", main = paste("DFFITS for event ", i),ylim=c(-15,15),ylab="DFFITS")
  abline(h=1, lty=2,lwd=2)
  abline(h=-1, lty=2,lwd=2)
  #print(sum(DFFITS[,i]>1)/2000)
}

plot(DFFITS[,1],type="l", ylim=c(-15,15),lwd=2,main = "DFFITS(t) for the August 1995 Event", ylab = "DFFITS")
abline(h=1, lty=2,lwd=2)
abline(h=-1, lty=2,lwd=2)

plot(DFFITS[,10],lwd=2,ylim=c(-15,15),type="l", main = "DFFITS(t) for the February 2020 Event", ylab = "DFFITS")
abline(h=1, lty=2,lwd=2)
abline(h=-1, lty=2,lwd=2)

#Mean |DFFITS(t)| for each observation
round(apply(abs(DFFITS),2,mean),3)

OUonFloodEM<-function(sigma=1,drift=1,times,endtime,Data){
  dt<-endtime/times; #want this to be around .5 for most situations to keep it stable enough
  objects<-ncol(Data) #number of flood events I have
  normalvec<-matrix(rnorm(objects*times,sd=1),nrow=objects,ncol=times); #this is randn from matlab code
  ouvalue<-matrix(0:0,nrow=objects,ncol=times) #null matrix of 0s
  
  # Defining the cluster sizes and the true clustering structure: No clusters here so just 1
  timevec<-seq(0,endtime,length=times)
  
  # loop to create approximation of Ornstein-Uhlenbeck process	# Using Euler-Maruyama approximation
  ouvalue[1:objects,1]=normalvec[,1] #So we don't start at 0.
  for (timeint in 2:times)
  {
    ouvalue[1:objects,timeint]<- ouvalue[1:objects,timeint-1]-drift*ouvalue[1:objects,timeint-1]*dt + sigma*normalvec[1:objects,timeint]*sqrt(dt) #normalvec[1:clussizes[1],timeint] n random values from normal(0, sigma)
  }
  
  samesizeOU<-matrix(0:0,nrow=nrow(Data),ncol=objects) #makes sure it is the same size as the vector it is being applied to
  for(i in 1:objects){
    samesizeOU[,i]<-approx(x=1:ncol(ouvalue),y = ouvalue[i,], n=nrow(Data))$y
  }
  # Defining and naming the observed curves
  obscurves<-Data+samesizeOU;
  out<-list(1:nrow(Data),obscurves,samesizeOU)
  names(out)<-c("timevec","obscurves","samesizeOU")
  return(out)
}

samplecurves<-function(xdata, ydata, prob_i, size=10){
  #generate 10 numbers independently.
  index<-sample(seq(1,size), size = size, prob = prob_i, replace = TRUE)
  output10x<-xdata[,index]
  output10y<-ydata[,index]
  out<-list(output10x,output10y, index)
  names(out)=c("xData", "yData", "index")
  return(out) #need y matrix
}

##Bootstrapping with means of Influence Measures
##Do the bootstrapping here with different Alphas
##Do I want perturbations. Yes
OUprocess_Metrics<-function(xData ,yData , basisList,PredVec, metric,L2=TRUE,prob_i,sigma=3, drift=.9,times=2000,endtime = 1000,BasisList){
  #Step 1: use Sim10 function to grab the bootstrapped events wanted
  tempdata<-samplecurves(xData, yData, prob_i = prob_i, size =ncol(xData)) #Samples 10 curves from the sample with replacement
  #Step 2: Apply the OH process to allow the response curves to have more variability
  OUyData<- OUonFloodEM(sigma=sigma,drift=drift,times=times, endtime = endtime,Data=tempdata$yData)
  newY<-matrix(NA,nrow=nrow(yData),ncol=ncol(yData))
  
  for(i in 1:ncol(yData)){
    newY[,i]<-supsmu(x=OUyData$timevec,y=OUyData$obscurves[,i])$y #generically smooth the OU data to look more like my data
  }
  
  #Step 3: Calculate "Full" Basis numbers needed for the iteration will use what was used for the original data
  #fullbasisNum<-L2bestEst(tempdata$xData, newY, basisList=BasisList)$`Best Basis` #chooses the best number of basis functions to use
  #fullbasisNum<-as.integer(names(fullbasisNum))
  fullbasisNum=BasisList
  
  #Step 4: Calculate Full model results getting B0, B1, Full prediction etc.
  #It's all saved under "FullResults"
  FullResults<-BetaOnly(tempdata$xData, newY,BasisNum=fullbasisNum)
  FullReconstructedY<-FullResults$B0t+FullResults$B1t*PredVec
  
  #Step 5: If L2 or L1 run the leave 1 out model and find the
  #L1 or L2 distance for each left out
  if(metric == "deltaP"){
    if(L2==TRUE){
      tempEachLeftOut<-Leave1OutFixedBasis(PredictorMat=tempdata$xData, ResponseMat=newY, PredictorVec=PredVec,  basisnumberused = fullbasisNum)
      tempL2<-c()
      for(j in 1:ncol(xdata)){
        squaredifference<-(FullReconstructedY - tempEachLeftOut$`Predicted Responses Without i`[,j])^2
        tempL2[j]<-sqrt(trapz(x=1:length(squaredifference), y = squaredifference))
      }
      L_distanceout<-tempL2
    }
    else{
      tempEachLeftOut<-Leave1OutFixedBasis(PredictorMat=tempdata$xData, ResponseMat=newY, PredictorVec=PredVec,  basislist = fullbasisNum)
      tempL1<-c()
      for(j in 1:ncol(xData)){
        absdifference<-abs(FullReconstructedY - tempEachLeftOut$`Predicted Responses Without i`[,j])
        tempL1[j]<-(trapz(x=1:length(absdifference), y = absdifference))
      }
      L_distanceout<-tempL1
    }
    
    #metric can be DFFBETAS0, DFFBETAS1, DFFITS, CooksD, deltaP
    #Does 1 iteration of the 
    output_metric =  L_distanceout
  }
  #Step 6  Calculate DFBETAs0 and DFBETAs1 for each event at each point return this matrix
  #6b Also return the mean abs difference of these for each event
  #Get Beta0t and Beta1t when each event is leftout
  if(metric == "DFBETAS0" || metric == "DFBETAS1" || metric == "DFFITS" || metric =="CooksD"){
    
    B0_i_<-matrix(NA,nrow=nrow(xData),ncol=ncol(xData))
    B1_i_<-matrix(NA,nrow=nrow(xData),ncol=ncol(xData))
    propB0<-c() #gives proportion of sample that is larger than cutoff
    propB1<-c()
    stderrB1_i_<-matrix(NA,nrow=nrow(xData),ncol=ncol(xData)) #std errort B0 and B1 each i left out
    stderrB0_i_<-matrix(NA,nrow=nrow(xData),ncol=ncol(xData))
    DFBETAS0<-matrix(NA,nrow=nrow(xData),ncol=ncol(xData))
    DFBETAS1<-matrix(NA,nrow=nrow(xData),ncol=ncol(xData))
    meanAbsB0<-c()
    meanAbsB1<-c()
    for(i in 1:ncol(xData)){
      #Optimizing number of basis functions for the 9 events in case that matters
      tempBetas<-BetaOnly(tempdata$xData[,-i],newY[,-i],BetaOnly=FALSE,BasisNum=fullbasisNum)
      B0_i_[,i]<-tempBetas$B0t
      B1_i_[,i]<-tempBetas$B1t
      stderrB0_i_[,i]<-tempBetas$stderr.B0t
      stderrB1_i_[,i]<-tempBetas$stderr.B1t
      DFBETAS0[,i]<-(FullResults$B0t-B0_i_[,i])/stderrB0_i_[,i]
      DFBETAS1[,i]<-(FullResults$B1t-B1_i_[,i])/stderrB1_i_[,i]
      meanAbsB0[i]<-mean(abs(DFBETAS0[,i]))
      meanAbsB1[i]<-mean(abs(DFBETAS1[,i]))
    }
    if(metric == "DFBETAS0"){
      output_metric = meanAbsB0
    }
    
    if(metric == "DFBETAS1"){
      output_metric = meanAbsB1
    }
    
  }
  
  #Step 7 Calculate the Hat Matrix Can probably use what is in step 3 (YES)
  #FullResults$AdjHatMatrix is the pointwise hat matrix
  #FullResults$MeanHat
  if(metric == "DFFITS"){
    sigHat<-c()
    for(i in 1:ncol(xData)){
      sigHat[i]<-sum(FullResults$AdjHatMatrix[i,]>=.4)/nrow(xData)
    }
    
    #Step 8  Calculate DFFITs for each event at each point return this matrix
    yHat_mini<-matrix(NA,nrow=nrow(xData),ncol=ncol(xData))
    for(i in 1:ncol(xData)){
      yHat_mini[,i]<-B0_i_[,i]+B1_i_[,i]*tempdata$xData[,i] #yHat for i when ith event is left out. B1_i and B0_i calculated before
    }
    DFFIT<-FullResults$FullyHat-yHat_mini
    #8b Also return the mean abs difference of these for each event
    MSE_mini_t<-matrix(NA,nrow=nrow(tempdata$xData),ncol=ncol(tempdata$xData))
    yHatsNoi<-matrix(NA,nrow=nrow(tempdata$xData),ncol=ncol(tempdata$xData)-1)
    n=1:ncol(tempdata$xData) #1:10
    for(i in n){
      k=1
      for(j in n[-i]){
        yHatsNoi[,k]<-B0_i_[,i]+B1_i_[,i]*tempdata$xData[,j] #goes through all but i using the betas with i left out
        k=k+1
      }
      MSE_mini_t[,i]<-apply((newY[,-i]-yHatsNoi)^2,1,mean) #Get an MSE across the other events when event i is left out
    }
    DFFITS<-matrix(NA,nrow=nrow(DFFIT),ncol=ncol(DFFIT))
    for(t in 1:nrow(DFFIT)){
      DFFITS[t,]<-DFFIT[t,]/sqrt(MSE_mini_t[t,]*FullResults$AdjHatMatrix[,t])
    }
    MeanAbsDFFITS<-apply(abs(DFFITS),2,mean)
    output_metric = MeanAbsDFFITS
    
  }
  
  #Step 9 Calculate Cooks D  and mean across each event and return cuttoff criteria based on F dist
  
  if(metric == "CooksD"){
    MSE_t<-apply((FullResults$Residuals)^2,1,mean)
    CooksD_ti<-matrix(0:0, nrow=nrow(tempdata$xData),ncol=ncol(tempdata$xData))
    for(i in 1:ncol(tempdata$xData)){
      #Get predicted value of each event using B0 and B1 calculated without i
      yHatj_<-matrix(NA,nrow=nrow(tempdata$xData),ncol=ncol(tempdata$xData))
      for(j in 1:ncol(tempdata$xData)){
        yHatj_[,j]<-B0_i_[,i]+B1_i_[,i]*tempdata$xData[,j] #B0_i_[,i]+B1_i_[,i] calculation when i is not included
      }
      for(t in 1:nrow(xData)){
        CooksD_ti[t,i]<-sum((as.numeric(FullResults$FullyHat[t,i])-yHatj_[t,i])^2)/((1+1)*MSE_t[t])
      }
    }
    
    MeanCooksD_i<-apply((CooksD_ti),2,mean)
    
    output_metric = MeanCooksD_i
  }
  
  out<-list(tempdata$xData,newY,fullbasisNum,tempdata$index,output_metric)
  names(out)<-c("xDataUsed","yDataUsed","fullbasisNum","IndicesUsed","output_metric")
  return(out)
}

meanDFB0<-apply(abs(DFBETAS011),2,mean)
meanDFB1<-apply(abs(DFBETAS111),2,mean)
meanDFFITS<-round(apply(abs(DFFITS),2,mean),3)
meanCD<-round(apply(abs(CooksD_ti),2,mean),3)

#Set Alpha and number of iterations
#Cooks D
CD_results_alpha0=c()
alpha = 0.5
#obs_mean_abscd<-apply(abs(meanDFB0),2,mean)
theta<-as.vector((1/meanCD)^alpha/sum((1/meanCD)^alpha))
nBoot=500
for(i in 1:nBoot){
  sigma = runif(1,3,5)
  drift = runif(1,0.5,1)
  test=OUprocess_Metrics(xData=FinalXtL1star, yData = FinalYtL1star,metric="CooksD",sigma = sigma, drift = drift, PredVec = Oct15CongHt$CongHt,prob_i=theta,times=2000,endtime = 1000,BasisList=11)
  CD_results_alpha0=c(CD_results_alpha0,test$output_metric)
}
quantile(CD_results_alpha0,c(.9,.95,.99))


#DFFITS
DFFIT_results_alphapoint1=c()
alpha = 0
theta<-as.vector((1/meanDFFITS)^alpha/sum((1/meanDFFITS)^alpha))
nBoot=10
for(i in 1:nBoot){
  sigma = runif(1,3,5)
  drift = runif(1,0.5,1)
  test=OUprocess_Metrics(xData=FinalXtL1star, yData = FinalYtL1star,metric="DFFITS",sigma = sigma, drift = drift, PredVec = Oct15CongHt$CongHt,prob_i=theta,times=2000,endtime = 1000,BasisList=11)
  DFFIT_results_alphapoint1=c(DFFIT_results_alphapoint1,test$output_metric)
}
quantile(DFFIT_results_alphapoint1,c(.9,.95,.99))
#DFFIT_results_alpha0

#DFBetas0 and 1
obs_mean_absdfbetas0<-apply(abs(DFBETAS011),2,mean)
theta<-as.vector((1/mean_abs_dfbetas0)^alpha/sum((1/mean_abs_dfbetas0)^alpha))
obs_mean_absdfbetas1<-apply(abs(DFBETAS111),2,mean)

DFBETA0boot=c()
alpha = 0
theta<-as.vector((1/obs_mean_absdfbetas0)^alpha/sum((1/obs_mean_absdfbetas0)^alpha))
nBoot=10
for(i in 1:nBoot){
  sigma = runif(1,3,5)
  drift = runif(1,0.5,1)
  test=OUprocess_Metrics(xData=FinalXtL1star, yData = FinalYtL1star,metric="DFBETAS0",sigma = sigma, drift = drift, PredVec = Oct15CongHt$CongHt,prob_i=theta,times=2000,endtime = 1000,BasisList=11)
  DFBETA0boot=c(DFBETA0boot,test$output_metric)
}
quantile(DFBETA0boot,c(.9,.95,.99))

DFBETA1boot=c()
alpha = 0
theta<-as.vector((1/obs_mean_absdfbetas1)^alpha/sum((1/obs_mean_absdfbetas1)^alpha))
nBoot=10
for(i in 1:nBoot){
  sigma = runif(1,3,5)
  drift = runif(1,0.5,1)
  test=OUprocess_Metrics(xData=FinalXtL1star, yData = FinalYtL1star,metric="DFBETAS1",sigma = sigma, drift = drift, PredVec = Oct15CongHt$CongHt,prob_i=theta,times=2000,endtime = 1000,BasisList=11)
  DFBETA1boot=c(DFBETA1boot,test$output_metric)
}
quantile(DFBETA1boot,c(.9,.95,.99))


##Application: Air and Water Temperature 
library(mapproj)
map(database= "world", ylim=c(20,70), xlim=c(-170,-30), col="grey80", fill=TRUE, projection="gilbert", orientation= c(90,0,225))
lon <- c(-68.204 ,-71.050 , -70.246 ,-73.181 ,-76.039 ,-74.418 ,
         -75.548 ,-76.671 ,-77.786 ,-79.924 ,-80.903 ,-81.465 ,
         -80.034 ,-81.807,-82.553 , -82.832 ,-85.880 ,-89.326 ,
         -91.338 ,-93.343 , -118.500 ,-120.754 ,-124.498 ,
         -124.105 ,-123.441 ,-177.361 ,-145.752 ,
         -131.625,-135.327 ,-124.184, -122.039 ,-162.327,-157.79,-97.215,
         -164.067  )  #longitude vector. needs -
lat <- c(44.392 ,42.355 , 43.656 ,41.174 ,38.220 , 39.357 ,
         35.796 ,34.717 ,34.213,32.781 ,32.035 ,30.675 ,
         26.613 ,26.132,27.858 ,27.978 ,30.213 ,30.326 ,
         29.450 , 29.768 ,34.008 ,35.169 ,42.739 ,
         46.904 , 48.125  ,28.215  ,60.558 ,
         55.331 ,59.450 ,41.746 ,38.056,55.062, 21.433, 26.061,
         67.575
)  #latitude vector
coord <- mapproject(lon, lat, proj="gilbert", orientation=c(90, 0, 225))  #convert points to projected lat/long
points(coord, pch=20, cex=c(rep(1,34),2), col=c(rep("red",34),"blue"))  #plot converted points


##air water graphs
AirTempMat<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/MainAirTemp.csv")
AirTempMat=AirTempMat[,-1]
WaterTempMat<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/MainWaterTemp.csv")
WaterTempMat=WaterTempMat[,-1]
datetimeplot<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/datetimeforgraphs.csv")
datetimeplot<-datetimeplot[,-1]
rockport_air<-read.csv(file="/Users/Ryan/Desktop/Research Material/Air_Water_Temp/RawDatasets/AirWaterOnly/rockportair.csv")
rockport_air=rockport_air[,-1]

#All air and water temp plots
plot(x=as_datetime(datetimeplot),y=AirTempMat[,1],ylim = c(0,35),col = "red",type = "l", lwd = 2,ylab = "Temperature (\u00B0C)", main = "Air and Water Temperature: Amerada Pass, LA",xlab = "2020 Date")
lines(x=as_datetime(datetimeplot),y=WaterTempMat[,1],col="blue",lwd=2,lty = 2)
legend( x="bottom",
        legend=c("Air Temp.", "Water Temp"), col = c("red","blue"), lty = c(1,2),lwd=c(2,2))


plot(x=as_datetime(datetimeplot),y=AirTempMat[,35],ylim = c(0,35),col = "red",type = "l", lwd = 2,ylab = "Temperature (\u00B0C)", main = "Air and Water Temperature: Westport, WA ",xlab = "2020 Date")
lines(x=as_datetime(datetimeplot),y=WaterTempMat[,35],col="blue",lwd=2,lty = 2)
legend( x="bottom",
        legend=c("Air Temp.", "Water Temp"), col = c("red","blue"), lty = c(1,2),lwd=c(2,2))

plot(as_datetime(datetimeplot),AirTempMat[,1],lwd=2,type ="l", ylab = "Temperature (\u00B0C)", xlab = "Date", ylim = c(-30,50),main="Air Temperature")
for(i in 2:ncol(AirTempMat)){
  lines(as_datetime(datetimeplot),AirTempMat[,i],lwd=2,col = i)
}

plot(as_datetime(datetimeplot),WaterTempMat[,1],lwd=2,type ="l", ylab = "Temperature (\u00B0C)", xlab = "Date", ylim = c(-30,50),main="Water Temperature")
for(i in 2:ncol(AirTempMat)){
  lines(as_datetime(datetimeplot),WaterTempMat[,i],lwd=2,col = i)
}

#Full model example with Rockport
FullTempResults<-PredictFRegressNormTestUpdatedForRUpdate(AirTempMat,WaterTempMat, rockport_air,Basis = "BSpline", nBasis = 21, Plot = TRUE)
par(mfrow=c(1,1))
plot(FullTempResults$PredictedResponse,type = "l", ylim = c(0,40), ylab = "Temp C.",col = "blue",main = "Rockport: AirTemp (red) vs Pred. WaterTemp (blue)")
lines(FullTempResults$Lower,col="green")
lines(FullTempResults$Upper,col="green")
#lines(rockport_air,col = "Red")


#Calculating Observed Influence Measures:
BetaOnlyBSpline<-function(xData,yData,BasisNum=21, BetaOnly=FALSE,lambda=10^(-1)){
  n<-nrow(xData)
  nPred<-ncol(xData) #Number of predictors (ie flood events)
  
  gaittime <- seq(1:n)
  gaitrange <- c(0,n)
  gaitfine = seq(0,n,.2)
  
  mygaitExp <- array(NA, dim = c(n,nPred,2))
  mygaitExp[1:n, ,] <- seq(1:n)
  
  for(i in 1:nPred){
    mygaitExp[,i, 1] <- xData[,i]
    mygaitExp[,i, 2] <- yData[,i]
  }
  
  gaitbasis <- create.bspline.basis(gaitrange, nbasis=BasisNum,norder=6) #original 20 leaving norder at 6.
  D2fdPar = fdPar(gaitbasis, lambda=lambda)
  gaitSmooth = smooth.basis(gaittime, mygaitExp, D2fdPar)
  betalist  = list(const=D2fdPar, cong=D2fdPar)
  
  gaitfd = gaitSmooth$fd
  
  names(gaitfd$fdnames) = c("Normalized time", "Event", "Height")
  gaitfd$fdnames[[3]] = c("Cong", "Cedar")
  
  congfd  = gaitfd[,1] #Predictor Stuff
  cedfd = gaitfd[,2] # Response Stuff This seems like a more detailed version of what we get from the other textbook.
  
  xfdlist   = list(const=rep(1,nPred), cong=congfd)
  
  gaitRegress= fRegress(cedfd, xfdlist, betalist)
  betaestlist = gaitRegress$betaestlist
  cedIntercept = predict(betaestlist$const$fd, gaitfine) #B0
  congCoef = predict(betaestlist$cong$fd, gaitfine) #Slope Term B1
  
  Intercept<-approx(x=1:length(cedIntercept), y=cedIntercept, n = n)$y
  Slope<-approx(x=1:length(congCoef), y=congCoef , n = n)$y
  
  if(BetaOnly==FALSE){
    cedhatfd = gaitRegress$yhatfd$fd #
    cedhatmat = eval.fd(gaittime, cedhatfd)  #1 to 2000 compare this to yHat.t.mat
    resmat. = mygaitExp[,,2] - cedhatmat #The 2 represents cedar (response)
    SigmaE = cov(t(resmat.))
    
    cedfinemat   = eval.fd(gaitfine, cedfd)
    cedmeanvec   = eval.fd(gaitfine, mean.fd(cedfd))
    cedhatfinemat= eval.fd(gaitfine, cedhatfd)
    resmat        = cedfinemat - cedhatfinemat #finer grid than the . needed for computation below
    ncurve        = dim(mygaitExp)[2] #Says the number of flood events
    resmat0 = cedfinemat - cedmeanvec %*% matrix(1,1,ncurve)
    SSE0 = apply((resmat0)^2, 1, sum)
    SSE1 = apply(resmat^2, 1, sum)
    ced.R2 = (SSE0-SSE1)/SSE0
    
    ylim2=c(0, max(congCoef, ced.R2)) #And this is what I had plotted before exactly.
    gaitbasismat = eval.basis(gaitfine, gaitbasis)
    y2cMap = gaitSmooth$y2cMap #Needed to get error
    
    cedhatfd = gaitRegress$yhatfd$fd #
    cedhatmat = eval.fd(gaittime, cedhatfd)  #functional yHat matrix
    resmat. = mygaitExp[,,2] - cedhatmat #The 2 represents cedar (response)
    SigmaE = cov(t(resmat.))
    
    fRegressList1 = fRegress(cedfd, xfdlist, betalist,
                             y2cMap=y2cMap, SigmaE=SigmaE)
    
    fRegressList2 = fRegress.stderr(fRegressList1, y2cMap, SigmaE)
    betastderrlist = fRegressList2$betastderrlist
    
    error.B0<-predict(betastderrlist[[1]], gaitfine)
    error.B1<-predict(betastderrlist[[2]], gaitfine)
    #This gives the standard errors for B0t and B1t to use for DFBETAS
    stderr.B0<-approx(x=1:length(error.B0), y=error.B0, n = n)$y
    stderr.B1<-approx(x=1:length(error.B1), y=error.B1, n = n)$y
    
    
    #Adjust to be nRow Ncol
    meanHat_i<-c()
    data1<-cbind(rep(1,nPred),t(xData)) #adds the column of 1's
    Hatiit<-matrix(NA,nrow=nPred,ncol=n)
    for(t in 2:(n+1)){
      X<-data1[,c(1,t)]
      Hat<-X%*%solve(t(X)%*%X)%*%t(X)
      Hatiit[,t-1]<-diag(Hat) #Each row is an event. The corresponding column is the Hii of that event at each timepoint
    }
    meanHat_i<-c(meanHat_i,apply(Hatiit,1,mean))#will return an average Hii for each of the 10 sampled curve
    
    out<-list(Intercept,Slope,stderr.B0,stderr.B1,resmat.,Hatiit,meanHat_i,cedhatmat)
    names(out)<-c("B0t","B1t","stderr.B0t","stderr.B1t","Residuals","AdjHatMatrix","MeanHat","FullyHat")
  }
  else {out<-list(Intercept,Slope)
  names(out)=c("B0t","B1t")
  }
  return(out)
}
FullTempResults<-PredictFRegressNormTestUpdatedForRUpdate(AirTempMat,WaterTempMat, rockport_air,Basis = "BSpline", nBasis = 21, Plot = TRUE)

FullResults<-BetaOnlyBSpline(AirTempMat,WaterTempMat,BasisNum=21)
B0_i_<-matrix(NA,nrow=nrow(AirTempMat),ncol=ncol(AirTempMat))
B1_i_<-matrix(NA,nrow=nrow(AirTempMat),ncol=ncol(AirTempMat))
stderrB1_i_2<-matrix(NA,nrow=nrow(AirTempMat),ncol=ncol(AirTempMat))
stderrB0_i_2<-matrix(NA,nrow=nrow(AirTempMat),ncol=ncol(AirTempMat))
BasisToUse<-c()
for(i in 1:ncol(AirTempMat)){
  BasisToUse[i]<-21
  tempBetas<-BetaOnlyBSpline(AirTempMat[,-i],WaterTempMat[,-i],BetaOnly=FALSE,BasisNum=BasisToUse[i])
  B0_i_[,i]<-tempBetas$B0t
  B1_i_[,i]<-tempBetas$B1t
  stderrB0_i_2[,i]<-tempBetas$stderr.B0t
  stderrB1_i_2[,i]<-tempBetas$stderr.B1t
}

DFBETAS011<-matrix(NA,nrow=nrow(AirTempMat),ncol=ncol(AirTempMat))
DFBETAS111<-matrix(NA,nrow=nrow(AirTempMat),ncol=ncol(AirTempMat))
DFBETA011<-matrix(NA,nrow=nrow(AirTempMat),ncol=ncol(AirTempMat))
DFBETA111<-matrix(NA,nrow=nrow(AirTempMat),ncol=ncol(AirTempMat))
for(i in 1:ncol(AirTempMat)){
  DFBETAS011[,i]<-(FullResults$B0t-B0_i_[,i])/stderrB0_i_2[,i]
  DFBETAS111[,i]<-(FullResults$B1t-B1_i_[,i])/stderrB1_i_2[,i]
  DFBETA011[,i]<-(FullResults$B0t-B0_i_[,i])
  DFBETA111[,i]<-(FullResults$B1t-B1_i_[,i])
}


##Cooks Distance
MSE_t<-apply((FullResults$Residuals)^2,1,mean) #Gives MSE at each t

##Need the predicted values of Yj when i is left out.
yHat_mini<-matrix(NA,nrow=nrow(AirTempMat),ncol=ncol(AirTempMat))
for(i in 1:ncol(AirTempMat)){
  yHat_mini[,i]<-B0_i_[,i]+B1_i_[,i]*AirTempMat[,i]
}

CooksD_ti<-matrix(0:0, nrow=nrow(AirTempMat),ncol=ncol(AirTempMat))
for(i in 1:ncol(AirTempMat)){
  #Get predicted value of each event using B0 and B1 calculated without i
  yHatj_<-matrix(NA,nrow=nrow(AirTempMat),ncol=ncol(AirTempMat))
  for(j in 1:ncol(AirTempMat)){
    yHatj_[,j]<-B0_i_[,i]+B1_i_[,i]*AirTempMat[,j] #B0_i_[,i]+B1_i_[,i] calculation of B0 and B1 when i is not included
  }
  for(t in 1:nrow(AirTempMat)){
    CooksD_ti[t,i]<-sum((as.numeric(FullResults$FullyHat[t,i])-yHatj_[t,i])^2)/((1+1)*MSE_t[t])
  }
}

##DFFITS
yHat_mini<-matrix(NA,nrow=nrow(AirTempMat),ncol=ncol(AirTempMat))
for(i in 1:ncol(AirTempMat)){
  yHat_mini[,i]<-B0_i_[,i]+B1_i_[,i]*AirTempMat[,i]
}
par(mfrow=c(1,1))
dev.new()
plot(yHat_mini[,1],type="l",main =i)
for(i in 1:ncol(AirTempMat)){
  plot(yHat_mini[,i],type="l",main =i)
  
  lines(FullResults$FullyHat[,i],col="red")
  legend("bottom", legend = c("yHati", "yHati(-i)"),
         lty = 1,col = c(1,2)) }

DFFIT<-FullResults$FullyHat-yHat_mini

#For MSE_mini minus i we want the true-predicted at each t
#for all 34 events when the 35th event is left out
MSE_mini_t<-matrix(NA,nrow=nrow(AirTempMat),ncol=ncol(AirTempMat))
yHatsNoi<-matrix(NA,nrow=nrow(AirTempMat),ncol=ncol(AirTempMat)-1)
n=1:ncol(AirTempMat)
for(i in n){
  k=1
  for(j in n[-i]){
    yHatsNoi[,k]<-B0_i_[,i]+B1_i_[,i]*AirTempMat[,j] #goes through all but i using the betas with i left out
    k=k+1
  }
  MSE_mini_t[,i]<-apply((WaterTempMat[,-i]-yHatsNoi)^2,1,mean) #Get an MSE across the other events when event i is left out
}
DFFITS<-matrix(NA,nrow=nrow(DFFIT),ncol=ncol(DFFIT))
for(t in 1:nrow(DFFIT)){
  DFFITS[t,]<-DFFIT[t,]/sqrt(MSE_mini_t[t,]*FullResults$AdjHatMatrix[,t])
}

DB0<-round(apply(abs(DFBETAS011),2,mean),3)
DB1<-round(apply(abs(DFBETAS111),2,mean),3)
CD<-round(apply(abs(CooksD_ti),2,mean),3)
meanDFFITS<-round(apply(abs(DFFITS),2,mean),3)

#Now bootstrapping null distributions to get percentile
leave1outBSpline<-function(PredictorMat, ResponseMat, PredictorVec, nBasis, lambda = 10^(-1)){
  AllPredictedResponse<-matrix(NA, ncol = ncol(ResponseMat), nrow = nrow(ResponseMat))
  AllPredictedLower<-matrix(NA, ncol = ncol(ResponseMat), nrow = nrow(ResponseMat))
  AllPredictedUpper<-matrix(NA, ncol = ncol(ResponseMat), nrow = nrow(ResponseMat))
  #tells it how many Fourier basis functions to use eachtime
  for(i in 1:ncol(PredictorMat)){
    PredictedResponse<-PredictFRegressNormTestUpdatedForRUpdate(PredictorMat[,-i], ResponseMat[,-i],PredictorVec , lambda = lambda,Basis="BSpline", nBasis = nBasis) #Gets predicted Response for the ith column left out
    PredictedResponseEstimate<-PredictedResponse$PredictedResponse
    PredictedResponseLower<-PredictedResponse$Lower
    PredictedResponseUpper<-PredictedResponse$Upper
    AllPredictedResponse[,i]<-PredictedResponseEstimate #stores the actual prediction if needed
    AllPredictedUpper[,i]<-PredictedResponseUpper #stores the actual prediction if needed
    AllPredictedLower[,i]<-PredictedResponseLower #stores the actual prediction if needed
  }
  out<-list(AllPredictedResponse, AllPredictedLower, AllPredictedUpper)
  names(out)<-c("Predicted Responses Without i", "Lower Bounds", "Upper Bounds")
  return(out)
}
OUprocess_MetricsBSpline<-function(xData ,yData , basisList,PredVec, metric,L2=TRUE,prob_i,sigma=3, drift=.9,times=2000,endtime = 1000,BasisList){
  #Step 1: use Sim10 function to grab the bootstrapped events wanted
  tempdata<-samplecurves(xData, yData, prob_i = prob_i, size =ncol(xData)) #Samples 10 curves from the sample with replacement
  #Step 2: Apply the OH process to allow the response curves to have more variability
  OUyData<- OUonFloodEM(sigma=sigma,drift=drift,times=times, endtime = endtime,Data=tempdata$yData)
  newY<-matrix(NA,nrow=nrow(yData),ncol=ncol(yData))
  
  for(i in 1:ncol(yData)){
    newY[,i]<-supsmu(x=OUyData$timevec,y=OUyData$obscurves[,i])$y #generically smooth the OU data to look more like my data
  }
  
  #Step 3: Calculate "Full" Basis numbers needed for the iteration will use what was used for the original data
  #fullbasisNum<-L2bestEst(tempdata$xData, newY, basisList=BasisList)$`Best Basis` #chooses the best number of basis functions to use
  #fullbasisNum<-as.integer(names(fullbasisNum))
  fullbasisNum=BasisList
  
  #Step 4: Calculate Full model results getting B0, B1, Full prediction etc.
  #It's all saved under "FullResults"
  FullResults<-BetaOnlyBSpline(tempdata$xData, newY,BasisNum=fullbasisNum)
  FullReconstructedY<-FullResults$B0t+FullResults$B1t*PredVec
  
  #Step 5: If L2 or L1 run the leave 1 out model and find the
  #L1 or L2 distance for each left out
  if(metric == "deltaP"){
    if(L2==TRUE){
      tempEachLeftOut<-leave1outBSpline(PredictorMat=tempdata$xData, ResponseMat=newY, PredictorVec=PredVec,  basisnumberused = fullbasisNum)
      tempL2<-c()
      for(j in 1:ncol(xdata)){
        squaredifference<-(FullReconstructedY - tempEachLeftOut$`Predicted Responses Without i`[,j])^2
        tempL2[j]<-sqrt(trapz(x=1:length(squaredifference), y = squaredifference))
      }
      L_distanceout<-tempL2
    }
    else{
      tempEachLeftOut<-leave1outBSpline(PredictorMat=tempdata$xData, ResponseMat=newY, PredictorVec=PredVec,  basislist = fullbasisNum)
      tempL1<-c()
      for(j in 1:ncol(xData)){
        absdifference<-abs(FullReconstructedY - tempEachLeftOut$`Predicted Responses Without i`[,j])
        tempL1[j]<-(trapz(x=1:length(absdifference), y = absdifference))
      }
      L_distanceout<-tempL1
    }
    output_metric =  L_distanceout
  }
  #Step 6  Calculate DFBETAs0 and DFBETAs1 for each event at each point return this matrix
  #6b Also return the mean abs difference of these for each event
  #Get Beta0t and Beta1t when each event is leftout
  if(metric == "DFBETAS0" || metric == "DFBETAS1" || metric == "DFFITS" || metric =="CooksD"){
    
    B0_i_<-matrix(NA,nrow=nrow(xData),ncol=ncol(xData))
    B1_i_<-matrix(NA,nrow=nrow(xData),ncol=ncol(xData))
    propB0<-c() #gives proportion of sample that is larger than cutoff
    propB1<-c()
    stderrB1_i_<-matrix(NA,nrow=nrow(xData),ncol=ncol(xData)) #std errort B0 and B1 each i left out
    stderrB0_i_<-matrix(NA,nrow=nrow(xData),ncol=ncol(xData))
    DFBETAS0<-matrix(NA,nrow=nrow(xData),ncol=ncol(xData))
    DFBETAS1<-matrix(NA,nrow=nrow(xData),ncol=ncol(xData))
    meanAbsB0<-c()
    meanAbsB1<-c()
    for(i in 1:ncol(xData)){
      #Optimizing number of basis functions for the 9 events in case that matters
      tempBetas<-BetaOnlyBSpline(tempdata$xData[,-i],newY[,-i],BetaOnly=FALSE,BasisNum=fullbasisNum) #Gets B0t and B1t with each i left out
      B0_i_[,i]<-tempBetas$B0t
      B1_i_[,i]<-tempBetas$B1t
      stderrB0_i_[,i]<-tempBetas$stderr.B0t
      stderrB1_i_[,i]<-tempBetas$stderr.B1t
      DFBETAS0[,i]<-(FullResults$B0t-B0_i_[,i])/stderrB0_i_[,i]
      DFBETAS1[,i]<-(FullResults$B1t-B1_i_[,i])/stderrB1_i_[,i]
      meanAbsB0[i]<-mean(abs(DFBETAS0[,i]))
      meanAbsB1[i]<-mean(abs(DFBETAS1[,i]))
    }
    if(metric == "DFBETAS0"){
      output_metric = meanAbsB0
    }
    
    if(metric == "DFBETAS1"){
      output_metric = meanAbsB1
    }
    
  }
  
  #Step 7 Calculate the Hat Matrix Can probably use what is in step 3 (YES)
  #FullResults$AdjHatMatrix is the pointwise hat matrix
  #FullResults$MeanHat
  if(metric == "DFFITS"){
    
    #Step 8  Calculate DFFITs for each event at each point return this matrix
    yHat_mini<-matrix(NA,nrow=nrow(xData),ncol=ncol(xData))
    for(i in 1:ncol(xData)){
      yHat_mini[,i]<-B0_i_[,i]+B1_i_[,i]*tempdata$xData[,i] #yHat for i when ith event is left out. B1_i and B0_i calculated before
    }
    DFFIT<-FullResults$FullyHat-yHat_mini
    #8b Also return the mean abs difference of these for each event
    MSE_mini_t<-matrix(NA,nrow=nrow(tempdata$xData),ncol=ncol(tempdata$xData))
    yHatsNoi<-matrix(NA,nrow=nrow(tempdata$xData),ncol=ncol(tempdata$xData)-1)
    n=1:ncol(tempdata$xData) #1:10
    for(i in n){
      k=1
      for(j in n[-i]){
        yHatsNoi[,k]<-B0_i_[,i]+B1_i_[,i]*tempdata$xData[,j] #goes through all but i using the betas with i left out
        k=k+1
      }
      MSE_mini_t[,i]<-apply((newY[,-i]-yHatsNoi)^2,1,mean) #Get an MSE across the other events when event i is left out
    }
    DFFITS<-matrix(NA,nrow=nrow(DFFIT),ncol=ncol(DFFIT))
    for(t in 1:nrow(DFFIT)){
      DFFITS[t,]<-DFFIT[t,]/sqrt(MSE_mini_t[t,]*FullResults$AdjHatMatrix[,t])
    }
    MeanAbsDFFITS<-apply(abs(DFFITS),2,mean)
    output_metric = MeanAbsDFFITS
    
  }
  
  #Step 9 Calculate Cooks D  and mean across each event and return cuttoff criteria based on F dist
  
  if(metric == "CooksD"){
    MSE_t<-apply((FullResults$Residuals)^2,1,mean)
    CooksD_ti<-matrix(0:0, nrow=nrow(tempdata$xData),ncol=ncol(tempdata$xData))
    for(i in 1:ncol(tempdata$xData)){
      #Get predicted value of each event using B0 and B1 calculated without i
      yHatj_<-matrix(NA,nrow=nrow(tempdata$xData),ncol=ncol(tempdata$xData))
      for(j in 1:ncol(tempdata$xData)){
        yHatj_[,j]<-B0_i_[,i]+B1_i_[,i]*tempdata$xData[,j] #B0_i_[,i]+B1_i_[,i] calculation when i is not included
      }
      for(t in 1:nrow(xData)){
        CooksD_ti[t,i]<-sum((as.numeric(FullResults$FullyHat[t,i])-yHatj_[t,i])^2)/((1+1)*MSE_t[t])
      }
    }
    
    MeanCooksD_i<-apply((CooksD_ti),2,mean)
    
    output_metric = MeanCooksD_i
  }
  
  out<-list(tempdata$xData,newY,fullbasisNum,tempdata$index,output_metric)
  names(out)<-c("xDataUsed","yDataUsed","fullbasisNum","IndicesUsed","output_metric")
  return(out)
}

#Gettting theta:
range_water<-c()
myrange<-function(x){max(x)-min(x)}
range_water<-apply(WaterTempMat,2,myrange)
mean(range_water)/2 #6.6
mean(range_water)/3 #4.4

#Simply adjust alpha and the metric to DFFITS, DFBETAS0, or DFBETAS1 to get all of the percentiles
#While also changing alpha and the theta calculation to include all observed measures of interest
cd_results_alpha0=c()
alpha = 0.5
#obs_mean_abscd<-apply(abs(meanDFB0),2,mean)
theta<-as.vector((1/CD)^alpha/sum((1/CD)^alpha))
for(i in 1:100){
  sigma = runif(1,4.4,6.6)
  drift = runif(1,0.5,1)
  test=OUprocess_MetricsBSpline(xData=AirTempMat, yData = WaterTempMat,metric="CooksD",sigma = sigma, drift = drift, PredVec = Rockport,prob_i=theta,times=1000,endtime = 500,BasisList=21)
  cd_results_alpha0=c(cd_results_alpha0,test$output_metric)
}
quantile(cd_results_alpha0,c(.9,.95,.99))
