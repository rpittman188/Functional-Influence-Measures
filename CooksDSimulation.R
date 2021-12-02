#This is the Code for Cook's D Simulation and does 1 iteration

##Can manually change things like BasisNum, the B0fun, B1fun, B1funout, functions, and length
library(fda)
.libPaths("~/Rlibs")
library(MASS)
library(corpcor) 
iter<-as.numeric(commandArgs()[8])

sin_a=function()sample(c(-3,-2,-1,0,1,2,3),1)
sin_k=function()sample(c(-300,-200,-100,100,200,300),1)
sin_d=function()sample(c(-100,-50,0,50,100),1)
sin_c=function()sample(c(-3,-2,-1,0,1,2,3),1)

cos_a=function()sample(c(-3,-2,-1,0,1,2,3),1)
cos_k=function()sample(c(-300,-200,-100,100,200,300),1)
cos_d=function()sample(c(-100,-50,0,50,100),1)
cos_c=function()sample(c(-3,-2,-1,0,1,2,3),1)

#used to generate X data
mu2t<-function(t){ 1/12*t+(sin_a()*sin(1/sin_k()*(t-sin_d()))+sin_c())*(cos_a()*cos(1/cos_k()*(t-cos_d()))+cos_c())}


OUonFloodEM<-function(sigma=1,drift=1,times,endtime,Data){
  dt<-endtime/times; #want this to be around .5 for most situations to keep it stable enough
  objects<-ncol(Data) #number of flood events I have
  normalvec<-matrix(rnorm(objects*times,sd=1),nrow=objects,ncol=times); #this is randn from matlab code
  ouvalue<-matrix(0:0,nrow=objects,ncol=times) #null matrix of 0s
  
  # Defining the cluster sizes and the true clustering structure: No clusters here so just 1
  timevec<-seq(0,endtime,length=times)
  
  # loop to create approximation of Ornstein-Uhlenbeck process	# Using Euler-Maruyama approximation
  ouvalue[1:objects,1]=sigma*normalvec[,1] #So we don't start at 0.
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

#Need this to be established for the above function to run
B0fun<-function(t) cos(t/200)+2
B1fun<-function(t) sin(t/200)+2
B1funout<-function(t,lambda) lambda*sin(t/200)+2
length = 1000
t=1:length #This is the length

#Generates Functional Data with an outlier
GenerateFunctionalDataOut<-function(N = 100,length=1000,lambda=1)
{
  NewXout<-matrix(NA, nrow = length,ncol = N)
  t<-1:length
  for(i in 1:N){
    NewXout[,i]<-mu2t(t)
  }
  #Runs the functions created at the top of the page 
  B0t<-B0fun(t)
  B1t<-B1fun(t)
  B1funout<-B1funout(t,lambda)
  #add an outlier
  outlier_obs<-sample(1:N,1)
  holdY<-NewXout ##just sets the dimensions
  for(i in 1:N){
    if(i != outlier_obs){
      holdY[,i]<-NewXout[,i]*B1t+B0t
    }
    else holdY[,i]<-NewXout[,i]*B1funout+B0t
  }
  
  NewYout<-holdY #just gives it same dimensions before editing
  run2<-OUonFloodEM(sigma=15,drift=.5,times= 2*length, endtime = length,Data=as.matrix(holdY))
  for(i in 1:N){
    NewYout[,i]<-supsmu(x=run2$timevec,y=run2$obscurves[,i])$y
  }
  out<-list(NewXout, NewYout,outlier_obs)
  names(out)<-c("Xdata","Ydata","OutlierNumber")
  return(out)
}

#Function to return some measures from the concurrent model
BetaOnlyBSpline<-function(xData,yData,BasisNum=15, BetaOnly=FALSE,lambda=10^(-1)){
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
    cedhatfd = gaitRegress$yhatfd$fd #THIS IS WHERE THE ERROR OCCURS unless using fda version 5.1.5.1
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

# temp1<-GenerateFunctionalDataOut(N = 20,length=1000,lambda=2)
# Xtmat<-temp1$Xdata
# Ytmat<-temp1$Ydata

##Calculates DFBETAS0 and DFBETAS1
CooksD_fun<-function(Xtmat,Ytmat,BasisNum = 15){
  FullResults<-BetaOnlyBSpline(Xtmat,Ytmat,BasisNum=BasisNum)
  #Then use the above function to get Betas with things left out.
  B0_i_<-matrix(NA,nrow=nrow(Xtmat),ncol=ncol(Xtmat))
  B1_i_<-matrix(NA,nrow=nrow(Xtmat),ncol=ncol(Xtmat))
  stderrB1_i_2<-matrix(NA,nrow=nrow(Xtmat),ncol=ncol(Xtmat))
  stderrB0_i_2<-matrix(NA,nrow=nrow(Xtmat),ncol=ncol(Xtmat))
  #BasisToUse<-c()
  for(i in 1:ncol(Xtmat)){
    tempBetas<-BetaOnlyBSpline(Xtmat[,-i],Ytmat[,-i],BasisNum=BasisNum) #calculate B0t and B1t without each i.
    B0_i_[,i]<-tempBetas$B0t
    B1_i_[,i]<-tempBetas$B1t
    stderrB0_i_2[,i]<-tempBetas$stderr.B0t
    stderrB1_i_2[,i]<-tempBetas$stderr.B1t
  }
  
  MSE_t<-apply((FullResults$Residuals)^2,1,mean) #Gives MSE at each t
  
  yHat_mini<-matrix(NA,nrow=nrow(Xtmat),ncol=ncol(Xtmat))
  for(i in 1:ncol(Xtmat)){
    yHat_mini[,i]<-B0_i_[,i]+B1_i_[,i]*Xtmat[,i]
  }
  
  CooksD_ti<-matrix(0:0, nrow=nrow(Xtmat),ncol=ncol(Xtmat))
  for(i in 1:ncol(Xtmat)){
    #Get predicted value of each event using B0 and B1 calculated without i
    yHatj_<-matrix(NA,nrow=nrow(Xtmat),ncol=ncol(Xtmat))
    for(j in 1:ncol(Xtmat)){
      yHatj_[,j]<-B0_i_[,i]+B1_i_[,i]*Xtmat[,j] #B0_i_[,i]+B1_i_[,i] calculation of B0 and B1 when i is not included
    }
    for(t in 1:nrow(Xtmat)){
      CooksD_ti[t,i]<-sum((as.numeric(FullResults$FullyHat[t,i])-yHatj_[t,i])^2)/((1+1)*MSE_t[t])
    }
  }

  out<-list(CooksD_ti)
  names(out)<-c("CooksD_ti")
  return(out)
} #Returns DFBETAS0t and DFBETAS1t for each of the N events. Each event is a column

absCooksD_bootsamedata<-function(alpha, xdata = sim_data$Xdata,ydata = sim_data$Ydata, outlier = sim_data$OutlierNumber, nboot = 100){
  N = ncol(xdata) #Say the number of observations we have
  #Now run the concurrent model and calculate mean |dffits| for each of the observation
  raw_CooksD<-CooksD_fun(Xtmat = xdata, Ytmat = ydata,BasisNum = 15)$CooksD_ti #pointwise DFBETAS0 for each observation
  mean_abs_CooksD<-apply(abs(raw_CooksD),2,mean)
  adjusted_obs_absCooksD<-mean_abs_CooksD[outlier]
  theta<-as.vector((1/mean_abs_CooksD)^alpha/sum((1/mean_abs_CooksD)^alpha))
  #Need to set up the bootstrap next SEE OrdinaryRegressionDistribution for Details
  #Set blank objects here
  boot_meanabsCooksD<-c()
  for(i in 1:nboot){
    index_samp<-sample(1:N,N,replace=TRUE,prob = theta) #resamples the needed data with the given probability 
    temp_CooksD<-CooksD_fun(Xtmat = xdata[,index_samp],Ytmat = ydata[,index_samp],BasisNum = 15)$CooksD_ti
    boot_meanabsCooksD<-c(boot_meanabsCooksD,apply(abs(temp_CooksD),2,mean))
  }
  percent_95<-quantile(boot_meanabsCooksD,0.95) #Gives the 95Th percentile of N*nboot influence measures
  percentile <- ecdf(boot_meanabsCooksD)
  adjusted_point_percentile<- percentile(mean_abs_CooksD[outlier])
  
  out<-list(percent_95,adjusted_point_percentile,mean_abs_CooksD) #returns the 95th percentile from the N*nboot |dffits| calculated
  names(out)<-c("Boot95","AdjPercentile","AllObsMeanAbsCooksD") #And returns the percentile of this list that the observed |dffits| falls into
  return(out)
}

##Can adjust this so I input more into the function like alpha, lambda etc
one_fulltable_bootstrap<-function(N = 100, nboot = 100, alpha, lambda, length = 1000){
  
  results_output95<-matrix(NA,nrow = length(lambda), ncol = length(alpha)) #Will store the 95%ile from the bootstrap results
  results_outputpercentile<-matrix(NA,nrow = length(lambda), ncol = length(alpha)) #Will store the percentile that the adjusted observation falls into 
  results_numberobs_over95<-matrix(NA,nrow = length(lambda), ncol = length(alpha)) #Will store the percentile that the adjusted observation falls into 
  
  
  colnames(results_output95)<-as.character(alpha)
  rownames(results_output95)<-as.character(lambda)
  colnames(results_outputpercentile)<-as.character(alpha)
  rownames(results_outputpercentile)<-as.character(lambda)
  colnames(results_numberobs_over95)<-as.character(alpha)
  rownames(results_numberobs_over95)<-as.character(lambda)
  
  ##Below will do 1 iteration of each cell
  i = 1 #starts filling in column 1
  for(l in lambda){
    #Generate new data for each lambda. then keep same for different alpha functions are defined at the top
    # B0fun(t<-function(t) cos(t/200)+2
    # B1fun<-function(t) sin(t/200)+2
    # B1funout<-function(t,lambda) l*sin(t/200)+2
    t=1:length
    sim_data<-GenerateFunctionalDataOut(N=N,length=length,lambda=l)
    j = 1 #start on row 1 each time
    for(a in alpha){
      tempout<-absCooksD_bootsamedata(xdata = sim_data$Xdata, ydata = sim_data$Ydata, outlier = sim_data$OutlierNumber, alpha = a,  nboot = nboot)
      results_output95[i,j]<- tempout$Boot95
      results_outputpercentile[i,j]<- tempout$AdjPercentile
      results_numberobs_over95[i,j] <-sum(ifelse(tempout$AllObsMeanAbsCooksD>=tempout$Boot95,1,0))
      #results_numberobs_over95[i,j]   count number from tempout$raw that are over tempout$boot95
      j=j+1 #move to the next column which is a new Alpha
    }
    i=i+1 #moves to the next column 
  }
  
  out<-list(results_output95,results_outputpercentile,results_numberobs_over95)
  names(out)<-c("the95percentile", "outlierpercentile","numberobs_over95")
  return(out)
}

##Run that function and repeat "iterations" time and get the mean of each cell
alpha<-c(0,.1,.3,.5)
#lambda<-c(.25,.5,.75,.9,1,1.1,1.25,1.5,2,3,4)
lambda = 1.75
#This will fill both the 95th table and the adjusted observations's percentile for each combination of alpha and lambda
#test<-one_fulltable_bootstrap(N=10, nboot = 10, alpha = alpha, lambda = lambda, length = 1000)

# repetitions = 100 #change as needed
# 
# sum_of_reps95<-matrix(0,nrow = length(lambda), ncol = length(alpha))
# sum_of_repspercentile<-matrix(0,nrow = length(lambda), ncol = length(alpha))

#This is what runs the sim, returning a matrix based on the size of alpha and lambda for one generation of data bootstrapped 100 times

###CHANGE N AND nboot BACK TO 100
hold<-one_fulltable_bootstrap(N=100, nboot = 100, alpha = alpha, lambda = lambda, length = length)
#hold$the95percentile
#hold$outlierpercentile
Adjobs_above_95<-ifelse(hold$outlierpercentile>=0.95,1,0) # is the adjusted observation over 95%ile yes/no

#Want to save avg_95 and avg_percentile
filename95=paste("~/InfluenceMeasures/data/N50CooksD175_avg_95_",iter,".RData",sep="")
x<-hold$the95percentile #the 95th percentile itself
save(x, file=filename95)

filenamepercentile=paste("~/InfluenceMeasures/data/N50CooksD175_avg_adjpercentile_",iter,".RData",sep="")
y<-hold$outlierpercentile #The outlier's percentile from this run
save(y, file=filenamepercentile)

filenameOver95=paste("~/InfluenceMeasures/data/N50CooksD175_avg_numberobs_over95_",iter,".RData",sep="")
z<-hold$numberobs_over95 #Number of all intitial observations over 95%ile
save(z, file=filenameOver95)

filenameAdjpercentileover95=paste("~/InfluenceMeasures/data/N50CooksD175_Above_95_",iter,".RData",sep="")
save(Adjobs_above_95, file=filenameAdjpercentileover95)

