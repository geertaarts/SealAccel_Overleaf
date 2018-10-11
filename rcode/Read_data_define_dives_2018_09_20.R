# Create dropbox token in windows  
  library(rdrop2)
  token <- drop_auth() # opens dropbox in internet explorer
  saveRDS(token, file = "d:/token.rds")
 
# Copy file from dropbox onto W; drop_download saves it into setwd   
  # Dropbox stuff on nemo (or windows)
  library(rdrop2)
  setwd("//home/aarts012/")
  setwd("wur/W/IMARES/IJmuiden/Afdeling/Projecten/Nemo/Aarts/versnelling/14477/")
  token<-readRDS("token.rds")
  drop_download(path = "/Collaborative projects/Accelerometer/Data Paul/14477/14477_accel.txt", dtoken=token, overwrite=TRUE)
  #drop_upload(file_name,path = "Seal_database/environmental/distance_to_haulout/done_north_sea", dtoken=token)

###################################################################3    

# Set working directory
  setwd("//home/aarts012/")
  setwd("wur/W/IMARES/IJmuiden/Afdeling/Projecten/Nemo/Aarts/versnelling/14477/")
  #setwd("E:/Dropbox/Collaborative projects/Accelerometer/Data Paul/14477")
  
# Read accelerometer data
  #header<-names(read.table("14477_accel.txt",nrows=1,header=TRUE))
  #PVAc<-read.table("14477_accel.txt",nrows=100000,header=TRUE) # first bit of data
  PVAc<-read.table("14477_accel.txt",header=TRUE) 
  #names(PVAc)<-header

  PVAc.all<-PVAc
  
# Loop through all the days
  for (DD in 1:length(unique(PVAc.all$DATE))){
    
    # Only select one day
    PVAc<-PVAc.all[PVAc.all$DATE==unique(PVAc.all$DATE)[DD],]
    
  
# Only select one day
  #PVAc<-PVAc[PVAc$DATE=="2017/04/25",]
  
# Get ddate of Accelerometer
  PVAc$ddate<-as.POSIXct(paste(PVAc$DATE,PVAc$X.UTC.),format="%Y/%m/%d %H:%M:%OS",tz="UTC") # cant convert to 3 digits, but last is 0, OS dependent?
  
# Calculate static component and standard deviation (expression variation within a window)
  library(zoo)
  numrows=12
  
  PVAc$meanX=rollapply(PVAc$X,numrows,function(x){mean(x,na.rm=TRUE)},fill=NA)
  PVAc$meanY=rollapply(PVAc$Y,numrows,function(x){mean(x,na.rm=TRUE)},fill=NA)
  PVAc$meanZ=rollapply(PVAc$Z,numrows,function(x){mean(x,na.rm=TRUE)},fill=NA)
  
  PVAc$sdX=rollapply(PVAc$X,numrows,function(x){sd(x,na.rm=TRUE)},fill=NA)
  PVAc$sdY=rollapply(PVAc$Y,numrows,function(x){sd(x,na.rm=TRUE)},fill=NA)
  PVAc$sdZ=rollapply(PVAc$Z,numrows,function(x){sd(x,na.rm=TRUE)},fill=NA)
  
  PVAc$sum_sd<-PVAc$sdX+PVAc$sdY+PVAc$sdZ
  
  
# Select only data with small variation over time  
  qq<-quantile(PVAc$sum_sd,0.05,na.rm=TRUE) # was 0.05
  qq<-quantile(PVAc$sum_sd,0.025,na.rm=TRUE) # was 0.05
  qq<-quantile(PVAc$sum_sd,0.01,na.rm=TRUE) # was 0.05
  PVAc.s<-PVAc[PVAc$sum_sd<qq & is.na(PVAc$sum_sd)==FALSE,]
  plot(PVAc.s$X,ylim=c(-1.2,1.2),col="red",pch=20,cex=0.3)
  points(PVAc.s$Y,col="blue",cex=0.3, pch=20)
  points(PVAc.s$Z,col="orange",cex=0.3, pch=20)

# New select for each X, Y, Z data points accross the range between -1 and 1, with the smallest standard devation; i.e. stable sections   
  store<-PVAc.s[0,]
  for (XYZ in c("meanX","meanY","meanZ")){
    for (i in seq(-1,0.9,by=0.1)){
        PVAc.ss<-PVAc.s[(PVAc.s[,XYZ]>i & PVAc.s[,XYZ]<(i+0.1)),]
      print(i)
    store<-rbind(store,PVAc.ss[PVAc.ss$sum_sd==min(PVAc.ss$sum_sd,na.rm=TRUE),][1,])
  }}

# Remove na's
  store<-store[is.na(store$DATE)==FALSE,]
  
# GEt the parameter estimates of offset and scale
  min.RSS20dp <- function(data, par) {
    with(store, sum((sqrt((par[1] + (par[2]+1) * meanX)^2 + (par[3] + (par[4]+1)*meanY)^2 + (par[5] + (par[6]+1)*meanZ)^2)- 1)^2))}
  
  result20dp <- optim(par=c(0,0,0,0,0,0), fn=min.RSS20dp,  data = OFF,method="L-BFGS-B",lower=rep(-0.5,6),upper=rep(0.5,6),control=list(maxit=20000,pgtol=0.00000000001),hessian = TRUE)
  result20dp
  library(HelpersMG)
  SE<-SEfromHessian(result20dp$hessian)
  SEfromHessian(result20dp$hessian)

# Note the SEs are very large, why???  

# Correct the mean X, Y and Z and see if Gstat is indeed closer to 1  
  pars<-result20dp$par
  print(pars)
  
  if(DD==1) pars.dat<-pars
    pars.dat<-rbind(pars.dat,pars)
  }
  
  pars.dat2<-pars.dat[-1,]
  row.names(pars.dat2)<-unique(PVAc.all$DATE)
  pars.dat2<-as.data.frame(pars.dat2)
  pars.dat2$ddate<-as.POSIXct(unique(PVAc.all$DATE),format="%Y/%m/%d",tz="UTC") # cant convert to 3 digits, but last is 0, OS dependent?
  
  
  store$Xn<--0.034 +(1+0.005)*store$meanX
  store$Yn<-0.035+(1+0.023)*store$meanY
  store$Zn<-0.262+(1+0.002)*store$meanZ
  
  store$Gstat_cor_old<-sqrt(store$Xn^2 + store$Yn^2 + store$Zn^2)
  
  store$Xn<-pars[1]+(1+pars[2])*store$meanX
  store$Yn<-pars[3]+(1+pars[4])*store$meanY
  store$Zn<-pars[5]+(1+pars[6])*store$meanZ
  
  store$Gstat_cor<-sqrt(store$Xn^2 + store$Yn^2 + store$Zn^2)
  store$Gstat<-sqrt(store$meanX^2 + store$meanY^2 + store$meanZ^2)

  store

# Save workspace
  save(pars.dat2,file="pv14477_pars.dat.rdata")
  drop_upload("pv14477_pars.dat.rdata",path = "/Collaborative projects/Accelerometer/Data Paul/14477/", dtoken=token, mode="overwrite")

  v=1
  if (v==1) the.main="offset X"
  if (v==2) the.main="scale X"
  if (v==3) the.main="offset Y"
  if (v==4) the.main="scale Y"
  if (v==5) the.main="offset Z"
  if (v==6) the.main="scale Z"
  par(mai=c(0.7,0.8,0.3,0.1))
  plot(pars.dat2$ddate,pars.dat2[,paste("V",v,sep="")],xlab="",ylab="parameter estimate",main=the.main)
  
  
  
# CLASIFY DIVES  
  
  # Read tdr dive depth data
    PVDEPTH<-readLines("14477_tdr..txt")
  
  # Get Dive depth data (Problem is that it also contains haul-out data inbetween the lines)
    qq.split<-strsplit(PVDEPTH,"\t")
    qTF<-which(substr(PVDEPTH,1,4)!="HAUL")
    dive_depths<-do.call("rbind", qq.split[qTF])
    colnames(dive_depths)<-dive_depths[1,]
    dive_depths<-as.data.frame(dive_depths[-1,])
    names(dive_depths)<-c("DATE_UTC","DEPTH_m","TEMP_C")
  
  # Process the dive data  
    str(dive_depths)
    dive_depths$ddate<-as.POSIXct(as.character(dive_depths$DATE_UTC),format="%Y/%m/%d %H:%M:%S",tz="UTC")
    #plot(dive_depths$ddate,pch=20,cex=0.4)
    dive_depths$DEPTH_m<-as.numeric(as.character(dive_depths$DEPTH_m))
    dive_depths$TEMP_C<-as.numeric(as.character(dive_depths$TEMP_C))
    dive_depths$DEPTH_m_1<-c(0,dive_depths$DEPTH_m[1:(nrow(dive_depths)-1)])
    dive_depths<-dive_depths[order(dive_depths$ddate),]
  
  # Read accelerometer depth data  
    accel_depth<-read.table("14477_accel_depth.txt",header=TRUE, fill=TRUE)
  
  # Get similar column names as dive_depth
    accel_depth$ddate<-as.POSIXct(paste(accel_depth$DATE, accel_depth$X.UTC.,sep=""),format="%Y/%m/%d %H:%M:%S",tz="UTC")
    #plot(accel_depth$ddate,pch=20,cex=0.4)
    names(accel_depth)[names(accel_depth)=="DEPTH"]<-"DEPTH_m"  
    accel_depth<-accel_depth[order(accel_depth$ddate),]
    accel_depth$DEPTH_m_1<-c(0,accel_depth$DEPTH_m[1:(nrow(accel_depth)-1)])
  
  # Use accelerometer depths
    dive_depths<-accel_depth
  
  # Get start and end times of the dives  
    dive_start<-dive_depths$ddate[dive_depths$DEPTH_m>1.5 & dive_depths$DEPTH_m_1<=1.5]
    dive_end<-dive_depths$ddate[dive_depths$DEPTH_m_1>1.5 & dive_depths$DEPTH_m<=1.5]
    dive.start.end<-data.frame(id=1:length(dive_start),dive_start,dive_end)
    dive.start.end$dive_start_nextdive<-c(dive_start[2:length(dive_start)],Inf)
    dive.start.end$dive_start_nextdive<-format(dive.start.end$dive_start_nextdive, tz="UTC", usetz=TRUE)
    dive.start.end$dive_start_nextdive<-as.POSIXct(dive.start.end$dive_start_nextdive, tz="UTC")
  
  # Only get dive data for specific period
    dive.start.end<-dive.start.end[dive.start.end$dive_start>(min(PVAc$ddate,na.rm=TRUE)-3600*24) & dive.start.end$dive_start<(max(PVAc$ddate,na.rm=TRUE)+3600*24),]
  
  # Get dive id
    for (i in 1:nrow(dive.start.end)){
      PVAc$dive.id[PVAc$ddate>dive.start.end$dive_start[i] & PVAc$ddate<=dive.start.end$dive_start_nextdive[i] ]<-dive.start.end$id[i]
    }
    unique(PVAc$dive.id)
  
  