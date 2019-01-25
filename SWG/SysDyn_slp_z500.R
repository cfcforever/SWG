library(ncdf4)
# nc = nc_open("/home/users/schen/nc/msl_1979-01-01to2018-07-31_NA_1.5x1.5.nc")
nc = nc_open("/Volumes/Data-ExFAT/DATA/nc/NA/msl_1979-01-01to2018-07-31_NA_1.5x1.5.nc")
data=ncvar_get(nc,"msl")
lon=ncvar_get(nc,'longitude')
lat=ncvar_get(nc,'latitude')
time=ncvar_get(nc,'time')

data = (data[,,seq(from=1,to=length(time),by=4)] + 
          data[,,seq(from=2,to=length(time),by=4)] + 
          data[,,seq(from=3,to=length(time),by=4)] + 
          data[,,seq(from=4,to=length(time),by=4)])/4

nx=dim(data)[2];ny=dim(data)[1]
nt=dim(data)[3]
dat = data*NA; dim(dat) <- c(nt,ny,nx)
for (i in 1:nt) dat[i,,] <- t(as.matrix(data[,,i]))
# two dimentions
dim(dat)=c(nt,nx*ny)
nc_close(nc)

# nc = nc_open("/home/users/schen/nc/z500_1979-01-01to2018-07-31_NA_1.5x1.5.nc")
nc = nc_open("/Volumes/Data-ExFAT/DATA/nc/NA/z500_1979-01-01to2018-07-31_NA_1.5x1.5.nc")
data=ncvar_get(nc,"z")
lon=ncvar_get(nc,'longitude')
lat=ncvar_get(nc,'latitude')
time=ncvar_get(nc,'time')

data = (data[,,seq(from=1,to=length(time),by=4)] + 
          data[,,seq(from=2,to=length(time),by=4)] + 
          data[,,seq(from=3,to=length(time),by=4)] + 
          data[,,seq(from=4,to=length(time),by=4)])/4

nx=dim(data)[2];ny=dim(data)[1]
nt=dim(data)[3]
dat1 = data*NA; dim(dat1) <- c(nt,ny,nx)
for (i in 1:nt) dat1[i,,] <- t(as.matrix(data[,,i]))
# two dimentions
dim(dat1)=c(nt,nx*ny)
nc_close(nc)

dat = cbind(dat, dat1)
rm(dat1, data)
# save(dat, file = "/Volumes/Data-ExFAT/dat.RData")

distance = dist(dat, upper = T, diag = T)
# save(distance, file = "/Volumes/Data-ExFAT/distance.RData")
distance = as.matrix(distance)
# save(distance, file = "/Volumes/Data-ExFAT/distance.RData")

#Quantile definition
quanti=0.98
npoints=nt
dim=numeric(npoints)
theta=numeric(npoints)

load(file = "/Volumes/Data-ExFAT/distance.RData")
d = distance
rm(distance)

l=1
for (l in 1:npoints){
  #computation of the distances from the iteration l
  #m=matrix(c(1,0,0,1,1,1,2,1,1),nrow=3,ncol=3)
  distance=d[l,]
  #distance=sqrt((x[l+1000,1]-x[,1])^2+(x[l+1000,2]-x[,2])^2+(x[l+1000,3]-x[,3])^2)
  #Taking the logarithm of the distance
  logdista=-log(distance)
  #Compute the thheshold corresponding to the quantile
  thresh=quantile(logdista, quanti,na.rm=TRUE)
  
  #Computation of theta
  Li<-which(logdista>thresh)
  #Length of each cluster
  Ti<-diff(Li)
  N=length(Ti)
  #The next formula is an estimator of theta Ferro (average cluster length)
  theta[l]=2*(sum(Ti-1))^2/(N*sum((Ti-1)*(Ti-2)))
  
  #Sort the exceedances  
  logdista=sort(logdista)
  #Find all the Peaks over Thresholds. Question, why over threshold?
  findidx=which(logdista>thresh)
  logextr=logdista[findidx[[1]]:(length(logdista)-1)]
  #The inverse of the dimension is just the average of the exceedances, why?
  dim[l]=1/mean(logextr-thresh)
  l=l+1
  print(l)
}
# plotting dim vs theta
plot(dim, theta)

# save(dim, file = "~/Documents/LSCE/Dynamic System/msl_z500_19790101-20180731_NA_1.5x1.5_dim.RData")
# save(theta, file = "~/Documents/LSCE/Dynamic System/msl_z500_19790101-20180731_NA_1.5x1.5_theta.RData")
