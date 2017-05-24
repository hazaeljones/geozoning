rm(list=ls())

source("srcZ.R") # source libraries and functions - params are in initParam.R

##########################################################################
# simulate map
seed=3

#seed=27
#map=genMap(DataObj=NULL,seed=seed,disp=FALSE)
# genMap calls genData and other functions from func-carteAlea.R

#resGene=genData(DataObj=NULL,nbPoints=nbPoints,seed=seed,typeMod=typeMod,Vpsill=Vpsill,Vrange=Vrange,Vmean=Vmean,Vang=Vang,Vanis=Vanis,bd=bd)


 #Generate random (x,y) values 
 set.seed(seed)
 nbPoints=500

 x=runif(nbPoints, min=0, max=1) 
 y=runif(nbPoints, min=0, max=1)

m=0
ps=10
r=0.25
nug=0
mod="Gau"
#mod="Lin"
 mod1=vgm(model=mod,range=r,psill=ps,nugget=nug);dmod1=variogramLine(mod1,maxdist=0.75,min=0)
plot(dmod1,ylim=c(0,ps*1.5))
lines(c(r,r),c(0,ps*1.5))


 #Generate z values according to Gaussian field
  modRM=RMgauss(var=ps,scale=r)
  modRM=RMgneiting(var=ps,scale=r)
  modRM=RMexp(var=ps,scale=r)
  
  modRMG = modRM+RMtrend(mean=m)+RMnugget(var=nug)
 	 
  map<-RFsimulate(modRMG,x,y)
  empvario=RFempiricalvariogram(data=map)
  plot(empvario,model=modRMG)


     
 #coordinates(tabData)=~x+y
 #vario=variogram(z~1,data=tabData,cutoff=0.75)
 #plot(vario,mod1,ylim=c(0,ps*1.5),xlim=c(0,0.75))
  #lines(dmod1,col="red")
 #lines(c(r,r),c(0,ps*1.5))
#lines(c(0,0.75),c(nug,nug))



x = runif(n=n, min=0,max=1)
y = runif(n=n, min=0, max=1)
s=0.25
# with vgm
X11()
 mod2=vgm(model="Gau",range=s,psill=1)
dmod2=variogramLine(mod2,maxdist=1,min=1e-3)
plot(dmod2,ylim=c(0,1))

# with RMmodel

mod1=RMgauss(scale=s);data=RFsimulate(model = mod1,x=x, y=y);empvario=RFempiricalvariogram(data=data)
plot(empvario,model=mod1,ylim=c(0,1),xlim=c(0,1))
points(dmod2,col="blue")


x11()
mod2=vgm(model="Gau",range=s,psill=1)
dmod2=variogramLine(mod2,maxdist=1,min=1e-3)
plot(dmod2,ylim=c(0,1))
mod1=RMgauss(scale=s);data=RFsimulate(model = mod1,x=x, y=y);empvario=RFempiricalvariogram(data=data)
points(empvario@centers,empvario@sd,col="blue")

vario=variogram(variable1~1,data=data,cutoff=0.75)
points(vario$dist,vario$gamma,col="red")