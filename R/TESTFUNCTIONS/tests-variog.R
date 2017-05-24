rm(list=ls())

source("srcZ.R") # source libraries and functions - params are in initParam.R

m=0
ps=10
nug=0

n=450
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
#points(dmod2,col="blue")


x11()
mod2=vgm(model="Gau",range=s,psill=1)
dmod2=variogramLine(mod2,maxdist=1,min=1e-3)
plot(dmod2,ylim=c(0,1))
mod1=RMgauss(scale=s);data=RFsimulate(model = mod1,x=x, y=y);empvario=RFempiricalvariogram(data=data)
points(empvario@centers,empvario@sd,col="blue")

vario=variogram(variable1~1,data=data,cutoff=0.75)
points(vario$dist,vario$gamma,col="red")

# different ranges (scale parameter)
# different nuggets
# different models
# Gauss gneiting exp
# different psills (var parameter)
pdf("variog.pdf")
oldpar=par(mfrow=c(9,6))

for(k in 1:3)
{
for(v in c(1,5))
{
for(s in c(0.1,0.25,0.5))
{
for(nug in c(0,1))
{

if(k==1) mod1=RMgauss(var=v,scale=s)+RMnugget(var=nug);simudata=RFsimulate(model = mod1,x=x, y=y);empvario=RFempiricalvariogram(data=simudata);mod="gauss"
if(k==2)  mod1=RMgneiting(var=v,scale=s)+RMnugget(var=nug);simudata=RFsimulate(model = mod1,x=x, y=y);empvario=RFempiricalvariogram(data=simudata);mod="gneiting"
if (k==3)  mod1=RMexp(var=v,scale=s)+RMnugget(var=nug);simudata=RFsimulate(model = mod1,x=x, y=y);empvario=RFempiricalvariogram(data=simudata);mod="exp"
plot(empvario,model=mod1,ylim=c(0,v))
text(0.5,0.5,paste("model=",mod,",var=",v,",scale=",s,",nugget=",nug,sep=""))
}
}
}
}
par(oldpar)
dev.off()

#vgm
s=0.1
for (ps in c(1,5))
{
x11()
modvgm="Gau"
mod2=vgm(model=modvgm,range=s,psill=ps);dmod2=variogramLine(mod2,maxdist=1,min=0)
plot(dmod2,ylim=c(0,ps),main=paste("vgm model=",modvgm,",psill=",ps,",range=",s,sep=""))
}