#####################################################################
testZ=function(qProb,map,pErr=0.9,simplitol=1e-3,optiCrit=2,disp=0)
#####################################################################
{
#
valRef= quantile(map$krigGrid,na.rm=TRUE,prob=qProb)
#initial zoning
ZK = initialZoning(qProb, map,pErr=pErr,simplitol=simplitol,optiCrit=optiCrit,disp=disp)
crit=ZK$resCrit
Z=ZK$resZ$zonePolygone
K=ZK$resZ
if(disp)
	dispZ(map$pas,map$krigMatTest,zonePolygone=Z,boundary=map$boundary,nbLvl=0)
return(list(Z=Z,K=K))
}