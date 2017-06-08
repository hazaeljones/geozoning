###############################
#' visuZ
#'
#' @details description, a paragraph
#' @param qProb xxxx
#' @param map xxxx
#' @param nq xxxx
#' @param zk xxxx
#' @param disp xxxx
#'
#' @return a ?
#' @importFrom grDevices contourLines
#'
#' @export
#'
#' @examples
#' # not run
visuZ=function(qProb,map,nq=3,zk,disp=TRUE)
###############################
{
# carte et zk var. globales
resC=correctionTree(qProb,map,SAVE=TRUE)
zk=resC$zk
criterion=resC$criterion
cost=resC$cost
costL=resC$costL
nz=resC$nz
best=searchNODcrit(qProb,length(zk),zk,criterion,cost,costL,nz)

if (nq==5)
{
    K=zk[[length(zk)]][[best$ind$q5[1]]]
   }
else if (nq==4)
{
     K=zk[[length(zk)]][[best$ind$q4[1]]]
   }
else if (nq==3)
{
    K=zk[[length(zk)]][[best$ind$q3[1]]]
   }
else if (nq==2)
{
    K=zk[[length(zk)]][[best$ind$q2[1]]]
   }
else
{
     K=zk[[length(zk)]][[best$ind$q1[1]]]
   }
 Z=K$zonePolygone
if(disp) dispZ(map$step,map$krigGrid,zonePolygone=Z,K=K,boundary=map$boundary,nbLvl=0)
return(K)
}

#################################################################################################################################

#---------------------------------------------------------------------------------------------------------------------------------#

#fonction qui affiche une représentation couleur des valeurs (pts sur un grille), superposée avec les isocontours

#entrée:step de répartition des valeurs sur la grille(numeric),vecteur des valeurs(numeric),nombre d'isocontours(numeric),
  #liste des polygones(list(list(numeric))),nombre de polygone(numeric)
#sortie:-  (affichage)


###############################################################################
#' dispZ
#'
#' @details description, a paragraph
#' @param step xxxx
#' @param matVal xxxx
#' @param nbLvl xxxx
#' @param zonePolygone xxxx
#' @param K xxxx
#' @param coulBreaks xxxx
#' @param texMain xxxx
#' @param boundary xxxx
#' @param id xxxx
#' @param valeurQuant xxxx
#' @param palCoul xxxx
#' @param indiceArg xxxx
#' @param mu xxxx
#'
#' @return a ?
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics text
#' @importFrom fields image.plot
#'
#' @export
#'
#' @examples
#' # not run
dispZ=function(step,matVal,nbLvl=20,zonePolygone=NULL,K=NULL,coulBreaks=0,texMain="",boundary=NULL,id=FALSE,valQ=NULL,palCoul=colorRampPalette(c("brown","yellow")),noXY=FALSE,indiceArg=0,mu=1,cex=1)
###############################################################################
{
  Z=zonePolygone
  nbPoly=length(Z) #0 if Z is NULL
  levelsList=NULL
  if (!is.null(Z))
     {
	vs=getZsize(Z)
	xsize=vs["x"]
  	ysize=vs["y"]
	}
	
  #génération des nbLvls isocontours dans un vecteur
  #soit les valeurs précises des iso sont données
  if(!is.null(valQ))
  {
    levelsList=valQ     
   
  }  else if(nbLvl!=0) #soit on a précisé combien on en voulait (par défaut 20),
  {
    levelsList=seq(min(matVal,na.rm=TRUE),max(matVal,na.rm=TRUE),length=nbLvl)
  }
  texMain=c(texMain,valQ)
  #si le nombre de paliers pour l'affichage est donné en entrée et valide
  if(length(coulBreaks)!=1){
    
    image.plot(main=texMain,seq(step,xsize,by=step), seq(step,ysize,by=step),matVal,
          breaks=coulBreaks,col=palCoul(length(coulBreaks)-1),legend.shrink=0.4,legend.width=0.6)   
  }
  else{
    #sinon on prend arbitrairement 20 couleurs de paliers uniformément entre la plus petite et la plus grande valeur des données en entrée
    #affichage de la carte des valeurs en entrée(nécessite une grille régulière...)
    if (!is.null(matVal))
       image.plot(seq(step,xsize-step,by=step), seq(step,ysize-step,by=step),matVal,col=palCoul(20),legend.shrink=0.4,legend.width=0.6,xlab="X",ylab="Y")
    else
	plot(0,0,xlim=c(0,xsize),ylim=c(0,ysize),type="n") #no data, only set axes
  }
    
    #affichage des contours si on en a généré
  if(!is.null(levelsList))
  {
       if (!is.null(matVal))
       	  {
		listLines=contourLines(seq(step,xsize,by=step), seq(step,ysize,by=step),
                matVal, levels = levelsList)
                for(i in (1:length(listLines)))
      		      {
        	      #affichage des isocontours par dessus la carte affichée
        	      lines(listLines[[i]])
		      }
   
           }
  }
  
 if(!is.null(boundary))
  {
    lines(boundary)
  }
 
 
  
  
  if(nbPoly!=0)
  {    
    
    for (j in (1:nbPoly))
    {
      #affichage des éventuels polygones de zones déja définis
      lines(Z[[j]]@polygons[[1]]@Polygons[[1]]@coords,type='l')
     
      pointLabel=Z[[j]]@polygons[[1]]@Polygons[[1]]@labpt
      # utiliser id Zone au lieu de num Zone
      if(id)
	jid = getZoneId(Z[[j]])
      else
	jid=j
	if(mu==2 & !is.null(K))
	{
	 vecmu=K$meanZone
     	 text(pointLabel[1],pointLabel[2],paste(jid,"(",round(vecmu[j],2),")",sep=""))
	 }
	 else if(mu==1)
	 text(pointLabel[1],pointLabel[2],jid,cex=cex)
	 
	 
      
      
    }
    
    if (indiceArg !=0)
    {
      lines(Z[[indiceArg]]@polygons[[1]]@Polygons[[1]]@coords,type='l',col='red')  
      lines(boundary)
    }
  }
  return()
}
#########################################################
#'
#' dispZmap
#' @details description, a paragraph
#' @param map xxxx
#' @param Z
#' @param m
#' @param valbp
#' @param scale
#' @param lev
#' @param palCoul
#' @param legend.width
#' # not run
dispZmap=function(map,Z=NULL,m=0,valbp=NULL,scale=NULL,lev=20,palCoul=colorRampPalette(c("brown","yellow")),legend.width=1)
#########################################################
{
step=map$step
xsize=map$xsize
ysize=map$ysize
matVal=map$krigGrid
bd=map$boundary
xlim=c(0,xsize);ylim=c(0,ysize)

par(fig=c(0,1,0.3,1))#c(x1,x2,y1,y2) bottom=0,0 - top=1,1
par(mar=c(0,0,5,0))#c(bottom, left, top, right))
image(seq(step,xsize-step,by=step), seq(step,ysize-step,by=step),matVal,col=palCoul(lev),asp=1,xlim=xlim,ylim=ylim,xlab="",ylab="",xaxt="n",yaxt="n")
image.plot( zlim=range(matVal,na.rm=TRUE), nlevel=20,legend.only=TRUE, vertical=FALSE,col=palCoul(20))
if (m!=0) title(paste("Best zoning for nL=",length(m)+1,"\nm=",paste(m,collapse=","),sep=""),cex.main=1.5)
#
dec=0.03
decy=0.05
xn=-0.2
yn=1.8
p1=polygon(list(x=c(xn,xn+dec,xn+dec,xn),y=c(yn,yn+decy,yn+2*decy,yn)))
p2=polygon(list(x=c(xn+dec,xn+dec,xn+2*dec,xn+dec),y=c(yn+decy,yn+2*decy,yn,yn+decy)),density=100,angle=45,col="black")
text(xn+dec,yn-decy,"N",font=2,cex=1)
#
# add scale
if(!is.null(scale))
{
met=100/scale
ys=0.78
xs=0.8
lines(c(xs,xs+met),c(ys,ys))
tau=0.02
lines(c(xs,xs),c(ys-tau,ys+tau))
lines(c(xs+met,xs+met),c(ys-tau,ys+tau))
text(xs,ys-3*tau,"0")
text(xs+met,ys-3*tau,"100m")
}
# add zones
par(fig=c(0,1,0.3,1))#c(x1,x2,y1,y2) bottom=0,0 - top=1,1
par(mar=c(0,0,5,0))#c(bottom, left, top, right))

if(!is.null(Z))
{
for (k in 1:length(Z))
{
lines(Z[[k]])
pointLabel=Z[[k]]@polygons[[1]]@Polygons[[1]]@labpt
jid=getZoneId(Z[[k]])
text(pointLabel[1],pointLabel[2],jid,cex=1.5)
}
}
#
if(!is.null(valbp))
{
par(fig=c(0.2,0.9,0,0.3),new=TRUE)
par(mar=c(2,0,3,0))#c(bottom, left, top, right)
boxplot(valbp,cex.axis=1,col=palCoul(length(valbp)))
title("Distribution of values within zones",cex.main=1)
}
return()
}

