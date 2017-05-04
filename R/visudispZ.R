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
#'
#' @export
#'
#' @examples
#' # not run
dispZ=function(step,matVal,nbLvl=20,zonePolygone=NULL,K=NULL,coulBreaks=0,texMain="",boundary=NULL,
               id=TRUE,valeurQuant=NULL,palCoul=colorRampPalette(c("brown","yellow")),indiceArg=0,mu=1)
###############################################################################
{
  nbPoly=length(zonePolygone) #0 if zonePolygone is NULL
  listeLevels=NULL
  #génération des nbLvls isocontours dans un vecteur
  #soit les valeurs précises des iso sont données
  if(!is.null(valeurQuant))
  {
    listeLevels=valeurQuant

  }  else if(nbLvl!=0) #soit on a précisé combien on en voulait (par défaut 20),
  {
    listeLevels=seq(min(matVal,na.rm=TRUE),max(matVal,na.rm=TRUE),length=nbLvl)
  }
  texMain=c(texMain,valeurQuant)
  #si le nombre de paliers pour l'affichage est donné en entrée et valide
  if(length(coulBreaks)!=1){

    image.plot(main=texMain,seq(1/step, 1 - 1/step, by=1/step), seq(1/step, 1 - 1/step, by=1/step),matVal,
          breaks=coulBreaks,col=palCoul(length(coulBreaks)-1),legend.shrink=0.4,legend.width=0.6)
  }
  else{
    #sinon on prend arbitrairement 20 couleurs de paliers uniformément entre la plus petite et la plus grande valeur des données en entrée
    #affichage de la carte des valeurs en entrée(nécessite une grille régulière...)
    if (!is.null(matVal))
       image.plot(seq(1/step, 1 - 1/step, by=1/step), seq(1/step, 1 - 1/step, by=1/step),matVal,col=palCoul(20),legend.shrink=0.4,legend.width=0.6,xlab="X",ylab="Y")
    else
	plot(0,0,xlim=c(0,1),ylim=c(0,1),type="n") #no data, only set axes
  }

    #affichage des contours si on en a généré
  if(!is.null(listeLevels))
  {
       if (!is.null(matVal))
       	  {
		listeLignes=contourLines(seq(1/step, 1 - 1/step, by=1/step),seq(1/step, 1 - 1/step, by=1/step),
                matVal, levels = listeLevels)
                for(i in (1:length(listeLignes)))
      		      {
        	      #affichage des isocontours par dessus la carte affichée
        	      lines(listeLignes[[i]])
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
      lines(zonePolygone[[j]]@polygons[[1]]@Polygons[[1]]@coords,type='l')

      pointLabel=zonePolygone[[j]]@polygons[[1]]@Polygons[[1]]@labpt
      # utiliser id Zone au lieu de num Zone
      if(id)
	jid = getZoneId(zonePolygone[[j]])
      else
	jid=j
	if(mu==2 & !is.null(K))
	{
	 vecmu=K$meanZone
     	 text(pointLabel[1],pointLabel[2],paste(jid,"(",round(vecmu[j],2),")",sep=""))
	 }
	 else if(mu==1)
	 text(pointLabel[1],pointLabel[2],jid)




    }

    if (indiceArg !=0)
    {
      lines(zonePolygone[[indiceArg]]@polygons[[1]]@Polygons[[1]]@coords,type='l',col='red')
      lines(boundary)
    }
  }
  return()
}
