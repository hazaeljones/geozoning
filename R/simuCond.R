################################################################################
#fonction qui effectue la simulation conditionnelle

#entree:1/step=difference en abscisse et en ordonnee entre deux points kriges sur la grille,taille x et y du cadre(numeric),
#coord et valeurs des points conditionnants(dataframe),nom du dossier du test courant(character),nombre de simulations a effectuer
#modele du variogramme,numero d'execution sur cette map
#sortie:structure contenant les valeurs et coordonnees de points pour differentes iterations de la simulation conditionnelle
######################################################################
#' simulCond
#'
#' @details description, a paragraph
#' @param step xxxx
#' @param tabAlea xxxx
#' @param dosstest xxxx
#' @param nbSimul xxxx
#' @param modelVar xxxx
#' @param bordure xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
simulCond=function(step,tabAlea,dosstest,nbSimul,modelVar,bordure)
######################################################################
{

  #on effectue les iterations de la simulation conditionnelle et on les sauvegarde dans un fichier RData
  #dans le cas ou les iterations ont deja ete effectuees et sauvegardees pour ce test,on les charge au lieu de resimuler (gain temps d'execution)
  if(file.exists(paste(getwd(),"/resultats/",dosstest,"/cond",step,1,1,nbSimul,".RData",sep="")))
  {
    load(paste(getwd(),"/resultats/",dosstest,"/cond",step,1,1,nbSimul,".RData",sep=""))
  }
  else
  {
    #grille pour la simulation conditionnelle(definit en quels points les valeurs seront simulees)
    #grille sur toute la map
    xempty=rep(seq(1/step, 1 -1/step, by=(1/step)),(1*step -1))
    yempty=as.vector( t( matrix( xempty,  1*step -1, (1*step -1)) ) )
    valCond=data.frame(matrix(rep(NA,length(xempty)*nbSimul),ncol=nbSimul))
    #restriction de la structure a remplir(grille a l'interieur de la bordure seulement)
    #on regarde quels points de la grille sont interieurs
    pointsInterieurs=point.in.polygon(xempty,yempty,bordure$x,bordure$y)

    #on creee les coordonnees x et y adequates
    xcond=xempty[pointsInterieurs==1]
    ycond=yempty[pointsInterieurs==1]

    #matrice des positions des points conditionnants
    dataVar=matrix(c(coordinates(tabAlea)[,1],coordinates(tabAlea)[,2]),ncol=2)
    colnames(dataVar)=c("data.x","data.y")
    #dataframe contenant les valeurs des points conditionnants et leur position
    dataCond=conventional2RFspDataFrame(data.frame(tabAlea)[[3]],dataVar)

    #on effectue la simulation conditionnelle a l'interieur du contour
    condPart=RFsimulate(model=modelVar,x=xcond,y=ycond,data=dataCond,n=nbSimul)#,spConform=FALSE)

    #on etend a la grille entiere en rajoutant des NA
    #valCond[pointsInterieurs==1]=condPart[[1]]
    valCond[pointsInterieurs==1,]=data.frame(condPart)[,1:nbSimul]
    cond=conventional2RFspDataFrame(as.matrix(valCond),matrix(c(xempty,yempty),nrow=length(xempty)),n=nbSimul)

    plot(condPart,dataCond)
    save(dataCond,cond,condPart,file=paste(getwd(),"/resultats/",dosstest,"/cond",step,1,1,nbSimul,".RData",sep=""))
  }

  return(cond)
}


#################################################
#fonction qui supprime les bords aberrants des valeurs moyennes de la simulation conditionnelle et les remplace par des valeurs krigees
#(on recupere les valeurs de krigTabAlea au bord et on les met a la palce des originales)
#entree:moyenne des valeurs de simulations conditionelles en chaque point(dataframe),step(numeric)
#sortie:moyCondTab modifiee(dataframe)
#################################################
#' krigeBord
#'
#' @details description, a paragraph
#' @param moyCondTab xxxx
#' @param step xxxx
#' @param krigTabAlea xxxx
#'
#' @return a ?
#'
#' @export
#'
#' @examples
#' # not run
krigeBord=function(moyCondTab,step,krigTabAlea)
#################################################
{
  #structure non spatiale
  krigTabAlea=data.frame(krigTabAlea)

  #on divise le data frame en deux: partie centrale(tabCentre) et contour(tabBord)
  tabCentre=subset(moyCondTab,x!=1/step & y!=1/step&x!=(1-1/step) & y!=(1-1/step))

  #contour que l'on veut eliminer(effets de bord)
  tabBord=subset(moyCondTab,x==1/step |y==1/step|x==(1-1/step)|y==(1-1/step))

  #contour par lequel on le remplace
  krigTabBord=subset(krigTabAlea,x==1/step |y==1/step|x==(1-1/step)|y==(1-1/step))[,1:3]

  #on reconstruit la dataframe avec la partie centrale inchangee et le contour krige
  moyCondTab=rbind(tabCentre,krigTabBord)

  #on remet dans l'ordre les valeurs du dataframe en vue de la transformation en matrice
  masque=order(moyCondTab[,"y"],moyCondTab[,"x"])
  moyCondTab=moyCondTab[masque,]

  return(moyCondTab)
}

###############################################################
#fonction d'affichage de differentes donnees sur un champ gaussien
#entree:matrices de: valeurs krigees,moyenne des simulations conditionnelles,ecart type et ecart-type normalise des simulations conditionnelles,
#taille du cadre en x et en y(numeric),1/step=distance entre 2 points sur la grille
#sortie:-(affichage)
#---------------------------------------------------------------------------------------------------------------------------------#
###############################################################
#' dispResCond
#'
#' @details description, a paragraph
#' @param krigMatTest xxxx
#' @param matMoyCond xxxx
#' @param matEcartCond xxxx
#' @param matEcartNorm xxxx
#' @param step xxxx
#'
#' @return a ?
#' @importFrom RandomFields conventional2RFspDataFrame
#'
#' @export
#'
#' @examples
#' # not run
dispResCond=function(krigMatTest,matMoyCond,matEcartCond,matEcartNorm,step)
###############################################################
{
  par(mfrow=c(3,2))
  palCoul=colorRampPalette(c("brown","yellow"))
  coulBreaks=seq(min(krigMatTest,na.rm=TRUE),max(krigMatTest,na.rm=TRUE),length.out=max(20,trunc(max(krigMatTest,na.rm=TRUE)-min(krigMatTest,na.rm=TRUE))))
  coulBreaks2=seq(min(matMoyCond,na.rm=TRUE),max(matMoyCond,na.rm=TRUE),length.out=max(20,trunc(max(matMoyCond,na.rm=TRUE)-min(matMoyCond,na.rm=TRUE))))

  image(main="Difference relative",x=seq(1/step, 1 - 1/step, by=1/step), y=seq(1/step, 1 - 1/step, by=1/step),abs((krigMatTest-matMoyCond)/krigMatTest))#,col=palCoul(200))  #on affiche sur une grille en couleur les variances normalisees en chaque point de la simulation conditionnelle
  dispZ(step,krigMatTest,nbLvl=20,coulBreaks=coulBreaks,texMain="Valeurs krigees")
  dispZ(step,matMoyCond,nbLvl=20,coulBreaks=coulBreaks2,texMain="Moyenne des simulations conditionnelles")
  image(main="Ecart-type des simulations conditionnelles",matEcartCond,x=seq(1/step, 1 - 1/step, by=1/step),y= seq(1/step, 1 - 1/step, by=1/step)) #on affiche sur une grille en couleur la valeur moyenne en chaque point de la simulation conditionnelle
  coulBreaks=seq(min(matEcartNorm),(2*sd(matEcartNorm)+mean(matEcartNorm)),length.out=20)
  image(main="Ecart-type/Moyenne",seq(1/step, 1 - 1/step, by=1/step), seq(1/step, 1 - 1/step, by=1/step),matEcartNorm,breaks=coulBreaks,col=palCoul(length(coulBreaks)-1))  #on affiche sur une grille en couleur les variances normalisees en chaque point de la simulation conditionnelle

  persp(krigMatTest)                    #affichage d'une vue en perspective du champ gaussien
  par(mfrow=c(1,1))
  return()
}
