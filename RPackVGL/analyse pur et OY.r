
#determiner le path du fichier actuel et le recuper 
dir <- dirname(rstudioapi::getSourceEditorContext()$path)
#marche pas hors de rstudio/ligne de commande? (https://stackoverflow.com/questions/47044068/get-the-path-of-current-script)

#dir <- choose.dir()
#"C:\\devel\\l-egume\\legume\\test"



source(paste(dir, "fonctions_analyses.r",sep="\\"))
source(paste(dir, "fonctions_mef.r",sep="\\"))




#lecture des fichiers de moyenne

#tabpur <- read.csv("H://Travail//simul//test variance BLW//save_tablePurmultN0.csv", header=T)
tabpur <- read.csv("H://Travail//simul//test variance BLW//save_tablePurmultN.csv", header=T)
#tabpur <- read.csv("H://Travail//simul//test variance BLW//save_tablePur1pN.csv", header=T)
#tabpur <- read.csv("H://Travail//simul//test variance BLW//save_tablePur1pN0.csv", header=T)



#a partir de ca peut recalculer les OY de chaque simuls
#lecture des dtoto mult...

#dtoto <- read.csv("C://inputs//inputs test variance BLW//dtoto//save_dtoto_allvar-multiN0.csv",header=T, sep=";")
dtoto <- read.csv("C://inputs//inputs test variance BLW//dtoto//save_dtoto_allvar-multiN+.csv",header=T, sep=";")
#dtoto <- read.csv("C://inputs//inputs test variance BLW//dtoto//save_dtoto1pN-allvar-new.csv",header=T, sep=",")
#dtoto <- read.csv("C://inputs//inputs test variance BLW//dtoto//save_dtoto1p0N-allvar-new.csv",header=T, sep=",")



#### ajouter data pur1 et pur2 pour chaque ligne
tabpur$trait <- as.character(tabpur$trait)
trait1 <- paste(dtoto$scenario1, dtoto$sd, dtoto$Mng)
trait2 <- paste(dtoto$scenario2, dtoto$sd, dtoto$Mng)


#?? pas ttes les lignes pour trait1? -> certains pas dans le tableau des purs!
#completer tabpur!

lsmanquantspur <- unique(trait1[which(! trait1 %in% tabpur$trait)])

#sauvegarde le dtoto
dtoto2 <- dtoto

#exrait les manquant
dtoto <- dtoto[which(! trait1 %in% tabpur$trait),]


dtoto$trait <- paste(dtoto$scenario1, dtoto$sd, dtoto$Mng)

res<-by(dtoto$Ytot, dtoto$trait, mean)
res <- data.frame(trait=names(res), Ytot=as.numeric(res))

ressup <- data.frame(t(data.frame(strsplit(as.character(res$trait), " "))))
row.names(ressup) <- 1:dim(ressup)[1]
names(ressup) <- c("scenario1","sd","Mng")

res <- cbind(ressup,res)


#ajout des PARi12 et QNupttot
x<-by(dtoto$QNupttot, dtoto$trait, mean)
res$QNupttot <- as.numeric(x)
x<-by(dtoto$Pari1, dtoto$trait, mean)
res$Pari1 <- as.numeric(x)
x<-by(dtoto$Pari2, dtoto$trait, mean)
res$Pari2 <- as.numeric(x)
res$Pari <- res$Pari1+res$Pari2

#met a jour tabpur
tabpur <- rbind(tabpur,res)


#relis le dtoto
dtoto <- dtoto2


res1 <- NULL
for(trait in trait1)
{
  res1 <- rbind(res1, tabpur[tabpur$trait==trait,c("Ytot","QNupttot","Pari")])
}
names(res1) <- c("Ytot_Pur1","QNupttot_Pur1","Pari_Pur1")

res2 <- NULL
for(trait in trait2)
{
  res2 <- rbind(res2, tabpur[tabpur$trait==trait,c("Ytot","QNupttot","Pari")])
}
names(res2) <- c("Ytot_Pur2","QNupttot_Pur2","Pari_Pur2")


dtoto <- cbind(dtoto, res1)
dtoto <- cbind(dtoto, res2)
names(dtoto)
dim(dtoto)


#### calcul OY et RYT

dtoto$RY1 <- dtoto$YEsp1 / dtoto$Ytot_Pur1
dtoto$RY2 <- dtoto$YEsp2 / dtoto$Ytot_Pur2
dtoto$RYT <- dtoto$RY1 + dtoto$RY2
dtoto$OY <- (dtoto$YEsp1+dtoto$YEsp2) - 0.5*(dtoto$Ytot_Pur1+dtoto$Ytot_Pur2)

#plot(OY~RYT, dtoto)


#### calcul indices Loreau... SE/CE



dtoto$deltaRY1 <- dtoto$RY1-dtoto$Semprop1
dtoto$deltaRY2 <- dtoto$RY2 - (1-dtoto$Semprop1)
#calcul de fonction Calc_CESE_diag qui est ok (pb avec transfertN??)
#pourquoi 1 dans transfertN?

vCE <- NULL
vSE <- NULL
for (idY in 1:dim(dtoto)[1])
{
  #idY <- 3
  
  ex_M <- as.numeric(dtoto[idY, c("Ytot_Pur1","Ytot_Pur2")]) #rdmt en pur!!
  ex_deltaRY <- as.numeric(dtoto[idY, c("deltaRY1","deltaRY2")])
  
  CE <- Calc_CEi(2,ex_M,ex_deltaRY)
  SE <- Calc_SEi(2,ex_M,ex_deltaRY)
  vCE <- rbind(vCE, CE)
  vSE <- rbind(vSE, SE)
}

vCE + 0.5*vSE #c'est bon!

dtoto$CE <- as.numeric(vCE)
dtoto$SE <- as.numeric(vSE)
#yes!


#write.csv(dtoto, "save_dtotomultiN-OY.csv", row.names=F)






###########################
#faire boxplot - les moyennes pas traitement!
# qqs graph

dtoto$trait <- paste(dtoto$scenario2, dtoto$sd)#, dtoto$Mng)

boxplot(Ytot~trait, dtoto, las=3, cex.axis=0.5, main=dtoto$Mng[1])
boxplot(OY~trait, dtoto, las=3, cex.axis=0.5, main=dtoto$Mng[1])
boxplot(SE~trait, dtoto, las=3, cex.axis=0.5, main=dtoto$Mng[1])
boxplot(CE~trait, dtoto, las=3, cex.axis=0.5, main=dtoto$Mng[1])


plot(dtoto$OY, dtoto$SE, main=dtoto$Mng[1])
plot(dtoto$OY, dtoto$CE, main=dtoto$Mng[1])
plot(dtoto$SE, dtoto$CE, main=dtoto$Mng[1])

plot(dtoto$Yprop1, dtoto$SE, main=dtoto$Mng[1])
plot(dtoto$Yprop1, dtoto$CE, main=dtoto$Mng[1])


#pour definir niveaux de facteur et leur ordre!
lsLum <- c(46,45,44,1,50,51,52)
lsN <- c(49,48,47,1,53,54,55)
lsall <- c(40,39,38,1,41,42,43)

Ntrait <-"N+"#"N-"#"N/2"# 
var_ <- "Ytot"#"CE"# "SE"#"OY"#"Yprop1"#

if (var_ == "Yprop1") {yrange <- c(0,1)
}else {yrange <-c(-100,400)}#c(-100,250)

if (var_ == "Ytot") {yrange <- c(0,1900)}


layout(matrix(1:3, 1,3))

#Light
sc_ <- "SD1-1"#"SD23-23"#"SD24-24"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsLum,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsLum)) #remise en ordre des niveaux de facteur!
x$var <- x[,c(var_)]
boxplot(var~scenario2, x, ylim=yrange, main=paste(Ntrait,"Light"), ylab=var_)
sc_ <- "SD24-24"#"SD1-1"#"SD23-23"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsLum,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsLum)) #remise en ordre des niveaux de facteur!
x$var <- x[,c(var_)]
boxplot(var~scenario2, x, add=T, border=2)


#Nitrogen
sc_ <- "SD1-1"#"SD26-26"#"SD25-25"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsN,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsN)) #remise en ordre des niveaux de facteur!
x$var <- x[,c(var_)]
boxplot(var~scenario2, x, ylim=yrange, main=paste(Ntrait, "Nitrogen"),ylab=var_)
sc_ <- "SD25-25"#"SD1-1"#"SD26-26"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsN,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsN)) #remise en ordre des niveaux de facteur!
x$var <- x[,c(var_)]
boxplot(var~scenario2, x, add=T, border=2)


#All
sc_ <-"SD1-1"#"SD4-4"# "SD3-3"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsall,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsall)) #remise en ordre des niveaux de facteur!
x$var <- x[,c(var_)]
boxplot(var~scenario2, x, ylim=yrange, main=paste(Ntrait, "All"), ylab=var_)
sc_ <-"SD4-4"#"SD1-1"# "SD3-3"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsall,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsall)) #remise en ordre des niveaux de facteur!
x$var <- x[,c(var_)]
boxplot(var~scenario2, x, add=T,  border=2)







#######################
#calcul initial moyenne des rendements par traitement  en pur a partir du dtoto

dtoto <- read.csv("H://Travail//simul//test variance BLW//save_dtotomultiPurN.csv",header=T)


#"Lusignan30IrrNN"#"Lusignan30Irr"

dtoto$trait <- paste(dtoto$scenario2, dtoto$sd, dtoto$Mng)

res<-by(dtoto$Ytot, dtoto$trait, mean)
res <- data.frame(trait=names(res), Ytot=as.numeric(res))

ressup <- data.frame(t(data.frame(strsplit(as.character(res$trait), " "))))
row.names(ressup) <- 1:dim(ressup)[1]
names(ressup) <- c("scenario1","sd","Mng")

res <- cbind(ressup,res)


#ajout des PARi12 et QNupttot
x<-by(dtoto$QNupttot, dtoto$trait, mean)
res$QNupttot <- as.numeric(x)
x<-by(dtoto$Pari1, dtoto$trait, mean)
res$Pari1 <- as.numeric(x)
x<-by(dtoto$Pari2, dtoto$trait, mean)
res$Pari2 <- as.numeric(x)
res$Pari <- res$Pari1+res$Pari2


#write.csv(res, "save_tablePur1pN0.csv", row.names=F)
#tabpur <- res

#boxplot(Ytot~trait, res,las=3)
boxplot(Ytot~trait, tabpur,las=3)
#rq dans presque tous les cas diversite a un effet positif sur le rendement
#effet de selection! pas les indiv moyens mais les elites qui pilotent
#dans l'absolu, fait gagne peu
#pour forts rendement peut faire perdre


#?? j'ai des simuls avec des N- dasn 1pN
#dtoto[dtoto$trait == "8 SD1-1",c("Ytot", "Mng")]

