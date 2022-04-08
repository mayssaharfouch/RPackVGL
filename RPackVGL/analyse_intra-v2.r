
library(ineq)

#determiner le path du fichier actuel et le recuper 
dir <- dirname(rstudioapi::getSourceEditorContext()$path)
#marche pas hors de rstudio/ligne de commande? (https://stackoverflow.com/questions/47044068/get-the-path-of-current-script)

#dir <- choose.dir()
#"C:\\devel\\l-egume\\legume\\test"



source(paste(dir, "fonctions_analyses.r",sep="\\"))
source(paste(dir, "fonctions_mef.r",sep="\\"))


dirlast <- paste(dir, "testjul",sep="\\")
dirlast <- paste(dirlast, "3parSD16",sep="\\")#"bea 553"#"3parSD16"#"neutre16"
dirlast <- paste(dirlast, "test30d",sep="\\")
#dirlast <- paste(dir, "test3\\sla",sep="\\")
#dirlast <- paste(dir, "histor",sep="\\")

#dirdsk <- choose.dir()
#"D:\\outputs\\SD-multi"
#dirlast <-dirdsk


setwd(dirlast)#(dir0)#
ls_files <- list.files(dirlast)#(dir0)#





#unzip the files
ls_zip <- ls_files[grepl('.zip', ls_files)]
for (file_ in ls_zip)
{unzip(file_, exdir=dirlast)}
ls_files <- list.files(dirlast)#reliste les fichier apres dezippage


#recupere la liste des toto file names du dossier de travail
ls_toto <- ls_files[grepl('toto', ls_files) & !grepl('.zip', ls_files)]
ls_paramSD <- ls_files[grepl('paramSD', ls_files)]



#creation du dataFrame dtoto et recup des info fichier



#11 col (avec sd)
cols_ <- strsplit(ls_toto, '_')
test_long <- as.numeric(lapply(cols_, length)) #pour separer selon nb de champs (avec sd)

dtoto <- as.data.frame(t(as.data.frame(cols_[test_long==11])))#as.data.frame(t(as.data.frame(strsplit(ls_toto, '_'))))#
row.names(dtoto) <- 1: length(dtoto[,1])
dtoto <- dtoto[,c(2,3,4,5,6,7,8,10)]
names(dtoto) <- c('usm','lsystem','mix','damier','scenario','Mng', 'seed','sd')
dtoto$name <- ls_toto[test_long==11]
dtoto$seed <- as.numeric(as.character(dtoto$seed))#substr(as.character(dtoto$seed), 1, 1)
dtoto$scenario <- substr(as.character(dtoto$scenario), 9, nchar(as.character(dtoto$scenario)))

#dtoto <- rbind(temp, dtoto) #merge des 2
dtoto$keysc <- paste(dtoto$scenario, dtoto$mix, dtoto$Mng, dtoto$sd)# ajout d'une cle unique par scenario
#dtoto$damier <- as.numeric(substr(as.character(dtoto$damier), 7, 7))

nomscenar <- as.data.frame(t(as.data.frame(strsplit(dtoto$scenario, "-"))))
names(nomscenar) <- c("scenario2", "scenario1")#inverse?
row.names(nomscenar) <- 1:length(nomscenar$scenario1)
dtoto$scenario1 <- as.character(nomscenar$scenario1)
dtoto$scenario2 <- as.character(nomscenar$scenario2)


#split de dtoto et stockage dans une liste de scenatios
dtoto$keysc <- paste(dtoto$keysc, dtoto$usm) #modif cle pour le cas histor
sp_dtoto <- split(dtoto, dtoto$keysc)

#pour recup des fichiers SD
ls_tabSDall <- vector("list", length(names(sp_dtoto)))
ls_MStotall <- vector("list", length(names(sp_dtoto)))

#pour recuperation des correlation MSindividuelles
ls_res_cor <- vector("list", length=length(names(sp_dtoto)))
names(ls_res_cor) <- names(sp_dtoto)



#didcols <- as.data.frame(residcols[seq(1,21,3), ])
#didcols$damier <- unique(as.character(dtoto$damier))
#write.csv(didcols, "didcols.csv", row.names=F)

#fichier d'id colones esp1 pour damier 8 (pour les cas ou bug/oubli dans les noms de colonnes)
#didcols <- read.csv("C:/devel/l-egume/legume/multisim/didcols.csv") 




#pour intra (old a la fin si pb)
DOYdeb<- 60
#DOYScoupe <- c(165,199,231,271,334) #Avignon - an 1
DOYScoupe <- c(187,229,282,334) #Lusignan - an 1

ls_parsd <- c("Len","Lfeuille","phyllochron", "Vmax2", "ELmax", "PPtreshh")#c("Len","Vmax2", "ELmax", "Lfeuille", "phyllochron", "PPtreshh") #a lire dans fichier sd idealement

#indices des voisins pour un damier 1*16 a 50/50
cote <- 16
nblignes <- 16
ls_idv <- def_indice_vois5050(cote, nblignes)
ls_idvois <- ls_idv[[1]]
ls_idKin <- ls_idv[[2]]
ls_idnonKin <- ls_idv[[3]]
# a generaliser!

for (key in names(sp_dtoto))#key <- names(sp_dtoto)[227]
{
  ls_toto_paquet <- sp_dtoto[[key]]$name
  
  #recuperation par paquet des fichiers de base (pas de stockage de l'ensemble des fichiers en memoire)
  dtoto <- build_dtoto(sp_dtoto, key, DOYdeb, DOYScoupe)
  
  
  #remise du dtoto locl dans sp_dtoto
  sp_dtoto[[key]] <- dtoto
  
  #lecture fichier SD et de tables MStot, puis calcul des indices de decile par esp et parametre
  ltoto <- read_ltoto(ls_toto_paquet)
  
  resdec <- data.frame(key)
  for (param_name in ls_parsd) #c("Len","Vmax2"))
  {
    #calcul decile /param
    resread <- read_lsSD_MStot(ltoto, ls_paramSD, param_name)
    
    #marche pour cas de resread simul unitaire
    sp_tabSD <- split(resread[["ls_tabSD"]][[1]], resread[["ls_tabSD"]][[1]]$name)
    MStot <- resread[["ls_MStot"]][[1]]
    res <- BuildResDecil(MStot, sp_tabSD)
    
    resdec <- cbind(resdec, res)

  }
  
  ls_tabSDall[[key]] <- resdec
  #ls_tabSDall[[key]] <- resread[["ls_tabSD"]]
  #ls_MStotall[[key]] <- resread[["ls_MStot"]]
  
  
  #calcul des correlations MSindividuelles
  ls_res_cor_i <- Calc_MSindiv_Corr(ltoto, ls_toto_paquet, ls_paramSD, lspar=ls_parsd)
  ls_res_cor[[key]] <- ls_res_cor_i[["tabCorMSindiv"]]
  tabindices <- ls_res_cor_i[["datIndices"]] #a sauvegarder eventuellement
  write.table(tabindices, paste("savetabindices ", key,".csv"), sep=";", col.names = T, row.names = F)
  
  
}


#reagrege dtoto
dtoto <- do.call("rbind", sp_dtoto)
row.names(dtoto) <- 1:length(dtoto$usm)

dresdec <- do.call("rbind", ls_tabSDall)
row.names(dresdec) <- 1:dim(dresdec)[1]

res_cor <- do.call("rbind",ls_res_cor)
row.names(res_cor) <- 1:dim(res_cor)[1]



#merge( dtoto,dresdec, by="key", all=T)

#pourquoi merge fonctionne pas?
dtoto <- cbind(dtoto, dresdec[,2:dim(dresdec)[2]])

dtoto <- cbind(dtoto, res_cor[,1:(dim(res_cor)[2]-1)])


#prevoir calcul autres variables! (racines, stressN, eau..., biomasse pa coupe)
#write.csv(dtoto, "save_dtoto1p0N.csv", row.names=F)

sp_dtoto <- split(dtoto, dtoto$keysc)
names(sp_dtoto)

#split par scenario
sc_dtoto <- split(dtoto, dtoto$scenario)
names(sc_dtoto)


#boxplot(Ytot~scenario, dtoto, main="Ytot")
#boxplot(YEsp1~scenario, dtoto, main="YEsp1")

#boxplot(Yprop1~sd,dtoto, ylim=c(0,1))
#boxplot(Ytot~sd,dtoto, ylim=c(0,2000))
#boxplot(dec5_Len_Fix0~sd,dtoto, ylim=c(0,100))
#boxplot(dec5_phyllochron_Fix0~sd,dtoto, ylim=c(0,100))


#split par scenario sd
sd_dtoto <- split(dtoto, dtoto$sd)

#write.table(dtoto, paste("save_dtoto4.csv"), sep=";", col.names = T, row.names = F)

#f_ <- file.choose() #"D:\\outputs\\SD-multi\\dtoto\\save_dtoto_allvar-1pN0.csv"
#dtoto <- read.table(f_, header=T, sep=";")


# visu rapide effet sd par scenario

for(sc_ in names(sc_dtoto))
{
  #sc_ <- "1-1"
  
  layout(matrix(1:4,1,4,byrow=T))
  boxplot(Ytot~sd, sc_dtoto[[sc_]], ylim=c(0,2000), main=sc_)
  boxplot(Yprop1~sd, sc_dtoto[[sc_]], ylim=c(0,1), main=sc_)
  boxplot(gini1~sd, sc_dtoto[[sc_]], ylim=c(0,1), main=sc_)
  boxplot(alive1~sd, sc_dtoto[[sc_]], ylim=c(0,500), main=sc_)
  
  # evolution des 2 pops pour les differents parametres
  layout(matrix(1:12,2,6,byrow=T))
  boxplot(dec5_Len_Fix0~sd, sc_dtoto[[sc_]], ylim=c(0,100), main=c(sc_, paste("Len", "Fix0")))
  boxplot(dec5_Lfeuille_Fix0~sd, sc_dtoto[[sc_]], ylim=c(0,100), main=c(sc_, paste("Lfeuille", "Fix0")))
  boxplot(dec5_Vmax2_Fix0~sd, sc_dtoto[[sc_]], ylim=c(0,100), main=c(sc_, paste("Vmax2", "Fix0")))
  boxplot(dec5_ELmax_Fix0~sd, sc_dtoto[[sc_]], ylim=c(0,100), main=c(sc_, paste("ELmax", "Fix0")))
  boxplot(dec5_phyllochron_Fix0~sd, sc_dtoto[[sc_]], ylim=c(0,100), main=c(sc_, paste("phyllochron", "Fix0")))
  boxplot(dec5_PPtreshh_Fix0~sd, sc_dtoto[[sc_]], ylim=c(0,100), main=c(sc_, paste("PPtreshh", "Fix0")))
  
  boxplot(dec5_Len_Fix1~sd, sc_dtoto[[sc_]], ylim=c(0,100), main=c(sc_, paste("Len", "Fix1")))
  boxplot(dec5_Lfeuille_Fix1~sd, sc_dtoto[[sc_]], ylim=c(0,100), main=c(sc_, paste("Lfeuille", "Fix1")))
  boxplot(dec5_Vmax2_Fix1~sd, sc_dtoto[[sc_]], ylim=c(0,100), main=c(sc_, paste("Vmax2", "Fix1")))
  boxplot(dec5_ELmax_Fix1~sd, sc_dtoto[[sc_]], ylim=c(0,100), main=c(sc_, paste("ELmax", "Fix1")))
  boxplot(dec5_phyllochron_Fix1~sd, sc_dtoto[[sc_]], ylim=c(0,100), main=c(sc_, paste("phyllochron", "Fix1")))
  boxplot(dec5_PPtreshh_Fix1~sd, sc_dtoto[[sc_]], ylim=c(0,100), main=c(sc_, paste("PPtreshh", "Fix1")))

}


#to do
# faire une fonction gini perso pour plus dependre de ineq!
# faire fonction unitaire pour // serveur





for(sc_ in names(sd_dtoto))
{
  #sc_ <- "SD1-1"#"SD10-10"#"SD9-9"#
  
  layout(matrix(1:4,1,4,byrow=T))
  boxplot(Ytot~scenario2, sd_dtoto[[sc_]], ylim=c(0,2000), main=sc_)
  boxplot(Yprop1~scenario2, sd_dtoto[[sc_]], ylim=c(0,1), main=sc_)
  boxplot(gini1~scenario2, sd_dtoto[[sc_]], ylim=c(0,1), main=sc_)
  boxplot(gini2~scenario2, sd_dtoto[[sc_]], ylim=c(0,1), main=sc_)
  boxplot(alive1~scenario2, sd_dtoto[[sc_]], ylim=c(0,500), main=sc_)
  
}


#pour definir niveaux de facteur et leur ordre!
lsLum <- c(46,45,44,1,50,51,52)
lsN <- c(49,48,47,1,53,54,55)
lsall <- c(40,39,38,1,41,42,43)

Ntrait <-"N-"#"N+"#"N/2"# 
sc_ <- "SD1-1"#"SD23-23"#"SD24-24"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsLum,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsLum)) #remise en ordre des niveaux de facteur!
boxplot(Yprop1~scenario2, x, ylim=c(0,1), main=paste(sc_, Ntrait))


sc_ <- "SD1-1"#"SD26-26"#"SD25-25"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsN,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsN)) #remise en ordre des niveaux de facteur!
boxplot(Yprop1~scenario2, x, ylim=c(0,1), main=paste(sc_, Ntrait))


sc_ <-"SD4-4"#"SD1-1"# "SD3-3"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsall,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsall)) #remise en ordre des niveaux de facteur!
boxplot(Yprop1~scenario2, x, ylim=c(0,1), main=paste(sc_, Ntrait))



#gini
sc_ <-"SD1-1"#"SD4-4"# "SD3-3"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsall,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsall)) #remise en ordre des niveaux de facteur!
boxplot(gini1~scenario2, x, ylim=c(0,1), main=paste(sc_, Ntrait))
boxplot(gini2~scenario2, x, ylim=c(0,1), main=paste(sc_, Ntrait))

#effet marque


#alive
sc_ <-"SD1-1"#"SD4-4"# "SD3-3"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsall,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsall)) #remise en ordre des niveaux de facteur!
boxplot(alive1~scenario2, x, ylim=c(0,400), main=paste(sc_, Ntrait))
#change peu 



#Cor_ParamAllNorm
sc_ <-"SD1-1"#"SD4-4"# "SD3-3"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsall,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsall)) #remise en ordre des niveaux de facteur!
boxplot(Cor_ParamAllNorm~scenario2, x, ylim=c(-1,1), main=paste(sc_, Ntrait))
#change peu 

Ntrait <-"N+"#"N/2"# "N-"#
sc_ <- "SD1-1"#"SD24-24"#"SD23-23"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsLum,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsLum)) #remise en ordre des niveaux de facteur!
boxplot(gini1~scenario2, x, ylim=c(0,1), main=paste(sc_, Ntrait))

#"SD26-26"#"SD25-25"#
Ntrait <-"N+"#"N/2"# "N-"#
sc_ <- "SD23-23"#"SD24-24"#"SD1-1"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsLum,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsLum)) #remise en ordre des niveaux de facteur!
boxplot(Fix0_Cor_ParamAllNorm~scenario2, x, ylim=c(-1,1), main=paste(sc_, Ntrait))
#ok a prendre
#tres domine -> correlation retends vers zero -> plus de selection?
dec5_phyllochron_Fix0
Fix0_Cor_ParamAllNorm
Fix0_Cor_phyllochron


Ntrait <-"N+"#"N/2"# "N-"#
sc_ <- "SD1-1"#"SD24-24"#"SD23-23"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsLum,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsLum)) #remise en ordre des niveaux de facteur!
boxplot(Cor_PARinonKin~scenario2, x, ylim=c(-1,1), main=paste(sc_, Ntrait))
#Ok a prendre!


Ntrait <-"N+"#"N/2"# "N-"#
sc_ <- "SD24-24"#"SD23-23"#"SD1-1"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsLum,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsLum)) #remise en ordre des niveaux de facteur!
boxplot(Fix1_Cor_diffMvoisNorm~scenario2, x, ylim=c(-1,1),main=paste(sc_, Ntrait))

#correlation genet diminue pour ecart interspe??

#dec5_phyllochron_Fix0
#Cor_ParamLightNorm
#Cor_ParamNNorm
#Fix1_Cor_diffMvoisNorm
#Fix1_Cor_MScumvois
#Fix1_Cor_ratioLight
#Fix1_Cor_ratioNupt
#Yprop1


Ntrait <-"N+"#"N/2"# "N-"#
sc_ <- "SD26-26"#"SD25-25"##"SD1-1"#
x <- dtoto[dtoto$sd== sc_ & dtoto$scenario2 %in% lsN,]
x$scenario2 <- factor(x$scenario2,levels=as.character(lsN)) #remise en ordre des niveaux de facteur!
boxplot(Fix1_Cor_ParamAllNorm~scenario2, x, ylim=c(-1,1), main=paste(sc_, Ntrait))
#ok a prendre
#tres domine -> correlation retends vers zero -> plus de selection?
dec5_phyllochron_Fix0
Fix0_Cor_ParamAllNorm
Fix0_Cor_phyllochron



factor(x$scenario2)
#
dtotosaveN <- dtoto
dtoto <- dtotosaveN


dtotosave <- dtoto
dtoto <- dtotosave


dtotosaveN2 <- dtoto
dtoto <- dtotosaveN2


dtotosave1pN <- dtoto



sp1 <- names(sp_tabSD)[1]
x1 <- sp_tabSD[[sp1]]


sp2 <- names(sp_tabSD)[2]
x2 <- sp_tabSD[[sp2]]

layout(matrix(1:12,2,6))
for (param_ in ls_parsd)
{
  hist(x1[,c(param_)], main=paste(param_,sp1))
  hist(x2[,c(param_)], main=paste(param_,sp2))
}


names(dtoto[250,])
dtoto[225:252,c("scenario" , "sd", "seed","YEsp1", "YEsp2", "Yprop1", "Yprop2")]
#calculs dans le dtoto sont bons!



dtoto$scenario

names(x)

hist(x1$Len)



########### test de lecture des fichiers de sortie spatiaux

couleurs <- colorRampPalette(c('white', 'red')) #default
couleursR <- colorRampPalette(c('white', 'red'))
couleursB <- colorRampPalette(c('white', 'blue'))
couleursV <- colorRampPalette(c('white', 'green'))
lscol10 <- colorRampPalette(c("blue", "red"))( 11 ) #palette de couleur des deciles

xx <- NULL
yy <- NULL
for(i in 1:16)
{
  xx <- c(xx, rep(i,16))
  yy <- c(yy, 16:1)#c(yy, 1:16)#
}



ls_tab <- ls_files[grepl('savetab', ls_files) & !grepl('.zip', ls_files)]

key <- names(sp_dtoto)[244]#[260]#
nomtab <- paste("savetabindices ", key,".csv")
tabindices <- read.table(nomtab, sep=";",header=T)
tabindices$x <- xx
tabindices$y <- yy
tabindices$phyllochron <- 1/tabindices$phyllochron

ls_resNorm <- calc_norm_par(tabindices ,lspar=c("Len"), plot_=T, main_ = "1")
ls_resNorm <- calc_norm_par(tabindices ,lspar=c("Len","Lfeuille","PPtreshh"), plot_=T, main_ = "3")
ls_resNorm <- calc_norm_par(tabindices ,lspar=c("phyllochron","Vmax2","ELmax"), plot_=T, main_ = "3")
ls_resNorm <- calc_norm_par(tabindices ,lspar=c("Len","Lfeuille","PPtreshh","phyllochron","Vmax2","ELmax"), plot_=T, main_ = "6")
ls_resNorm <- calc_norm_par(tabindices ,lspar=c("ParamAllNorm"), plot_=T, main_ = "AllNorm")

#pb avec phyllochron tjrs plus grand??? transfo en 1/x?? -> oui
# pb avec variance faible avec prametre qui ont effet inverse dans visualisation



#examiner correlation entre variables par pairs
vars <- c("MStot_fin", "ParamAllNorm","ParamLightNorm","ParamNNorm","ratioLight","ratioNupt","PARinonKin","NuptakenonKin")

pairs(tabindices[,vars], lower.panel = panel.smooth, upper.panel = panel.cor,gap=0, row1attop=FALSE, main=paste(key,"Fix1"))


range_ <- seq(0,16,0.5)
hist(tabindices[,"PARinonKin"], breaks=range_,col="grey")
hist(tabindices[tabindices$name=="Fix1","PARinonKin"], breaks=range_,add=T,col="red")
#hist(tabindices[tabindices$name=="Fix0","PARinonKin"], breaks=range_,add=T,col="white")



#test figure pour voir les correlation
#calcul des moyennes de correlation par traitement

var_ <- "Cor_PARinonKin" 
x <- dtoto[,c("scenario","sd",var_)]
x$trait <- paste(x$scenario,x$sd)

res <- by(x[,var_], x$trait, mean)
df <- data.frame(trait=names(res), cor=as.numeric(res))

noms_trait <- as.data.frame(t(as.data.frame(strsplit(as.character(df$trait), " "))))
names(noms_trait) <- c("scenario","sd")
df$scenario <- as.character(noms_trait$scenario)
df$sd <- as.character(noms_trait$sd)


#faire des plotde taille et couleur normalise!
#dans quel ordre??




########### boucle test des analyses indiv



#test plot correlation voisin / paramètres / coupe1???


#pour calcul voisins d'ordre 1 (All / kin / noKin)
cote <- 16
nblignes <- 16

ls_idvois <- vector("list", length=(cote*nblignes))
names(ls_idvois) <- 1:(cote*nblignes)
ls_idKin <- vector("list", length=(cote*nblignes))
names(ls_idKin) <- 1:(cote*nblignes)
ls_idnonKin <- vector("list", length=(cote*nblignes))
names(ls_idnonKin) <- 1:(cote*nblignes)

for (i in 1:(cote*nblignes))
{
  idvois <- ls_idvois_ordre1(i-1, cote, nblignes) +1 # appel avec i-1 (pour comme nump python) # ajout 1 a sortie pour rang R
  ls_idvois[[i]] <- idvois 
  ls_idKin[[i]] <- idvois[c(1,3,6,8)]
  ls_idnonKin[[i]] <- idvois[c(2,4,5,7)]
}


#coordonnes des pts pour 16*16
#maintenant sorti dans fichier paramSD du model!
#matrix(1:(16*16),16,16)

#xx <- NULL
#yy <- NULL

#for(i in 1:16)
#{
#  xx <- c(xx, rep(i,16))
#  yy <- c(yy, 16:1)#c(yy, 1:16)#
#}
#dfcoord <- data.frame(x=xx, y=yy, MStot_fin, MSnorm=MStot_fin/max(MStot_fin),nump=(1:256)-1)
#Esp <- substring(names(dat)[3:dim(dat)[2]],1,4)
#dfcoord$Esp <- Esp

ls_res_cor <- vector("list", length=length(names(sp_dtoto)))
names(ls_res_cor) <- names(sp_dtoto)


for (key in names(sp_dtoto))
{
  
  
  
  #layout(matrix(1:5,1,5))
  #key <- names(sp_dtoto)[260]#[330]#[16]#[3]#[31]#[19]#
  ls_toto_paquet <- sp_dtoto[[key]]$name
  ltoto <- read_ltoto(ls_toto_paquet)
  #names(ltoto[[ls_toto_paquet]])
  #ltoto[[ls_toto_paquet]]$V1
  dat <- ltoto[[ls_toto_paquet]]
  nb <- dim(dat)[2]-2
  
  lspar <-  c("Len","Lfeuille","phyllochron", "Vmax2", "ELmax", "PPtreshh")
  param_name <- "phyllochron"#"Len"#ls_par[1]
  resread <- read_lsSD_MStot(ltoto, ls_paramSD, param_name)
  sp_tabSD <- split(resread[["ls_tabSD"]][[1]], resread[["ls_tabSD"]][[1]]$name)
  MStot <- resread[["ls_MStot"]][[1]]
  #res <- BuildResDecil(MStot, sp_tabSD)
  
  #c(187,229,282,334) 
  MStot_ini <- as.numeric(MStot[60,])#30
  MStot_coupe1 <- as.numeric(MStot[127,])#65
  MStot_coupe2 <- as.numeric(MStot[169,])#100
  MStot_coupe3 <- as.numeric(MStot[222,])#150
  MStot_fin <- as.numeric(MStot[dim(MStot)[1],])
  #MStot_coupe4 <- as.numeric(MStot[200,])
  #MStot_coupe5 <- as.numeric(MStot[250,])
  #hist(MStot_fin, main=key)
  #hist(MStot_coupe1, main=key)
  
  #temptab <- resread[["ls_tabSD"]][[1]][, c("nump","name","x","y","retard","Len","Vmax2","ELmax","Lfeuille","phyllochron","PPtreshh")]
  temptab <- resread[["ls_tabSD"]][[1]][, c("nump","name","retard","Len","Vmax2","ELmax","Lfeuille","phyllochron","PPtreshh")]
  
  
  #ordonne dans l'ordre des nump!!
  temptab <- temptab[order(temptab$nump),]
  
  
  #retard <- temptab$retard
  Val_param <- temptab[,c(param_name)]#temptab$phyllochron
  #Val_param <- temptab$Len
  #hist(Val_param, main=key)
  
  #cacuk= de la valeur normalisee des parametres (multi-trait)
  temptab$phyllochron[temptab$phyllochron<8] <- 8 #pour les valeur <0 mise a 10-10!
  temptab$phyllochron <- 1/(temptab$phyllochron)
  temptab$PPtreshh <- 24-temptab$PPtreshh
  ParamAllNorm <- calc_norm_par(temptab[,lspar] ,lspar, plot_=F)$mean_norm_par
  temptab$ParamAllNorm <- ParamAllNorm 
  
  lightPar <- c("Len","Lfeuille","phyllochron")
  ParamLightNorm <- calc_norm_par(temptab[,lightPar] ,lightPar, plot_=F)$mean_norm_par
  temptab$ParamLightNorm <- ParamLightNorm
  NPar <- c("Vmax2", "ELmax", "PPtreshh")
  ParamNNorm <- calc_norm_par(temptab[,NPar] ,NPar, plot_=F)$mean_norm_par
  temptab$ParamNNorm <- ParamNNorm
  
  
  
  #calcul des moyenne des voisins
  x <- temptab[,c(lspar,"ParamAllNorm","ParamLightNorm","ParamNNorm")]
  #transforme param phyllochrone et PPtreshh pour avoir effet positif pour valeur croissante
  #x$phyllochron <- 1/x$phyllochron
  #x$PPtreshh <- 24-x$PPtreshh
  
  resN <- calc_neighb_param(x,c(lspar,"ParamAllNorm","ParamLightNorm","ParamNNorm"), ls_idvois, ls_idKin, ls_idnonKin)
  temptab <- cbind(temptab, resN)
  #caluler les difference pour sp1 et sp2
  temptab$diffMvoisNorm <- temptab$ParamAllNormMvois - temptab$ParamAllNorm
  temptab$diffMKinNorm <- temptab$ParamAllNormMKin - temptab$ParamAllNorm
  temptab$diffMnonKinNorm <- temptab$ParamAllNormMnonKin - temptab$ParamAllNorm
  temptab$diffMvoisLightNorm <- temptab$ParamLightNormMvois - temptab$ParamLightNorm
  temptab$diffMvoisNNorm <- temptab$ParamNNormMvois - temptab$ParamNNorm
  
  
  
  #recup PARiPlante et N uptake plante et faire cumul
  PARi <- dat[dat$V1=='PARiPlante',3:(3+nb-1)] #
  for (i in 1:nb) {PARi[,i] <- cumsum(PARi[,i])}
  Nuptake <- dat[dat$V1=='Nuptake_sol',3:(3+nb-1)] #sans fixation!!!
  for (i in 1:nb) {Nuptake[,i] <- cumsum(Nuptake[,i])}
  PARi_fin <- as.numeric(PARi[dim(PARi)[1],])
  Nuptake_fin <- as.numeric(Nuptake[dim(Nuptake)[1],])
  
  #plot(PARi_fin, MStot_fin)
  #plot(Nuptake_fin , MStot_fin)
  
  
  
  #df <- data.frame(MStot_ini, MStot_fin, MStot_coupe1, MStot_coupe2, MStot_coupe3, MStot_coupe4, MStot_coupe5)
  
  
  
  #calcul du cumul de biomasse, note moyenne, uptake des voisins
  MScumvois <- NULL
  MScumKin <- NULL
  MScumnonKin <- NULL
  PARivois <- NULL
  PARiKin <- NULL
  PARinonKin <- NULL
  Nuptakevois <- NULL
  NuptakeKin <- NULL
  NuptakenonKin <- NULL
  for (i in 1:(cote*nblignes))
  {
    #ALL ordre 1
    MSvois <- sum(MStot_fin[ls_idvois[[i]]])
    MScumvois <- cbind(MScumvois, MSvois)
    PARivois <- cbind(PARivois, sum(PARi_fin[ls_idvois[[i]]]) )
    Nuptakevois <- cbind(Nuptakevois, sum(Nuptake_fin[ls_idvois[[i]]]) )
    
    #Kin/NonKin
    MScumKin <- cbind(MScumKin, sum(MStot_fin[ls_idKin[[i]]]))
    MScumnonKin <- cbind(MScumnonKin, sum(MStot_fin[ls_idnonKin[[i]]]))
    PARiKin <- cbind(PARiKin, sum(PARi_fin[ls_idKin[[i]]]) )
    NuptakeKin <- cbind(NuptakeKin, sum(Nuptake_fin[ls_idKin[[i]]]) )
    PARinonKin <- cbind(PARinonKin, sum(PARi_fin[ls_idnonKin[[i]]]) )
    NuptakenonKin <- cbind(NuptakenonKin, sum(Nuptake_fin[ls_idnonKin[[i]]]) )
    
  }
  MScumvois <- as.numeric(MScumvois)
  PARivois <- as.numeric(PARivois)
  Nuptakevois <- as.numeric(Nuptakevois)
  MScumKin <- as.numeric(MScumKin)
  MScumnonKin <- as.numeric(MScumnonKin)
  PARiKin <- as.numeric(PARiKin)
  NuptakeKin <- as.numeric(PARiKin)
  PARinonKin <- as.numeric(PARinonKin)
  NuptakenonKin <- as.numeric(PARinonKin)

  
  dfMS <- data.frame(nump=temptab$nump, MStot_fin, MStot_ini, MStot_coupe1,MStot_coupe2,MStot_coupe3,PARi=PARi_fin, Nuptake=Nuptake_fin, MScumvois, MScumKin, MScumnonKin, PARivois, PARiKin, PARinonKin, Nuptakevois, NuptakeKin, NuptakenonKin)
  #ratio de capture des ressources avec voisins
  dfMS$ratioLight <- PARi_fin/PARivois
  dfMS$ratioNupt <- Nuptake_fin/Nuptakevois
  #
  #EcardPotentiel <-  MStot_fin/mean(MStot_fin) - ParamAllNorm #pas tres logique en multitrait / simple trait
  
  
  temptab <- merge(temptab, dfMS, by="nump")
  
  #correlation MSindiv avec valeur des parametres / valeur des voisins / ecart des voisins / ressources / ressources des voisins
  #subx <- temptab[,5:dim(temptab)[2]]#avec x,y
  subx <- temptab[,3:dim(temptab)[2]]#sans x,y
  rescor <- as.data.frame(cor(subx))
  valcorAll <- rescor$MStot_fin
  #barplot(valcorAll, names.arg =row.names(rescor), las=3,cex.names=0.6,main=key)
  
  #faire un data.frame de ca
  res <- data.frame(t(valcorAll))
  names(res) <- paste("Cor_", row.names(rescor),sep="")
  
  #Corr par espece
  s_temp <- split(temptab, temptab$name)
  sp <- names(s_temp)[1]#"Fix0"
  #subx <- s_temp[[sp]][,5:dim(temptab)[2]]#avec x,y
  subx <- s_temp[[sp]][,3:dim(temptab)[2]]#sans x,y
  rescor <- as.data.frame(cor(subx))
  valcorSp1 <- rescor$MStot_fin
  res1 <- data.frame(t(valcorSp1))
  names(res1) <- paste(sp,"_Cor_", row.names(rescor),sep="")
  
  sp <- names(s_temp)[2]#"Fix1"
  #subx <- s_temp[[sp]][,5:dim(temptab)[2]]#avec x,y
  subx <- s_temp[[sp]][,3:dim(temptab)[2]]#sans x,y
  rescor <- as.data.frame(cor(subx))
  valcorSp2 <- rescor$MStot_fin
  res2 <- data.frame(t(valcorSp2))
  names(res2) <- paste(sp,"_Cor_", row.names(rescor),sep="")
  ##barplot(valcorSp2, names.arg =row.names(rescor), las=3,cex.names=0.6,main=key)
  #plot(subx$PPtreshh, subx$MStot_fin)
  #plot(temptab$PARinonKin, temptab$MStot_fin)
  
  res <- cbind(res,res1,res2)
  #res$key <- key
  res <- cbind(dtoto[dtoto$key==key,], res)
  write.table(res, paste("savedic ", key,".csv"), sep=";", col.names = T, row.names = F)
  
  ls_res_cor[[key]] <- res
  
  
}

#res_cor <- do.call("rbind", ls_res_cor)
#write.table(res_cor, paste("save_dtoto3.csv"), sep=";", col.names = T, row.names = F)

#dtoto1 <- res_cor
#dtoto2 <- res_cor
#dtoto3 <- res_cor
#dtoto <- dtoto2

  #!! gere seulement 1 paramtre
  df <- data.frame(Val_param,ParaMvois,ParaMKin,ParaMnonKin,MStot_fin,MScumvois,MScumnonKin)
                   #PARivois, MScumnonKin,MScumvois, MStot_ini, MStot_fin, MStot_coupe1, MStot_coupe2, MStot_coupe3, MStot_coupe4, MStot_coupe5)
  
    
  #hist(MScumvois, main=key)
  #hist(ParaMvois, main=key)
  #hist(MScumKin)
  
  #plot(MStot_coupe1, MStot_fin)
  #plot(Val_param, MStot_fin, main=key)
  
  #plot(MScumvois, MStot_fin, main=key)
  #mod_ <- lm( MStot_fin~MScumvois)
  #abline(mod_)
  #summary(mod_)
  
  #plot(MScumvois, MStot_coupe1, main=key)
  #mod_ <- lm( MStot_coupe1~MScumvois)
  #abline(mod_)
  #summary(mod_)
  
  # !!
  # pb decalage de 1 avec les nump!
  #l'appeler avec id-1 et ajouter 1 au sorties!
  #??MStot coupe 1 presque egal MStot fin?? -> pas bon id
  

  
  #pairs(df, lower.panel = panel.smooth, upper.panel = panel.cor,gap=0, row1attop=FALSE)
  
  

  #faire un graph visu spatiale des MStot
  
  
  
  df <- data.frame(retard,Val_param,ParaMvois,PARivois, MScumnonKin,MScumvois, MStot_ini, MStot_fin, MStot_coupe1, MStot_coupe2, MStot_coupe3)
  pairs(df, lower.panel = panel.smooth, upper.panel = panel.cor,gap=0, row1attop=FALSE, main=key)
  df$Esp <- Esp

  
  
  lscols <- couleursR(101)
  cols_ <- col100((dfcoord$MSnorm-0.00001)*100, lscols)
  cols_[is.na(cols_)] <- lscols[101] #cas superieur a 1
  
  
  plot(dfcoord$x, dfcoord$y, cex=4*dfcoord$MSnorm,pch=16,col=cols_, main=key)
  points(dfcoord$x, dfcoord$y, cex=4*(1/Val_param)/max(1/Val_param), col="blue")
  #points(dfcoord$x, dfcoord$y, cex=4*(Val_param)/max(Val_param), col="blue")

  #plot en separant les esp
  sp_dfcoord <- split(dfcoord, dfcoord$Esp)
  
  lscols <- couleursR(101)
  cols_ <- col100((sp_dfcoord[["Fix0"]]$MSnorm-0.00001)*100, lscols)
  cols_[is.na(cols_)] <- lscols[101] #cas superieur a 1
  plot(sp_dfcoord[["Fix0"]]$x, sp_dfcoord[["Fix0"]]$y, cex=4*sp_dfcoord[["Fix0"]]$MSnorm,pch=16,col=cols_, main=key)
  
  lscols <- couleursV(101)
  cols_ <- col100((sp_dfcoord[["Fix1"]]$MSnorm-0.00001)*100, lscols)
  cols_[is.na(cols_)] <- lscols[101] #cas superieur a 1
  points(sp_dfcoord[["Fix1"]]$x, sp_dfcoord[["Fix1"]]$y, cex=4*sp_dfcoord[["Fix1"]]$MSnorm,pch=16,col=cols_)
  
  #pairs en separant les especes
  sp_df <- split(df, df$Esp)
  pairs(sp_df[["Fix0"]][1:13], lower.panel = panel.smooth, upper.panel = panel.cor,gap=0, row1attop=FALSE, main=paste(key,"Fix0"))
  pairs(sp_df[["Fix1"]][1:13], lower.panel = panel.smooth, upper.panel = panel.cor,gap=0, row1attop=FALSE, main=paste(key,"Fix1"))
  
  
  #plot Area
  sp <- names(sp_tabSD)[1]
  res <- Build_EvolProportions(MStot, sp_tabSD, sp)
  don <- res[,2:12] #t et les deciles
  titre <- paste(sp, key)#esps,
  My_AreaPlot(don, titre=titre, xlab="t", ylab=paste("decile ", param_name), lscol=rev(lscol10))
  
  
  sp <- names(sp_tabSD)[2]
  res <- Build_EvolProportions(MStot, sp_tabSD, sp)
  don <- res[,2:12] #t et les deciles
  titre <- paste(sp, key)#esps,
  My_AreaPlot(don, titre=titre, xlab="t", ylab=paste("decile ", param_name), lscol=rev(lscol10))
  
  
  
}


plot(MStot_fin, PARi_fin/(PARi_fin+PARivois))
plot(MStot_fin, PARi_fin/(PARi_fin+PARivois))
plot(MStot_fin, PARi_fin)

plot(Val_param/ParaMvois, PARi_fin/(PARi_fin+PARivois))


ratioValparam <- Val_param/ParaMvois
EcartValparam <- Val_param-ParaMvois
ratiolight <- PARi_fin/PARivois

RatioPotentielRealise <- dfcoord$MSnorm / ((1/Val_param)/max(1/Val_param))
EcardPotentielRealise <- dfcoord$MSnorm - ((1/Val_param)/max(1/Val_param))

df <- data.frame(RatioPotentielRealise,EcardPotentielRealise,ParaMvois, ratioValparam, EcartValparam,ratiolight)
pairs(df, lower.panel = panel.smooth, upper.panel = panel.cor,gap=0, row1attop=FALSE, main=key)



#to do kin/non kin parametre = a faire pour tous les parmetre (ou une liste)
#calculer valeur prametre moyen normalise!




#test stat spatailes
#calcul du Ripley's
library(spatstat)

?ppp #objet pour stat spatial
mypattern <- ppp(xx,yy,c(0.5,16.5),c(0.5,16.5))
plot(mypattern)
#ajout des donnees
df <- data.frame(Val_param, MScumvois,  MStot_fin)
marks(mypattern) <- MStot_fin#df
plot(Smooth(mypattern))
plot.ppp(mypattern, use.marks = TRUE, which.marks = "MStot_fin")
plot.ppp(mypattern, use.marks = TRUE, which.marks = "Val_param")


plot(Kest(mypattern)) #Ripley's K-function

#seulement poisition des points!!

#existe: mark-weighted K-function -> ca qui faut!
#https://www.routledgehandbooks.com/doi/10.1201/b16195-4
#http://www.stat.ucla.edu/~frederic/222/S19/SpatStatIntro.pdf
#Kmark a apriori!
K <- Kmark(mypattern)
plot(K,main=key)
summary(K)

#?Kmark
?pcf
#spatial pair correlation
spp <- pcf(K)
plot(spp, main=key)
spp[["pcf"]]
spp[["r"]]
#traduit quoi? le fait qu'il y a tjts 1 points ou valeur des points?


#
#X_ <- simdat
#p <- pcf(X_)
#plot(p)
#plot(X_)

#Homogene / autocorrele?




#temptab a seulement decile d'un parametre!

names(temptab)



#plot des param d'entree
lspar <- c("Len","Lfeuille","phyllochron", "Vmax2", "ELmax", "PPtreshh")
x <- temptab[,lspar]

#transforme param phyllochrone et PPtreshh pour avoir effet positif
x$phyllochron <- 1/x$phyllochron
x$PPtreshh <- 24-x$PPtreshh

#ajout des coord (a reprendre de temptab )
x$x <- dfcoord$x
x$y <- dfcoord$y


calc_norm_par <- function(tabpar,lspar, plot_=F)
{
  # Calcul les valeur normalisee (par la moyenne) des parametresSD et la moyenne des valeurs normalisee
  # avec plot_ a True et les coord x,y,  fait un graph de visu
  
  nbpar <- length(lspar)
  ls_resNorm <- vector("list", length=(nbpar+1))
  names(ls_resNorm) <- c(lspar, "mean_norm_par")
  
  ls_col_ <- 1:nbpar #a passer en argument eventuellement
                       
  
  Val_par <- tabpar[,c(lspar[1])]
  #normalise par la moyenne (de ce qui est donne en entree: communaute ou population)
  norm_par <- Val_par/mean(Val_par) 
  ls_resNorm[[lspar[1]]] <- norm_par
  mean_norm_par <- norm_par
  
  if (plot_ == T)
  {plot(tabpar$x, tabpar$y, cex=1.5*norm_par,col="blue", main="params")}
  
  for (i in 2:nbpar)
  {
    Val_par <- tabpar[,c(lspar[i])]
    norm_par <- Val_par/mean(Val_par)
    ls_resNorm[[lspar[i]]] <- norm_par
    mean_norm_par <- mean_norm_par+norm_par
    
    if (plot_ == T)
    {points(tabpar$x, tabpar$y, cex=1.5*norm_par,col=ls_col_[i])}
    
  }
  ls_resNorm[["mean_norm_par"]] <- mean_norm_par/nbpar
  names(ls_resNorm)[1:nbpar] <- paste(lspar,"Norm", sep="")
  ls_resNorm <- as.data.frame(ls_resNorm)
  ls_resNorm
}

ls_resNorm <- calc_norm_par(x ,lspar, plot_=F)
ParamNorm <- calc_norm_par(temptab[,lspar] ,lspar, plot_=F)$mean_norm_par




#ecart de parmetre des voisin pour tous les parametres

lspar <- c("Len","Lfeuille","phyllochron", "Vmax2", "ELmax", "PPtreshh")
x <- temptab[,lspar]

#transforme param phyllochrone et PPtreshh pour avoir effet positif
x$phyllochron <- 1/x$phyllochron
x$PPtreshh <- 24-x$PPtreshh



#faire une boucle pour tous les param + moyenne de all

calc_neighb_param <- function(tabpar,lspar, ls_idvois, ls_idKin, ls_idnonKin)
{
  #calculate average parameter value of order 1 neighbours, with kin and non kin in a binary mixture
  #can be used for any vector and any list of lspar (not only parameters)
  
  nbpar <- length(lspar)
  ls_res <- vector("list", length=nbpar)
  names(ls_res) <- c(lspar)
  
  for (param_name in lspar)
  {
    
    #param_name <- lspar[1]
    Val_param <- tabpar[,c(param_name)]
    
    
    ParaMvois <- NULL
    ParaMKin <- NULL
    ParaMnonKin <- NULL
    for (i in 1:(cote*nblignes))
    {
      #ALL ordre 1
      ParaMvois <- cbind(ParaMvois, mean(Val_param[ls_idvois[[i]]]))
      #Kin/NonKin
      ParaMKin <- cbind(ParaMKin, mean(Val_param[ls_idKin[[i]]]) )
      ParaMnonKin <- cbind(ParaMnonKin, mean(Val_param[ls_idnonKin[[i]]]) )
    }
    ParaMvois <- as.numeric(ParaMvois)
    ParaMKin <- as.numeric(ParaMKin)
    ParaMnonKin <- as.numeric(ParaMnonKin)
    
    res <- data.frame(ParaMvois, ParaMKin, ParaMnonKin)
    names(res) <- paste(param_name, c("Mvois", "MKin", "MnonKin"), sep="")
    
    
    ls_res[[param_name]] <- res
  
  }
  
  res <- as.data.frame(ls_res)
  names(res) <- as.character(as.data.frame(t(as.data.frame(strsplit(names(res),"\\."))))$V2)
  
  res
}
#pourrait prevoir de mettre ls_idvois, ls_idKin, ls_idnonKin dans la table d'entree au prelablable
#calc_neighb_param(x,lspar, ls_idvois, ls_idKin, ls_idnonKin)
#avoir une option mean/sum?
#en theorie marche aussi pour ordre 2 ou> si donne les bonnes listes d'ID!


#mettre en application



#recuperer les coordonnee x et y dans les fichiers!!!
#puis relancer les simuls



Calc_MSindiv_Corr <- function(ltoto, ls_toto_paquet, ls_paramSD, lspar=c("Len","Lfeuille","phyllochron", "Vmax2", "ELmax", "PPtreshh"))
{
  ## fonction pour mettre en forme valeurs de parametres normalise, valeur d'effet de voisinnage, et calculer les correlations entre indices
  
  #key <- names(sp_dtoto)[260]#[330]#[16]#[3]#[31]#[19]#
  #ls_toto_paquet <- sp_dtoto[[key]]$name
  #ltoto <- read_ltoto(ls_toto_paquet)
  #names(ltoto[[ls_toto_paquet]])
  #ltoto[[ls_toto_paquet]]$V1
  dat <- ltoto[[ls_toto_paquet]]
  nb <- dim(dat)[2]-2
  
  #lspar <-  c("Len","Lfeuille","phyllochron", "Vmax2", "ELmax", "PPtreshh")
  param_name <- "phyllochron"#"Len"#ls_par[1] #pour exemple, pas utilise
  resread <- read_lsSD_MStot(ltoto, ls_paramSD, param_name)
  sp_tabSD <- split(resread[["ls_tabSD"]][[1]], resread[["ls_tabSD"]][[1]]$name)
  MStot <- resread[["ls_MStot"]][[1]]
  #res <- BuildResDecil(MStot, sp_tabSD)
  
  #c(187,229,282,334) #dates de coupes fixes
  MStot_ini <- as.numeric(MStot[60,])#30
  MStot_coupe1 <- as.numeric(MStot[127,])#65
  MStot_coupe2 <- as.numeric(MStot[169,])#100
  MStot_coupe3 <- as.numeric(MStot[222,])#150
  MStot_fin <- as.numeric(MStot[dim(MStot)[1],])
  #MStot_coupe4 <- as.numeric(MStot[200,])
  #MStot_coupe5 <- as.numeric(MStot[250,])
  #hist(MStot_fin, main=key)
  #hist(MStot_coupe1, main=key)
  
  #temptab <- resread[["ls_tabSD"]][[1]][, c("nump","name","x","y","retard","Len","Vmax2","ELmax","Lfeuille","phyllochron","PPtreshh")]
  temptab <- resread[["ls_tabSD"]][[1]][, c("nump","name","retard","Len","Vmax2","ELmax","Lfeuille","phyllochron","PPtreshh")]
  
  
  #ordonne dans l'ordre des nump!!
  temptab <- temptab[order(temptab$nump),]
  
  
  #Val_param <- temptab[,c(param_name)]#temptab$phyllochron
  #Val_param <- temptab$Len
  #hist(Val_param, main=key)
  
  #calcul de la valeur normalisee des parametres (multi-trait)
  temptab$phyllochron[temptab$phyllochron<8] <- 8 #pour les valeur <0 mise a 10-10!
  temptab$phyllochron <- 1/(temptab$phyllochron)
  temptab$PPtreshh <- 24-temptab$PPtreshh
  ParamAllNorm <- calc_norm_par(temptab[,lspar] ,lspar, plot_=F)$mean_norm_par
  temptab$ParamAllNorm <- ParamAllNorm 
  
  #agrege par Light / N (specifique papier beatrice)
  lightPar <- c("Len","Lfeuille","phyllochron")
  ParamLightNorm <- calc_norm_par(temptab[,lightPar] ,lightPar, plot_=F)$mean_norm_par
  temptab$ParamLightNorm <- ParamLightNorm
  NPar <- c("Vmax2", "ELmax", "PPtreshh")
  ParamNNorm <- calc_norm_par(temptab[,NPar] ,NPar, plot_=F)$mean_norm_par
  temptab$ParamNNorm <- ParamNNorm
  
  
  
  #calcul des moyenne des voisins
  x <- temptab[,c(lspar,"ParamAllNorm","ParamLightNorm","ParamNNorm")]
  #transforme param phyllochrone et PPtreshh pour avoir effet positif pour valeur croissante

  
  resN <- calc_neighb_param(x,c(lspar,"ParamAllNorm","ParamLightNorm","ParamNNorm"), ls_idvois, ls_idKin, ls_idnonKin)
  temptab <- cbind(temptab, resN)
  #caluler les difference pour sp1 et sp2
  temptab$diffMvoisNorm <- temptab$ParamAllNormMvois - temptab$ParamAllNorm
  temptab$diffMKinNorm <- temptab$ParamAllNormMKin - temptab$ParamAllNorm
  temptab$diffMnonKinNorm <- temptab$ParamAllNormMnonKin - temptab$ParamAllNorm
  temptab$diffMvoisLightNorm <- temptab$ParamLightNormMvois - temptab$ParamLightNorm
  temptab$diffMvoisNNorm <- temptab$ParamNNormMvois - temptab$ParamNNorm
  

  #recup PARiPlante et N uptake plante et faire cumul
  PARi <- dat[dat$V1=='PARiPlante',3:(3+nb-1)] #
  for (i in 1:nb) {PARi[,i] <- cumsum(PARi[,i])}
  Nuptake <- dat[dat$V1=='Nuptake_sol',3:(3+nb-1)] #sans fixation!!!
  for (i in 1:nb) {Nuptake[,i] <- cumsum(Nuptake[,i])}
  PARi_fin <- as.numeric(PARi[dim(PARi)[1],])
  Nuptake_fin <- as.numeric(Nuptake[dim(Nuptake)[1],])
  
  
  #calcul du cumul de biomasse, note moyenne, uptake des voisins
  MScumvois <- NULL
  MScumKin <- NULL
  MScumnonKin <- NULL
  PARivois <- NULL
  PARiKin <- NULL
  PARinonKin <- NULL
  Nuptakevois <- NULL
  NuptakeKin <- NULL
  NuptakenonKin <- NULL
  for (i in 1:(cote*nblignes))
  {
    #ALL ordre 1
    MSvois <- sum(MStot_fin[ls_idvois[[i]]])
    MScumvois <- cbind(MScumvois, MSvois)
    PARivois <- cbind(PARivois, sum(PARi_fin[ls_idvois[[i]]]) )
    Nuptakevois <- cbind(Nuptakevois, sum(Nuptake_fin[ls_idvois[[i]]]) )
    
    #Kin/NonKin
    MScumKin <- cbind(MScumKin, sum(MStot_fin[ls_idKin[[i]]]))
    MScumnonKin <- cbind(MScumnonKin, sum(MStot_fin[ls_idnonKin[[i]]]))
    PARiKin <- cbind(PARiKin, sum(PARi_fin[ls_idKin[[i]]]) )
    NuptakeKin <- cbind(NuptakeKin, sum(Nuptake_fin[ls_idKin[[i]]]) )
    PARinonKin <- cbind(PARinonKin, sum(PARi_fin[ls_idnonKin[[i]]]) )
    NuptakenonKin <- cbind(NuptakenonKin, sum(Nuptake_fin[ls_idnonKin[[i]]]) )
    
  }
  MScumvois <- as.numeric(MScumvois)
  PARivois <- as.numeric(PARivois)
  Nuptakevois <- as.numeric(Nuptakevois)
  MScumKin <- as.numeric(MScumKin)
  MScumnonKin <- as.numeric(MScumnonKin)
  PARiKin <- as.numeric(PARiKin)
  NuptakeKin <- as.numeric(PARiKin)
  PARinonKin <- as.numeric(PARinonKin)
  NuptakenonKin <- as.numeric(PARinonKin)
  
  
  dfMS <- data.frame(nump=temptab$nump, MStot_fin, MStot_ini, MStot_coupe1,MStot_coupe2,MStot_coupe3,PARi=PARi_fin, Nuptake=Nuptake_fin, MScumvois, MScumKin, MScumnonKin, PARivois, PARiKin, PARinonKin, Nuptakevois, NuptakeKin, NuptakenonKin)
  #ratio de capture des ressources avec voisins
  dfMS$ratioLight <- PARi_fin/PARivois
  dfMS$ratioNupt <- Nuptake_fin/Nuptakevois
  #EcardPotentiel <-  MStot_fin/mean(MStot_fin) - ParamAllNorm #pas tres logique en multitrait / simple trait
  
  
  temptab <- merge(temptab, dfMS, by="nump")
  
  
  #correlation MSindiv avec valeur des parametres / valeur des voisins / ecart des voisins / ressources / ressources des voisins
  #subx <- temptab[,5:dim(temptab)[2]]#new: avec x,y
  subx <- temptab[,3:dim(temptab)[2]]#old: sans x,y
  rescor <- as.data.frame(cor(subx))
  valcorAll <- rescor$MStot_fin
  #barplot(valcorAll, names.arg =row.names(rescor), las=3,cex.names=0.6,main=key)
  
  #faire un data.frame de ca
  res <- data.frame(t(valcorAll))
  names(res) <- paste("Cor_", row.names(rescor),sep="")
  
  #Corr par espece
  s_temp <- split(temptab, temptab$name)
  sp <- names(s_temp)[1]#"Fix0"
  #subx <- s_temp[[sp]][,5:dim(temptab)[2]]#new: avec x,y
  subx <- s_temp[[sp]][,3:dim(temptab)[2]]#old: sans x,y
  rescor <- as.data.frame(cor(subx))
  valcorSp1 <- rescor$MStot_fin
  res1 <- data.frame(t(valcorSp1))
  names(res1) <- paste(sp,"_Cor_", row.names(rescor),sep="")
  
  sp <- names(s_temp)[2]#"Fix1"
  #subx <- s_temp[[sp]][,5:dim(temptab)[2]]#new: avec x,y
  subx <- s_temp[[sp]][,3:dim(temptab)[2]]#old: sans x,y
  rescor <- as.data.frame(cor(subx))
  valcorSp2 <- rescor$MStot_fin
  res2 <- data.frame(t(valcorSp2))
  names(res2) <- paste(sp,"_Cor_", row.names(rescor),sep="")
  ##barplot(valcorSp2, names.arg =row.names(rescor), las=3,cex.names=0.6,main=key)

  res <- cbind(res,res1,res2)
  res$key <- key

  res
  
  #renvoie aussi du tableau des donnees : temptab
  ls_resOK <- list(res, temptab)
  names(ls_resOK) <- c("tabCorMSindiv", "datIndices")
  ls_resOK
  
}
#distinguer 2 fonctions? voir 3?: mef et calcul des correlations?

ls_res_cor_i <- Calc_MSindiv_Corr(ltoto, ls_toto_paquet, ls_paramSD, lspar=c("Len","Lfeuille","phyllochron", "Vmax2", "ELmax", "PPtreshh"))
ls_res_cor_i[["tabCorMSindiv"]]
ls_res_cor_i[["datIndices"]]




