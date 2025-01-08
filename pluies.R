
library(geoR)
library(fields)
# library(akima)
# library(RandomFields)

dta <- read.table('/home/metienne/Cours/ENSAI/spatial/docs/Cours_liliane/geostat/sic_obs.dat',
                  col.names = c('id','x','y','pluies'),sep=',') |> 
    select(x,y, pluies) 
xmin <- min(dta$x)
ymin <- min(dta$y)

dta <- dta |>  
    mutate(x = (x-xmin)/1000, y = (y- ymin)/1000)

dta_full <-  read.table('/home/metienne/Cours/ENSAI/spatial/docs/Cours_liliane/geostat/sic_full.dat',
                        col.names = c('id','x','y','pluies'),sep=',',skip=6) |> 
    select(x,y, pluies) |> 
    mutate(x = (x-xmin)/1000, y = (y- ymin)/1000)

bords <- read.table('/home/metienne/Cours/ENSAI/spatial/docs/Cours_liliane/geostat/bords.txt',col.names = c('x','y'))|> 
    mutate(x = (x-xmin)/1000, y = (y- ymin)/1000)

elevation <- read.table('/home/metienne/Cours/ENSAI/spatial/docs/Cours_liliane/geostat/surfdem.grd',skip=6)

taille = dta_full$pluies/100
symbols(dta_full$x,dta_full$y,circles=taille,bg='gray',inches=FALSE,asp=1,
     xlab= 'x (km)',ylab ='y (km)',main='Pluies suisses après Tchernobyl',
     ylim=range (bords$y),xlim=range (bords$x))
lines(bords$x,bords$y)

xllcorner  =   -185556.375/1000 -xmin
yllcorner  =   -127261.5234375/1000 -ymin
cellsize    =  1009.975/1000

xelev = seq(xllcorner,cellsize*375+xllcorner,cellsize)
yelev = seq(yllcorner,cellsize*252+yllcorner,cellsize)

image.plot(xelev,yelev,t(as.matrix(elevation)),xlab= 'x (km)',ylab ='y (km)',
               col = terrain.colors(64),asp=1,main ='Altitude (m)',
               ylim=range (bords$y),xlim=range(bords$x))              
lines(bords$x,bords$y)

#------------
# variogrammemin(elevation)
#------------
# 1. Variogramme empirique
#-------------------------

geodata = as.geodata(dta)
m.d = 150                  # distance maximale
interv = seq(0,150,20)    # intervalles
p.m = 10                     # nombre minimal de paire

vario.c = variog(geodata,max.dist=m.d,pairs.min=p.m,uvec=interv,op="cloud")
plot(vario.c,main = "Nuée variographique")
vario.b = variog(geodata,max.dist=m.d,pairs.min=p.m,uvec=interv,messages.screen=FALSE)
plot(vario.b,main = "Variogramme empirique")

vario.bc = variog(geodata,max.dist=m.d,pairs.min=p.m,uvec=interv,bin.cloud=TRUE,messages.screen=FALSE)
plot(vario.bc,main = "Box-plot sur le variogramme empirique",bin.cloud=TRUE)

# 2. Ajustement du variogramme
#-----------------------------
c.m = "exponential"
i.c = c(15000,60)
varioest = variofit(vario.b,cov.model = c.m,minimisation.function = "nls",ini.cov.pars=i.c,fix.nugget=F,fix.kappa=TRUE,max.dist=vario.b$max.dist)
titre = paste("modele ",c.m,", portee =",round(varioest$cov.pars[2]*100)/100,", palier =",round(varioest$cov.pars[1]*100)/100)
plot(vario.b,main=titre)
lines(varioest)

   
# 3. Krigeage
#------------
Kfull = krige.conv(geodata,loc=cbind(donneesfull$x,donneesfull$y),krige=krige.control(type.krige="ok",obj.model=varioest))
RMSEP = sqrt(mean((Kfull$predict-donneesfull$pluies)^2))

grx = seq(-30,318,2)
gry = seq(-20,202,2)
grille = expand.grid(grx,gry)
K = krige.conv(geodata,loc=grille,krige=krige.control(type.krige="ok",obj.model=varioest))

fichier = paste(FicFig,3,sep="")
imprime(impression,fichier,form)
Z = matrix(K$predict,nrow=length(grx),ncol=length(gry),byrow=F)
titre = paste("Krigeage avec un modele",c.m)
image.plot(grx,gry,Z,xlab= 'x (km)',ylab ='y (km)',
           col=couleur,main =titre,zlim=c(0,600),
               ylim=range (bords$y),asp=1,xlim=range(bords$x))  
points(donnees$x,donnees$y,pch=20)            
lines(bords$x,bords$y)
ferme(impression)

fichier = paste(FicFig,4,sep="")
imprime(impression,fichier,form)
titre = "Ecart-type"
s=apply(cbind(K$krige.var,rep(0,length(K$krige.var))),1,max)
Sigma = sqrt(matrix(s,nrow=length(grx),ncol=length(gry),byrow=F))
titre = "Ecart-type de krigeage"
image.plot(grx,gry,Sigma,xlab= 'x (km)',ylab ='y (km)',
               col = tim.colors(64),asp=1,main =titre,
               ylim=range (bords$y),xlim=range(bords$x))              
lines(bords$x,bords$y)
ferme(impression)

# Anisotropie
#------------
vpsiA = seq(0,pi,l=13)
vpsiB = 1:10
mRMSEP = matrix(0,length(vpsiA),length(vpsiB))
for (ipsiA in 1:13){
    for (ipsiB in 1:10){
        paraniso = c(vpsiA[ipsiA],vpsiB[ipsiB])
        Kfull = krige.conv(geodata,loc=cbind(donneesfull$x,donneesfull$y),krige=krige.control(type.krige="ok",obj.model=varioest,aniso.pars=paraniso))
        mRMSEP[ipsiA,ipsiB] = sqrt(mean((Kfull$predict-donneesfull$pluies)^2))
    }
}

print(round(mRMSEP,1))
paraniso = c(pi/6,3)
Kaniso = krige.conv(geodata,loc=grille,krige=krige.control(type.krige="ok",obj.model=varioest,aniso.pars=paraniso))

fichier = paste(FicFig,5,sep="")
imprime(impression,fichier,form)
Z = matrix(Kaniso$predict,nrow=length(grx),ncol=length(gry),byrow=F)
titre = paste("Krigeage avec un modele",c.m,'anisotrope')
image.plot(grx,gry,Z,xlab= 'x (km)',ylab ='y (km)',zlim=c(0,600),
               col = couleur,asp=1,main =titre,
               ylim=range (bords$y),xlim=range(bords$x))     
points(donnees$x,donnees$y,pch=20)         
lines(bords$x,bords$y)
ferme(impression)

titre = "Ecart-type"
s=apply(cbind(K$krige.var,rep(0,length(K$krige.var))),1,max)
Sigma = sqrt(matrix(s,nrow=length(grx),ncol=length(gry),byrow=F))
titre = "Ecart-type de krigeage"
filled.contour(grx,gry,Sigma,levels=seq(0,200,20),main=titre,col =palette.perso(seq(0,200,20),0,200), plot.axes={ axis(1); axis(2); points(donnees$x,donnees$y,pch=20)})

# Avec l'elevation
#-----------------
paraniso = c(pi/6,3)
grElev = expand.grid(xelev,yelev)
elev = as.vector(t(as.matrix(elevation)))

elevobs = interpp(grElev[,1],grElev[,2],elev,donnees$x,donnees$y)
elevgr = interpp(grElev[,1],grElev[,2],elev,grille[,1],grille[,2])

sousind = seq(1,dim(grElev)[1],2)
Kelevaniso = krige.conv(geodata,loc=grille,
        krige=krige.control(type.krige="ok",obj.model=varioest,trend.d = ~ elevobs$z, 
        trend.l = ~elevgr$z,aniso.pars=paraniso))


fichier = paste(FicFig,6,sep="")
imprime(impression,fichier,form)
Z = matrix(Kelevaniso$predict,nrow=length(grx),ncol=length(gry),byrow=F)
titre = paste("Krigeage avec un modele",c.m,"avec l'altitude")
image.plot(grx,gry,Z,xlab= 'x (km)',ylab ='y (km)',
               col = couleur,asp=1,main =titre,zlim=c(0,600),
               ylim=range (bords$y),xlim=range(bords$x))              
lines(bords$x,bords$y)
points(donnees$x,donnees$y,pch=20)
ferme(impression)

V = variog(geodata,trend.d=~elevobs$z)
ml = likfit(geodata, trend=~elevobs$z,ini.cov.pars=i.c)
KC <- krige.control(trend.d=~elevobs$z, trend.l=~elevgr$z, obj=ml)
Keleviso = krige.conv(geodata,loc=grille,krige=KC)


fichier = paste(FicFig,"6b",sep="")
imprime(impression,fichier,form)
Z = matrix(Keleviso$predict,nrow=length(grx),ncol=length(gry),byrow=F)
titre = paste("Krigeage avec un modele",c.m,"avec l'altitude")
image.plot(grx,gry,Z,xlab= 'x (km)',ylab ='y (km)',
           col = couleur,asp=1,main =titre,zlim=c(0,600),
           ylim=range (bords$y),xlim=range(bords$x))              
lines(bords$x,bords$y)
points(donnees$x,donnees$y,pch=20)
ferme(impression)

Z = matrix(Keleviso$predict-K$predict,nrow=length(grx),ncol=length(gry),byrow=F)
image.plot(grx,gry,Z,xlab= 'x (km)',ylab ='y (km)',
               col = couleur,asp=1,main =titre,zlim=c(-65,25),
               ylim=range (bords$y),xlim=range(bords$x))              
lines(bords$x,bords$y)
points(donnees$x,donnees$y,pch=20)

par(mfrow = c(1,2))
Z = matrix(Kelev$predict,nrow=length(grx),ncol=length(gry),byrow=F)
titre = paste("Krigeage avec un modele",c.m,"avec l'altitude")
image.plot(grx,gry,Z,xlab= 'x (km)',ylab ='y (km)',
               col = couleur,asp=1,main =titre,zlim=c(13,571),
               ylim=range (bords$y),xlim=range(bords$x))              
lines(bords$x,bords$y)
points(donnees$x,donnees$y,pch=20)
Z = matrix(Kaniso$predict,nrow=length(grx),ncol=length(gry),byrow=F)
titre = paste("Krigeage avec un modele",c.m,'anisotrope')
image.plot(grx,gry,Z,xlab= 'x (km)',ylab ='y (km)',zlim=c(13,571),
               col = couleur,asp=1,main =titre,
               ylim=range (bords$y),xlim=range(bords$x))     
points(donnees$x,donnees$y,pch=20)         
lines(bords$x,bords$y)
par(mfrow = c(1,1))

elevfull = interpp(grElev[,1],grElev[,2],elev,donneesfull$x,donneesfull$y)
Kelevfull = krige.conv(geodata,loc=cbind(donneesfull$x,donneesfull$y),
            krige=krige.control(type.krige="ok",obj.model=varioest,trend.d = ~ elevobs$z, trend.l = ~elevfull$z,aniso.pars=paraniso))
RMSEPelev = sqrt(mean((Kelevfull$predict-donneesfull$pluies)^2))

# Simulations conditionnelles
#----------------------------

grx10 = seq(-30,320,5)
gry10 = seq(-20,210,5)
grille10 = expand.grid(grx10,gry10)
s.out = output.control(n.predictive = 100, quant=0.9, thres=300)
KS =  krige.conv(geodata, loc = grille10,
                 krige = krige.control(type.krige="ok",obj.model=varioest),
                 output = s.out)
simus=KS$simulations
simus[simus<0]=0

fichier = paste(FicFig,7,sep="")
imprime(impression,fichier,form)
image.plot(grx10,gry10,matrix(simus[,1],length(grx10),ncol=length(gry10),byrow=F), 
           col=couleur,main="Simulation conditionnelle",zlim=c(0,600),asp=1,
           xlab= 'x (km)',ylab ='y (km)')
lines(bords$x,bords$y)
points(donnees$x,donnees$y,pch=20)
ferme(impression)

fichier = paste(FicFig,8,sep="")
imprime(impression,fichier,form)
image.plot(grx10,gry10,matrix((1 - KS$prob),length(grx10),ncol=length(gry10),byrow=F),
           col=couleur, main="Proba(Y > 300)",zlim=c(0,1),xlab= 'x (km)',
           ylab ='y (km)',asp=1)
lines(bords$x,bords$y)
points(donnees$x,donnees$y,pch=20)
ferme(impression)

grx2 = seq(-30,320,2)
gry2 = seq(-20,210,2)
grille2 = expand.grid(grx2,gry2)
s.out = output.control(n.predictive = 100, thres=300)
KS =  krige.conv(geodata, loc = grille2,
                 krige = krige.control(type.krige="ok",obj.model=varioest),
                 output = s.out)
simus=KS$simulations
simus[simus<0]=0

fichier = paste(FicFig,"7b",sep="")
imprime(impression,fichier,form)
image.plot(grx2,gry2,matrix(simus[,1],length(grx2),ncol=length(gry2),byrow=F), 
           col=couleur,main="Simulation conditionnelle",zlim=c(0,600),asp=1,
           xlab= 'x (km)',ylab ='y (km)')
lines(bords$x,bords$y)
points(donnees$x,donnees$y,pch=20)
ferme(impression)

fichier = paste(FicFig,"8b",sep="")
imprime(impression,fichier,form)
image.plot(grx2,gry2,matrix((1 - KS$prob),length(grx2),ncol=length(gry2),byrow=F),
           col=couleur, main="Proba(Y > 300)",zlim=c(0,1),xlab= 'x (km)',
           ylab ='y (km)',asp=1)
lines(bords$x,bords$y)
points(donnees$x,donnees$y,pch=20)
ferme(impression)

# En prenant l'altitude en compte
#--------------------------------
altitude = elevobs$z
reg = lm(geodata$data~altitude)
georesid = geodata
georesid$data = residuals(reg)
vario.r = variog(georesid,max.dist=m.d,pairs.min=p.m,uvec=interv,messages.screen=FALSE)
plot(vario.r,main = "Variogramme empirique")

varioest.r = variofit(vario.r,cov.model = c.m,minimisation.function = "nls",ini.cov.pars=i.c,fix.nugget=F,fix.kappa=TRUE,max.dist=vario.b$max.dist)
titre = paste("modele ",c.m,", portee =",round(varioest.r$cov.pars[2]*100)/100,", palier =",round(varioest$cov.pars[1]*100)/100)
plot(vario.r,main=titre)
lines(varioest.r)

K.r = krige.conv(georesid,loc=grille,krige=krige.control(type.krige="ok",obj.model=varioest.r))
nouveau = data.frame(altitude=elevgr$z)
pred = predict(reg,nouveau)
Kelev = K.r$predict+pred
fichier = paste(FicFig,"6b",sep="")
imprime(impression,fichier,form)
Z = matrix(Kelev,nrow=length(grx),ncol=length(gry),byrow=F)
titre = paste("Krigeage avec un modele",c.m,"avec l'altitude")
image.plot(grx,gry,Z,xlab= 'x (km)',ylab ='y (km)',
           col = couleur,asp=1,main =titre,zlim=c(0,600),
           ylim=range (bords$y),xlim=range(bords$x))              
lines(bords$x,bords$y)
points(donnees$x,donnees$y,pch=20)
ferme(impression)

