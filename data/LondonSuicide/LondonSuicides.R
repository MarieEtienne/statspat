rm(list=ls())
graphics.off()

library(spdep)
library(spatialreg)

# Lecture des donnees
#--------------------
load("LondonSuicides.RData")
london.gen <- readShapePoly("LDNSuicides.shp")
plot(london.gen)
n = length(y)
Y = y/E
brks = seq(0.5,1.75,0.25) # brks = seq(50,250,20)
paletteBisque = colorRampPalette(c("blue", "white","red"))
couleur = paletteBisque(5)[findInterval(Y,brks)]
titre = "London Suicides"
plot(london.gen,col=couleur,main=titre)
legend('bottomleft', legend=leglabs(brks,'<','>'),fill=paletteBisque(5), 
       bty="n",cex=0.8)

# Etude de la correlation sur les prix au m2
#-------------------------------------------
coords <- coordinates(london.gen)

# systemes de voisinage
#---------------------
voisinage1 <- tri2nb(coords) # par triangulation
titre = 'Voisinages par triangulation'
plot(london.gen,main=titre)
plot(voisinage1,coords,add=TRUE,col='red')

ppv <- knearneigh(coords, k=4)
voisinage2 = knn2nb(ppv)
titre = 'Voisinages par 4-plus proches voisins'
plot(london.gen,main=titre)
plot(voisinage2,coords,add=TRUE,col='red')

voisinage3 <- dnearneigh(coords, d1=0, d2=10000 )
titre = 'Voisinages par plus petites distances (<10000)'
plot(london.gen,main=titre)
plot(voisinage3,coords,add=TRUE,col='red')

voisinage4 <- poly2nb(london.gen)
titre = 'Voisinages par frontieres'
plot(london.gen,main=titre)
plot(voisinage4,coords,add=TRUE,col='red')

voisinage = voisinage4
poids.vois <- nb2listw(voisinage,style="W")
mat.vois <- nb2mat(voisinage,style ="W")

# Donnees binaires
#-----------------
seuil = 1.3#140
fort= Y > seuil
faible = Y < seuil
couleur = rep("red",n)
couleur[faible] = "white"
titre = paste("London Suicides >",seuil)
plot(london.gen,col=couleur,main=titre)
legend('bottomleft', legend=leglabs(c(50,seuil,250),'<','>'),fill=c("white","red"), 
       bty="n")
Fort = rep("fort",n)
Fort[faible] = "faible"
Fort = as.factor(Fort)
print(joincount.mc(Fort,poids.vois,nsim=100))
print(joincount.test(Fort,poids.vois))

# Test de Moran et Geary
#-----------------------
moran(Y,poids.vois,n,Szero(poids.vois))
geary(Y,poids.vois,n,n-1,Szero(poids.vois))

moran.test(Y,randomisation=FALSE,poids.vois)
geary.test(Y,randomisation=FALSE,poids.vois)
moran.test(Y,poids.vois)
geary.test(Y,poids.vois)
moran.mc(Y,poids.vois,nsim=100)
geary.mc(Y,poids.vois,nsim=100)

print(sp.correlogram(voisinage,y,order=5,method="I",zero.policy=T))
plot(sp.correlogram(voisinage,y,order=5,method="I",zero.policy=T),main='London suicides correlogramme spatial')

#====================
# Regression spatiale
#====================

brks = seq(-1.50,2.50,0.5)
paletteBisque = colorRampPalette(c("blue", "white","ed"))
couleur = paletteBisque(8)[findInterval(x1,brks)]
titre = "London detresse sociale"
plot(london.gen,col=couleur,main=titre)
legend('bottomleft', legend=leglabs(brks,'<','>'),fill=paletteBisque(8), 
       bty="n",cex=0.8)

brks = seq(-2,2.50,0.5)
paletteBisque = colorRampPalette(c("blue", "white","red"))
couleur = paletteBisque(9)[findInterval(x2,brks)]
titre = "London fragmentation sociale"
plot(london.gen,col=couleur,main=titre)
legend('bottomleft', legend=leglabs(brks,'<','>'),fill=paletteBisque(9), 
       bty="n",cex=0.8)

london.df = data.frame(suicides = y, Esuicides = log(Y), detresse = x1, fragmentation = x2)
     
# 1. Modele lineaire standard
#----------------------------
london.lm = lm(Esuicides ~ x1+x2, london.df) 
summary(london.lm)
anova(london.lm)
lm.morantest(london.lm,poids.vois)

paletteRes = colorRampPalette(c('darkblue',"white", "darkred"))
brks= seq(-0.25,0.25,0.1)

res0 <- residuals(london.lm)
I0 = moran(res0,poids.vois,n,Szero(poids.vois))
moran.test(res0,poids.vois) 
couleur = paletteRes(5)[findInterval(res0,brks)]
titre = paste("Residus modele lineaire, I = ",round(I0$I,2))
plot(london.gen,col=couleur,main=titre)
legend('bottomleft', legend=leglabs(brks,'<','>'),fill=paletteRes(5), 
       bty="n",cex=0.8)

# 2. Modele lineaire avec corr?lation spatiale
#---------------------------------------------
prixm2.gls = gls(prix ~ popu, Pop,
                 correlation = corExp(form = ~ coords[,1]+coords[,2],nugget=TRUE))
summary(prixm2.gls)
anova(prixm2.gls)

res1 <- residuals(prixm2.gls)
I1 = moran(res1,poids.vois,length(prixm2$prix.m2),Szero(poids.vois))
moran.test(res1,poids.vois)

couleur = coulRes[findInterval(res1,brksRes)]
titre = paste("Residus modele lineaire g?n?ralis?, I = ",round(I1$I,2))
traceQuartiers(titre,couleur)
legend('topleft', legend=leglabs(brksRes,'<','>',':'),fill=coulRes, 
       bty="n",cex=0.8)

# 3. Modele SAR sur les erreurs
#------------------------------
london.errorsarlm <- errorsarlm(Esuicides ~ x1+x2,london.df,poids.vois) 
summary(london.errorsarlm)

res2 = residuals(london.errorsarlm)
moran.test(res2,poids.vois)
I2 = moran(res2,poids.vois,n,Szero(poids.vois))
titre = paste("Residus modele SAR sur les erreurs, I = ",round(I2$I,2))

couleur = coulRes[findInterval(res2,brksRes)]
traceQuartiers(titre,couleur)
legend('topleft', legend=leglabs(brksRes,'<','>',':'),fill=coulRes, 
       bty="n",cex=0.8)


# 4. Modele SAR a retard
#-----------------------
london.lagsarlm <- lagsarlm(Esuicides ~ x1+x2,london.df,poids.vois)
summary(london.lagsarlm)

res3 = residuals(prixm2.lagsarlm)
moran.test(res3,poids.vois.sym)
I3 = moran(res3,poids.vois,length(prixm2$prix.m2),Szero(poids.vois))
titre = paste("Residus modele SAR ? retard, I = ",round(I3$I,2))

couleur = coulRes[findInterval(res3,brksRes)]
traceQuartiers(titre,couleur)
legend('topleft', legend=leglabs(brksRes,'<','>',':'),fill=coulRes, 
       bty="n",cex=0.8)

# 5. Modele mixte SAR a retard
#-----------------------------
prixm2.mixte <- lagsarlm(prix ~ popu,Pop,poids.vois.sym, type="mixed",
                                method="Matrix", 
                                control=list(tol.opt = .Machine$double.eps^2), quiet=FALSE) 
summary(prixm2.mixte)

res4 = residuals(prixm2.mixte)
moran.test(res4,poids.vois.sym)
I4 = moran(res4,poids.vois,length(prixm2$prix.m2),Szero(poids.vois))
titre = paste("Residus modele SAR mixte, I = ",round(I4$I,2))

couleur = coulRes[findInterval(res4,brksRes)]
traceQuartiers(titre,couleur)
legend('topleft', legend=leglabs(brksRes,'<','>',':'),fill=coulRes, 
       bty="n",cex=0.8)

# Comparaison des modeles
#------------------------
print(anova(prixm2.mixte,prixm2.errorsarlm,prixm2.lm))
print(anova(prixm2.mixte,prixm2.lagsarlm,prixm2.lm))

# tests du multiplicateur de Lagrange
#------------------------------------
lm.LMtests(prixm2.lm, poids.vois.sym, test=c("LMerr", "LMlag", "RLMerr",
       "RLMlag", "SARMA"))
lm.LMtests(prixm2.lm, poids.vois.sym)
lm.LMtests(residuals(prixm2.lm), poids.vois.sym)

