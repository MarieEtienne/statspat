### R code for Chapter 6
### Last update: 18/08/2014

###################################################
### Set working directory and load packages
###################################################
remove(list=ls())
my.dir <- paste(getwd(),"/",sep="")

require(INLA)
inla.setOption(scale.model.default=FALSE)

require(gstat)
require(geoR)
require(fields)
require(maptools)
require(lattice)
require(spdep)


###################################################
### Code for Section 6.1.2
###################################################
# You neeed a folder called "London Suicides" inside your working directory (my.dir)
# with the data downloaded from
# https://sites.google.com/a/r-inla.org/stbook/datasets

library(maptools)
library(spdep)

load(paste(my.dir,"London Suicides/LondonSuicides.RData",sep=""))
london.gen <- readShapePoly(paste(my.dir,"London Suicides/LDNSuicides.shp",sep=""))

temp <- poly2nb(london.gen)
nb2INLA("LDN.graph", temp)
LDN.adj <- paste(getwd(),"/LDN.graph",sep="")

H <- inla.read.graph(filename="LDN.graph")
image(inla.graph2matrix(H),xlab="",ylab="")

formula <- y ~ 1 + f(ID, model="bym",graph=LDN.adj, scale.model=TRUE,
          hyper=list(prec.unstruct=list(prior="loggamma",param=c(1,0.001)), prec.spatial=list(prior="loggamma",param=c(1,0.001))))

names <- sort(london.gen$NAME) 
data.suicides <- data.frame(NAME=names, y=y, E=E, x1=x1, x2=x2)
Nareas <- length(data.suicides[,1])

# The order of the areas needs to be the same between the data and the spatial polygon object obtained importing the shapefile, so we re-order the data
data.boroughs <- attr(london.gen, "data")
order <- match(data.boroughs$NAME,data.suicides$NAME)
data.suicides <- data.suicides[order,]
data.suicides$ID <- seq(1,Nareas)
# Include in the london.gen dataframe the ID variable  1:32 (because SP_Id is 0:31)
attr(london.gen, "data") <- merge(data.boroughs,data.suicides,by="NAME")

mod.suicides <- inla(formula,family="poisson",
              data=data.suicides,E=E,
              control.compute=list(dic=TRUE))

round(mod.suicides$summary.fixed,3) 
round(head(mod.suicides$summary.random$ID),3) #partial output

exp.b0.mean <- inla.emarginal(exp,mod.suicides$marginals.fixed[[1]])
exp.b0.mean
exp.b0.95CI <- inla.qmarginal(c(0.025,0.975), inla.tmarginal(exp,mod.suicides$marginals.fixed[[1]]))
exp.b0.95CI

# *** Code for Figure 6.6 left
csi <- mod.suicides$marginals.random$ID[1:Nareas]
zeta <- lapply(csi,function(x) inla.emarginal(exp,x))
# Define the cutoff for zeta
zeta.cutoff <- c(0.6, 0.9, 1.0, 1.1, 1.8)
# Transform zeta in categorical variable
cat.zeta <- cut(unlist(zeta),breaks=zeta.cutoff,include.lowest=TRUE)
# Create a dataframe with all the information needed for the map
maps.cat.zeta <- data.frame(ID=data.suicides$ID, cat.zeta=cat.zeta)
# Add the categorized zeta to the spatial polygon
data.boroughs <- attr(london.gen, "data")
attr(london.gen, "data") <- merge(data.boroughs, maps.cat.zeta, by="ID")

trellis.par.set(axis.line=list(col=NA))
spplot(obj=london.gen, zcol= "cat.zeta", col.regions=gray(seq(0.9,0.1,length=4)), asp=1)
# ***

# *** Code for Figure 6.6 right
a <- 0
prob.csi <- lapply(csi, function(x) {1 - inla.pmarginal(a, x)})
prob.csi.cutoff <- c(0,0.2,0.8,1)
cat.prob.csi <- cut(unlist(prob.csi),breaks=prob.csi.cutoff, include.lowest=TRUE)
maps.cat.prob.csi <- data.frame(ID=data.suicides$ID, cat.prob.csi=cat.prob.csi)
data.boroughs <- attr(london.gen, "data")
attr(london.gen, "data") <- merge(data.boroughs, maps.cat.prob.csi, by="ID")

spplot(obj=london.gen, zcol= "cat.prob.csi", col.regions=gray(seq(0.9,0.1,length=3)))
# ***

mat.marg <- matrix(NA, nrow=Nareas, ncol=100000)
m <- mod.suicides$marginals.random$ID
for (i in 1:Nareas){
  # Remember that the first Nareas values of the random effects
  # are u+v, while u values are stored in the Nareas+1 to 2*Nareas elements.
  u <- m[[Nareas+i]]
  mat.marg[i,] <- inla.rmarginal(100000, u)
}
var.u <- apply(mat.marg, 2, var)

var.v <- inla.rmarginal(100000,inla.tmarginal(function(x) 1/x,
          mod.suicides$marginals.hyper$"Precision for ID (iid component)"))
perc.var.u <- mean(var.u/(var.u+var.v))
perc.var.u

marg.hyper <- inla.hyperpar.sample(100000,mod.suicides)
perc.var.u1 <- mean(marg.hyper[,1] / (marg.hyper[,1]+marg.hyper[,2]))
perc.var.u1

###################################################
### Code for Section 6.2 (run the code for Section 6.1.2 first) 
###################################################
formula.eco.reg <- y ~ 1 + x1 + x2 + f(ID,model="bym", graph=LDN.adj)

mod.eco.reg <- inla(formula.eco.reg,family="poisson",
                    data=data.suicides,E=E,
                    control.compute=list(dic=TRUE))

mod.eco.reg$dic$dic

# *** Code for Figure 6.7 left
csi.eco.reg <- mod.eco.reg$marginals.random$ID[1:Nareas]
zeta.eco.reg <- lapply(csi.eco.reg,function(x) inla.emarginal(exp,x))
# Define the cutoff for zeta
zeta.cutoff <- c(0.6, 0.9, 1.0, 1.1, 1.8)
# Transform zeta in categorical variable
cat.zeta.eco.reg <- cut(unlist(zeta.eco.reg),breaks=zeta.cutoff,include.lowest=TRUE)
# Create a dataframe with all the information needed for the map
maps.cat.zeta.eco.reg <- data.frame(ID=data.suicides$ID, cat.zeta.eco.reg=cat.zeta.eco.reg)
# Add the categorized zeta to the spatial polygon dataframe object london.gen
data.boroughs <- attr(london.gen, "data")
attr(london.gen, "data") <- merge(data.boroughs, maps.cat.zeta.eco.reg, by="ID")

spplot(obj=london.gen, zcol= "cat.zeta.eco.reg", col.regions=gray(seq(0.9,0.1,length=4)))
# ***

# *** Code for Figure 6.7 right
a <- 0
prob.csi.eco.reg <- lapply(csi.eco.reg, function(x) {1 - inla.pmarginal(a, x)})
prob.csi.cutoff <- c(0,0.2,0.8,1)
cat.prob.csi.eco.reg <- cut(unlist(prob.csi.eco.reg),breaks=prob.csi.cutoff,include.lowest=TRUE)
maps.cat.prob.csi.eco.reg <- data.frame(ID=data.suicides$ID, cat.prob.csi.eco.reg=cat.prob.csi.eco.reg)
data.boroughs <- attr(london.gen, "data")
attr(london.gen, "data") <- merge(data.boroughs, maps.cat.prob.csi.eco.reg, by="ID")

spplot(obj=london.gen, zcol= "cat.prob.csi.eco.reg", col.regions=gray(seq(0.9,0.1,length=3)))
# ***

marg.hyper <- inla.hyperpar.sample(100000,mod.eco.reg)
perc.var.u1.ecoreg <- mean(marg.hyper[,1]/(marg.hyper[,1]+marg.hyper[,2]))
perc.var.u1.ecoreg

###################################################
### Code for Section 6.3.1
###################################################
# You neeed a folder called "Brain Cancer Navarra" inside your working directory (my.dir)
# with the data downloaded from
# https://sites.google.com/a/r-inla.org/stbook/datasets

load(paste(my.dir,"Brain Cancer Navarra/Navarre.RData",sep=""))
names(brainnav)
navarra.graph <- inla.read.graph(paste(my.dir,"Brain Cancer Navarra/Navarra.graph",sep=""))

data.navarra <- data.frame(ZBS=brainnav$ZBS,Y=brainnav$OBSERVED,E=brainnav$EXPECTED)

formula.zip <- Y ~ 1 + f(ZBS,  model="bym", graph=navarra.graph, 
                hyper=list(prec.unstruct=list(prior="gaussian",param=c(0,1)),
                           prec.spatial=list(prior="gaussian",param=c(0,1))))

mod.zip1 <- inla(formula.zip,family="zeroinflatedpoisson1",
                 data=data.navarra, offset = log(E),
                 control.predictor=list(compute=TRUE))

round(mod.zip1$summary.hyperpar,3)

mod.zip0 <- inla(formula.zip,family="zeroinflatedpoisson0",
                 data=data.navarra, E = E,
                 control.predictor=list(compute=TRUE))

round(mod.zip0$summary.hyperpar,3)

# *** Code for Figure 6.8
Nareas <- nrow(data.navarra)
# Extract the random effects
zeta.navarra0 <- data.frame(zeta=unlist(lapply(mod.zip0$marginals.random$ZBS[1:Nareas],function(x)inla.emarginal(exp,x))))
zeta.navarra1 <- data.frame(zeta=unlist(lapply(mod.zip1$marginals.random$ZBS[1:Nareas],function(x)inla.emarginal(exp,x))))
# Create factor variables
RR.cutoff<- c(0.5, 0.8, 0.95,  1.05,  1.2, 2.5)
RR.navarra0 <- cut(zeta.navarra0$zeta,breaks=RR.cutoff,include.lowest=TRUE)
RR.navarra1 <- cut(zeta.navarra1$zeta,breaks=RR.cutoff,include.lowest=TRUE)

results <- data.frame(ZBS=data.navarra$ZBS,RR.navarra0, RR.navarra1)
data.navarra.shp <- attr(brainnav, "data")
attr(brainnav, "data") <- merge(data.navarra.shp, results, by="ZBS")

trellis.par.set(axis.line=list(col=NA))
spplot(obj=brainnav, zcol="RR.navarra1", col.regions=gray(4.5:0.5/5),main="")
# ***

###################################################
### Code for Section 6.3.2
###################################################
# You neeed a folder called "SDO Piemonte" inside your working directory (my.dir)
# with the data downloaded from
# https://sites.google.com/a/r-inla.org/stbook/datasets

data <- read.csv(paste(my.dir,"SDO Piemonte/dataResp.csv",sep=""))

# Import the shapefile
torino <- readShapePoly(paste(my.dir,"SDO Piemonte/Shapefile/Torino1.shp",sep=""))
# Adjacency matrix
torino.adj <- paste(my.dir,"SDO Piemonte/torino.graph",sep="")

formula.inla <- RESP ~ PM10 + f(ID.MUNICIPALITY, model="besag",graph=torino.adj) 

mod.zib <- inla(formula.inla,family= "zeroinflated.binomial.0", 
                Ntrials=TotalPop,data=data,
                control.compute=list(dic=TRUE))
mod.zib$dic$dic
round(mod.zib$summary.hyperpar,3)
mod.zib$summary.fixed[2,]*10

# *** Code for Figure 6.9
zeta.piemonte <- data.frame(zeta1=unlist(lapply(mod.zib$marginals.random$ID.MUNICIPALITY[1:315], function(x)inla.emarginal(exp,x))),
                            ISTAT_CODE=data$MUNICIPALITY)
# Create factor variables
RR.cutoff.zeta <- c(0.3, 1,  1.2,  5.5)
RR.piemonte <- cut(zeta.piemonte$zeta1,breaks=RR.cutoff.zeta,include.lowest=TRUE)
results <- data.frame(ISTAT_CODE=data$MUNICIPALITY,RR.piemonte) 

data.communes <- attr(torino, "data")
attr(torino, "data") <- merge(data.communes, results, by="ISTAT_CODE")
spplot(obj=torino, zcol="RR.piemonte", col.regions=gray(2.5:0.5/3),main="")
# ***

# Include a random slope in the model
formula.inla1 <- RESP ~ f(ID.PM10,  PM10, model = "iid", graph=torino.adj, constr=TRUE,
                        hyper=list(prec = list(prior="gaussian",param=c(0,1)))) +                                                                                                  
                  f(ID.MUNICIPALITY, model="besag",graph=torino.adj) 

mod.zib1 <- inla(formula.inla1, family= "zeroinflated.binomial.0", 
                 Ntrials=TotalPop, data=data,
                 control.compute=list(dic=TRUE))

mod.zib1$summary.hyperpar
mod.zib1$dic$dic

###################################################
### Code for Section 6.7
###################################################
data(SPDEtoy)
dim(SPDEtoy)
head(SPDEtoy, n=3)

# *** Code for Figure 6.11
cutoff <- quantile(SPDEtoy$y, 0:5/5)
y.factor <- cut(SPDEtoy$y,breaks=cutoff,include.lowest=TRUE,right=TRUE)

point.cex <- seq(0.5,2,length=5)
y.cex <- point.cex[as.numeric(y.factor)]

plot(SPDEtoy[,1:2], xlim=c(0,1.4), cex=y.cex, xlab="", ylab="",xaxt="n", pch=21, bg=grey(0.9))
axis(side=1, at=seq(0,1,l=6),label=seq(0,1,l=6))
legend("topright",  pch=21, pt.bg=rep(grey(0.9),5), pt.cex=point.cex, legend=levels(y.factor), bty="n")
# ***

###################################################
### Code for Section 6.7.1
###################################################
coords <- as.matrix(SPDEtoy[,1:2])
mesh0 <- inla.mesh.2d(loc=coords, max.edge=0.1)
mesh1 <- inla.mesh.2d(loc=coords, max.edge=c(0.1, 0.1))
mesh2 <- inla.mesh.2d(loc=coords, max.edge=c(0.1, 0.2))

# *** Code for Figure 6.12
# I need mesh4 for setting xlim for all meshes
mesh4 <- inla.mesh.2d(loc=coords, max.edge=c(0.1, 0.2),offset=c(0.2,0.4))

plot(mesh0, main = "", asp=1, xlim=summary(mesh4)$xlim)
points(coords, pch=21, bg=1, col="white", cex=1.8)

plot(mesh1, main="", asp=1, xlim=summary(mesh4)$xlim)
points(coords, pch=21, bg=1, col="white", cex=1.8)

plot(mesh2, main="", asp=1, xlim=summary(mesh4)$xlim)
points(coords, pch=21, bg=1, col="white", cex=1.8)
# ***

mesh3 <- inla.mesh.2d(loc=coords, max.edge=c(0.1, 0.2), offset=c(0.4,0.1))
mesh4 <- inla.mesh.2d(loc=coords, max.edge=c(0.1, 0.2), offset=c(0.1,0.4))

# *** Code for Figure 6.13
plot(mesh3, main="",xlim=summary(mesh4)$xlim,asp=1)
points(coords, pch=21, bg=1, col="white", cex=1.8)

plot(mesh4, main="",xlim=summary(mesh4)$xlim,asp=1)
points(coords, pch=21, bg=1, col="white", cex=1.8)
# ***

domain <- matrix(cbind(c(0,1,1,0.7,0), c(0,0,0.7,1,1)),ncol=2)
mesh5 <- inla.mesh.2d(loc.domain=domain, max.edge=c(0.04, 0.2), cutoff=0.015, offset = c(0.1, 0.4))
mesh6 <- inla.mesh.2d(loc.domain=domain, max.edge=c(0.04, 0.2), cutoff=0.05, offset = c(0.1, 0.4))

# *** Code for Figure 6.14
plot(mesh5, main="",xlim=summary(mesh4)$xlim,asp=1)
points(coords, pch=21, bg=1, col="white", cex=1.8)
points(domain, pch=23, bg="grey", col=1, cex=2, lwd=2.5)

plot(mesh6, main="",xlim=summary(mesh4)$xlim,asp=1)
points(coords, pch=21, bg=1, col="white", cex=1.8)
points(domain, pch=23, bg="grey", col=1, cex=2, lwd=2.5)
# ***

# Obtain the number of vertices
inla.spde2.matern(mesh0,alpha=2)$n.spde

vertices <- c()
for(i in 0:6){
  vertices[i+1] <- inla.spde2.matern(get(paste("mesh",i,sep="")),alpha=2)$n.spde
}
vertices

bnd <- inla.nonconvex.hull(as.matrix(coords),convex=0.07)
mesh7 <- inla.mesh.2d(loc=coords, boundary=bnd, max.edge=c(0.04, 0.2), cutoff=0.05, offset = c(0.1, 0.4))

# *** Code for Figure 6.15
plot(mesh7, main="",xlim=summary(mesh4)$xlim,asp=1)
points(coords, pch=21, bg=1, col="white", cex=1.8)
# ***

###################################################
### Code for Section 6.7.2
###################################################
A.est1 <- inla.spde.make.A(mesh=mesh1, loc=coords)
dim(A.est1)
A.est6 <- inla.spde.make.A(mesh=mesh6, loc=coords)
dim(A.est6)

table(apply(A.est6,1,nnzero))
table(apply(A.est6,1,sum))
table(apply(A.est6,2,sum) > 0) 
table(apply(A.est1,2,sum))

###################################################
### Code for Section 6.7.3
###################################################
spde <- inla.spde2.matern(mesh=mesh6, alpha=2)
formula <- y ~ -1 + intercept + f(spatial.field, model=spde)

output6 <- inla(formula,
                data = list(y=SPDEtoy$y, intercept=rep(1,spde$n.spde),spatial.field=1:spde$n.spde),
                control.predictor=list(A=A.est6,compute=TRUE))

round(output6$summary.fixed,3)
round(output6$summary.hyperpar[1,],3)

inla.emarginal(function(x) 1/x, output6$marginals.hyper[[1]])

output6.field <- inla.spde2.result(inla=output6, name="spatial.field", spde=spde, do.transf=TRUE)

inla.emarginal(function(x) x, output6.field$marginals.kappa[[1]])
inla.emarginal(function(x) x, output6.field$marginals.variance.nominal[[1]])
inla.emarginal(function(x) x, output6.field$marginals.range.nominal[[1]])

inla.hpdmarginal(0.95, output6.field$marginals.kappa[[1]])
inla.hpdmarginal(0.95, output6.field$marginals.variance.nominal[[1]])
inla.hpdmarginal(0.95, output6.field$marginals.range.nominal[[1]])

###################################################
### Code for Section 6.8
###################################################
s.index <- inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde)
names(s.index)

stack.est <- inla.stack(data=list(y=SPDEtoy$y),
    A=list(A.est6),
    effects=list(c(s.index, list(intercept=1))), 
    tag="est") #Estimation

output6.stack <- inla(formula,
   data=inla.stack.data(stack.est, spde=spde),
   family="gaussian",
   control.predictor=list(A=inla.stack.A(stack.est), 
   compute=TRUE))             


# *** Code for Figure 6.16
lres <- lrf <- list()
n.meshes <- 8
    
for (k in 1:n.meshes) {
      cat(k)
      A <- inla.spde.make.A(get(paste('mesh', k-1, sep='')), loc=coords)
      spde <- inla.spde2.matern(get(paste('mesh', k-1, sep='')), alpha=2)
      
      formula <- y ~ -1 + intercept + f(spatial.field, model=spde)
      s.index <- inla.spde.make.index(name="spatial.field", 
                                      n.spde=spde$n.spde)
      stack.est <- inla.stack(data=list(y=SPDEtoy$y),
                              A=list(A),
                              effects=list(c(s.index, list(intercept=1))),
                              tag="est")
      
      lres[[k]] <- inla(formula,
                        data=inla.stack.data(stack.est, spde=spde),
                        control.predictor=list(A=inla.stack.A(stack.est)))
      
      lrf[[k]] <- inla.spde2.result(lres[[k]], "spatial.field",
                                    spde, do.transf=TRUE)
}
    
# True parameter values
beta0 <- 10; sigma2e <- 0.3; sigma2x <- 5; kappa <- 7; nu <- 1
# Plot setting
lty.vec=c(1,2,3,2,3,2,3)
lwd.vec=c(2,1,1,2,2,3,3)
    
# Plot for beta0
xrange <- range(sapply(lres, function(x) range(x$marginals.fix[[1]][,1])))
yrange <- range(sapply(lres, function(x) range(x$marginals.fix[[1]][,2])))

plot(lres[[1]]$marginals.fix[[1]], type='l', xlim=xrange, ylim=yrange, xlab=expression(b[0]),ylab="")
for (k in 2:n.meshes)
    lines(lres[[k]]$marginals.fix[[1]], lty=lty.vec[[k-1]],lwd=lwd.vec[[k-1]])
abline(v=beta0, lty=2, lwd=2)
legend('topright', c(paste('mesh', 0:7, sep='')),lty=c(1,lty.vec), lwd=c(1,lwd.vec),  bty='n')
    
# Plot for the variance sigma2e
s2e.marg <- lapply(lres, function(m) inla.tmarginal(function(x) 1/x, m$marginals.hy[[1]]))
xrange <- range(sapply(s2e.marg, function(x) range(x[,1])))
yrange <- range(sapply(s2e.marg, function(x) range(x[,2])))

plot.default(s2e.marg[[1]], type='l', xlim=xrange, ylim=yrange, xlab=expression(sigma[e]^2),ylab="")
for (k in 2:n.meshes)
      lines(s2e.marg[[k]], lty=lty.vec[[k-1]],lwd=lwd.vec[[k-1]])
abline(v=sigma2e, lty=2, lwd=2)
legend('topright', c(paste('mesh', 0:7, sep='')),lty=c(1,lty.vec), lwd=c(1,lwd.vec),  bty='n')
  
# Plot for the variance sigma2
xrange <- range(sapply(lrf, function(r) range(r$marginals.variance.nominal[[1]][,1])))
yrange <- range(sapply(lrf, function(r) range(r$marginals.variance.nominal[[1]][,2])))

plot(lrf[[1]]$marginals.variance.nominal[[1]], type='l', xlim=xrange, ylim=yrange, xlab=expression(sigma^2),ylab="")
for (k in 2:n.meshes)
  lines(lrf[[k]]$marginals.variance.nominal[[1]], lty=lty.vec[[k-1]],lwd=lwd.vec[[k-1]])
abline(v=sigma2x, lty=2, lwd=2)
legend('topright', c(paste('mesh', 0:7, sep='')),lty=c(1,lty.vec), lwd=c(1,lwd.vec),  bty='n')
    
#Plot for the range
xrange <- range(sapply(lrf, function(r) range(r$marginals.range.nominal[[1]][,1])))
yrange <- range(sapply(lrf, function(r) range(r$marginals.range.nominal[[1]][,2])))

plot(lrf[[1]]$marginals.range.nominal[[1]], type='l', xlim=xrange, ylim=yrange, xlab=expression(r),ylab="")
for (k in 2:n.meshes)
      lines(lrf[[k]]$marginals.range.nominal[[1]], lty=lty.vec[[k-1]],lwd=lwd.vec[[k-1]])
abline(v=sqrt(8)/kappa, lty=2, lwd=2)
legend('topright', c(paste('mesh', 0:7, sep='')),lty=c(1,lty.vec), lwd=c(1,lwd.vec),  bty='n')
# ***

###################################################
### Code for Section 6.8.1
###################################################
grid.x <- 50
grid.y <- 50
pred.grid <- expand.grid(x = seq(0, 1, length.out = grid.x), y = seq(0, 1, length.out = grid.y))
dim(pred.grid)

A.pred6 <- inla.spde.make.A(mesh=mesh6, loc=as.matrix(pred.grid))
dim(A.pred6)

spde <- inla.spde2.matern(mesh=mesh6, alpha=2)
s.index <- inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde)

formula <- y ~ -1 + intercept + f(spatial.field, model=spde)

stack.est <- inla.stack(data=list(y=SPDEtoy$y),
         A=list(A.est6),
         effects=list(c(s.index, list(intercept=1))),
         tag="est")

stack.pred.latent <- inla.stack(data=list(xi=NA), 
         A=list(A.pred6),
         effects=list(s.index),
         tag="pred.latent") 

stack.pred.response <- inla.stack(data=list(y=NA), 
         A=list(A.pred6),
         effects=list(c(s.index, list(intercept=1))),
         tag="pred.response") 

join.stack <- inla.stack(stack.est, stack.pred.latent, stack.pred.response)

join.output <- inla(formula,
         data=inla.stack.data(join.stack),
         control.predictor=list(A=inla.stack.A(join.stack), compute=TRUE))

index.pred.latent <- inla.stack.index(join.stack, tag="pred.latent")$data
index.pred.response <- inla.stack.index(join.stack, tag="pred.response")$data

round(head(join.output$summary.linear.predictor[index.pred.latent,1:5],n=3),3)
round(head(join.output$summary.fitted.values[index.pred.latent,1:5],n=3),3)

post.mean.pred.latent <- join.output$summary.linear.predictor[index.pred.latent,"mean"]
post.sd.pred.latent <- join.output$summary.linear.predictor[index.pred.latent,"sd"]
post.mean.pred.response <- join.output$summary.fitted.values[index.pred.response,"mean"]
post.sd.pred.response <- join.output$summary.fitted.values[index.pred.response,"sd"]

# *** Code for Figure 6.17
levelplot(matrix(post.mean.pred.latent,grid.x,grid.y),
          col.regions=gray(seq(0.85,0,length=50)),
          xlab="", ylab="", scales=list(draw=FALSE))

levelplot(matrix(post.sd.pred.latent,grid.x,grid.y),
          col.regions=gray(seq(0.85,0,length=50)),
          xlab="", ylab="", scales=list(draw=FALSE),range=c(0,2),at=seq(0,2,by=0.2))

levelplot(matrix(post.mean.pred.response,grid.x,grid.y),
          col.regions=gray(50:0/50),
          xlab="", ylab="", scales=list(draw=FALSE))

levelplot(matrix(post.sd.pred.response,grid.x,grid.y),
          col.regions=gray(50:0/50),
          xlab="", ylab="", scales=list(draw=FALSE),range=c(0,2),at=seq(0,2,by=0.2))
# ***

###################################################
### Code for Section 6.9.1
###################################################
limits <- cbind(c(0,1,1,0,0), c(0,0,1,1,0))
n.loc <- 100

set.seed(1653)
locations <- cbind(s1=sample(1:n.loc / n.loc), s2=sample(1:n.loc / n.loc))

mesh <- inla.mesh.2d(loc=locations,loc.domain= limits, cutoff=0.03, max.edge=c(0.07,.12)) 

# *** Code for Figure 6.18
plot(mesh,main="",asp=1)
points(locations ,cex=1.2, pch=19)
# ***

range0 <- 0.3
sigma0 <- 1
kappa0 <- sqrt(8)/range0
tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)

spde_stat <- inla.spde2.matern(mesh,
                  B.tau=matrix(c(log(tau0), -1, +1),nrow=1,ncol=3),
                  B.kappa=matrix(c(log(kappa0), 0, -1),nrow=1,ncol=3), 
                  theta.prior.mean=c(0, 0),
                  theta.prior.prec=c(0.1, 0.1))

Q_stat <- inla.spde2.precision(spde=spde_stat, theta=c(0,0))
sample_stat <- as.vector(inla.qsample(n=1, Q=Q_stat, seed=1434))
length(sample_stat)

spde_stat_v2 <- inla.spde2.matern(mesh,
                            B.tau=matrix(c(0, 1, 0),nrow=1,ncol=3),
                            B.kappa=matrix(c(0, 0, 1),nrow=1,ncol=3),
                            theta.prior.mean=c(0, 0),
                            theta.prior.prec=c(0.1, 0.1))
Q_stat_v2 <- inla.spde2.precision(spde_stat_v2, 
                                  theta=c(log(tau0),log(kappa0)))
sample_stat_v2 <- as.vector(inla.qsample(n=1, Q=Q_stat_v2, seed=1434))

# Check that the simulated values of the GMRF are the same:
sum(sample_stat - sample_stat_v2) 

# Simulate the covariate values
set.seed(344)
covariate <- rnorm(n.loc,mean=0,sd=1)
# Compute the observation matrix
A <- inla.spde.make.A(mesh, loc=locations)
dim(A)
# Simulate the observations y
set.seed(545)
y <- 10 + 3*covariate + as.vector(A %*%  sample_stat) + rnorm(n.loc,mean=0,sd=sqrt(0.25))

mesh.index <- inla.spde.make.index(name="field", n.spde=spde_stat$n.spde)
    
stack.est <- inla.stack(data=list(y=y),
                       A=list(A,1),
                       effects=list(c(mesh.index,list(intercept=1)), list(x=covariate)),
                       tag="est")

formula <- y ~ -1 + intercept + x + f(field,model=spde_stat)

output_stat <- inla(formula,
                    data=inla.stack.data(stack.est,spde=spde_stat),
                    family="normal",
                    control.predictor=list(A=inla.stack.A(stack.est),compute=TRUE))

spde.result <- inla.spde2.result(inla=output_stat,name="field",spde=spde_stat)

# *** Code for Figure 6.19
plot(spde.result$marginals.range.nominal[[1]],t="l",xlab="r",ylab="")
abline(v=range0)

plot(spde.result$marginals.variance.nominal[[1]],t="l",xlab=expression(sigma^2),ylab="")
abline(v=sigma0)

plot(inla.smarginal(output_stat$marginals.fixed[[1]]),t="l",xlab=expression(b[0]),ylab="")
abline(v=10)

plot(inla.smarginal(output_stat$marginals.fixed[[2]]),t="l",xlab=expression(beta[1]),ylab="")
abline(v=3)
# ***

###################################################
### Code for Section 6.10
###################################################
remove(list=ls())

library(geoR)
data(SIC)

# Data preparation  
sic.all$coords[,1] <- (sic.all$coords[,1] - apply(sic.borders,2,mean)[1])
sic.all$coords[,2] <- (sic.all$coords[,2] - apply(sic.borders,2,mean)[2])
sic.borders <- apply(sic.borders, 2, scale, scale=F)

# Find in sic.all the 100 locations given in sic.100 and used for estimation  
index.est <- which(as.numeric(rownames(sic.all$coords)) %in% as.numeric(rownames(sic.100$coords)))

# Prepare data for estimation (100 locations)
est.coord <- sic.all$coords[index.est,]
est.data <-  sqrt(sic.all$data[index.est])
est.elevation <- sic.all$altitude[index.est]/1000

# Prepare data for validation (367 locations)
val.coord <- sic.all$coords[-index.est,]
val.data <- sqrt(sic.all$data[-index.est])
val.elevation <- sic.all$altitude[-index.est]/1000

# *** Code for Figure 6.20
color.map <- gray(c(1, .8, .5, 0))
sic.cutoff <- quantile(sqrt(sic.all$data))

sic.factor.est <- cut(est.data,breaks=sic.cutoff,include.lowest=TRUE,right=TRUE)
q4col.est <- color.map[as.numeric(sic.factor.est)]
sic.factor.val <- cut(val.data,breaks=sic.cutoff,include.lowest=TRUE,right=TRUE)
q4col.val <- color.map[as.numeric(sic.factor.val)]

plot(sic.borders, type="l",axes=F, xlab="",ylab="",asp=1)
points(est.coord[,"V2"], est.coord[,"V3"], bg=q4col.est, pch=21, cex=1.2, lwd=2)
points(val.coord[,"V2"], val.coord[,"V3"], bg=q4col.val, pch=24, cex=1.2)
legend(x=-200, y=106, pch=rep(21,4), pt.cex=rep(1.2,4),
       pt.bg=rev(color.map), legend=rev(levels(sic.factor.val)), bty="n")
legend(x=-205, y=106, pch=rep(24,4), pt.cex=rep(1.2,4),
       pt.bg=rev(color.map), legend=c("","","",""), bty="n")
# ***

library(gstat)
data(sic97)

x.res <- demstd@grid@cells.dim[1]  #376
y.res <- demstd@grid@cells.dim[2]  #253
pred.elevation <- as.matrix(demstd@data)
elevation.grid <- matrix(pred.elevation, nrow=x.res,ncol=y.res,byrow=F)
# Set the right orientation for the grid
elevation.grid <- elevation.grid[,y.res:1]/1000 
dim(elevation.grid)

seq.x.grid <- seq(from=demstd@coords[1,1],to=demstd@coords[2,1],length=x.res)/1000 
seq.y.grid <- seq(from=demstd@coords[1,2],to=demstd@coords[2,2],length=y.res)/1000 
pred.grid <- as.matrix(expand.grid(x=seq.x.grid,y=seq.y.grid))

# *** Code for Figure 6.21
library(fields)
image.plot(seq.x.grid,seq.y.grid,elevation.grid, col=grey(50:0/50),cex=2,
           xlab="W-E (km)", ylab="N-S (km)",axis.args = list(cex.axis = 2),legend.mar=7)
lines(sic.borders,lwd=2)
# ***

# Model fitting and spatial prediction at the validation stations 
Swiss.mesh <- inla.mesh.2d(loc.domain=sic.borders, max.edge=c(35,100))

# *** Code for Figure 6.22
plot(Swiss.mesh,main="",asp=1)
points(est.coord,pch=19, col=1, cex=1)
lines(sic.borders,lwd=4,col=1)
# ***

Swiss.spde <- inla.spde2.matern(mesh=Swiss.mesh,alpha=2)

A.est <- inla.spde.make.A(mesh=Swiss.mesh, loc=est.coord)
A.val <- inla.spde.make.A(mesh=Swiss.mesh, loc=val.coord)

s.index <- inla.spde.make.index(name="spatial.field", n.spde=Swiss.spde$n.spde) 

stack.est <- inla.stack(data = list(rain=est.data),
              A = list(A.est, 1),
              effects = list(c(s.index, list(Intercept=1)), list(Elevation=est.elevation)),
              tag="est")

stack.val <- inla.stack(data = list(rain=NA),
              A = list(A.val,1),
              effects = list(c(s.index, list(Intercept=1)), list(Elevation=val.elevation)),
              tag="val")

join.stack <- inla.stack(stack.est, stack.val)

formula <- rain ~ -1 + Intercept + Elevation + f(spatial.field, model=spde)

Swiss.output <-  inla(formula, 
                  data=inla.stack.data(join.stack, spde=Swiss.spde),family="gaussian",
                  control.predictor=list(A=inla.stack.A(join.stack), compute=TRUE),
                  control.compute=list(cpo=TRUE, dic=TRUE))

# Extract the parameter estimates for the model WITH Elevation
fixed.out <- round(Swiss.output$summary.fixed[,1:5],3)

sigma2e_marg <- inla.tmarginal(function(x) 1/x, Swiss.output$marginals.hyperpar[[1]])
sigma2e_m1 <- inla.emarginal(function(x) x, sigma2e_marg)
sigma2e_m2 <- inla.emarginal(function(x) x^2, sigma2e_marg)
sigma2e_stdev <- sqrt(sigma2e_m2 - sigma2e_m1^2)
sigma2e_quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), sigma2e_marg)

mod.field <- inla.spde2.result(Swiss.output, name="spatial.field", Swiss.spde)

var.nom.marg <- mod.field$marginals.variance.nominal[[1]]
var.nom.m1 <- inla.emarginal(function(x) x, var.nom.marg)
var.nom.m2 <- inla.emarginal(function(x) x^2, var.nom.marg)
var.nom.stdev <- sqrt(var.nom.m2 - var.nom.m1^2)
var.nom.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)

range.nom.marg <- mod.field$marginals.range.nominal[[1]]
range.nom.m1 <- inla.emarginal(function(x) x, range.nom.marg)
range.nom.m2 <- inla.emarginal(function(x) x^2, range.nom.marg)
range.nom.stdev <- sqrt(range.nom.m2 - range.nom.m1^2)
range.nom.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)

# Implement the model WITHOUT elevation and then extract the parameters
stack.est.noelev <- inla.stack(data = list(rain=est.data),
                    A = list(A.est),
                    effects = list(c(s.index, list(Intercept=1))),
                    tag="est")

stack.val.noelev <- inla.stack(data=list(rain=NA),
                    A=list(A.val),
                    effects = list(c(s.index, list(Intercept=1))),
                    tag="val")

join.stack.noelev <- inla.stack(stack.est.noelev, stack.val.noelev)

formula <- rain ~ -1 + Intercept +  f(spatial.field, model=spde)

Swiss.output.noelev <-  inla(formula, 
              data=inla.stack.data(join.stack.noelev, spde=Swiss.spde), family="gaussian",
              control.predictor=list(A=inla.stack.A(join.stack.noelev), compute=TRUE),
              control.compute=list(cpo=TRUE, dic=TRUE))

fixed.out.noelev <- round(Swiss.output.noelev$summary.fixed[,1:5],3)

sigma2e_marg.noelev <- inla.tmarginal(function(x) 1/x,Swiss.output.noelev$marginals.hyperpar[[1]])
sigma2e_m1.noelev <- inla.emarginal(function(x) x, sigma2e_marg.noelev)
sigma2e_m2.noelev <- inla.emarginal(function(x) x^2, sigma2e_marg.noelev)
sigma2e_stdev.noelev <- sqrt(sigma2e_m2.noelev - sigma2e_m1.noelev^2)
sigma2e_quantiles.noelev <- inla.qmarginal(c(0.025, 0.5, 0.975), sigma2e_marg.noelev)

mod.field.noelev <- inla.spde2.result(Swiss.output.noelev, name="spatial.field", Swiss.spde)

var.nom.marg.noelev <- mod.field.noelev$marginals.variance.nominal[[1]]
var.nom.m1.noelev <- inla.emarginal(function(x) x, var.nom.marg.noelev)
var.nom.m2.noelev <- inla.emarginal(function(x) x^2, var.nom.marg.noelev)
var.nom.stdev.noelev <- sqrt(var.nom.m2.noelev - var.nom.m1.noelev^2)
var.nom.quantiles.noelev <- inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg.noelev)

range.nom.marg.noelev <- mod.field.noelev$marginals.range.nominal[[1]]
range.nom.m1.noelev <- inla.emarginal(function(x) x, range.nom.marg.noelev)
range.nom.m2.noelev <- inla.emarginal(function(x) x^2, range.nom.marg.noelev)
range.nom.stdev.noelev <- sqrt(range.nom.m2.noelev - range.nom.m1.noelev^2)
range.nom.quantiles.noelev <- inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg.noelev)

index.val <- inla.stack.index(join.stack,"val")$data

post.mean.val <- Swiss.output.noelev$summary.linear.predictor[index.val,"mean"]
post.sd.val <- Swiss.output.noelev$summary.linear.predictor[index.val,"sd"]

res <- val.data - post.mean.val
RMSE <- sqrt(mean(res^2)) #RMSE
RMSE
correl <- cor(val.data, post.mean.val)
correl

# Spatial prediction at the grid locations 
A.pred <- inla.spde.make.A(mesh=Swiss.mesh)

stack.pred <- inla.stack(data = list(rain=NA),
              A = list(A.pred),
              effects = list(c(s.index, list(Intercept=1))),
              tag="pred")

join.stack.noelev <- inla.stack(stack.est.noelev, stack.pred)

Swiss.output.noelev <-  inla(formula, 
                data=inla.stack.data(join.stack.noelev, spde=Swiss.spde),
                family="gaussian",
                control.predictor=list(A=inla.stack.A(join.stack.noelev), compute=TRUE))

index.pred <- inla.stack.index(join.stack.noelev,"pred")$data
post.mean.pred <- Swiss.output.noelev$summary.linear.predictor[index.pred, "mean"]
post.sd.pred <- Swiss.output.noelev$summary.linear.predictor[index.pred, "sd"]

proj.grid <- inla.mesh.projector(Swiss.mesh, xlim=range(pred.grid[,1]), ylim=range(pred.grid[,2]), dims=c(x.res,y.res))
post.mean.pred.grid <- inla.mesh.project(proj.grid, post.mean.pred)
post.sd.pred.grid <- inla.mesh.project(proj.grid, post.sd.pred)

# *** Code for Figure 6.23
image.plot(seq.x.grid,seq.y.grid,post.mean.pred.grid,
      xlab="W-E (km)", ylab="N-S (km)", col=grey(20:0/20), axes=F,cex=2,legend.mar=7)
contour(seq.x.grid,seq.y.grid,post.mean.pred.grid, add=T, lwd=2, labcex=1)
lines(sic.borders,lwd=4)

image.plot(seq.x.grid,seq.y.grid,post.sd.pred.grid,
      xlab="W-E (km)", ylab="N-S (km)",col=grey(20:0/20), axes=F,cex=2,legend.mar=7)
contour(seq.x.grid,seq.y.grid,post.sd.pred.grid, add=T, lwd=2, labcex=1)
lines(sic.borders,lwd=3)
# ***

###################################################
### Code for Section 6.11
###################################################
remove(list=ls())

library(geoR)
data(gambia)

coords <-  as.matrix(gambia[,1:2])/1000 #in km

# Create an index at the village level
ind <- paste("x",coords[,1], "y", coords[,2], sep="")
# Detect non duplicated villages
which.nodupl <- which(!duplicated(ind))
village.index <- c(NA, length=nrow(gambia))
village.index[1 : (which.nodupl[length(which.nodupl)]-1)] <-
rep(1:64,times=as.numeric(diff(which.nodupl)))
village.index[which.nodupl[length(which.nodupl)] : nrow(gambia)] <- 65
gambia$village.index <- village.index

bnd <- inla.nonconvex.hull(coords,convex=-0.1)
gambia.mesh <- inla.mesh.2d(boundary = bnd,offset=c(30, 60), max.edge=c(20,40))

# *** Code for Figure 6.24
plot(gambia.mesh,main="",asp=1)
points(coords,pch=21,bg="white",cex=1.5,lwd=1.5)
# ***

gambia.spde <- inla.spde2.matern(mesh=gambia.mesh, alpha=2)
A.est <- inla.spde.make.A(mesh=gambia.mesh, loc=coords)
s.index <- inla.spde.make.index(name="spatial.field",n.spde=gambia.spde$n.spde)

gambia.stack.est <- inla.stack(data=list(y=gambia$pos),
    A=list(A.est, 1, 1, 1, 1, 1, 1),
    effects=
      list(c(s.index, list(Intercept=1)),
      list(age=gambia$age/365),
      list(treated=gambia$treated),
      list(netuse=gambia$netuse),
      list(green=gambia$green),
      list(phc=gambia$phc),
      list(village.index=gambia$village.index)),
    tag="est")

formula <- y ~ -1 + Intercept + treated + netuse + age + green + phc + 
          f(spatial.field, model=gambia.spde) + f(village.index, model="iid")

gambia.output <- inla(formula,
      data=inla.stack.data(gambia.stack.est, spde=gambia.spde),
      family="binomial",Ntrials=1,
      control.predictor=list(A=inla.stack.A(gambia.stack.est), compute=TRUE),
      control.compute=list(dic=TRUE))

fixed.out <- round(gambia.output$summary.fixed[,1:5],3)

# Precision for village.index
sigma2u_marg <- inla.tmarginal(function(x) 1/x, gambia.output$marginals.hyperpar[[3]])
sigma2u_m1 <- inla.emarginal(function(x) x, sigma2u_marg)
sigma2u_m2 <- inla.emarginal(function(x) x^2, sigma2u_marg)
sigma2u_stdev <- sqrt(sigma2u_m2 - sigma2u_m1^2)
sigma2u_quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), sigma2u_marg)

mod.field <- inla.spde2.result(gambia.output, name="spatial.field", gambia.spde)

var.nom.marg <- mod.field$marginals.variance.nominal[[1]]
var.nom.m1 <- inla.emarginal(function(x) x, var.nom.marg)
var.nom.m2 <- inla.emarginal(function(x) x^2, var.nom.marg)
var.nom.stdev <- sqrt(var.nom.m2 - var.nom.m1^2)
var.nom.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), var.nom.marg)

range.nom.marg <- mod.field$marginals.range.nominal[[1]]
range.nom.m1 <- inla.emarginal(function(x) x, range.nom.marg)
range.nom.m2 <- inla.emarginal(function(x) x^2, range.nom.marg)
range.nom.stdev <- sqrt(range.nom.m2 - range.nom.m1^2)
range.nom.quantiles <- inla.qmarginal(c(0.025, 0.5, 0.975), range.nom.marg)

inla.emarginal(function(x) exp(x)/(1+exp(x)),gambia.output$marginals.fixed[["netuse"]])
inla.emarginal(inla.link.invlogit, gambia.output$marginals.fixed[["netuse"]])

# Model WITHOUT the spatial effect
gambia.stack.est.noGF <- inla.stack(data=list(y=gambia$pos),
  A=list(1, 1, 1, 1, 1, 1, 1),
  effects=
    list(list(Intercept=rep(1,nrow(gambia))),
    list(age=gambia$age/365),
    list(treated=gambia$treated),
    list(netuse=gambia$netuse),
    list(green=gambia$green),
    list(phc=gambia$phc),
    list(village.index=gambia$village.index)),
    tag="est")

formula <- y ~ -1 + Intercept + treated + netuse + age + green + phc + f(village.index,model="iid")

gambia.output.noGF <- inla(formula,
      data=inla.stack.data(gambia.stack.est.noGF, spde=gambia.spde),
      family="binomial",Ntrials=1,
      control.predictor=list(A=inla.stack.A(gambia.stack.est.noGF), compute=TRUE),
      control.compute=list(dic=TRUE))

gambia.output.noGF$dic$dic
gambia.output$dic$dic

###################################################
### Code for Section 6.12.1
###################################################
limits <- cbind(c(0,1,1,0,0), c(0,0,1,1,0))
n <- 100
set.seed(1653)
locations <- cbind(s1=sample(1 : n / n), s2=sample(1:n / n))

mesh <- inla.mesh.2d(locations,
                     loc.domain= limits, 
                     cutoff=0.03,
                     max.edge=c(0.07,.12))

# Simulate the covariate values
set.seed(344)
covariate <- rnorm(n)
# Compute the observation matrix
A <- inla.spde.make.A(mesh, loc=locations)

spde_nstat <- inla.spde2.matern(mesh, 
                  B.tau=matrix(cbind(0, 1, 0, sin(pi*mesh$loc[,1])),ncol=4),
                  B.kappa=matrix(c(0, 0, 1, 0),nrow=1,ncol=4), 
                  theta.prior.mean=c(0, 0, 0),
                  theta.prior.prec=c(0.1, 0.1, 0.1))
theta <- c(-1, 2, -1) 

Q_nstat <- inla.spde2.precision(spde=spde_nstat, theta=theta)
sample_nstat <- as.vector(inla.qsample(n=1, Q=Q_nstat, seed=1))

set.seed(5514)
y_nstat <- 10 + 3*covariate + as.vector(A %*% sample_nstat) + rnorm(n,mean=0,sd=sqrt(0.25))

mesh.index_nstat <- inla.spde.make.index(name="field", n.spde=spde_nstat$n.spde)

stack.est_nstat <- inla.stack(data=list(y=y_nstat),
                  A=list(A,1),
                  effects=list(c(mesh.index_nstat,list(intercept=1)),
                              list(x=covariate)),
                  tag="est")

formula_nstat <- y ~ -1+ intercept + x + f(field,model=spde_nstat)

output_nstat <- inla(formula_nstat,
                  data=inla.stack.data(stack.est_nstat,spde=spde_nstat),
                  family="normal",
                  control.predictor=list(A=inla.stack.A(stack.est_nstat),compute=TRUE))

spde_nstat_prior2 <- inla.spde2.matern(mesh, 
                    B.tau=matrix(cbind(0, 1, 0, sin(pi*mesh$loc[,1])),ncol=4),
                    B.kappa=matrix(c(0, 0, 1, 0),nrow=1,ncol=4),
                    theta.prior.mean=c(0, 2, 0),
                    theta.prior.prec=c(1,1,1))

stack.est_nstat_prior2 <- inla.stack(data=list(y=y_nstat),
                             A=list(A,1),
                             effects=list(c(mesh.index_nstat,list(intercept=1)), list(x=covariate)),
                             tag="est")

formula_nstat_prior2 <- y ~ -1+ intercept + x + f(field,model=spde_nstat_prior2)

output_nstat_prior2 <- inla(formula_nstat_prior2,
                         data=inla.stack.data(stack.est_nstat_prior2,spde=spde_nstat),
                         family="normal",
                         control.predictor=list(A=inla.stack.A(stack.est_nstat_prior2),compute=TRUE))

# *** Code for Figure 6.25
# theta1
myrange.y= c(0,
             max(max(output_nstat$marginals.hyperpar[[2]][,2]),
             max(output_nstat_prior2$marginals.hyperpar[[2]][,2])))
plot(inla.smarginal(output_nstat$marginals.hyperpar[[2]]),t="l",xlab=expression(theta[1]),ylab="",ylim=myrange.y)
lines(inla.smarginal(output_nstat_prior2$marginals.hyperpar[[2]]),lty=2)
abline(v=theta[1])

# theta2
myrange.y= c(0,
             max(max(output_nstat$marginals.hyperpar[[3]][,2]),
             max(output_nstat_prior2$marginals.hyperpar[[3]][,2])))
plot(inla.smarginal(output_nstat$marginals.hyperpar[[3]]),t="l",xlab=expression(theta[2]),ylab="",ylim=myrange.y)
lines(inla.smarginal(output_nstat_prior2$marginals.hyperpar[[3]]),lty=2)
abline(v=theta[2]) 

# theta3
myrange.y= c(0,
             max(max(output_nstat$marginals.hyperpar[[4]][,2]),
             max(output_nstat_prior2$marginals.hyperpar[[4]][,2])))
plot(inla.smarginal(output_nstat$marginals.hyperpar[[4]]),t="l",xlab=expression(theta[3]),ylab="",ylim=myrange.y)
lines(inla.smarginal(output_nstat_prior2$marginals.hyperpar[[4]]),lty=2)
abline(v=theta[3])
# ***

