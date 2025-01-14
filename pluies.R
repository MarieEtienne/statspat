library(tidyverse)
library(ggforce)
library(geoR)
library(fields)
library(gstat)
library(sf)
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


dta_full_sf <- st_as_sf(dta_full, coords = c("x", "y"))
dta_sf <- st_as_sf(dta, coords = c("x", "y"))

bords <- read.table('/home/metienne/Cours/ENSAI/spatial/docs/Cours_liliane/geostat/bords.txt',col.names = c('x','y'))|> 
    mutate(x = (x-xmin)/1000, y = (y- ymin)/1000)

elevation <- read.table('/home/metienne/Cours/ENSAI/spatial/docs/Cours_liliane/geostat/surfdem.grd',skip=6)

taille = dta_full$pluies/100

dta |> ggplot() + geom_circle(aes(x0=x, y0=y, r=pluies/100, fill = pluies)) +
    scale_fill_viridis_c(name = "Pluies", option = "mako", direction = -1) +
    geom_path(data = bords, aes(x=x, y=y))


dta_full |> ggplot() + geom_circle(aes(x0=x, y0=y, r=pluies/100, fill = pluies)) +
    scale_fill_viridis_c(name = "Pluies", option = "mako", direction = -1) +
    geom_path(data = bords, aes(x=x, y=y))

### remettre dans le même référentiel
xllcorner  =   (-185556.375 -xmin)/1000
yllcorner  =   (-127261.5234375 -ymin)/1000
cellsize    =  1009.975/1000

xelev = seq(xllcorner,cellsize*375+xllcorner,cellsize)
yelev = seq(yllcorner,cellsize*252+yllcorner,cellsize)
elev <- expand.grid(xelev,yelev) |> mutate(elevation = as.numeric(t(as.matrix(elevation))))
elev |> ggplot() + 
    geom_point(aes(x=Var1, y = Var2, col = elevation)) + 
    geom_path(data = bords, aes(x=x, y=y), linewidth = 1) +
    theme_minimal() +
    theme(legend.position = "right") +
    scale_color_viridis_c(name = "Z", option = "mako", transform= "log10")

#------------
# variogrammemin(elevation)
#------------
# 1. Variogramme empirique
#-------------------------

m.d = 200                  # distance maximale essayer 150 et 500
interv = seq(0,m.d,20)    # intervalles
p.m = 20                     # nombre minimal de paire regarder l'impact

dta_sf <- st_as_sf(dta, coords = c("x", "y"))
vario.cloud = variogram(pluies~1, data = dta_sf, cloud = TRUE)
vario.cloud |> ggplot(aes(x=dist, y =gamma)) + geom_point() + ggtitle("Nuée variographique") 


# 2. Ajustement du variogramme
#-----------------------------
vario.b = variogram(pluies~1, data = dta_sf, alpha=c(0,45,90,135))
ggplot(data = vario.b, aes(x=dist, y =gamma)) + facet_wrap(~dir.hor) + geom_point()

vario.iso = variogram(pluies~1, data = dta_sf)
v.fit = fit.variogram(vario.iso, vgm(model =  "Exp", 15000, 60))
vario_fit_exp <- variogramLine(v.fit, maxdist = m.d)

ggplot(data = vario.iso, aes(x=dist, y =gamma))  + geom_point() + 
    geom_line(data= vario_fit_exp , aes(x=dist, y = gamma))



   
# 3. Krigeage
#------------
## transform in sf object
dta_full_sf <- st_as_sf(dta_full_sf)
Kfull <- krige(formula = pluies ~ 1, locations = dta_sf, newdata = dta_full_sf,  model = v.fit)
Kfull |> st_join(dta_full_sf) |> 
    mutate( error = pluies- var1.pred) |> ggplot() + geom_sf(aes(col = error)) +
    theme_minimal() +
    theme(legend.position = "right") +
    scale_color_viridis_c(name = "Température", option = "plasma")
RMSEP = Kfull |> st_join(dta_full_sf) |>  mutate( error2 = (pluies- var1.pred)^2) |> 
    summarise(rmsep = sum(error2)) |> pull(rmsep)

# krigeage sur une grille
grid <- st_make_grid(dta_full_sf, n= c(50, 50))
dta |> ggplot()  +
    scale_fill_viridis_c(name = "Pluies", option = "mako", direction = -1)  +
    geom_sf(data= grid) +
    geom_path(data = bords, aes(x=x, y=y))+ 
    geom_circle(aes(x0=x, y0=y, r=pluies/100, fill = pluies))
Kfull <- krige(formula = pluies ~ 1, locations = dta_sf, newdata = grid,  model = v.fit)
Kfull |> ggplot() + geom_sf(aes(fill = var1.pred)) 


## Tester différentes formes de variagoram

# Explorer l'anisotropie
#----------------- 

# aJOUT DE LA COVARIABLE elevation
#-----------------

# Simulations conditionnelles
#----------------------------
