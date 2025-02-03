library(tidyverse)
library(ggforce)
library(gstat)
library(sf)


dta_full_sf <- st_read(dsn = "swiss/swiss_rain_full.shp") 
dta_sf <- st_read(dsn = "swiss/swiss_rain.shp")
bords_sf <- st_read(dsn = "swiss/bords.shp") 
elev_sf <- st_read(dsn = "swiss/elev.shp") 

dta_sf |> st_coordinates() |> bind_cols(dta_sf) |> ggplot() + 
    geom_sf(data = bords_sf, linewidth = 0.1) +
    geom_circle(aes(x0=X, y0=Y, r=pluies*10, fill = pluies)) +
    scale_fill_viridis_c(name = "Pluies", option = "mako", direction = -1) 


dta_full_sf |> st_coordinates() |> bind_cols(dta_full_sf) |> 
    ggplot() + 
    geom_sf(data = bords_sf, linewidth = 0.1) +
    geom_circle(aes(x0=X, y0=Y, r=pluies*10, fill = pluies)) +
    scale_fill_viridis_c(name = "Pluies", option = "mako", direction = -1) 

elev_sf |> ggplot() + 
    geom_sf(aes(col = elevation)) + 
    geom_sf(data = bords_sf, fill = NA, linewidth = 1.5) + 
    theme_minimal() +
    theme(legend.position = "right") +
    scale_color_viridis_c(name = "Z", option = "mako", transform= "log10")

#------------
# 0. Nuees variographiques à la main
#-------------------------

nrow_dta <- nrow(dta_sf)
ecart2_list <- lapply(1:nrow_dta, 
                      function(i){
    xy <- dta_sf|> st_coordinates()
    x0 = xy[i,1]
    y0 = xy[i,2]
    pluie0 = dta_sf$pluies[i]
    
    dta_sf |> rownames_to_column('Id')  |> 
        mutate(delta_x = xy[,1] - x0, 
                  delta_y = xy[,2]- y0,
                  ecart = pluies - pluie0) |> 
        filter(Id > i) |> 
        mutate( i = i, j = Id, h = sqrt(delta_x^2 +delta_y^2), ecart2 = ecart^2/2) |> 
        select(i,j, h, ecart2)
})
nuees_dta <- do.call(rbind, ecart2_list)
nuees_dta |> st_drop_geometry() |> ggplot()  + geom_point(aes(x=h, y = ecart2)) + xlim(c(0, 120000))


#------------
# 1. Variogramme empirique
#-------------------------

m.d = 120000
vario.cloud = variogram(pluies~1, data = dta_sf, cloud = TRUE)
vario.cloud |> ggplot(aes(x=dist, y =gamma)) + geom_point() + ggtitle("Nuée variographique") 


# 2. Ajustement du variogramme
#-----------------------------
vario.b = variogram(pluies~1, data = dta_sf, alpha=c(0,45,90,135))
ggplot(data = vario.b, aes(x=dist, y =gamma)) + facet_wrap(~dir.hor) + geom_point()


## choix d'un modèle exponentiel
vario.iso = variogram(pluies~1, data = dta_sf)
v.fit = fit.variogram(vario.iso, vgm(model =  "Exp", psill= 15000, range = 50000))
vario_fit_exp <- variogramLine(v.fit, maxdist = m.d)

ggplot(data = vario.iso, aes(x=dist, y =gamma))  + geom_point() + 
    geom_line(data= vario_fit_exp , aes(x=dist, y = gamma))

   
# 3. Krigeage
#------------
## transform in sf object
Kfull <- krige(formula = pluies ~ 1, locations = dta_sf, 
               newdata = dta_full_sf,  model = v.fit)
Kfull |> st_join(dta_full_sf) |> 
    mutate( error = pluies- var1.pred) |> 
    ggplot() + geom_sf(aes(col = error), size =4) +
    theme_minimal() +
    theme(legend.position = "right")  +
    scale_color_viridis_c(name = "Erreur", option = "plasma")

RMSEP = Kfull |> st_join(dta_full_sf) |>  mutate( error2 = (pluies- var1.pred)^2) |> 
    summarise(rmsep = sqrt(sum(error2)/n())) |> pull(rmsep)

# 4. Krigeage sur une grille
#------------ 

grid <- st_make_grid(bords_sf, n= c(50, 50))
grid |>  st_intersection(bords_sf) |> ggplot() + geom_sf() 

grid |> ggplot() + geom_sf()

dta_sf |>  st_coordinates() |> bind_cols(dta_sf)  |> ggplot()  +
    scale_fill_viridis_c(name = "Pluies", option = "mako", direction = -1)  +
    geom_sf(data = bords_sf, linewidth = 0.1) +
    geom_sf(data= grid) +
    geom_circle(aes(x0=X, y0=Y, r=pluies*10, fill = pluies))


Kfull <- krige(formula = pluies ~ 1, locations = dta_sf, newdata = grid,  model = v.fit)


Kfull |> ggplot() + geom_sf(aes(fill = var1.pred)) + 
    geom_sf(data = dta_sf, aes(col = pluies), size = 3) +
    scale_color_viridis_c(name = "Z", option = "mako", transform= "log10")+
    scale_fill_viridis_c( option = "mako", transform= "log10") 

## représenter la variance de prédiction
Kfull |> ggplot() + 
    geom_sf(aes(fill = var1.var)) +  
    geom_sf(data = dta_sf, col = "red") +
    scale_color_viridis_c( option = "mako", transform= "log10") 


# 5. A votre tour
#------------
## 5.1 Tester différentes formes de variogram


# 5.2 Ajout de la covariable
#-----------------

res <- do.call('rbind', 
               lapply(
                   split(dta_sf, 1:nrow(dta_sf)),
                   function(x) {
                       st_join(x, elev_sf, join = st_nearest_feature)
                   }))


