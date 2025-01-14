### simulation of a spatial gaussian random fields with exponential covariance 

library(fields) # Pour la fonction rdist
# library(sp)    # Pour les objets spatiaux
library(ggplot2)
library(gstat)


# 1. Simulation d'un champ spatial
set.seed(123)

# Définir les coordonnées spatiales
n <- 500  # Nombre de points
coords <- data.frame(x=runif(n, 0, 10), y = runif(n, 0, 10))

# Définir les paramètres de la covariance exponentielle
sigma2 <- 2.5    # Variance
phi <- 2       # Portée spatiale
nugget <- 0  # Effet pépite

# Matrice de distances
dist_matrix <- rdist(coords)

# Matrice de covariance
cov_matrix <- sigma2 * exp(-dist_matrix / phi) + diag(nugget, n)

# Simulation des valeurs spatiales
z <- t(chol(cov_matrix)) %*% rnorm(n)
dta <- coords |> add_column(z=z)

# Map
ggplot(dta) + geom_point(aes(x=x, y= y , col=z), size= 2) + 
  theme_minimal() +
  theme(legend.position = "right") +
  scale_color_viridis_c(name = "z", option = "plasma", limits=c(-3.5, 3.5)) 
  
  
## Nuees variagraphiques
## calcul de la nuée
dta_vario <- variogram(z~1, loc=~x+y, data=dta, cloud=TRUE)
ggplot(data.frame(dta_vario),aes(x=dist,y=gamma))+ geom_point() 


## calcul de la nuée
dta_vario <- variogram(z~1, loc=~x+y, data=dta )
ggplot(data.frame(dta_vario),aes(x=dist,y=gamma))+ geom_point() + ylim(c(0, 2.5))
