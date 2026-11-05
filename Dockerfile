
FROM rocker/geospatial:4.6
RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update 
RUN R -e "install.packages('remotes')"
RUN R -e "install.packages('purrr')" # map function
ENV R_CRAN_WEB="https://cran.rstudio.com/"
RUN R -e "install.packages(c('geostatsp', 'ggpubr', 'imager', 'RefManageR', 'wesanderson', 'plotly'))" # GET function
RUN R -e "install.packages('INLA',repos=c(getOption('repos'),INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)"
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
  mercurial gsl-bin libgsl-dev \
  libc6-i386 \
  libglpk40 libglpk-dev 
RUN R -e "install.packages(c('maps','fields'))"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/SpatialExtremes/SpatialExtremes_2.1-0.tar.gz', repos = NULL, type = 'source')"
RUN R -e "install.packages('igraph')"
RUN R -e "install.packages('ggthemes')"


