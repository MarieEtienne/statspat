
FROM rocker/geospatial:4.4
RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update \
 && apt-get install -y pandoc 
RUN R -e "install.packages('remotes')"
RUN R -e "install.packages('purrr')" # map function
ENV R_CRAN_WEB="https://cran.rstudio.com/"
RUN R -e "install.packages(c('geostatsp', 'parallel', 'ggpubr', 'imager', 'RefManageR', 'wesanderson', 'plotly'))" # GET function
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
  mercurial gdal-bin libgdal-dev gsl-bin libgsl-dev \
  libc6-i386

