FROM rocker/geospatial:4.6
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      mercurial gsl-bin libgsl-dev libc6-i386 libglpk40 libglpk-dev \
      ghostscript fontconfig \
 && apt-get clean && rm -rf /var/lib/apt/lists/*

# --- LaTeX pour les versions PDF (beamer) : complète la TeX Live de rocker ---
RUN tlmgr install \
      xetex beamer pgf translator tcolorbox environ fontspec unicode-math \
      lm lm-math tex-gyre xcolor etoolbox iftex amsmath tools caption \
      bookmark xurl booktabs fancyvrb framed upquote microtype parskip \
      footnotehyper \
 && fmtutil-sys --all

# --- R ---
ENV R_CRAN_WEB="https://cran.rstudio.com/"
RUN R -e "install.packages(c('remotes', 'purrr'))"
RUN R -e "install.packages(c('geostatsp', 'ggpubr', 'imager', 'RefManageR', 'wesanderson', 'plotly'))"
RUN R -e "install.packages('INLA',repos=c(getOption('repos'),INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)"
RUN R -e "install.packages(c('maps', 'fields', 'igraph', 'ggthemes'))"
RUN R -e "install.packages(c('gstat', 'ggforce', 'usmap', 'inlabru'))"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/SpatialExtremes/SpatialExtremes_2.1-0.tar.gz', repos = NULL, type = 'source')"
