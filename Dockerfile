FROM --platform=linux/amd64 rocker/shiny:4.6.1 

# system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    pkg-config \
    curl \
    zlib1g-dev \
    libuv1 \
    libuv1-dev \
    libnetcdf-dev \
    libhdf5-dev \
    libmagick++-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    openssl \
    libxml2-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libglpk40 \
    libglpk-dev \
    && rm -rf /var/lib/apt/lists/*
# CRAN packages
RUN install2.r --error --skipinstalled \
    Amelia \
    broom \
    bsicons \
    bslib \
    callr \
    curl \
    dplyr \
    DT \
    ggnewscale \
    ggplot2 \
    ggtext \
    ggupset \
    httr \
    knitr \
    later \
    minpack.lm \
    pheatmap \
    plotly \
    pracma \
    processx \
    purrr \
    RColorBrewer \
    readr \
    readxl \
    shiny \
    shinybusy \
    shinycssloaders \
    shinydashboard \
    shinyjs \
    stringr \
    tibble \
    tidyr \
    tidyverse \
    VIM

# Bioconductor installer
RUN install2.r --error --skipinstalled BiocManager remotes

# Bioconductor packages
RUN R -e "BiocManager::install(c( \
    'Biobase', \
    'BiocParallel', \
    'Biostrings', \
    'IRanges', \
    'GenomicAlignments', \
    'circlize', \
    'clusterProfiler', \
    'ComplexHeatmap', \
    'enrichplot', \
    'GO.db', \
    'impute', \
    'limma' \
), ask = FALSE, update = FALSE)"

#Uniprot
RUN R -e "install.packages('UniprotR')"

# DEP 
RUN R -e "remotes::install_github( \
  'arnesmits/DEP', \
  ref = 'b425d8d0db67b15df4b8bcf87729ef0bf5800256', \
  dependencies = TRUE \
)" && \
R -e "library(DEP)"

# TPP and NPARC (if not already pulled as dependencies)
RUN R -e "BiocManager::install(c('TPP','NPARC'), ask = FALSE, update = FALSE)"

# Additional packages discovered during testing
RUN install2.r --error --skipinstalled \
    ggridges

# Copy app
RUN rm -rf /srv/shiny-server/*
COPY App/Core /srv/shiny-server/

RUN chown -R shiny:shiny /srv/shiny-server && \
    chmod -R u+rwX /srv/shiny-server


EXPOSE 3838

USER shiny

CMD ["/usr/bin/shiny-server"]