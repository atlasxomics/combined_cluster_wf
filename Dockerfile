FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:5303-main
RUN apt-get update -y
RUN apt-get install -y gdebi-core 
RUN apt install -y aptitude
RUN aptitude install -y libjpeg-dev
RUN apt-get update -y

# Install R
RUN apt-get update -y && \
    apt-get install -y \
        r-base \
        r-base-dev \
        apt-transport-https \
        build-essential \
        gfortran \
        libhdf5-dev \
        libatlas-base-dev \
        libbz2-dev \        
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libgit2-dev \
        libgsl-dev \
        libicu-dev \
        liblzma-dev \
        libpango-1.0-0 \
        libpangocairo-1.0-0 \
        libpcre3-dev \
        libssl-dev \
        libtcl8.6 \
        libtiff5 \
        libtk8.6 \
        libxml2-dev \
        libxt-dev \
        libx11-dev \
        libtiff-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libglpk-dev \
        locales \
        make \
        pandoc \
        tzdata \
        vim \
        wget \
        zlib1g-dev \
        r-cran-rjava

RUN echo "alias ll='ls -l --color=auto'" >> .bashrc

# Fix systemd conflict with timedatectl
RUN echo "TZ=$( cat /etc/timezone )" >> /etc/R/Renviron.site

# Have to install devtools, cairo like this; see https://stackoverflow.com/questions/20923209
RUN apt-get install -y r-cran-devtools libcairo2-dev

# Install packages
RUN R -e "install.packages(c('Cairo', 'Matrix', 'knitr', 'patchwork', 'BiocManager', 'gridExtra', 'dplyr', 'tibble', 'hdf5r', 'stringer', 'rjson', 'rmarkdown', 'purrr', 'harmony', 'pheatmap', 'RColorBrewer', 'ggrepel'))"
RUN R -e "devtools::install_github('GreenleafLab/ArchR', ref='master', repos = BiocManager::repositories())"
RUN R -e "devtools::install_github('GreenleafLab/chromVARmotifs')"
RUN R -e "library('ArchR'); ArchR::installExtraPackages()"

# Upgrade R to version 4.3.0
RUN wget https://cran.r-project.org/src/base/R-4/R-4.3.0.tar.gz
RUN tar zxvf R-4.3.0.tar.gz
RUN cd R-4.3.0 && ./configure --enable-R-shlib
RUN cd R-4.3.0 && make && make install


# Install more R packages
RUN R -e "install.packages(c('pkgconfig', 'munsell', 'zip', 'zoo', 'xtable', 'listenv', 'lazyeval', 'bit64', 'rJava', 'labeling'), repos = 'http://cran.us.r-project.org')"
RUN R -e "ArchR::installExtraPackages()"
RUN R -e "BiocManager::install(version = '3.17',ask = FALSE)"
RUN R -e "BiocManager::install('EnhancedVolcano')"
RUN R -e "BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')"
RUN R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')"

# numpy needed to be install before macs2 v-2.2.6
RUN python3 -m pip install numpy
RUN python3 -m pip install macs2==2.2.6


RUN apt-get update -y
RUN apt-get install -y libmagick++-dev
RUN apt-get install -y libgdal-dev
RUN R -e "install.packages(c('Seurat'), dependencies = TRUE, repos = 'http://cran.us.r-project.org')"

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
RUN python3 -m pip install --upgrade latch
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root