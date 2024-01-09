FROM rocker/r-ver:4.1.3

WORKDIR /home

RUN apt update && apt upgrade -y

# R dependencies
RUN apt install -y libcurl4-openssl-dev libssl-dev libpng-dev libxml2-dev libfontconfig1-dev libhdf5-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libtiff5-dev libjpeg-dev libcairo2-dev pandoc libxt6

RUN Rscript -e "local({r <- getOption('repos'); \
       r['CRAN'] <- 'https://cran.biotools.fr'; \
       options(repos=r);\
}); \
install.packages('remotes'); \
install.packages('BiocManager'); \
remotes::install_version('Seurat', version='4.3.0'); \
remotes::install_version('SeuratObject', version='4.1.3'); \
remotes::install_version('dplyr', version='1.1.4'); \
remotes::install_version('kableExtra', version='1.3.4'); \
remotes::install_version('RcppTOML', version='0.2.2'); \
remotes::install_version('stringr', version='1.5.0'); \
remotes::install_version('optparse', version='1.7.3'); \
remotes::install_version('cowplot', version='1.1.1'); \
remotes::install_version('tibble', version='3.2.1'); \
remotes::install_version('ggside', version='0.2.2'); \
remotes::install_version('ROCR', version='1.0.11'); \
remotes::install_version('KernSmooth', version='2.23.21'); \
remotes::install_version('fields', version='14.1'); \
remotes::install_version('viridis', version='0.6.3'); \
remotes::install_version('viridisLite', version='0.4.2'); \
remotes::install_version('patchwork', version='1.1.2'); \
remotes::install_version('spam', version='2.9.1'); \
remotes::install_version('scales', version='1.3.0'); \
remotes::install_version('DT', version='0.27'); \
remotes::install_version('Matrix', version='1.5.4.1'); \
remotes::install_version('tidyr', version='1.3.0'); \
remotes::install_version('RColorBrewer', version='1.1.3'); \
remotes::install_version('ggplot2', version='3.4.4'); \
remotes::install_version('reticulate', version='1.28'); \
remotes::install_version('knitr', version='1.42'); \
BiocManager::install('biomaRt'); \
remotes::install_github('mojaveazure/seurat-disk'); \
remotes::install_github('ekernf01/DoubletFinder', force = T); \
remotes::install_version('scCustomize', version='1.1.1')"

# python dependencies
RUN apt install -y git libncurses-dev libgdbm-dev libz-dev tk-dev libsqlite3-dev libreadline-dev liblzma-dev libffi-dev libbz2-dev libgdbm-compat-dev

# manually build python to ensure both version and shared libraries for reticulate
RUN cd /usr/bin
RUN git clone --depth 1 --branch v3.10.13 https://github.com/python/cpython
RUN cd cpython
RUN ./configure --enable-shared --enable-optimizations --prefix /usr/local LDFLAGS="-Wl,-rpath /usr/local/lib"
RUN make
RUN make install

RUN pip install pandas numpy meld scanpy matplotlib scprep scikit-learn leidenalg

CMD ["bash"]