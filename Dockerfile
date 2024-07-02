FROM rocker/rstudio:4.1.3

WORKDIR /home

RUN apt update && apt upgrade -y

# R dependencies
RUN apt install -y libcurl4-openssl-dev libssl-dev libpng-dev libxml2-dev libfontconfig1-dev libhdf5-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libtiff5-dev libjpeg-dev libcairo2-dev pandoc libxt6

# Install Rstudio for Ubuntu focal (20.04)
RUN apt install -y libnss3 libatk1.0-0 libatk-bridge2.0-0 libcups2 libdrm2 libgtk-3-0 libgbm1 libasound2 curl libcurl4-gnutls-dev hdf5-tools

# RUN wget https://download1.rstudio.org/electron/focal/amd64/rstudio-2023.12.0-369-amd64.deb
# RUN apt install -y ./rstudio-2023.12.0-369-amd64.deb
# RUN rm ./rstudio-2023.12.0-369-amd64.deb
# RUN echo "alias rstudio='rstudio --no-sandbox'" > ~/.bash_aliases

RUN Rscript -e "local({r <- getOption('repos'); \
       r['CRAN'] <- 'https://cran.biotools.fr'; \
       options(repos=r);\
}); \
install.packages('remotes'); \
install.packages('BiocManager'); \
remotes::install_github('mojaveazure/seurat-disk'); \
remotes::install_github('ekernf01/DoubletFinder', force = T); \
remotes::install_github('cancerbits/DElegate');\
remotes::install_version('scCustomize', upgrade='never', version='1.1.1');\
remotes::install_version('Seurat', upgrade='never', version='4.3.0'); \
remotes::install_version('SeuratObject', upgrade='never', version='4.1.3'); \
remotes::install_version('dplyr', upgrade='never', version='1.1.4'); \
remotes::install_version('kableExtra', upgrade='never', version='1.3.4'); \
remotes::install_version('RcppTOML', upgrade='never', version='0.2.2'); \
remotes::install_version('stringr', upgrade='never', version='1.5.0'); \
remotes::install_version('optparse', upgrade='never', version='1.7.3'); \
remotes::install_version('cowplot', upgrade='never', version='1.1.1'); \
remotes::install_version('tibble', upgrade='never', version='3.2.1'); \
remotes::install_version('ggside', upgrade='never', version='0.2.2'); \
remotes::install_version('ROCR', upgrade='never', version='1.0.11'); \
remotes::install_version('KernSmooth', upgrade='never', version='2.23.21'); \
remotes::install_version('fields', upgrade='never', version='14.1'); \
remotes::install_version('viridis', upgrade='never', version='0.6.3'); \
remotes::install_version('viridisLite', upgrade='never', version='0.4.2'); \
remotes::install_version('patchwork', upgrade='never', version='1.1.2'); \
remotes::install_version('spam', upgrade='never', version='2.9.1'); \
remotes::install_version('scales', upgrade='never', version='1.3.0'); \
remotes::install_version('DT', upgrade='never', version='0.27'); \
remotes::install_version('Matrix', upgrade='never', version='1.6.3'); \
remotes::install_version('tidyr', upgrade='never', version='1.3.0'); \
remotes::install_version('RColorBrewer', upgrade='never', version='1.1.3'); \
remotes::install_version('ggplot2', upgrade='never', version='3.4.4'); \
remotes::install_version('reticulate', upgrade='never', version='1.35'); \
remotes::install_version('gprofiler2', upgrade='never', version='0.2.1'); \
remotes::install_version('data.table', upgrade='never', version='1.14.8'); \
remotes::install_version('knitr', upgrade='never', version='1.43'); \
remotes::install_version('rmdformats', upgrade='never', version='1.0.3'); \
remotes::install_version('R.utils', upgrade='never', version='2.12.2'); \
remotes::install_version('gghighlight', upgrade='never', version='0.4.0'); \
BiocManager::install('DESeq2', force=TRUE); \
BiocManager::install('fgsea', force=TRUE); \
BiocManager::install('msigdbr', force=TRUE);"
# BiocManager::install('biomaRt'); \


# python dependencies
RUN apt install -y git libncurses-dev libgdbm-dev libz-dev tk-dev libsqlite3-dev libreadline-dev liblzma-dev libffi-dev libbz2-dev libgdbm-compat-dev

# manually build python to ensure both version and shared libraries for reticulate
RUN git clone --depth 1 --branch v3.10.13 https://github.com/python/cpython /usr/bin/cpython && \
cd /usr/bin/cpython && ./configure --enable-shared --enable-optimizations --prefix /usr/local LDFLAGS="-Wl,-rpath /usr/local/lib" && \
make -C /usr/bin/cpython && make install -C /usr/bin/cpython

RUN ln -s /usr/local/bin/python3.10 /usr/bin/python
RUN ln -s /usr/local/bin/pip3 /usr/bin/pip

RUN pip install pandas numpy meld scanpy matplotlib scprep scikit-learn leidenalg


RUN chgrp -R 998 /usr /var /home && chmod -R g=u /usr /var /home

CMD ["/bin/bash"]

# To run the image connecting to the project's directory:
# docker run -it --mount type=bind,source=/path/to/your/project,target=/home -e SCDP_PATH_ROOT="/home/" scdatapipeline
# Add the following arguments to connect to the host's X server (might be a security vulnerability, use only if you want to use Rstudio):
# -v /tmp/.X11-unix:/tmp/.X11-unix  -e DISPLAY="${DISPLAY}" --network="host"