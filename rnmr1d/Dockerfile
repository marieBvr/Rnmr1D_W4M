FROM r-base:latest

MAINTAINER "Daniel Jacob" daniel.jacob@u-bordeaux.fr

# Install dependencies
RUN install.r Rcpp doParallel

# Install Rnmr1D package
RUN mkdir /home/docker/src
ADD . /home/docker/src
RUN R CMD INSTALL /home/docker/src
RUN cp /usr/local/lib/R/site-library/Rnmr1D/exec/Rnmr1D /usr/local/bin/ \
    && chmod 755 /usr/local/bin/Rnmr1D

# Volume as the root directory that will contain the output data
# the data will be stored within the directory /data/<DATAID>
RUN mkdir -p /data
VOLUME /data

# Volume as the temporary directory that should contain the input data
# i.e. the raw bruker directories, and the INI parameter file
RUN mkdir -p /srv
VOLUME /srv
WORKDIR /srv

# Rnmr1D entrypoint
ENTRYPOINT ["/usr/bin/Rscript","/usr/local/bin/Rnmr1D"]
