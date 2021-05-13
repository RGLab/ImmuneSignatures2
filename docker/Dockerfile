ARG VERSION=4.0.2
FROM rocker/r-ver:${VERSION}


RUN apt-get update && \
    apt-get install -y --no-install-recommends apt-utils && \
    apt-get install -y --no-install-recommends dos2unix nano procps libxml2-dev xvfb sudo && \
    apt-get install -y --no-install-recommends  \
    curl libcurl4-openssl-dev \
    libgd-dev libcairo2 libcairo2-dev libxt-dev pandoc \
    libssl-dev openssl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

## Set a default user. Available via runtime flag `--user docker`
## Add user to 'staff' group, granting them write privileges to /usr/local/lib/R/site.library
## User should also have & own a home directory.
RUN useradd docker \
    && mkdir /home/docker \
    && mkdir /home/docker/R_Sandbox \
    && chown -R docker /home/docker \
    && adduser docker staff \
    && adduser docker sudo \
    && echo 'docker ALL=(root) NOPASSWD:ALL' >> /etc/sudoers \
    && echo 'docker ALL=(root) NOPASSWD:ALL' >> /etc/sudoers.d/docker \
    && chmod 0440 /etc/sudoers.d/docker \
    && chmod -R 775 /home/docker

# Pass in your github pat when building so ImmuneSignatures2 private repo can install 
ENV GITHUB_PAT=${GITHUB_PAT}

# Use renv to install R packages
ENV RENV_VERSION 0.11.0
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

## do not use apt-get install r-cran-* to install R packages on this stack
## see rocker/r-ver (see https://hub.docker.com/r/rocker/r-ver/)
# COPY renv.lock /home/docker/
# WORKDIR /home/docker
# RUN dos2unix install.R && \
#     R --no-site-file -f install.R &&  \
#     chmod -R 755 /usr/local/lib/R/site-library

CMD ["su", "docker", "-c", "/bin/bash"]
