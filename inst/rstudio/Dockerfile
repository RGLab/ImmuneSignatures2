ARG version
FROM immunespace/rstudio:${version}
ARG github_pat

COPY install.R /
RUN dos2unix /install.R; Rscript --vanilla install.R ${github_pat}

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
