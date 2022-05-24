# ImmuneSignatures2 Docker Image

This directory contains sources for building the docker image used in the ImmuneSpace servers for running all ImmuneSignatures2 code. It can also be used to mimic the server environment if you want to run ImmuneSignatures2 code on your own machine.

## Quickstart

Option 1: Pull from Dockerhub

```
docker pull rglab/immunesignatures2
docker run --rm -ti rglab/immunesignatures2:latest
```

Note that data is not included with the docker image, so you will have to download it from  https://www.immunespace.org/is2.url

Option 2: build yourself. 

1. Build the image. 

```
docker build -t immunesignatures2:latest .
docker run --rm -ti rglab/immunesignatures2:latest
```

## How it's set up

Based on a combination of the recommended [LabKey sandboxed docker setup](https://github.com/LabKey/docker-rstudio/tree/develop/images/labkey/rsandbox-ver) and [renv docker setup](https://rstudio.github.io/renv/articles/docker.html). Additional changes and customizations are noted in comments in the scripts. 

```
docker
├── Dockerfile
├── README.md
├── install.R
└── renv.lock
```

## Installation on ImmuneSpace server

For ImmuneSpace developers: See Notion documentation for [server setup](https://www.notion.so/rglab/ImmSig2-server-setup-446532d9490841aab1bab51ecb3c8f10) and [updating docker image with new R packages](https://www.notion.so/rglab/Updating-ImmSig-docker-image-with-new-packages-252883ef06b04f06954c4d889d943b80)


